import commands as cmd
import argparse
import ROOT
from array import array
from math import *

#Some useful functions

def DeltaPhi(phi1,phi2):
    PHI=phi1-phi2
    if PHI >= pi:
        PHI -= 2*pi
    elif PHI < -1*pi:
        PHI += 2*pi
    return PHI

def MT(Lpt,MET,Lphi,METphi):
    return sqrt(2*Lpt*MET*(1-cos(DeltaPhi(METphi,Lphi))))

#Global parser

def GlobalAnalysisParser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--DIR', help='Directory to read input')
    parser.add_argument('--XS', help='Sample cross section')
    parser.add_argument('--LUMI', help='Luminosity [pb-1]')
    parser.add_argument('--MLL', help='Minimum mll')
    parser.add_argument('--MET', help='Minimum met')
    parser.add_argument('--OUT', help='ROOT output file')
    parser.add_argument('--ANA', help='Analysis selector')
    parser.add_argument('--DELPHES', help='Delphes libraries location')
    return parser.parse_args()

args = GlobalAnalysisParser()

###########Exceptions need to be implemented requiring the parsing of mandatory arguments

#Delphes Initialization

def DelphesInit(DelphesPath=args.DELPHES, DelphesTreeName="Delphes", SamplesDir=args.DIR, JetBranch="", MetBranch="", MuonBranch="", ElectronBranch=""):
    ROOT.gSystem.AddDynamicPath(DelphesPath)
    ROOT.gSystem.Load("libDelphes");
    try:
        ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
        ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
        print "Delphes classes imported"
    except:
        pass

    TreeName=DelphesTreeName
    DataChain=ROOT.TChain(TreeName)
    DataChain.Add(SamplesDir+"*.root")

    # Create object of class ExRootTreeReader
    treeReader = ROOT.ExRootTreeReader(DataChain)

    BranchDictionary = {}
    
    if len(JetBranch)!=0:
        branchJet = treeReader.UseBranch(JetBranch)
        BranchDictionary["Jet"] = branchJet
    if len(MetBranch)!=0:
        branchMET = treeReader.UseBranch(MetBranch)
        BranchDictionary["Met"] = branchMET
    if len(MuonBranch)!=0:
        branchMuon = treeReader.UseBranch(MuonBranch)
        BranchDictionary["Muon"] = branchMuon
    if len(ElectronBranch)!=0:
        branchElectron = treeReader.UseBranch(ElectronBranch)
        BranchDictionary["Electron"] = branchElectron

    NumberOfEventsToCheck=DataChain.GetEntries()

    print "--------------------------------------> ", NumberOfEventsToCheck, "events have been loaded and are ready for analysis!!!!!!!!!!!"

    return NumberOfEventsToCheck, treeReader, BranchDictionary

#ROOT file saver

def BasketFile(FileName=args.OUT, OpenOption="recreate"):
    try:
        if args.OUT is None: raise NameError('No file name was declared, please do!')
    except NameError:
        raise
    return ROOT.TFile(FileName, OpenOption)
    

#Lumi weighting function

def LumiWeighter(TotalEvts, Evts, Lumi=float(args.LUMI), xs=float(args.XS)):
    Weight=(Lumi*xs/TotalEvts)
    return Evts*Weight, sqrt(Evts)*Weight, Weight

#Weight histogram

def WeightHisto(Weight,Histo):
    Histo.Sumw2()
    Histo.Scale(Weight)
    print "Bin range,", "Value"
    for i in xrange(Histo.GetNbinsX()):
        print str(Histo.GetBinLowEdge(i+1))+"-"+str(Histo.GetBinLowEdge(i+1)+Histo.GetBinWidth(i+1))+",", Histo.GetBinContent(i+1)

#Shortcuts functions and variables

MuonMass = 0.10566

ElectronMass = 0.000511

def Np(ABranch): return ABranch.GetEntries()

def PT(ABranch,index): return ABranch.At(index).PT

def ETA(ABranch,index): return ABranch.At(index).Eta

def PHI(ABranch,index): return ABranch.At(index).Phi

def GetParticle(ABranch,index,mass):
    Particle = ROOT.TLorentzVector(0,0,0,0)
    Particle.SetPtEtaPhiM(PT(ABranch,index),ETA(ABranch,index),PHI(ABranch,index),mass)
    return Particle

def Veto(ABranch,PTCut=0,ETACut=5,ETAgap=[5.0,-5.0]):
    IN=Np(ABranch)
    Nparticles=0
    NoParticlesInTheEvent = False
    if IN > 0:
        for i in xrange(IN):
            if PT(ABranch,i)>PTCut and abs(ETA(ABranch,i))<ETACut:
                if abs(ETA(ABranch,i))<ETAgap[0] or abs(ETA(ABranch,i))>ETAgap[1]:
                    Nparticles+=1
    if Nparticles == 0:
        NoParticlesInTheEvent = True
    return NoParticlesInTheEvent


#Information Functions

def DisplayFinalInfo(Evtsi, Evtsf, Histo=None, MyFile=None):
    E, DeltaE, W = LumiWeighter(Evtsi, Evtsf)
    
    print "Raw number of events before full selection:", Evtsi
    print "Raw number of events that passed full selection:", Evtsf, "+-", sqrt(Evtsf)
    print "Used cross section", args.XS, "pb"
    print "Weighted number of events that passed full selection:", E, "+-", DeltaE

    if Histo is None: return

    WeightHisto(W,Histo)

    if args.OUT is not None and MyFile is not None:
        Histo.Write()
        MyFile.Close()
