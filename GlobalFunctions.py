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
    parser.add_argument('--OUT', help='ROOT output file')
    parser.add_argument('--ANA', help='Analysis selector')
    parser.add_argument('--DELPHES', help='Delphes libraries location')
    return parser.parse_args()

args = GlobalAnalysisParser()

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

def LumiWeighter(Lumi=float(args.LUMI), xs=float(args.XS), TotalEvts, Evts):
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

def Nl(ABranch): return ABranch.GetEntries()

def PT(ABranch,index): return ABranch.At(index).PT

def ETA(ABranch,index): return ABranch.At(index).Eta

def PHI(ABranch,index): return ABranch.At(index).Phi

def GetParticle(ABranch,index,mass):
    Particle = ROOT.TLorentzVector(0,0,0,0)
    Particle.SetPtEtaPhiM(PT(ABranch,index),ETA(ABranch,index),PHI(ABranch,index),mass)

