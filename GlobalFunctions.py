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
    numberOfEntries = treeReader.GetEntries()

    if len(JetBranch)!=0: branchJet = treeReader.UseBranch(JetBranch)
    if len(MetBranch)!=0: branchMET = treeReader.UseBranch(MetBranch)
    if len(MuonBranch)!=0: branchMuon = treeReader.UseBranch(MuonBranch)
    if len(ElectronBranch)!=0: branchElectron = treeReader.UseBranch(ElectronBranch)

    NumberOfEventsToCheck=DataChain.GetEntries()

    print "--------------------------------------> ", NumberOfEventsToCheck, "events have been loaded and are ready for analysis!!!!!!!!!!!"

    return NumberOfEventsToCheck

def LumiWeighter(Lumi=float(args.LUMI), xs=float(args.XS), TotalEvts,Evts):
    Weight=(Lumi*xs/TotalEvts)
    return Evts*Weight, sqrt(Evts)*Weight
