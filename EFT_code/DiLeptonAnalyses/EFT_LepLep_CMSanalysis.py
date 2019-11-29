import commands as cmd
import argparse
import ROOT
from array import array
from math import *

parser = argparse.ArgumentParser()
parser.add_argument('--DIR', help='Directory to read input')
parser.add_argument('--XS', help='Sample cross section')
parser.add_argument('--MLL', help='Minimum mll')
parser.add_argument('--OUT', help='ROOT output file')
parser.add_argument('--ANA', help='Analysis selector')
parser.add_argument('--DELPHES', help='Delphes libraries location')
args = parser.parse_args()

def DeltaPhi(phi1,phi2):
    PHI=phi1-phi2
    if PHI >= pi:
        PHI -= 2*pi
    elif PHI < -1*pi:
        PHI += 2*pi
    return PHI

def MT(Lpt,MET,Lphi,METphi):
    return sqrt(2*Lpt*MET*(1-cos(DeltaPhi(METphi,Lphi))))

#Delphes_Path="/cms/Jose/EFT_MG_PRL_version/MG5_aMC_v2_6_1/Delphes/"
#Delphes_Path="/home/joser/Dropbox/Vandy/EFT/ElectronNeutrino/Delphes-3.4.1/"
#ROOT.gSystem.AddDynamicPath(Delphes_Path)
ROOT.gSystem.AddDynamicPath(args.DELPHES)
#ROOT.gROOT.ProcessLine("#include <math.h>")
ROOT.gSystem.Load("libDelphes");

try:
  ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
  ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
  print "Delphes classes imported"
except:
  pass

TreeName="Delphes"
DataChain=ROOT.TChain(TreeName)
DataChain.Add(args.DIR+"*.root")

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(DataChain)
numberOfEntries = treeReader.GetEntries()

#branchJet = treeReader.UseBranch("Jet")
branchMuon = treeReader.UseBranch("Muon")
branchElectron = treeReader.UseBranch("Electron")
#branchMET = treeReader.UseBranch("MissingET")

#ListOfEventNumber=[]

NumberOfEventsToCheck=DataChain.GetEntries()

print "--------------------------------------> Going to analyze", NumberOfEventsToCheck, "events!!!!!!!!!!!"

CountingEvents=0

MinValBin1=600
MinValBin2=900
MaxValBin2=1300
Val1Counter=0
Val2Counter=0

if args.ANA=="Electron":

    #Electron selection: pT(e)> 35 GeV, |eta(e)|<1.44, 1.57<|eta(e)|<2.5

    Lumi=137000. #pb-1
    
    if args.OUT is not None:
        RootFile = ROOT.TFile(args.OUT,"recreate")
        BinArray=[1370,1430,1490,1550,1610,1680,1750,1820,1890,1970,2050,2130,2210,2290,2370,2530,2760,4000]
        MLLhisto = ROOT.TH1F("Mll","Mll",len(BinArray)-1,array('d',BinArray))
        #MLLhisto = ROOT.TH1F("Mll","Mll",688,60,3500)

    for entry in xrange(NumberOfEventsToCheck):
        treeReader.ReadEntry(entry)
        #Selection from PAS-EXO-09-019
        #print "Processing event:", entry, "----------------------------------------------------------------------------------------------"
        if branchElectron.GetEntries() < 2: continue
        NE=0
        for i in xrange(branchElectron.GetEntries()):
            if branchElectron.At(i).PT>35 and abs(branchElectron.At(i).Eta)<2.5:
                if abs(branchElectron.At(i).Eta)<1.44 or abs(branchElectron.At(i).Eta)>1.57:
                    NE+=1
        #print NE
        if NE != 2: continue
        Electron1 = ROOT.TLorentzVector(0,0,0,0)
        Electron1.SetPtEtaPhiM(branchElectron.At(0).PT,branchElectron.At(0).Eta,branchElectron.At(0).Phi,0.000511)
        Electron2 = ROOT.TLorentzVector(0,0,0,0)
        Electron2.SetPtEtaPhiM(branchElectron.At(1).PT,branchElectron.At(1).Eta,branchElectron.At(1).Phi,0.000511)
        DiElectron = Electron1 + Electron2
        EventMLL=DiElectron.M()
        #print EventMLL
        if EventMLL>=MinValBin1 and EventMLL<MinValBin2:
            Val1Counter+=1
        if EventMLL>=MinValBin2 and EventMLL<MaxValBin2:
            Val2Counter+=1
        if EventMLL>float(args.MLL):
            CountingEvents=CountingEvents+1
            if args.OUT is not None: MLLhisto.Fill(EventMLL)

if args.ANA=="Muon":

    #Muon selesction: pt(mu)>53, |eta(mu)|<2.4, Mmumu>150

    Lumi=140000. #pb-1
    
    if args.OUT is not None:
        RootFile = ROOT.TFile(args.OUT,"recreate")
        BinArray=[1410,1530,1660,1790,1940,2100,2280,2480,2680,2900,3150,4000]
        MLLhisto = ROOT.TH1F("Mll","Mll",len(BinArray)-1,array('d',BinArray))
        #MLLhisto = ROOT.TH1F("Mll","Mll",670,150,3500)
        
    for entry in xrange(NumberOfEventsToCheck):
        treeReader.ReadEntry(entry)
        #Selection from PAS-EXO-09-019
        #print "Processing event:", entry, "----------------------------------------------------------------------------------------------"
        if branchMuon.GetEntries() < 2: continue
        NE=0
        for i in xrange(branchMuon.GetEntries()):
            if branchMuon.At(i).PT>53 and abs(branchMuon.At(i).Eta)<2.4:
                NE+=1
        if NE != 2: continue
        Muon1 = ROOT.TLorentzVector(0,0,0,0)
        Muon1.SetPtEtaPhiM(branchMuon.At(0).PT,branchMuon.At(0).Eta,branchMuon.At(0).Phi,0.10566)
        Muon2 = ROOT.TLorentzVector(0,0,0,0)
        Muon2.SetPtEtaPhiM(branchMuon.At(1).PT,branchMuon.At(1).Eta,branchMuon.At(1).Phi,0.10566)
        DiMuon = Muon1 + Muon2
        EventMLL=DiMuon.M()
        if EventMLL<150: continue
        if EventMLL>=MinValBin1 and EventMLL<MinValBin2:
            Val1Counter+=1
        if EventMLL>=MinValBin2 and EventMLL<MaxValBin2:
            Val2Counter+=1
        if EventMLL>float(args.MLL):
            CountingEvents=CountingEvents+1
            if args.OUT is not None: MLLhisto.Fill(EventMLL)
            
print "Number of events before full selection:", NumberOfEventsToCheck
print "Number of events that passed full selection:", CountingEvents, "+-", sqrt(CountingEvents)
#TotalXS=10150 #pb
TotalXS=float(args.XS)
#Lumi=35900. #pb-1
Weight=(Lumi*TotalXS/NumberOfEventsToCheck)
print "Used cross section", TotalXS, "pb"
print "Weighted number of events that passed full selection:", CountingEvents*Weight, "+-", sqrt(CountingEvents)*Weight
print "Bin content belong 600-900:", Val1Counter*Weight
print "Bin content belong 900-1300:", Val2Counter*Weight
if args.OUT is not None:
    MLLhisto.Sumw2()
    MLLhisto.Scale(Weight)
    print "Bin range,", "Value"
    for i in xrange(MLLhisto.GetNbinsX()):
        print str(MLLhisto.GetBinLowEdge(i+1))+"-"+str(MLLhisto.GetBinLowEdge(i+1)+MLLhisto.GetBinWidth(i+1))+",", MLLhisto.GetBinContent(i+1)
    MLLhisto.Write()
    RootFile.Close()

