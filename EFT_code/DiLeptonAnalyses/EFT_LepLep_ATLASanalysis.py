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

MinValBin1=602.522
MinValBin2=893.489
MaxValBin2=1200.68
Val1Counter=0
Val2Counter=0

if args.ANA=="Electron":

    #Electron selection: pT(e)> 30 GeV, |eta(e)|<1.37, 1.52<|eta(e)|<2.4, mee>225 GeV
    if args.OUT is not None:
        RootFile = ROOT.TFile(args.OUT,"recreate")
        BinArray=[602.522,622.634,643.417,664.894,687.087,710.022,733.722,758.213,783.521,809.674,836.7,864.629,893.489,923.313,954.133,985.981,1018.89,1052.9,1088.05,1124.36,1161.9,1200.68,1240.76,1282.17,1324.97,1369.2,1414.9,1462.13,1510.93,1561.36,1613.48,1667.34,1722.99,1780.5,1839.94,1901.35,1964.82,2030.4,2098.17,2168.21,2240.58,2315.37,2392.65,2472.52,2555.05,2640.34,2728.47,2819.54,2913.66,3010.91,3111.41,3215.27,3322.59,3433.5,3548.1,3666.54,3788.92,3915.39,4046.09,4181.14,4320.7,4464.92,4613.96,4767.97,4927.12,5091.58,5261.54,5437.16,5618.65]
        MLLhisto = ROOT.TH1F("MT","MT",len(BinArray)-1,array('d',BinArray))
        #MLLhisto = ROOT.TH1F("Mll","Mll",688,60,3500)

    for entry in xrange(NumberOfEventsToCheck):
        treeReader.ReadEntry(entry)
        #Selection from 1903.06248
        #print "Processing event:", entry, "----------------------------------------------------------------------------------------------"
        if branchElectron.GetEntries() < 2: continue
        NE=0
        for i in xrange(branchElectron.GetEntries()):
            if branchElectron.At(i).PT>60 and abs(branchElectron.At(i).Eta)<2.4:
                if abs(branchElectron.At(i).Eta)<1.37 or abs(branchElectron.At(i).Eta)>1.52:
                    NE+=1
        #print NE
        if NE != 2: continue
        Electron1 = ROOT.TLorentzVector(0,0,0,0)
        Electron1.SetPtEtaPhiM(branchElectron.At(0).PT,branchElectron.At(0).Eta,branchElectron.At(0).Phi,0.000511)
        Electron2 = ROOT.TLorentzVector(0,0,0,0)
        Electron2.SetPtEtaPhiM(branchElectron.At(1).PT,branchElectron.At(1).Eta,branchElectron.At(1).Phi,0.000511)
        DiElectron = Electron1 + Electron2
        EventMLL=DiElectron.M()
        if EventMLL<225: continue
        #print EventMLL
        if EventMLL>=MinValBin1 and EventMLL<MinValBin2:
            Val1Counter+=1
        if EventMLL>=MinValBin2 and EventMLL<MaxValBin2:
            Val2Counter+=1
        if EventMLL>float(args.MLL):
            CountingEvents=CountingEvents+1
            if args.OUT is not None: MLLhisto.Fill(EventMLL)

if args.ANA=="Muon":

    #Muon selesction: pt(mu)>30, |eta(mu)|<2.5, Mmumu>225
    
    if args.OUT is not None:
        RootFile = ROOT.TFile(args.OUT,"recreate")
        BinArray=[602.522,622.634,643.417,664.894,687.087,710.022,733.722,758.213,783.521,809.674,836.7,864.629,893.489,923.313,954.133,985.981,1018.89,1052.9,1088.05,1124.36,1161.9,1200.68,1240.76,1282.17,1324.97,1369.2,1414.9,1462.13,1510.93,1561.36,1613.48,1667.34,1722.99,1780.5,1839.94,1901.35,1964.82,2030.4,2098.17,2168.21,2240.58,2315.37,2392.65,2472.52,2555.05,2640.34,2728.47,2819.54,2913.66,3010.91,3111.41,3215.27,3322.59,3433.5,3548.1,3666.54,3788.92,3915.39,4046.09,4181.14,4320.7,4464.92,4613.96,4767.97,4927.12,5091.58,5261.54,5437.16,5618.65]
        MLLhisto = ROOT.TH1F("MT","MT",len(BinArray)-1,array('d',BinArray))
        #MThisto = ROOT.TH1F("MT","MT",95,250,4000)
        #MLLhisto = ROOT.TH1F("Mll","Mll",670,150,3500)
        
    for entry in xrange(NumberOfEventsToCheck):
        treeReader.ReadEntry(entry)
        #Selection from 1903.06248
        #print "Processing event:", entry, "----------------------------------------------------------------------------------------------"
        if branchMuon.GetEntries() < 2: continue
        NE=0
        for i in xrange(branchMuon.GetEntries()):
            if branchMuon.At(i).PT>60 and abs(branchMuon.At(i).Eta)<2.5:
                NE+=1
        if NE != 2: continue
        Muon1 = ROOT.TLorentzVector(0,0,0,0)
        Muon1.SetPtEtaPhiM(branchMuon.At(0).PT,branchMuon.At(0).Eta,branchMuon.At(0).Phi,0.10566)
        Muon2 = ROOT.TLorentzVector(0,0,0,0)
        Muon2.SetPtEtaPhiM(branchMuon.At(1).PT,branchMuon.At(1).Eta,branchMuon.At(1).Phi,0.10566)
        DiMuon = Muon1 + Muon2
        EventMLL=DiMuon.M()
        if EventMLL<225: continue
        if EventMLL>=MinValBin1 and EventMLL<MinValBin2:
            Val1Counter+=1
        if EventMLL>=MinValBin2 and EventMLL<MaxValBin2:
            Val2Counter+=1
        if EventMLL>float(args.MLL):
            CountingEvents=CountingEvents+1
            if args.OUT is not None: MLLhisto.Fill(EventMLL)
            
print "Number of events before full selection:", NumberOfEventsToCheck
print "Number of events that passed full selection:", CountingEvents, "+-", sqrt(CountingEvents)
#TotalXS=139000 #pb
TotalXS=float(args.XS)
Lumi=139000. #pb-1
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
