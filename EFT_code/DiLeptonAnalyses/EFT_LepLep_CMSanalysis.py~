import commands as cmd
import argparse
import ROOT
from array import array
from math import *

parser = argparse.ArgumentParser()
parser.add_argument('--DIR', help='Directory to read input')
parser.add_argument('--XS', help='Sample cross section')
#parser.add_argument('--MT', help='Minimum MT')
parser.add_argument('--OUT', help='ROOT output file')
parser.add_argument('--ANA', help='Analysis selector')
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

Delphes_Path="/cms/Jose/EFT_MG_PRL_version/MG5_aMC_v2_6_1/Delphes/"
#Delphes_Path="/home/joser/Dropbox/Vandy/EFT/ElectronNeutrino/Delphes-3.4.1/"
ROOT.gSystem.AddDynamicPath(Delphes_Path)
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

#branchElectron = treeReader.UseBranch("Electron")
#branchJet = treeReader.UseBranch("Jet")
branchMuon = treeReader.UseBranch("Muon")
branchElectron = treeReader.UseBranch("Electron")
branchMET = treeReader.UseBranch("MissingET")

#ListOfEventNumber=[]

NumberOfEventsToCheck=DataChain.GetEntries()

print "--------------------------------------> Going to analyze", NumberOfEventsToCheck, "events!!!!!!!!!!!"

CountingEvents=0

if args.ANA=="Electron":

    #Electron selection: MET>65, exactly one tight muon, pt electron > 65, evt veto if more lep with medium and 20 GeV pt, MT>130
    
    if args.OUT is not None:
        RootFile = ROOT.TFile(args.OUT,"recreate")
        #BinArray=[130,138.93,148.473,158.672,169.572,181.22,193.668,206.972,221.189,236.383,252.621,269.974,288.519,308.338,329.518,352.154,376.344,402.195,429.823,459.349,490.902,524.623,560.661,599.174,640.332,684.318,731.325,781.561,835.248,892.623,953.939,1019.47,1089.5,1164.34,1244.32,1329.79,1421.14,1518.76,1623.09,1734.58,1853.73,1981.07,2117.15,2262.58,2418,2584.1,2761.61,2951.31,3154.04,3370.7,3602.24,3849.68,4114.12,4396.73,4698.75,5021.52,5366.46,5735.09,6129.05,6550.06,7000.]
        #MThisto = ROOT.TH1F("MT","MT",len(BinArray)-1,array('d',BinArray))
        MThisto = ROOT.TH1F("MT","MT",94,250,4000)

    for entry in xrange(NumberOfEventsToCheck):
        treeReader.ReadEntry(entry)
        #Selection from 1706.04786
        #print "Processing event:", entry, "----------------------------------------------------------------------------------------------"
        if branchMET.At(0).MET<=150: continue
        if branchElectron.GetEntries() < 1: continue
        NE=0
        for i in xrange(branchElectron.GetEntries()):
            if branchElectron.At(i).PT>25:
                NE+=1
        if NE != 1: continue
        if abs(branchElectron.At(0).Eta)>=2.5: continue
        if abs(branchElectron.At(0).Eta)>1.44 and abs(branchElectron.At(0).Eta)<1.56: continue
        if branchElectron.At(0).PT<=135: continue
        if (branchElectron.At(0).PT/branchMET.At(0).MET)>=1.5 or (branchElectron.At(0).PT/branchMET.At(0).MET)<=0.4: continue
        if abs(DeltaPhi(branchElectron.At(0).Phi,branchMET.At(0).Phi))<=2.5: continue
        EventMT=MT(branchElectron.At(0).PT,branchMET.At(0).MET,branchElectron.At(0).Phi,branchMET.At(0).Phi)
        if EventMT>250.:
            CountingEvents=CountingEvents+1
            MThisto.Fill(EventMT)

if args.ANA=="Muon":

    #Muon selesction: MET>55, exactly one tight muon, pt muon > 55, evt veto if more lep with medium and 20 GeV pt, MT>110
    
    if args.OUT is not None:
        RootFile = ROOT.TFile(args.OUT,"recreate")
        BinArray=[1040.,1100.,1200.,1300.,1400.,1500.,1650.,1800.,1950.,2150.,2350.,2550.,4000.]
        #MThisto = ROOT.TH1F("MT","MT",len(BinArray)-1,array('d',BinArray))
        #MThisto = ROOT.TH1F("MT","MT",95,250,4000)
        
    for entry in xrange(NumberOfEventsToCheck):
        treeReader.ReadEntry(entry)
        #Selection from 1706.04786
        #print "Processing event:", entry, "----------------------------------------------------------------------------------------------"
        if branchMET.At(0).MET<=50: continue
        if branchMuon.GetEntries() < 1: continue
        NE=0
        for i in xrange(branchMuon.GetEntries()):
            if branchMuon.At(i).PT>25:
                NE+=1
        if NE != 1: continue
        if abs(branchMuon.At(0).Eta)>=2.4: continue
        #if abs(branchMuon.At(0).Eta)>1.37 and abs(branchMuon.At(0).Eta)<1.52: continue
        if branchMuon.At(0).PT<=53: continue
        if (branchMuon.At(0).PT/branchMET.At(0).MET)>=1.5 or (branchMuon.At(0).PT/branchMET.At(0).MET)<=0.4: continue
        if abs(DeltaPhi(branchMuon.At(0).Phi,branchMET.At(0).Phi))<=2.5: continue
        EventMT=MT(branchMuon.At(0).PT,branchMET.At(0).MET,branchMuon.At(0).Phi,branchMET.At(0).Phi)
        if EventMT>100.:
            CountingEvents=CountingEvents+1
            MThisto.Fill(EventMT)
            
print "Number of events before full selection:", NumberOfEventsToCheck
print "Number of events that passed full selection:", CountingEvents, "+-", sqrt(CountingEvents)
#TotalXS=10150 #pb
TotalXS=float(args.XS)
Lumi=35900. #pb-1
Weight=(Lumi*TotalXS/NumberOfEventsToCheck)
print "Used cross section", TotalXS, "pb"
print "Weighted number of events that passed full selection:", CountingEvents*Weight, "+-", sqrt(CountingEvents)*Weight
if args.OUT is not None:
    MThisto.Sumw2()
    MThisto.Scale(Weight)
    print "Bin range,", "Value"
    for i in xrange(MThisto.GetNbinsX()):
        print str(MThisto.GetBinLowEdge(i+1))+"-"+str(MThisto.GetBinLowEdge(i+1)+MThisto.GetBinWidth(i+1))+",", MThisto.GetBinContent(i+1)
    MThisto.Write()
    RootFile.Close()

