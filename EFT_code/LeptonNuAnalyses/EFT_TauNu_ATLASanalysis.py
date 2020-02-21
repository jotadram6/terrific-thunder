import commands as cmd
import argparse
import ROOT
from array import array
from math import *

parser = argparse.ArgumentParser()
parser.add_argument('--DIR', help='Directory to read input')
parser.add_argument('--XS', help='Sample cross section')
parser.add_argument('--MT', help='Minimum MT')
parser.add_argument('--OUT', help='ROOT output file')
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

#Delphes_Path="/home/joser/Pheno_Studies/MG5_aMC_v2_6_1/Delphes/"
#Delphes_Path="/home/joser/Dropbox/Vandy/EFT/ElectronNeutrino/Delphes-3.4.1/"
Delphes_Path="/cms/Jose/EFT_MG_PRL_version/MG5_aMC_v2_6_1/Delphes/"
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
branchTau = treeReader.UseBranch("Jet")
branchMET = treeReader.UseBranch("MissingET")

#ListOfEventNumber=[]

NumberOfEventsToCheck=DataChain.GetEntries()

print "--------------------------------------> Going to analyze", NumberOfEventsToCheck, "events!!!!!!!!!!!"

CountingEvents=0

if args.OUT is not None:
    RootFile = ROOT.TFile(args.OUT,"recreate")
    BinArray=[501,563,632,709,797,894,1004,1128,1266,1422,1597,1793,2013,2260,2538,2850,3200] #Double_t xAxis2[13] = {1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3200,3800,5000};
    MThisto = ROOT.TH1F("MT","MT",16,array('d',BinArray))
    #MThisto = ROOT.TH1F("MT","MT",30,float(args.MT),2000.)

for entry in xrange(NumberOfEventsToCheck):
    treeReader.ReadEntry(entry)
    #Selection from 1801.06992
    #print "Processing event:", entry, "----------------------------------------------------------------------------------------------"
    if branchMET.At(0).MET<=150: continue
    if branchTau.GetEntries() < 1: continue
    TauJet=-1
    NTau=0
    for i in xrange(branchTau.GetEntries()):
        #print "Jet:", i, "Tautagged:", branchTau.At(i).TauTag
        if branchTau.At(i).TauTag:
            NTau+=1
            TauJet=i
    #print NTau, TauJet
    #if branchTau.GetEntries() != 1: continue
    if NTau != 1: continue
    if abs(branchTau.At(TauJet).Eta)>=2.4: continue
    if abs(branchTau.At(TauJet).Eta)>1.37 and abs(branchTau.At(TauJet).Eta)<1.52: continue
    if branchTau.At(TauJet).PT<=50: continue
    #if branchTau.At(0).SumPt>=5: continue
    if abs(DeltaPhi(branchTau.At(TauJet).Phi,branchMET.At(0).Phi))<=2.4: continue
    if (branchTau.At(TauJet).PT/branchMET.At(0).MET)>0.7 and (branchTau.At(TauJet).PT/branchMET.At(0).MET)<1.3:
        EventMT=MT(branchTau.At(TauJet).PT,branchMET.At(0).MET,branchTau.At(TauJet).Phi,branchMET.At(0).Phi)
        if EventMT>float(args.MT):
            CountingEvents=CountingEvents+1
            MThisto.Fill(EventMT)



print "Number of events before full selection:", NumberOfEventsToCheck
print "Number of events that passed full selection:", CountingEvents, "+-", sqrt(CountingEvents)
#TotalXS=10150 #pb
TotalXS=float(args.XS)
Lumi=36100. #pb-1
Weight=(Lumi*TotalXS/NumberOfEventsToCheck)
print "Weighted number of events that passed full selection:", CountingEvents*Weight, "+-", sqrt(CountingEvents)*Weight
if args.OUT is not None:
    MThisto.Sumw2()
    MThisto.Scale(Weight)
    print "Bin range,", "Value"
    for i in xrange(MThisto.GetNbinsX()):
        print str(MThisto.GetBinLowEdge(i+1))+"-"+str(MThisto.GetBinLowEdge(i+1)+MThisto.GetBinWidth(i+1))+",", MThisto.GetBinContent(i+1)
    MThisto.Write()
    RootFile.Close()
