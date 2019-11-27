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
Delphes_Path="/home/joser/Dropbox/Vandy/EFT/ElectronNeutrino/Delphes-3.4.1/"
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
    #BinArray=[480.,560.,640.,720.,800.,880.,960.,1040.,1120.,1200.,1280.,1360.,1440.,1520.,1600.,1680.,1760.,1840.,1920.,2000.,2080.,2160.,2240.,2320.,2400.,2480.,2560.,2640.,2720.,2800.,2880.,2960.,3040.,3120.,3200.,3280.,3360.,3440.,3520.,3600.,3680.,3760.,3840.,3920.,4000.,4080.] #Double_t xAxis2[13] = {1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3200,3800,5000};
    #print "Numbers of bins:", len(BinArray)
    #MThisto = ROOT.TH1F("MT","MT",len(BinArray)+1,array('d',BinArray))
    #MThisto = ROOT.TH1F("MT","MT",len(BinArray)+1,array('d',BinArray))
    MThisto = ROOT.TH1F("MT","MT",45,480.,4080.)
    #MThisto = ROOT.TH1F("MT","MT",30,float(args.MT),2000.)

for entry in xrange(NumberOfEventsToCheck):
    treeReader.ReadEntry(entry)
    #Selection from 1801.06992
    #print "Processing event:", entry, "----------------------------------------------------------------------------------------------"
    if branchMET.At(0).MET<=200: continue
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
    #if abs(branchTau.At(TauJet).Eta)>1.37 and abs(branchTau.At(TauJet).Eta)<1.52: continue
    if branchTau.At(TauJet).PT<=80: continue
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
