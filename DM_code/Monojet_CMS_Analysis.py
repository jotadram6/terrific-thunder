import commands as cmd
import argparse
import ROOT
from array import array
from math import *

parser = argparse.ArgumentParser()
parser.add_argument('--DIR', help='Directory to read input')
parser.add_argument('--XS', help='Sample cross section')
parser.add_argument('--MET', help='Minimum met')
parser.add_argument('--OUT', help='ROOT output file')
#parser.add_argument('--ANA', help='Analysis selector')
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

branchJet = treeReader.UseBranch("Jet")
branchMuon = treeReader.UseBranch("Muon")
branchElectron = treeReader.UseBranch("Electron")
branchMET = treeReader.UseBranch("MissingET")

#ListOfEventNumber=[]

NumberOfEventsToCheck=DataChain.GetEntries()

print "--------------------------------------> Going to analyze", NumberOfEventsToCheck, "events!!!!!!!!!!!"

CountingEvents=0

#################################
#Monojet CMS analysis 1703.01651#
#################################

MinValBin1=200
MaxValBin1=350
Val1Counter=0

#Selection: MET>200 GeV, pT(j1)> 100 GeV, |eta(j1)|<2.5, HT>110 GeV, No lep with pT>10 or tau with pT>18, |eta(e)|<2.5, |eta(mu)|<2.4, No b-tagged jets with pT > 15 with Medium WP, |DeltaPhi(j1-4,MET)|>0.5
if args.OUT is not None:
    RootFile = ROOT.TFile(args.OUT,"recreate")
    #BinArray=[130,138.93,148.473,158.672,169.572,181.22,193.668,206.972,221.189,236.383,252.621,269.974,288.519,308.338,329.518,352.154,376.344,402.195,429.823,459.349,490.902,524.623,560.661,599.174,640.332,684.318,731.325,781.561,835.248,892.623,953.939,1019.47,1089.5,1164.34,1244.32,1329.79,1421.14,1518.76,1623.09,1734.58,1853.73,1981.07,2117.15,2262.58,2418,2584.1,2761.61,2951.31,3154.04,3370.7,3602.24,3849.68,4114.12,4396.73,4698.75,5021.52,5366.46,5735.09,6129.05,6550.06,7000.]
    #METhisto = ROOT.TH1F("MET","MET",len(BinArray)-1,array('d',BinArray))
    METhisto = ROOT.TH1F("MET","MET",20,200,2000)

for entry in xrange(NumberOfEventsToCheck):
    treeReader.ReadEntry(entry)
    if branchMET.At(0).MET<=200: continue
    if branchJet.GetEntries() < 1: continue
    if branchJet.At(0).PT<100 or abs(branchJet.At(0).Eta)>2.5: continue
    NTau=0
    HT=0

    Dp1=DeltaPhi(branchJet.At(0).Phi,branchMET.At(0).Phi)
    if branchJet.GetEntries()>=2: Dp2=DeltaPhi(branchJet.At(1).Phi,branchMET.At(0).Phi)
    if branchJet.GetEntries()>=3: Dp3=DeltaPhi(branchJet.At(2).Phi,branchMET.At(0).Phi)
    if branchJet.GetEntries()>=4: Dp4=DeltaPhi(branchJet.At(3).Phi,branchMET.At(0).Phi)        

    for i in xrange(branchJet.GetEntries()):
        if branchJet.At(i).PT>20:
            HT+=branchJet.At(i).PT
        #print "Jet:", i, "Tautagged:", branchTau.At(i).TauTag
        if branchJet.At(i).TauTag and branchJet.At(i).PT>15:
            NTau+=1
    if NTau != 0: continue
    if branchElectron.GetEntries() > 0: continue
    NE=0
    for i in xrange(branchElectron.GetEntries()):
        if branchElectron.At(i).PT>10 and abs(branchElectron.At(i).Eta)<2.5:
            NE+=1
    if NE != 0: continue
    if branchMuon.GetEntries() > 0: continue
    NMu=0
    for i in xrange(branchMuon.GetEntries()):
        if branchMuon.At(i).PT>10 and abs(branchMuon.At(i).Eta)<2.4:
            NMu+=1
    if NMu != 0: continue

    if branchMET.At(0).MET>=MinValBin1 and branchMET.At(0).MET<MaxValBin1:
        Val1Counter+=1
    if branchMET.At(0).MET>float(args.MET):
        CountingEvents=CountingEvents+1
        if args.OUT is not None: METhisto.Fill(branchMET.At(0).MET)


print "Number of events before full selection:", NumberOfEventsToCheck
print "Number of events that passed full selection:", CountingEvents, "+-", sqrt(CountingEvents)
#TotalXS=10150 #pb
TotalXS=float(args.XS)
Lumi=12900. #pb-1
Weight=(Lumi*TotalXS/NumberOfEventsToCheck)
print "Used cross section", TotalXS, "pb"
print "Weighted number of events that passed full selection:", CountingEvents*Weight, "+-", sqrt(CountingEvents)*Weight
print "Bin content belong 200-350:", Val1Counter*Weight
if args.OUT is not None:
    METhisto.Sumw2()
    METhisto.Scale(Weight)
    print "Bin range,", "Value"
    for i in xrange(METhisto.GetNbinsX()):
        print str(METhisto.GetBinLowEdge(i+1))+"-"+str(METhisto.GetBinLowEdge(i+1)+METhisto.GetBinWidth(i+1))+",", METhisto.GetBinContent(i+1)
    METhisto.Write()
    RootFile.Close()

