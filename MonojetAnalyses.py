from GlobalFunctions import *

if args.ANA=="CMS Monojet":

    CountingEvents=0
    
    InitialEvts, TreeReader, Branches = DelphesInit(JetBranch="Jet", MetBranch="MissingET", MuonBranch="Muon", ElectronBranch="Electron")

    RootFile = BasketFile()
    METhisto = ROOT.TH1F("MET","MET",20,200,2000)

    for entry in xrange(InitialEvts):
        TreeReader.ReadEntry(entry)

        NJets=Np(Branches["Jet"])
        
        if Branches["Met"].At(0).MET<=200: continue
        if NJets < 1: continue
        if PT(Branches["Jet"],0)<100 or abs(ETA(Branches["Jet"],0))>2.5: continue

        Dp1=DeltaPhi(PHI(Branches["Jet"],0),PHI(Branches["Met"],0))
        if NJets>=2: Dp2=DeltaPhi(PHI(Branches["Jet"],1),PHI(Branches["Met"],0))
        if NJets>=3: Dp3=DeltaPhi(PHI(Branches["Jet"],2),PHI(Branches["Met"],0))
        if NJets>=4: Dp4=DeltaPhi(PHI(Branches["Jet"],3),PHI(Branches["Met"],0))

        #Tau veto and recalculation of HT
        NTau=0
        HT=0
        for i in xrange(NJets):
            if PT(Branches["Jet"],i)>20:
                HT+=PT(Branches["Jet"],i)
                #print "Jet:", i, "Tautagged:", branchTau.At(i).TauTag
            if Branches["Jet"].At(i).TauTag and PT(Branches["Jet"],i)>15:
                NTau+=1
        if NTau != 0: continue

        #Electron veto
        if not Veto(Branches["Electron"],PTCut=10,ETACut=2.5,ETAgap=[1.44,1.57]): continue

        #Muon veto
        if not Veto(Branches["Muon"],PTCut=10,ETACut=2.4): continue

        if Branches["Met"].At(0).MET>float(args.MET):
            CountingEvents=CountingEvents+1
            METhisto.Fill(Branches["Met"].At(0).MET)

    DisplayFinalInfo(InitialEvts, CountingEvents, Histo=METhisto, MyFile=RootFile)
