from GlobalFunctions import *

if args.ANA=="CMS Dijet":

    #Based on 1806.00843: 
    CountingEvents=0
    
    InitialEvts, TreeReader, Branches = DelphesInit(JetBranch="Jet", MetBranch="", MuonBranch="", ElectronBranch="")

    RootFile = BasketFile()
    DijetHisto = ROOT.TH1F("mjj","mjj",20,200,2000)

    for entry in xrange(InitialEvts):
        TreeReader.ReadEntry(entry)

        NJets=Np(Branches["Jet"])
        
        if NJets < 2: continue
        if PT(Branches["Jet"],0)<500 or abs(ETA(Branches["Jet"],0))>2.5: continue
        if abs(ETA(Branches["Jet"],1))>2.5: continue

        #Fat jets
        FJet1=ROOT.TLorentzVector()
        FJet2=ROOT.TLorentzVector()
        FJet1.SetPtEtaPhiM(PT(Branches["Jet"],0),ETA(Branches["Jet"],0),PHI(Branches["Jet"],0),M(Branches["Jet"],0))
        FJet2.SetPtEtaPhiM(PT(Branches["Jet"],1),ETA(Branches["Jet"],1),PHI(Branches["Jet"],1),M(Branches["Jet"],1))

        Jeti=ROOT.TLorentzVector()
        
        HT=0
        for i in xrange(NJets):
            if PT(Branches["Jet"],i)>30 and abs(ETA(Branches["Jet"],i))<2.5:
                HT+=PT(Branches["Jet"],i)
                if i!=0 and i!=1:
                    Jeti.SetPtEtaPhiM(PT(Branches["Jet"],i),ETA(Branches["Jet"],i),PHI(Branches["Jet"],i),M(Branches["Jet"],i))
                    dj1=DeltaR(Branches["Jet"],0,Branches["Jet"],i)
                    dj2=DeltaR(Branches["Jet"],1,Branches["Jet"],i)
                    if dj1<1.1 and dj2>1.1:
                        FJet1=FJet1+Jeti
                    if dj2<1.1 and dj1>1.1:
                        FJet2=FJet2+Jeti
                    if dj1<1.1 and dj2<1.1 and dj1<dj2:
                        FJet1=FJet1+Jeti
                    if dj1<1.1 and dj2<1.1 and dj2<dj1:
                        FJet2=FJet2+Jeti

        if HT<900: continue
        if abs(FJet1.Eta()-FJet2.Eta())>1.3: continue

        dijet=FJet1+FJet2
        
        CountingEvents=CountingEvents+1
        DijetHisto.Fill(diijet.M())

    DisplayFinalInfo(InitialEvts, CountingEvents, Histo=DijetHisto, MyFile=RootFile)
