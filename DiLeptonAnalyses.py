from GlobalFunctions import *

FlagAnalysisDone=False

if args.ANA=="ATLAS Electron":

    #Selection from 1903.06248
    #Electron selection: pT(e)> 30 GeV, |eta(e)|<1.37, 1.52<|eta(e)|<2.4, mee>225 GeV

    CountingEvents=0
    
    InitialEvts, TreeReader, Branches = DelphesInit(JetBranch="", MetBranch="MissingET", MuonBranch="", ElectronBranch="Electron")

    RootFile = BasketFile()

    BinArray=[602.522,622.634,643.417,664.894,687.087,710.022,733.722,758.213,783.521,809.674,836.7,864.629,893.489,923.313,954.133,985.981,1018.89,1052.9,1088.05,1124.36,1161.9,1200.68,1240.76,1282.17,1324.97,1369.2,1414.9,1462.13,1510.93,1561.36,1613.48,1667.34,1722.99,1780.5,1839.94,1901.35,1964.82,2030.4,2098.17,2168.21,2240.58,2315.37,2392.65,2472.52,2555.05,2640.34,2728.47,2819.54,2913.66,3010.91,3111.41,3215.27,3322.59,3433.5,3548.1,3666.54,3788.92,3915.39,4046.09,4181.14,4320.7,4464.92,4613.96,4767.97,4927.12,5091.58,5261.54,5437.16,5618.65]
    MLLhisto = ROOT.TH1F("MLL","MLL",len(BinArray)-1,array('d',BinArray))

    for entry in xrange(InitialEvts):
        TreeReader.ReadEntry(entry)
        INE=Np(Branches["Electron"])
        if INE < 2: continue
        NE=0
        for i in xrange(INE):
            if PT(Branches["Electron"],i)>60 and abs(ETA(Branches["Electron"],i))<2.4:
                if abs(ETA(Branches["Electron"],i))<1.37 or abs(ETA(Branches["Electron"],i))>1.52:
                    NE+=1

        if NE != 2: continue
        Electron1 = GetParticle(Branches["Electron"],0,ElectronMass)
        Electron2 = GetParticle(Branches["Electron"],1,ElectronMass)
        DiElectron = Electron1 + Electron2
        EventMLL=DiElectron.M()
        if EventMLL<225: continue
        elif EventMLL>float(args.MLL):
            CountingEvents=CountingEvents+1
            MLLhisto.Fill(EventMLL)


    DisplayFinalInfo(InitialEvts, CountingEvents, Histo=MLLhisto, MyFile=RootFile)
    FlagAnalysisDone=True
    
if args.ANA=="ATLAS Muon":

    #Selection from 1903.06248
    #Muon selesction: pt(mu)>30, |eta(mu)|<2.5, Mmumu>225

    CountingEvents=0
    
    InitialEvts, TreeReader, Branches = DelphesInit(JetBranch="", MetBranch="", MuonBranch="Muon", ElectronBranch="")

    RootFile = BasketFile()

    BinArray=[602.522,622.634,643.417,664.894,687.087,710.022,733.722,758.213,783.521,809.674,836.7,864.629,893.489,923.313,954.133,985.981,1018.89,1052.9,1088.05,1124.36,1161.9,1200.68,1240.76,1282.17,1324.97,1369.2,1414.9,1462.13,1510.93,1561.36,1613.48,1667.34,1722.99,1780.5,1839.94,1901.35,1964.82,2030.4,2098.17,2168.21,2240.58,2315.37,2392.65,2472.52,2555.05,2640.34,2728.47,2819.54,2913.66,3010.91,3111.41,3215.27,3322.59,3433.5,3548.1,3666.54,3788.92,3915.39,4046.09,4181.14,4320.7,4464.92,4613.96,4767.97,4927.12,5091.58,5261.54,5437.16,5618.65]
    MLLhisto = ROOT.TH1F("MT","MT",len(BinArray)-1,array('d',BinArray))

    for entry in xrange(InitialEvts):
        TreeReader.ReadEntry(entry)
        INM=Np(Branches["Muon"])
        if INM < 2: continue
        NE=0
        for i in xrange(INM):
            if PT(Branches["Muon"],i)>60 and abs(ETA(Branches["Muon"],i))<2.5:
                NE+=1
        if NE != 2: continue
        Muon1 = GetParticle(Branches["Muon"],0,MuonMass)
        Muon2 = GetParticle(Branches["Muon"],1,MuonMass)
        DiMuon = Muon1 + Muon2
        EventMLL=DiMuon.M()
        if EventMLL<225: continue
        elif EventMLL>float(args.MLL):
            CountingEvents=CountingEvents+1
            MLLhisto.Fill(EventMLL)


    DisplayFinalInfo(InitialEvts, CountingEvents, Histo=MLLhisto, MyFile=RootFile)
    FlagAnalysisDone=True

if args.ANA=="ATLAS Tau":

    #Selection from 1709.07242
    #Electron selection: pT(tau)> 65 GeV with Ntau>=2, opposite charges for the two leading taus, |DeltaPhi(tau1,tau2)|>2.7

    CountingEvents=0
    
    InitialEvts, TreeReader, Branches = DelphesInit(JetBranch="Jet", MetBranch="MissingET", MuonBranch="Muon", ElectronBranch="Electron")

    RootFile = BasketFile()

    BinArray=[150,170,180,190,200,210,220,240,260,300,350,400,450,500,600,700,800,1000]
    MLLhisto = ROOT.TH1F("MLL","MLL",len(BinArray)-1,array('d',BinArray))
    #MToThisto = ROOT.TH1F("MToT","MToT",len(BinArray)-1,array('d',BinArray))

    for entry in xrange(InitialEvts):
        TreeReader.ReadEntry(entry)
        TausList=GetTaus(Branches["Jet"],PTcut=65)
        #if not Veto(Branches["Electron"],PTCut=15,ETACut=2.47,ETAgap=[1.37,1.52]): continue
        #if not Veto(Branches["Muon"],PTCut=7,ETACut=2.5): continue
        if len(TausList)<2: continue
        if PT(Branches["Jet"],TausList[0])<130: continue
        if abs(ETA(Branches["Jet"],TausList[0]))>2.5: continue
        if abs(ETA(Branches["Jet"],TausList[1]))>2.5: continue
        if abs(ETA(Branches["Jet"],TausList[0]))>1.37 and abs(ETA(Branches["Jet"],TausList[0]))<1.52: continue
        if abs(ETA(Branches["Jet"],TausList[1]))>1.37 and abs(ETA(Branches["Jet"],TausList[1]))<1.52: continue
        if CH(Branches["Jet"],TausList[0])*CH(Branches["Jet"],TausList[1])>=0: continue
        if abs(DeltaPhi(PHI(Branches["Jet"],TausList[0]),PHI(Branches["Jet"],TausList[1])))<2.7: continue
        PT1M=PT(Branches["Jet"],TausList[0])
        PT1phi=PHI(Branches["Jet"],TausList[0])
        PT2M=PT(Branches["Jet"],TausList[1])
        PT2phi=PHI(Branches["Jet"],TausList[1])
        METM=Branches["Met"].At(0).MET
        METphi=PHI(Branches["Met"],0)
        EventMLL=TotalTransMass(PT1M,PT2M,METM,
                                PT1M*cos(PT1phi),PT1M*sin(PT1phi),
                                PT2M*cos(PT2phi),PT2M*sin(PT2phi),
                                METM*cos(METphi),METM*sin(METphi))
        #Mu1LV=ROOT.TLorentzVector()
        #Mu2LV=ROOT.TLorentzVector()
        #METLV=ROOT.TLorentzVector()
        #Mu1LV.SetPtEtaPhiM(PT(Branches["Jet"],TausList[0]),ETA(Branches["Jet"],TausList[0]),PHI(Branches["Jet"],TausList[0]),MuonMass)
        #Mu2LV.SetPtEtaPhiM(PT(Branches["Jet"],TausList[1]),ETA(Branches["Jet"],TausList[1]),PHI(Branches["Jet"],TausList[1]),MuonMass)
        #METLV.SetPxPyPzE(METM*cos(METphi),METM*sin(METphi),0.0,METM)
        #FullMTVector=Mu1LV+Mu2LV+METLV
        if EventMLL<150: continue
        elif EventMLL>float(args.MLL):
            CountingEvents=CountingEvents+1
            MLLhisto.Fill(EventMLL)
            #MToThisto.Fill(FullMTVector.M())

    DisplayFinalInfo(InitialEvts, CountingEvents, Histo=MLLhisto, MyFile=RootFile)        
    #DisplayFinalInfo(InitialEvts, CountingEvents, Histo=MToThisto, MyFile=RootFile)
    FlagAnalysisDone=True
    
if args.ANA=="CMS Electron":

    #Selection from PAS-EXO-09-019
    #Electron selection: pT(e)> 35 GeV, |eta(e)|<1.44, 1.57<|eta(e)|<2.5

    CountingEvents=0
    
    InitialEvts, TreeReader, Branches = DelphesInit(JetBranch="", MetBranch="", MuonBranch="", ElectronBranch="Electron")

    RootFile = BasketFile()

    BinArray=[600,900,1370,1430,1490,1550,1610,1680,1750,1820,1890,1970,2050,2130,2210,2290,2370,2530,2760,4000]
    MLLhisto = ROOT.TH1F("Mll","Mll",len(BinArray)-1,array('d',BinArray))

    for entry in xrange(InitialEvts):
        TreeReader.ReadEntry(entry)
        INE=Np(Branches["Electron"])
        if INE < 2: continue
        NE=0
        for i in xrange(INE):
            if PT(Branches["Electron"],i)>35 and abs(ETA(Branches["Electron"],i))<2.5:
                if abs(ETA(Branches["Electron"],i))<1.44 or abs(ETA(Branches["Electron"],i))>1.57:
                    NE+=1

        if NE != 2: continue
        Electron1 = GetParticle(Branches["Electron"],0,ElectronMass)
        Electron2 = GetParticle(Branches["Electron"],1,ElectronMass)
        DiElectron = Electron1 + Electron2
        EventMLL=DiElectron.M()
        if EventMLL>float(args.MLL):
            CountingEvents=CountingEvents+1
            MLLhisto.Fill(EventMLL)


    DisplayFinalInfo(InitialEvts, CountingEvents, Histo=MLLhisto, MyFile=RootFile)
    FlagAnalysisDone=True

if args.ANA=="CMS Muon":

    #Selection from PAS-EXO-09-019
    #Muon selesction: pt(mu)>53, |eta(mu)|<2.4, Mmumu>150

    CountingEvents=0
    
    InitialEvts, TreeReader, Branches = DelphesInit(JetBranch="", MetBranch="", MuonBranch="Muon", ElectronBranch="")

    RootFile = BasketFile()

    BinArray=[600,900,1410,1530,1660,1790,1940,2100,2280,2480,2680,2900,3150,4000]
    MLLhisto = ROOT.TH1F("Mll","Mll",len(BinArray)-1,array('d',BinArray))

    for entry in xrange(InitialEvts):
        TreeReader.ReadEntry(entry)
        INM=Np(Branches["Muon"])
        if INM < 2: continue
        NE=0
        for i in xrange(INM):
            if PT(Branches["Muon"],i)>53 and abs(ETA(Branches["Muon"],i))<2.4:
                NE+=1
        if NE != 2: continue
        Muon1 = GetParticle(Branches["Muon"],0,MuonMass)
        Muon2 = GetParticle(Branches["Muon"],1,MuonMass)
        DiMuon = Muon1 + Muon2
        EventMLL=DiMuon.M()
        if EventMLL<150: continue
        elif EventMLL>float(args.MLL):
            CountingEvents=CountingEvents+1
            MLLhisto.Fill(EventMLL)


    DisplayFinalInfo(InitialEvts, CountingEvents, Histo=MLLhisto, MyFile=RootFile)
    FlagAnalysisDone=True

if args.ANA=="CMS Tau":

    print "Not implemented yet"
    FlagAnalysisDone=True
    
if args.ANA=="CMS Muon AFB":

    #Selection from 
    #Muon selesction: pt(mu)>53, |eta(mu)|<2.4, Mmumu>150

    CountingMuEvents=0
    CountingElEvents=0
    
    InitialEvts, TreeReader, Branches = DelphesInit(JetBranch="", MetBranch="", MuonBranch="Muon", ElectronBranch="Electron")

    RootFile = BasketFile()

    YBinArray=[0.0,1.0,1.25,1.5,2.4,5.0]
    FBBinArray=[-1.0,0.0,1.0]
    MuAFBhisto = ROOT.TH2F("MuAFB","MuAFB",len(YBinArray)-1,array('d',YBinArray),len(FBBinArray)-1,array('d',FBBinArray))
    ElAFBhisto = ROOT.TH2F("ElAFB","ElAFB",len(YBinArray)-1,array('d',YBinArray),len(FBBinArray)-1,array('d',FBBinArray))
    Mllpeakmin=86
    Mllpeakmax=96

    for entry in xrange(InitialEvts):
        TreeReader.ReadEntry(entry)
        INM=Np(Branches["Muon"])
        if INM < 2: continue
        NE=0
        for i in xrange(INM):
            if PT(Branches["Muon"],i)>20 and abs(ETA(Branches["Muon"],i))<2.4:
                NE+=1
        if NE != 2: continue
        Muon1Ch=CH(Branches["Muon"],0)
        Muon2Ch=CH(Branches["Muon"],1)
        if Muon1Ch*Muon2Ch>=0: continue
        Muon1 = GetParticle(Branches["Muon"],0,MuonMass)
        Muon2 = GetParticle(Branches["Muon"],1,MuonMass)
        DiMuon = Muon1 + Muon2
        EventMLL=DiMuon.M()
        if EventMLL<86 or EventMLL>96: continue
        Q2=DiMuon.Mag2()
        QT2=DiMuon.Pt()**2
        if Muon1Ch<0:
            P1m=(Muon1.Energy()-Muon1.Pz())/sqrt(2)
            P1p=(Muon1.Energy()+Muon1.Pz())/sqrt(2)
            P2m=(Muon2.Energy()-Muon2.Pz())/sqrt(2)
            P2p=(Muon2.Energy()+Muon2.Pz())/sqrt(2)
        else:
            P1m=(Muon2.Energy()-Muon2.Pz())/sqrt(2)
            P1p=(Muon2.Energy()+Muon2.Pz())/sqrt(2)
            P2m=(Muon1.Energy()-Muon1.Pz())/sqrt(2)
            P2p=(Muon1.Energy()+Muon1.Pz())/sqrt(2)
        CosThetaNumerator=2*((P1p*P2m)-(P1m*P2p))
        CosThetaDenominator=sqrt(Q2*(Q2+QT2))
        CosTheta=cos(CosThetaNumerator/CosThetaDenominator)
        TransformedCosTheta=(abs(DiMuon.Pz())/DiMuon.Pz())*CosTheta
        DiMuonRapidity=0.5*log((DiMuon.Energy()+DiMuon.Pz())/(DiMuon.Energy()-DiMuon.Pz()))
        #print CosTheta, TransformedCosTheta, DiMuonRapidity
        #Forwadr=+1, Backward=-1
        if TransformedCosTheta>0:
            MuAFBhisto.Fill(abs(DiMuonRapidity),0.5)
        else:
            MuAFBhisto.Fill(abs(DiMuonRapidity),-0.5)
        CountingMuEvents=CountingMuEvents+1

    for entry in xrange(InitialEvts):
        TreeReader.ReadEntry(entry)
        INM=Np(Branches["Electron"])
        if INM < 2: continue
        NE=0
        for i in xrange(INM):
            if PT(Branches["Electron"],i)>20 and abs(ETA(Branches["Electron"],i))<2.4:
                NE+=1
        if NE != 2: continue
        Electron1Ch=CH(Branches["Electron"],0)
        Electron2Ch=CH(Branches["Electron"],1)
        if Electron1Ch*Electron2Ch>=0: continue
        Electron1 = GetParticle(Branches["Electron"],0,ElectronMass)
        Electron2 = GetParticle(Branches["Electron"],1,ElectronMass)
        DiElectron = Electron1 + Electron2
        EventMLL=DiElectron.M()
        if EventMLL<86 or EventMLL>96: continue
        Q2=DiElectron.Mag2()
        QT2=DiElectron.Pt()**2
        if Electron1Ch<0:
            P1m=(Electron1.Energy()-Electron1.Pz())/sqrt(2)
            P1p=(Electron1.Energy()+Electron1.Pz())/sqrt(2)
            P2m=(Electron2.Energy()-Electron2.Pz())/sqrt(2)
            P2p=(Electron2.Energy()+Electron2.Pz())/sqrt(2)
        else:
            P1m=(Electron2.Energy()-Electron2.Pz())/sqrt(2)
            P1p=(Electron2.Energy()+Electron2.Pz())/sqrt(2)
            P2m=(Electron1.Energy()-Electron1.Pz())/sqrt(2)
            P2p=(Electron1.Energy()+Electron1.Pz())/sqrt(2)
        CosThetaNumerator=2*((P1p*P2m)-(P1m*P2p))
        CosThetaDenominator=sqrt(Q2*(Q2+QT2))
        CosTheta=cos(CosThetaNumerator/CosThetaDenominator)
        TransformedCosTheta=(abs(DiElectron.Pz())/DiElectron.Pz())*CosTheta
        DiElectronRapidity=0.5*log((DiElectron.Energy()+DiElectron.Pz())/(DiElectron.Energy()-DiElectron.Pz()))
        #print CosTheta, TransformedCosTheta, DiElectronRapidity
        #Forwadr=+1, Backward=-1
        if TransformedCosTheta>0:
            ElAFBhisto.Fill(abs(DiElectronRapidity),0.5)
        else:
            ElAFBhisto.Fill(abs(DiElectronRapidity),-0.5)
        CountingElEvents=CountingElEvents+1

    TotalHisto=MuAFBhisto.Clone("TotalHisto")
    TotalHisto.Sumw2()
    TotalHisto.Add(ElAFBhisto)
        
    DisplayFinalInfo(InitialEvts, CountingMuEvents, Histo=[MuAFBhisto,ElAFBhisto,TotalHisto], MyFile=RootFile)
    FlagAnalysisDone=True

#NO ANALYSIS APPLIED
if not FlagAnalysisDone: print "No analysis applied!"
