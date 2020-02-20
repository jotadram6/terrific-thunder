from GlobalFunctions import *

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


if args.ANA=="ATLAS Tau":

    #Selection from 1709.07242
    #Electron selection: pT(tau)> 65 GeV with Ntau>=2, opposite charges for the two leading taus, |DeltaPhi(tau1,tau2)|>2.7

    CountingEvents=0
    
    InitialEvts, TreeReader, Branches = DelphesInit(JetBranch="Jet", MetBranch="", MuonBranch="", ElectronBranch="")

    RootFile = BasketFile()

    BinArray=[150,170,180,190,200,210,220,240,260,300,350,400,450,500,600,700,800,1000]
    MLLhisto = ROOT.TH1F("MLL","MLL",len(BinArray)-1,array('d',BinArray))

    for entry in xrange(InitialEvts):
        TreeReader.ReadEntry(entry)
        TausList=GetTaus(Branches["Jet"],PTcut=65)
        if len(TausList)<2: continue
        if CH(Branches["Jet"],TausList[0])*CH(Branches["Jet"],TausList[1])>=0: continue
        if abs(DeltaPhi(PHI(Branches["Jet"],TausList[0]),PHI(Branches["Jet"],TausList[1])))<2.7: continue
        PT1M=PT(Branches["Jet"],TausList[0])
        PT1phi=PHIBranches["Jet"],TausList[0])
        PT2M=PT(Branches["Jet"],TausList[1])
        PT2phi=PHIBranches["Jet"],TausList[1])
        METM=Branches["MissingET"].At(0).MET
        METphi=PHI(Branches["MissingET"],0)
        EventMLL=TotalTransMass(PT1M,PT2M,METM,
                                PT1M*cos(PT1phi),PT1M*sin(PT1phi),
                                PT2M*cos(PT2phi),PT2M*sin(PT2phi),
                                METM*cos(METphi),METM*sin(METphi))
        if EventMLL<150: continue
        elif EventMLL>float(args.MLL):
            CountingEvents=CountingEvents+1
            MLLhisto.Fill(EventMLL)


    DisplayFinalInfo(InitialEvts, CountingEvents, Histo=MLLhisto, MyFile=RootFile)        

    
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


if args.ANA=="CMS Tau":

    print "Not implemented yet"

