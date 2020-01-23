from GlobalFunctions import *

if args.ANA=="ATLAS Electron":

    #Selection from 1903.06248
    #Electron selection: pT(e)> 30 GeV, |eta(e)|<1.37, 1.52<|eta(e)|<2.4, mee>225 GeV

    CountingEvents=0
    
    InitialEvts, TreeReader, Branches = DelphesInit(JetBranch="", MetBranch="", MuonBranch="", ElectronBranch="Electron")

    RootFile = BasketFile()

    BinArray=[602.522,622.634,643.417,664.894,687.087,710.022,733.722,758.213,783.521,809.674,836.7,864.629,893.489,923.313,954.133,985.981,1018.89,1052.9,1088.05,1124.36,1161.9,1200.68,1240.76,1282.17,1324.97,1369.2,1414.9,1462.13,1510.93,1561.36,1613.48,1667.34,1722.99,1780.5,1839.94,1901.35,1964.82,2030.4,2098.17,2168.21,2240.58,2315.37,2392.65,2472.52,2555.05,2640.34,2728.47,2819.54,2913.66,3010.91,3111.41,3215.27,3322.59,3433.5,3548.1,3666.54,3788.92,3915.39,4046.09,4181.14,4320.7,4464.92,4613.96,4767.97,4927.12,5091.58,5261.54,5437.16,5618.65]
    MLLhisto = ROOT.TH1F("MLL","MLL",len(BinArray)-1,array('d',BinArray))

    for entry in xrange(InitialEvts):
        TreeReader.ReadEntry(entry)
        INE=Nl(Branches["Electron"])
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
        INM=Nl(Branches["Muon"])
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
