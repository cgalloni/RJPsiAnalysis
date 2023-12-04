import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *



electronTrgSelector = cms.EDProducer("ElectronTriggerSelector",
                                 electronCollection = cms.InputTag("slimmedElectrons"), #same collection as in NanoAOD                                                           
                                 bits = cms.InputTag("TriggerResults","","HLT"),
                                 prescales = cms.InputTag("patTrigger"),
                                 objects = cms.InputTag("slimmedPatTrigger"),
                                 vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 
                                 ##for the output trigger matched collection
                                 maxdR_matching = cms.double(0.1),
                                 
                                 ## for the output selected collection (tag + all compatible in dZ)
                                 dzForCleaning_wrtTrgElectron = cms.double(1.),

                                 ptMin = cms.double(25),
                                 absEtaMax = cms.double(2.4),
                                 # keeps only electrons with at soft Quality flag
                                 softElectronsOnly = cms.bool(False)
                             )

countTrgElectrons = cms.EDFilter("PATCandViewCountFilter",
                            minNumber = cms.uint32(0),
                            maxNumber = cms.uint32(999999),
                            src = cms.InputTag("electronTrgSelector", "trgElectrons")
                         )


electronBParkTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("electronTrgSelector:SelectedElectrons"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Electron"),
    doc  = cms.string("slimmedElectrons for BPark after basic selection"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the electrons
    variables = cms.PSet(CandVars,
        ptErr   = Var("bestTrack().ptError()", float, doc = "ptError of the electron track", precision=6),
        xErr = Var("bestTrack().covariance(0,0)", float, doc = "xError of the electron track", precision=10),
        yErr = Var("bestTrack().covariance(1,1)", float, doc = "xError of the electron track", precision=10),
        zErr = Var("bestTrack().covariance(2,2)", float, doc = "xError of the electron track", precision=10),
        xyErr = Var("bestTrack().covariance(0,1)", float, doc = "xError of the electron track", precision=10),
        yzErr = Var("bestTrack().covariance(0,2)", float, doc = "xError of the electron track", precision=10),
        xzErr = Var("bestTrack().covariance(1,2)", float, doc = "xError of the electron track", precision=10),
        ## All the following properties will be calculated later with the PV we select.
        #G: dz = Var("dB('PVDZ')",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
        #G: dzErr = Var("abs(edB('PVDZ'))",float,doc="dz uncertainty, in cm",precision=6),
        #G: dxy = Var("dB('PV2D')",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
        #G: dxyErr = Var("edB('PV2D')",float,doc="dxy uncertainty, in cm",precision=6),
        #G: vx = Var("vx()",float,doc="x coordinate of vertex position, in cm",precision=6),
        #G: vy = Var("vy()",float,doc="y coordinate of vertex position, in cm",precision=6),
        #G: vz = Var("vz()",float,doc="z coordinate of vertex position, in cm",precision=6),
        #G: ip3d = Var("abs(dB('PV3D'))",float,doc="3D impact parameter wrt first PV, in cm",precision=10),
        #G: sip3d = Var("abs(dB('PV3D')/edB('PV3D'))",float,doc="3D impact parameter significance wrt first PV",precision=10),
#        segmentComp   = Var("segmentCompatibility()", float, doc = "electron segment compatibility", precision=14), # keep higher precision since people have cuts with 3 digits on this
#        nStations = Var("numberOfMatchedStations", int, doc = "number of matched stations with default arbitration (segment & track)"),
        #nTrackerLayers = Var("innerTrack().hitPattern().trackerLayersWithMeasurement()", int, doc = "number of layers in the tracker"),
#        pfRelIso03_chg = Var("pfIsolationR03().sumChargedHadronPt/pt",float,doc="PF relative isolation dR=0.3, charged component"),
#        tightCharge = Var("?(electronBestTrack().ptError()/electronBestTrack().pt() < 0.2)?2:0",int,doc="Tight charge criterion using pterr/pt of electronBestTrack (0:fail, 2:pass)"),
        #isPFcand = Var("isPFElectron",bool,doc="electron is PF candidate"),
                         

 
#        mediumPromptId = Var("passed('CutBasedIdMediumPrompt')",bool,doc="cut-based ID, medium prompt WP"),
     #        softMvaId = Var("passed('SoftMvaId')",bool,doc="soft MVA ID"),
#        highPtId = Var("?passed('CutBasedIdGlobalHighPt')?2:passed('CutBasedIdTrkHighPt')","uint8",doc="high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)"),
        #pfIsoId = Var("passed('PFIsoVeryLoose')+passed('PFIsoLoose')+passed('PFIsoMedium')+passed('PFIsoTight')+passed('PFIsoVeryTight')+passed('PFIsoVeryVeryTight')","uint8",doc="PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)"),
        #tkIsoId = Var("?passed('TkIsoTight')?2:passed('TkIsoLoose')","uint8",doc="TkIso ID (1=TkIsoLoose, 2=TkIsoTight)"),
#        mvaId = Var("passed('MvaLoose')+passed('MvaMedium')+passed('MvaTight')","uint8",doc="Mva ID from miniAOD selector (1=MvaLoose, 2=MvaMedium, 3=MvaTight)"),
#        miniIsoId = Var("passed('MiniIsoLoose')+passed('MiniIsoMedium')+passed('MiniIsoTight')+passed('MiniIsoVeryTight')","uint8",doc="MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)"),
#        multiIsoId = Var("?passed('MultiIsoMedium')?2:passed('MultiIsoLoose')","uint8",doc="MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium)"),
        #triggerIdLoose = Var("passed('TriggerIdLoose')",bool,doc="TriggerIdLoose ID"),
#        inTimeElectron = Var("passed('InTimeElectron')",bool,doc="inTimeElectron ID"),

        isTriggering = Var("userInt('isTriggering')", int,doc="flag the reco electron is also triggering"),

    ),
)


electronsBParkMCMatchForTable = cms.EDProducer("MCMatcher",            # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src         = electronBParkTable.src,                         # final reco collection
    matched     = cms.InputTag("finalGenParticlesBPark"),     # final mc-truth particle collection
                                               mcPdgId     = cms.vint32(11),                             # one or more PDG ID (13 = mu); absolute values (see below)
    checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
    mcStatus    = cms.vint32(1),                              # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR   = cms.double(0.03),                           # Minimum deltaR for the match
    maxDPtRel   = cms.double(0.5),                            # Minimum deltaPt/Pt for the match
    resolveAmbiguities    = cms.bool(True),                   # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True),                   # False = just match input in order; True = pick lowest deltaR pair first
)

electronTriggerMatchedTable = electronBParkTable.clone(
    src = cms.InputTag("electronTrgSelector:trgElectrons"),
    name = cms.string("TriggerElectron"),
    doc  = cms.string("reco ele matched to triggering ele"),
    variables = cms.PSet(CandVars,
        vx = Var("vx()",float,doc="x coordinate of vertex position, in cm",precision=6),
        vy = Var("vy()",float,doc="y coordinate of vertex position, in cm",precision=6),
        vz = Var("vz()",float,doc="z coordinate of vertex position, in cm",precision=6),
        trgElectronMatchingEle32_index = Var("userInt('trgElectronMatchingEle32_index')", int,doc="index in trigger ele collection"),
   )
)

#electronBParkSequence = cms.Sequence(electronTrgSelector * countTrgElectrons)
electronBParkSequence = cms.Sequence(electronTrgSelector)
electronBParkMC = cms.Sequence(electronBParkSequence )#+ electronsBParkMCMatchForTable + selectedElectronsMCMatchEmbedded + electronBParkMCTable)
#electronBParkMC = cms.Sequence(electronBParkSequence + electronsBParkMCMatchForTable + electronMCMatchSequence + electronBParkMCTable)
electronBParkTables = cms.Sequence(electronBParkTable)
#electronTriggerMatchedTables = cms.Sequence(electronTriggerMatchedTable)
