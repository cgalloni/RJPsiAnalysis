import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *



trgTable = cms.EDProducer( "TrgBitTableProducer",
                          hltresults = cms.InputTag("TriggerResults::HLT"),
                          l1results  = cms.InputTag("gtStage2Digis::RECO"),
                          #add interesting paths
                          paths      = cms.vstring(
                              "HLT_DoubleMu4_JpsiTrk_Displaced",
                              "HLT_Dimuon0_Jpsi3p5_Muon2",
                              "HLT_DoubleMu4_LowMassNonResonantTrk_Displaced",
                              "HLT_DoubleMu4_PsiPrimeTrk_Displaced",
                              "HLT_DoubleMu4_3_Jpsi_v",
                              "HLT_Dimuon0_Jpsi_NoVertexing",
                              "HLT_DoubleMu4_3_Jpsi_Displaced_v",              
                              "HLT_Ele32_WPTight_Gsf",
                              "HLT_PFMET120_PFMHT120_IDTight",
                              "HLT_Mu8_v",
                              "HLT_Dimuon0_Jpsi_L1_4R_0er1p5R",
                              "HLT_IsoMu24_v",
),

                           #add interesting seeds
                           seeds     = cms.vstring(
                               "L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9",
                               "L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9",
                               "L1_SingleMu5",
                               "L1_SingleMu7",
                               "L1_SingleMu22",
                               "L1_DoubleMu0er1p6_dEta_Max1p8_OS",
                               "L1_DoubleMu_10_0_dEta_Max1p8",
                               "L1_DoubleMu0er1p4_dEta_Max1p8_OS", 
                               "L1_DoubleMu_11_4", 
                               "L1_DoubleMu_12_5", 
                               "L1_DoubleMu_13_6", 
                               "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", 
                               "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4",  
                               "L1_DoubleMu4p5_SQ_OS_dR_Max1p2",  
                               "L1_DoubleMu4_SQ_OS_dR_Max1p2", 
                           ),
                            
)

trgTables = cms.Sequence(trgTable)



