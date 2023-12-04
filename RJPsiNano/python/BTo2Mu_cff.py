import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var, CandVars
from PhysicsTools.RJPsiNano.common_cff import RJpsiCandVars, ufloat, uint, ubool
from PhysicsTools.RJPsiNano.rjpsi_common_cff import JpsiMuonPairs, BuilderDefaultCfg, TableDefaultVariables, TableDefault, Final2MuonsTableVariables
from PhysicsTools.RJPsiNano.primaryVertices_cff import *

BTo2MuCfg = BuilderDefaultCfg.clone()
#BTo2MuCfg.dileptons             = cms.InputTag('JpsiMuonPairs')
#BTo2MuCfg.leptonTransientTracks = JpsiMuonPairs.transientTracksSrc

BTo2Mu = cms.EDProducer(
    'BTo2MuBuilder',
    BTo2MuCfg,
    muons = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    muonSelection = cms.string('pt > 4.0 && isLooseMuon '),
    srcGen = cms.InputTag("prunedGenParticles"),
    dimuons_fitter = cms.InputTag('JpsiMuonPairs','dimuonFitter')
)

#BTo2MuTableVariables = Final3MuonsTableVariables.clone()
BTo2MuTableVariables = Final2MuonsTableVariables.clone()      
#BTo2MuTableVariables = TableDefaultVariables.clone()   

BTo2MuTable = TableDefault.clone()
BTo2MuTable.src       = cms.InputTag("BTo2Mu")
BTo2MuTable.cut       = cms.string("")
BTo2MuTable.name      = cms.string("BTo2Mu")
BTo2MuTable.doc       = cms.string("BTo2Mu Variable")
BTo2MuTable.singleton = cms.bool(False)
BTo2MuTable.extension = cms.bool(False)
BTo2MuTable.variables = BTo2MuTableVariables

BTo2MuSequence = cms.Sequence(
    (JpsiMuonPairs * pvSelector * BTo2Mu)
)

CountBTo2Mu = cms.EDFilter(
    "PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BTo2Mu")
)
