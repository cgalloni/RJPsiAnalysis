#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"
constexpr bool debug = false;

class DiMuonBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<pat::Muon> MuonCollection;
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit DiMuonBuilder(const edm::ParameterSet &cfg):
    mu1_selection_{cfg.getParameter<std::string>("muon1Selection")},
    mu2_selection_{cfg.getParameter<std::string>("muon2Selection")},
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
    src_{consumes<MuonCollection>( cfg.getParameter<edm::InputTag>("src") )},
    ttracks_src_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("transientTracksSrc") )},
    vertexSrc_{consumes<reco::VertexCollection> ( cfg.getParameter<edm::InputTag>("vertexCollection"))},
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )}{
       produces<pat::CompositeCandidateCollection>("muonPairsForB");
       produces<TransientTrackCollection>("dimuonTransientTracks");
    }

  ~DiMuonBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  int  getPVIdx(const reco::VertexCollection*,const reco::TransientTrack&) const;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
  const StringCutObjectSelector<pat::Muon> mu1_selection_; // cut on leading muon
  const StringCutObjectSelector<pat::Muon> mu2_selection_; // cut on sub-leading muon
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-muon before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-muon after the SV fit
  const edm::EDGetTokenT<MuonCollection> src_;
  const edm::EDGetTokenT<TransientTrackCollection> ttracks_src_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;
};

void DiMuonBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &) const {

  //input
  edm::Handle<MuonCollection> muons;
  evt.getByToken(src_, muons);
  
  edm::Handle<TransientTrackCollection> ttracks;
  evt.getByToken(ttracks_src_, ttracks);

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);

  edm::Handle<reco::VertexCollection> thePrimaryVerticesHandle;
  evt.getByToken(vertexSrc_, thePrimaryVerticesHandle);
  const reco::VertexCollection* vertices = thePrimaryVerticesHandle.product();


  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_value(new pat::CompositeCandidateCollection());
  std::unique_ptr<TransientTrackCollection> dimuon_tt(new TransientTrackCollection);
  if(debug) std::cout<<"DIMUON builder "<<std::endl;
  for(size_t mu1_idx = 0; mu1_idx < muons->size(); ++mu1_idx) {
    edm::Ptr<pat::Muon> mu1_ptr(muons, mu1_idx);
    if(!mu1_selection_(*mu1_ptr)) continue; 
    int isMuonFromJpsi_dimuon0_1 = mu1_ptr->userInt("isMuonFromJpsi_dimuon0Trg");
    int isMuonFromJpsi_jpsiTrk_1 = mu1_ptr->userInt("isMuonFromJpsi_jpsiTrkTrg");
    int isMuonFromJpsi_jpsiTrk_PsiPrime_1 = mu1_ptr->userInt("isMuonFromJpsi_jpsiTrk_PsiPrimeTrg");
    int isMuonFromJpsi_jpsiTrk_NonResonant_1 = mu1_ptr->userInt("isMuonFromJpsi_jpsiTrk_NonResonantTrg");
    int isDimuon0Trg1 = mu1_ptr->userInt("isDimuon0Trg");
    int isJpsiTrkTrg1 = mu1_ptr->userInt("isJpsiTrkTrg");
    int isJpsiTrk_PsiPrimeTrg1 = mu1_ptr->userInt("isJpsiTrk_PsiPrimeTrg");
    int isJpsiTrk_NonResonantTrg1 = mu1_ptr->userInt("isJpsiTrk_NonResonantTrg");
    int isDimuon0_jpsi_Trg1 = mu1_ptr->userInt("isDimuon0_jpsi_Trg");
    int isDimuon43_jpsi_displaced_Trg1 = mu1_ptr->userInt("isDimuon43_jpsi_displaced_Trg");

    if(debug) std::cout<< "mu1 "<<mu1_ptr->pt()<<" isMuonFromJpsi_jpsiTrk_1 "<<isMuonFromJpsi_jpsiTrk_1<<" isJpsiTrkTrg1 "<<isJpsiTrkTrg1<<std::endl;

    for(size_t mu2_idx = mu1_idx + 1; mu2_idx < muons->size(); ++mu2_idx) {
      edm::Ptr<pat::Muon> mu2_ptr(muons, mu2_idx);
      if(!mu2_selection_(*mu2_ptr)) continue;
      // Form pairs only with triggered muons
      int isMuonFromJpsi_dimuon0_2 = mu2_ptr->userInt("isMuonFromJpsi_dimuon0Trg");
      int isMuonFromJpsi_jpsiTrk_2 = mu2_ptr->userInt("isMuonFromJpsi_jpsiTrkTrg");
      int isMuonFromJpsi_jpsiTrk_PsiPrime_2 = mu2_ptr->userInt("isMuonFromJpsi_jpsiTrk_PsiPrimeTrg");
      int isMuonFromJpsi_jpsiTrk_NonResonant_2 = mu2_ptr->userInt("isMuonFromJpsi_jpsiTrk_NonResonantTrg");
      int isDimuon0Trg2 = mu2_ptr->userInt("isDimuon0Trg");
      int isJpsiTrkTrg2 = mu2_ptr->userInt("isJpsiTrkTrg");
      int isJpsiTrk_PsiPrimeTrg2 = mu2_ptr->userInt("isJpsiTrk_PsiPrimeTrg");
      int isJpsiTrk_NonResonantTrg2 = mu2_ptr->userInt("isJpsiTrk_NonResonantTrg");
      int isDimuon0_jpsi_Trg2 = mu2_ptr->userInt("isDimuon0_jpsi_Trg");
      int isDimuon43_jpsi_displaced_Trg2 = mu2_ptr->userInt("isDimuon43_jpsi_displaced_Trg");

      int dimuon0_trigger = (isDimuon0Trg1 && isMuonFromJpsi_dimuon0_1) && (isDimuon0Trg2 && isMuonFromJpsi_dimuon0_2);
      int jpsitrk_trigger = (isJpsiTrkTrg1 && isMuonFromJpsi_jpsiTrk_1) && (isJpsiTrkTrg2 && isMuonFromJpsi_jpsiTrk_2);
      int jpsitrk_PsiPrime_trigger = (isJpsiTrk_PsiPrimeTrg1 && isMuonFromJpsi_jpsiTrk_PsiPrime_1) && (isJpsiTrk_PsiPrimeTrg2 && isMuonFromJpsi_jpsiTrk_PsiPrime_2);
      int jpsitrk_NonResonant_trigger = (isJpsiTrk_NonResonantTrg1 && isMuonFromJpsi_jpsiTrk_NonResonant_1) && (isJpsiTrk_NonResonantTrg2 && isMuonFromJpsi_jpsiTrk_NonResonant_2);
      int dimuon0_jpsi_trigger = isDimuon0_jpsi_Trg1 && isDimuon0_jpsi_Trg2 ;
      int dimuon43_jpsi_displaced_trigger = isDimuon43_jpsi_displaced_Trg1 && isDimuon43_jpsi_displaced_Trg2;


      // if(isJpsiTrk_NonResonantTrg1) std::cout << "DimuonBuilder::isJpsiTrk_NonResonantTrg1" << std::endl;
      // if(isJpsiTrk_NonResonantTrg2) std::cout << "DimuonBuilder::isJpsiTrk_NonResonantTrg2" << std::endl;
      // if(isMuonFromJpsi_jpsiTrk_NonResonant_1) std::cout << "DimuonBuilder::isMuonFromJpsi_jpsiTrk_NonResonant_1" << std::endl;
      // if(isMuonFromJpsi_jpsiTrk_NonResonant_2) std::cout << "DimuonBuilder::isMuonFromJpsi_jpsiTrk_NonResonant_2" << std::endl;

      if(debug) std::cout<< "mu2 "<<mu2_ptr->pt()<<" isMuonFromJpsi_jpsiTrk_2 "<<isMuonFromJpsi_jpsiTrk_2<<" isJpsiTrkTrg2 "<<isJpsiTrkTrg2<<std::endl;
      
      //Trig: HLT Eff if(!jpsitrk_trigger && !dimuon0_trigger && !jpsitrk_PsiPrime_trigger && !jpsitrk_NonResonant_trigger) continue;
      //  std::cout << "++++DimuonBuilder::jpsitrk_NonResonant_trigger" << std::endl;
      
      pat::CompositeCandidate muon_pair;
      muon_pair.setP4(mu1_ptr->p4() + mu2_ptr->p4());
      muon_pair.setCharge(mu1_ptr->charge() + mu2_ptr->charge());
      muon_pair.addUserFloat("muons12_deltaR", reco::deltaR(*mu1_ptr, *mu2_ptr));
      // Put the muon passing the corresponding selection
      muon_pair.addUserInt("mu1_idx", mu1_idx );
      muon_pair.addUserInt("mu2_idx", mu2_idx );
      muon_pair.addUserInt("isJpsiTrkTrg", jpsitrk_trigger);
      muon_pair.addUserInt("isJpsiTrk_PsiPrimeTrg", jpsitrk_PsiPrime_trigger);
      muon_pair.addUserInt("isJpsiTrk_NonResonantTrg", jpsitrk_NonResonant_trigger);
      muon_pair.addUserInt("isDimuon0Trg", dimuon0_trigger);
      muon_pair.addUserInt("isDimuon0_jpsi_Trg", dimuon0_jpsi_trigger);
      muon_pair.addUserInt("isDimuon43_jpsi_displaced_Trg", dimuon43_jpsi_displaced_trigger);
      // Use UserCands as they should not use memory but keep the Ptr itself
      muon_pair.addUserCand("mu1", mu1_ptr );
      muon_pair.addUserCand("mu2", mu2_ptr );
      if(debug) std::cout<<"l1 "<<mu1_ptr->pt()<<" l2 "<<mu2_ptr->pt()<<" mass "<<muon_pair.mass()<<" deltaR"<<reco::deltaR(*mu1_ptr, *mu2_ptr)<<" dz "<<mu1_ptr->bestTrack()->dz()-mu2_ptr->bestTrack()->dz()<<std::endl;
      if( !pre_vtx_selection_(muon_pair) ) {
	if(debug) std::cout<<"pre vtx selection dies"<<std::endl;

	continue; // before making the SV, cut on the info we have
      }
      KinVtxFitter fitter(
        {ttracks->at(mu1_idx), ttracks->at(mu2_idx)},
        {mu1_ptr->mass(), mu2_ptr->mass()},
        {LEP_SIGMA, LEP_SIGMA} //some small sigma for the particle mass
        );

      muon_pair.addUserFloat("sv_success", fitter.success()); // float??
      muon_pair.addUserFloat("sv_prob", fitter.prob());
      muon_pair.addUserFloat("sv_chi2", fitter.chi2());
      if( !post_vtx_selection_(muon_pair) ) continue;
      //if(debug) std::cout << "vx_: " << fitter.fitted_candidate() << std::endl;
      muon_pair.addUserFloat("sv_position", fitter.fitted_vtx().x()); // float??
      muon_pair.addUserFloat("sv_ndof", fitter.dof()); // float??
      auto fit_p4 = fitter.fitted_p4();
      muon_pair.addUserFloat("fitted_mass", fitter.success() ? fitter.fitted_candidate().mass() : -1);
      //std::cout << "Dimuon mass: " << fitter.fitted_candidate().mass() << std::endl;
      muon_pair.addUserFloat("fitted_massErr", fitter.success() ? sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6)) : -1);
      muon_pair.addUserFloat("vtx_x",fitter.fitted_vtx().x());
      muon_pair.addUserFloat("vtx_y",fitter.fitted_vtx().y());
      muon_pair.addUserFloat("vtx_z",fitter.fitted_vtx().z());
      muon_pair.addUserFloat("vtx_ex",sqrt(fitter.fitted_vtx_uncertainty().cxx()));
      muon_pair.addUserFloat("vtx_ey",sqrt(fitter.fitted_vtx_uncertainty().cyy()));
      muon_pair.addUserFloat("vtx_ez",sqrt(fitter.fitted_vtx_uncertainty().czz()));
      muon_pair.addUserFloat("sv_chi2", fitter.chi2());
      muon_pair.addUserFloat("fitted_pt"  , fit_p4.pt());
      muon_pair.addUserFloat("fitted_eta" , fit_p4.eta());
      muon_pair.addUserFloat("fitted_phi" , fit_p4.phi());
      muon_pair.addUserFloat("fitted_x",fit_p4.x());
      muon_pair.addUserFloat("fitted_y",fit_p4.y());
      muon_pair.addUserFloat("fitted_z",fit_p4.z());
      muon_pair.addUserFloat("jpsi_err00",fitter.fitted_candidate().kinematicParametersError().matrix()(0,0));
      muon_pair.addUserFloat("jpsi_err11",fitter.fitted_candidate().kinematicParametersError().matrix()(1,1));
      muon_pair.addUserFloat("jpsi_err22",fitter.fitted_candidate().kinematicParametersError().matrix()(2,2));
      muon_pair.addUserFloat("jpsi_err01",fitter.fitted_candidate().kinematicParametersError().matrix()(0,1));
      muon_pair.addUserFloat("jpsi_err02",fitter.fitted_candidate().kinematicParametersError().matrix()(0,2));
      muon_pair.addUserFloat("jpsi_err12",fitter.fitted_candidate().kinematicParametersError().matrix()(1,2));
      
      
      muon_pair.addUserFloat(
			       "fitted_cos_theta_2D",
			       cos_theta_2D(fitter, *beamspot, fit_p4)
			       );
      auto lxy = l_xy(fitter, *beamspot);
      muon_pair.addUserFloat("l_xy", lxy.value());
      muon_pair.addUserFloat("l_xy_unc", lxy.error());

      muon_pair.addUserFloat("vtx_x",muon_pair.vx());
      muon_pair.addUserFloat("vtx_y",muon_pair.vy());
      muon_pair.addUserFloat("vtx_z",muon_pair.vz());

      muon_pair.addUserFloat("vtx_ex",sqrt(fitter.fitted_vtx_uncertainty().cxx()));
      muon_pair.addUserFloat("vtx_ey",sqrt(fitter.fitted_vtx_uncertainty().cyy()));
      muon_pair.addUserFloat("vtx_ez",sqrt(fitter.fitted_vtx_uncertainty().czz()));

      muon_pair.addUserFloat("fitted_mu1_pt" , fitter.daughter_p4(0).pt());
      muon_pair.addUserFloat("fitted_mu1_eta", fitter.daughter_p4(0).eta());
      muon_pair.addUserFloat("fitted_mu1_phi", fitter.daughter_p4(0).phi());
      muon_pair.addUserFloat("fitted_mu2_pt" , fitter.daughter_p4(1).pt());
      muon_pair.addUserFloat("fitted_mu2_eta", fitter.daughter_p4(1).eta());
      muon_pair.addUserFloat("fitted_mu2_phi", fitter.daughter_p4(1).phi());




      muon_pair.addUserFloat(
                               "cos_theta_2D",
                               cos_theta_2D(fitter, *beamspot, muon_pair.p4())
                               );
     
      // if needed, add here more stuff

      // cut on the SV info
      // const reco::TransientTrack& fitted_candidate_ttrk()
      if(!fitter.fitted_candidate_ttrk().isValid()) continue;
      const reco::TransientTrack& dimuonTT = fitter.fitted_candidate_ttrk();
      int pvIdx = getPVIdx(vertices, dimuonTT);
      

      dimuon_tt->emplace_back(fitter.fitted_candidate_ttrk());
      //ret_value->push_back(muon_pair);
      ret_value->emplace_back(muon_pair);
      ret_value->back().addUserInt("muonpair_fromdimuon0", dimuon0_trigger);
      ret_value->back().addUserInt("muonpair_fromdimuon0_jpsi", dimuon0_jpsi_trigger);
      ret_value->back().addUserInt("muonpair_fromdimuon43_jpsi_displaced", dimuon43_jpsi_displaced_trigger);
      ret_value->back().addUserInt("muonpair_fromjpsitrk", jpsitrk_trigger);
      ret_value->back().addUserInt("muonpair_fromjpsitrk_PsiPrime", jpsitrk_PsiPrime_trigger);
      ret_value->back().addUserInt("muonpair_fromjpsitrk_NonResonant", jpsitrk_NonResonant_trigger);
      ret_value->back().addUserInt("pvIdx", pvIdx);
      
    }
  }
  
  evt.put(std::move(ret_value), "muonPairsForB");
  evt.put(std::move(dimuon_tt), "dimuonTransientTracks");
}
int DiMuonBuilder::getPVIdx(const reco::VertexCollection* vertices,const reco::TransientTrack& dimuonTT) const
{
    double dzMin = 1000000.;
    reco::Vertex bestVertex;
    int pvIdx = 0;
    //const reco::VertexCollection* vertices = thePrimaryVerticesHandle.product();
    for(size_t i = 0; i < vertices->size() ; i++)
    {
      reco::Vertex primVertex = vertices->at(i);
      //std::cout << "prim vertex z: " << primVertex->z() << std::endl;
      if (abs(dzMin) > abs(dimuonTT.track().dz(primVertex.position())))
      {
        bestVertex = primVertex;
        pvIdx = i;
        //bestVertex = primVertex;
        dzMin = dimuonTT.track().dz(primVertex.position());
      }
    }
    if(debug) std::cout<< "Best vertex x: " << bestVertex.x() << std::endl;
    if(debug) std::cout<< "Best vertex id: " << pvIdx << std::endl;
  return pvIdx;
}

DEFINE_FWK_MODULE(DiMuonBuilder);
