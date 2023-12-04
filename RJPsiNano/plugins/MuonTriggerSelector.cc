// class to produce 2 pat::MuonCollections
// one matched to the Park triggers
// another fitered wrt the Park triggers

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"

#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <TLorentzVector.h>
#include "helper.h"

using namespace std;

constexpr bool debug = true;
constexpr bool debugTrg = false;

class MuonTriggerSelector : public edm::EDProducer {
public:
  explicit MuonTriggerSelector(const edm::ParameterSet& iConfig);

  ~MuonTriggerSelector() override{};

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  edm::EDGetTokenT<std::vector<pat::Muon>> muonSrc_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;
  edm::EDGetTokenT<GlobalAlgBlkBxCollection> l1results_;

  //for trigger match
  const double maxdR_;

  //for filter wrt trigger
  const double dzTrg_cleaning_;  // selects primary vertex

  const double ptMin_;        // min pT in all muons for B candidates
  const double absEtaMax_;    //max eta ""
  const bool softMuonsOnly_;  //cuts muons without soft ID

  l1t::L1TGlobalUtil* fGtUtil;
};

MuonTriggerSelector::MuonTriggerSelector(const edm::ParameterSet& iConfig)
    : muonSrc_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection"))),
      triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
      triggerObjects_(
          consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("objects"))),
      triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
      vertexSrc_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"))),
      l1results_(consumes<GlobalAlgBlkBxCollection>(iConfig.getParameter<edm::InputTag>("l1results"))),
      maxdR_(iConfig.getParameter<double>("maxdR_matching")),
      dzTrg_cleaning_(iConfig.getParameter<double>("dzForCleaning_wrtTrgMuon")),
      ptMin_(iConfig.getParameter<double>("ptMin")),
      absEtaMax_(iConfig.getParameter<double>("absEtaMax")),
      softMuonsOnly_(iConfig.getParameter<bool>("softMuonsOnly")) {
  // produce 2 collections: trgMuons (tags) and SelectedMuons (probes & tags if survive preselection cuts)
  produces<pat::MuonCollection>("trgMuons");
  produces<pat::MuonCollection>("SelectedMuons");
  produces<TransientTrackCollection>("SelectedTransientMuons");
  fGtUtil = new l1t::L1TGlobalUtil(iConfig,
                                   consumesCollector(),
                                   *this,
                                   iConfig.getParameter<edm::InputTag>("l1results"),
                                   iConfig.getParameter<edm::InputTag>("l1results"));
}

void MuonTriggerSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  edm::Handle<reco::VertexCollection> vertexHandle;
  iEvent.getByToken(vertexSrc_, vertexHandle);
  //const reco::Vertex & PV = vertexHandle->front();

  //// L1 information
  edm::Handle<GlobalAlgBlkBxCollection> l1results;
  iEvent.getByToken(l1results_, l1results);

  if (debug)
    std::cout << " MuonTriggerSelector::produce " << std::endl;

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  const edm::TriggerNames& names = iEvent.triggerNames(*triggerBits);

  std::vector<pat::TriggerObjectStandAlone> triggeringMuons;

  //taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2016#Trigger
  edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  if (debug)
    std::cout << "\n TRIGGER OBJECTS " << std::endl;

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  edm::Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(muonSrc_, muons);

  std::unique_ptr<pat::MuonCollection> trgmuons_out(new pat::MuonCollection);
  std::unique_ptr<pat::MuonCollection> muons_out(new pat::MuonCollection);
  std::unique_ptr<TransientTrackCollection> trans_muons_out(new TransientTrackCollection);

  //Getting HLT, L1 struct
  std::map<std::string, std::vector<std::string>> hlt_l1_map{
      {"HLT_Dimuon0_Jpsi3p5_Muon2_v",
	  {"L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9",
	      "L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"}},
      {"HLT_Mu8_v", {"L1_SingleMu5","L1_SingleMu7"}},
      {"HLT_IsoMu24_v", {"L1_SingleMu22"}},
      {"HLT_DoubleMu4_JpsiTrk_Displaced_v",
	  {"L1_DoubleMu_10_0_dEta_Max1p8",
	    "L1_DoubleMu0er1p6_dEta_Max1p8_OS",
	    "L1_DoubleMu0er1p4_dEta_Max1p8_OS",
	    "L1_DoubleMu_11_4",
	    "L1_DoubleMu_12_5",
	    "L1_DoubleMu_13_6",
	    "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4",
	    "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", 
	    "L1_DoubleMu4p5_SQ_OS_dR_Max1p2",
	    "L1_DoubleMu4_SQ_OS_dR_Max1p2"}},
      {"HLT_DoubleMu4_3_Jpsi_v", {"L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4"}},
      {"HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v",{}},
	//{"L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4",
        //"L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4",
        //"L1_DoubleMu4p5_SQ_OS_dR_Max1p2",
        //"L1_DoubleMu4_SQ_OS_dR_Max1p2"}},
	{"HLT_Dimuon0_Jpsi_NoVertexing_v",{}},
		 //{"L1_DoubleMu0_SQ_OS", "L1_DoubleMu0_SQ"}},
	{"HLT_Dimuon0_Jpsi_v",{}},    
	//    {"L1_DoubleMu0_SQ_OS", "L1_DoubleMu0_SQ"}},
	{"HLT_DoubleMu4_PsiPrimeTrk_Displaced_v",{}},
	//      {"L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4",
        //"L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4",
        //"L1_DoubleMu4p5_SQ_OS_dR_Max1p2",
        //"L1_DoubleMu4_SQ_OS_dR_Max1p2"}},
       {"HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v",{}},
	//{"L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4",
        //"L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4",
        //"L1_DoubleMu4p5_SQ_OS_dR_Max1p2",
        //"L1_DoubleMu4_SQ_OS_dR_Max1p2"}},
      {"HLT_DoubleMu4_3_Jpsi_Displaced_v",
       {"L1_DoubleMu_10_0_dEta_Max1p8",
        "L1_DoubleMu0er1p6_dEta_Max1p8_OS",
        "L1_DoubleMu0er1p4_dEta_Max1p8_OS",
	"L1_DoubleMu_11_4",
	"L1_DoubleMu_12_5", 
	"L1_DoubleMu_13_6", 
	"L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", 
	"L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"}},
	};
 

  // Getting the indexes of the HLT paths
  unsigned int index_dimuon01 = names.triggerIndex("HLT_Dimuon0_Jpsi3p5_Muon2_v5");
  unsigned int index_dimuon02 = names.triggerIndex("HLT_Dimuon0_Jpsi3p5_Muon2_v6");
  unsigned int index_jpsiTrk1 = names.triggerIndex("HLT_DoubleMu4_JpsiTrk_Displaced_v14");
  unsigned int index_jpsiTrk2 = names.triggerIndex("HLT_DoubleMu4_JpsiTrk_Displaced_v15");
  unsigned int index_jpsiTrk3 = names.triggerIndex("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v14");
  unsigned int index_jpsiTrk4 = names.triggerIndex("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v15");
  unsigned int index_jpsiTrk5 = names.triggerIndex("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v14");
  unsigned int index_jpsiTrk6 = names.triggerIndex("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15");
  unsigned int index_doubleMu1 = names.triggerIndex("HLT_DoubleMu4_3_Jpsi_v1");
  unsigned int index_doubleMu2 = names.triggerIndex("HLT_DoubleMu4_3_Jpsi_v2");

  unsigned int index_dimuon0_jpsi_1 = names.triggerIndex("HLT_Dimuon0_Jpsi_NoVertexing_v7");
  unsigned int index_dimuon0_jpsi_2 = names.triggerIndex("HLT_Dimuon0_Jpsi_NoVertexing_v8");
  unsigned int index_dimuon43_jpsi_displaced = names.triggerIndex("HLT_DoubleMu4_3_Jpsi_Displaced_v7");

  unsigned int index_dimuon0_jpsi_4R_1 = names.triggerIndex("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v6");
  unsigned int index_dimuon0_jpsi_4R_2 = names.triggerIndex("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v7");

  unsigned int index_muon_1 = names.triggerIndex("HLT_Mu8_v12");
  unsigned int index_muon_2 = names.triggerIndex("HLT_Mu8_v11");

  bool pass_dimuon01_path = false;
  bool pass_dimuon02_path = false;
  bool pass_jpsiTrk1_path = false;
  bool pass_jpsiTrk2_path = false;
  bool pass_jpsiTrk3_path = false;
  bool pass_jpsiTrk4_path = false;
  bool pass_jpsiTrk5_path = false;
  bool pass_jpsiTrk6_path = false;
  bool pass_doubleMu1_path = false;
  bool pass_doubleMu2_path = false;
  bool pass_dimuon0_jpsi_1_path = false;
  bool pass_dimuon0_jpsi_2_path = false;
  bool pass_dimuon43_jpsi_displaced_path = false;

  bool pass_dimuon0_jpsi_4R_1_path = false;
  bool pass_dimuon0_jpsi_4R_2_path = false;

  bool pass_muon_1_path = false;
  bool pass_muon_2_path = false;

  bool pass_iso_muon_path = false;

  for (int v=0;v<17;v++){
    string name = "HLT_Dimuon0_Jpsi3p5_Muon2_v"+std::to_string(v);
    unsigned int index = names.triggerIndex(name);                                                                                                                
    if( index < triggerBits->size() ) {                                                                                                                               
      if( triggerBits->accept( index ) ) pass_dimuon01_path = true;                                                                                                     
    }
  } 
  for (int v=0;v<17;v++){
    string name = "HLT_DoubleMu4_JpsiTrk_Displaced_v"+std::to_string(v);
    unsigned int index = names.triggerIndex(name);
    if( index < triggerBits->size() ) {
      if( triggerBits->accept( index ) ) pass_jpsiTrk1_path = true;
    }
  }
  for (int v=0;v<17;v++){
    string name = "HLT_DoubleMu4_PsiPrimeTrk_Displaced_v"+std::to_string(v);                                                                                           
    unsigned int index = names.triggerIndex(name);
    if( index < triggerBits->size() ) {
      if( triggerBits->accept( index ) ) pass_jpsiTrk3_path = true;
    }
  }
  for (int v=0;v<17;v++){
    string name = "HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v"+std::to_string(v);
    unsigned int index = names.triggerIndex(name);
    if( index < triggerBits->size() ) {
      if( triggerBits->accept( index ) ) pass_jpsiTrk5_path = true;
    }
  }
  for (int v=0;v<17;v++){
    string name = "HLT_DoubleMu4_3_Jpsi_v"+std::to_string(v);
    unsigned int index = names.triggerIndex(name);
    if( index < triggerBits->size() ) {
      if( triggerBits->accept( index ) ) {
	pass_doubleMu1_path = true;
	std::cout << "HLT name " << name << std::endl; 
	std::cout << "- fired" << std::endl;        
      }
    }
  }
  
  for (int v=0;v<17;v++){
    string name = "HLT_Dimuon0_Jpsi_NoVertexing_v"+std::to_string(v);
    unsigned int index = names.triggerIndex(name);
    if( index < triggerBits->size() ) {
      if( triggerBits->accept( index ) ) pass_dimuon0_jpsi_1_path = true;
    }
  }

  for (int v=0;v<17;v++){
    string name = "HLT_DoubleMu4_3_Jpsi_Displaced_v"+std::to_string(v);
    unsigned int index = names.triggerIndex(name);
    if( index < triggerBits->size() ) {
      std::cout << "HLT name " << name << std::endl;
      if( triggerBits->accept( index ) ) {
	pass_dimuon43_jpsi_displaced_path = true;
	std::cout << "- fired" << std::endl;
      }
    }
  }
  
  for (int v=0;v<17;v++){
    string name = "HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v"+std::to_string(v);
    unsigned int index = names.triggerIndex(name);
    if( index < triggerBits->size() ) {
      if( triggerBits->accept( index ) ) pass_dimuon0_jpsi_4R_1_path = true;
    }
  }
  for (int v=0;v<17;v++){
    string name = "HLT_Mu8_v"+std::to_string(v);
    unsigned int index = names.triggerIndex(name);
    if( index < triggerBits->size() ) {


      if( triggerBits->accept( index ) ) pass_muon_1_path = true;
      std::cout << "HLT name " << name << std::endl;      
      std::cout << "- fired" << std::endl;
    }
  }
  for (int v=0;v<17;v++){
    string name = "HLT_IsoMu24_v"+std::to_string(v);
    unsigned int index = names.triggerIndex(name);
    std::cout << "HLT name " << name << std::endl;
    if( index < triggerBits->size() ) {
      std::cout << "- fired" << std::endl;
      if( triggerBits->accept( index ) ) pass_iso_muon_path = true;
    }
  }
       
  /*
  if (index_dimuon01 != triggerBits->size()) {
    pass_dimuon01_path = triggerBits->accept(index_dimuon01);
  }
  if (index_dimuon02 != triggerBits->size()) {
    pass_dimuon02_path = triggerBits->accept(index_dimuon02);
  }
  if (index_jpsiTrk1 != triggerBits->size())
    pass_jpsiTrk1_path = triggerBits->accept(index_jpsiTrk1);
  if (index_jpsiTrk2 != triggerBits->size())
    pass_jpsiTrk2_path = triggerBits->accept(index_jpsiTrk2);
  if (index_jpsiTrk3 != triggerBits->size())
    pass_jpsiTrk3_path = triggerBits->accept(index_jpsiTrk3);
  if (index_jpsiTrk4 != triggerBits->size())
    pass_jpsiTrk4_path = triggerBits->accept(index_jpsiTrk4);
  if (index_jpsiTrk5 != triggerBits->size())
    pass_jpsiTrk5_path = triggerBits->accept(index_jpsiTrk5);
  if (index_jpsiTrk6 != triggerBits->size())
    pass_jpsiTrk6_path = triggerBits->accept(index_jpsiTrk6);

  if (index_doubleMu1 != triggerBits->size())
    pass_doubleMu1_path = triggerBits->accept(index_doubleMu1);
  if (index_doubleMu2 != triggerBits->size())
    pass_doubleMu2_path = triggerBits->accept(index_doubleMu2);

  if (index_dimuon0_jpsi_1 != triggerBits->size())
    pass_dimuon0_jpsi_1_path = triggerBits->accept(index_dimuon0_jpsi_1);
  if (index_dimuon0_jpsi_2 != triggerBits->size())
    pass_dimuon0_jpsi_2_path = triggerBits->accept(index_dimuon0_jpsi_2);
  if (index_dimuon43_jpsi_displaced != triggerBits->size())
    pass_dimuon43_jpsi_displaced_path = triggerBits->accept(index_dimuon43_jpsi_displaced);

  if (index_dimuon0_jpsi_4R_1 != triggerBits->size())
    pass_dimuon0_jpsi_4R_1_path = triggerBits->accept(index_dimuon0_jpsi_4R_1);
  if (index_dimuon0_jpsi_4R_2 != triggerBits->size())
    pass_dimuon0_jpsi_4R_2_path = triggerBits->accept(index_dimuon0_jpsi_4R_2);

  if (index_muon_1 != triggerBits->size())
    pass_muon_1_path = triggerBits->accept(index_muon_1);
  if (index_muon_2 != triggerBits->size())
    pass_muon_2_path = triggerBits->accept(index_muon_2);
  */
  /*
  for( std::vector<std::string>::const_iterator it = hltPathsOfInterest.begin();
       it != hltPathsOfInterest.end(); ++it ) {
    int fired = 0;
    unsigned int index = names.triggerIndex(*it);
    if( index < triggerResults->size() ) {
      if( triggerResults->accept( index ) ) fired = 1;
    } else {
      edm::LogInfo("PromptAna_TriggerInfo") << "Requested HLT path \"" << (*it) << "\" does not exist";
    }
    hltresults->push_back( fired );
    */

  //Getting HLT, psHLT struct
  std::map<std::string, std::vector<int>> hlt_psHLT_map;
  //Initialize ps to 0 for all hlt paths
  for (auto& [key_hlt, value_ps] : hlt_l1_map) {
    //std::cout << "Inititalize hlt_psHLT_map key_hlt" << key_hlt << std::endl;
    hlt_psHLT_map[key_hlt].push_back(0);
  }
  //std::cout<< "After map creation hlt_psHLT_map " << std::endl;
  //for (auto& [key_hlt, value_ps] : hlt_psHLT_map) {
    //std::cout<< "key_hlt "<< key_hlt;
  //  for (auto& ps : value_ps) {
      //std::cout<< " ps "<< ps;
  //  }
    //std::cout<< std::endl;
  //}

  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    //for (unsigned int it = 0; it < hlt_l1_map.size(); it++){
    for (auto& [key_hlt, value_l1] : hlt_l1_map) {
      if (names.triggerName(i).find(key_hlt) != std::string::npos) {
        //std::cout << "key_hlt "<< key_hlt<< " Trig name " << names.triggerName(i) ;
        if (triggerBits->accept(i)) {
          hlt_psHLT_map[key_hlt].push_back(triggerPrescales->getPrescaleForIndex(i));
          //  std::cout << " passed  prescale " <<triggerPrescales->getPrescaleForIndex(i)<<std::endl;
        } else
          hlt_psHLT_map[key_hlt].push_back(0);
        //std::cout << "not passed  prescale " << 0 <<std::endl;
      }
    }
  }

  //L1
  std::map<std::string, std::vector<int>> hlt_psL1_map;
  //Initialize ps to 0 for all hlt paths
  for (auto& [key_hlt, value_ps] : hlt_l1_map) {
    //    std::cout << "Inititalize hlt_psL1_map key_hlt" << key_hlt << std::endl;
    hlt_psL1_map[key_hlt].push_back(0);
  }
  if (l1results.isValid()) {
    fGtUtil->retrieveL1(iEvent, iSetup, l1results_);
    const std::vector<std::pair<std::string, bool>> finalDecisions = fGtUtil->decisionsFinal();
    const std::vector<std::pair<std::string, int>> prescales = fGtUtil->prescales();
    for (unsigned int i = 0; i < finalDecisions.size(); ++i) {
      std::string name = (finalDecisions.at(i)).first;
      if (name == "NULL")
        continue;
      for (auto& [key_hlt, value_l1] : hlt_l1_map) {
	//std::cout << " Name "<<name << " Key_hlt "<< key_hlt<< std::endl ;
        // for (unsigned int it = 0; it < l1Table_.size(); it++){
        for (auto& it_l1 : value_l1) {
	  //std::cout << " it_l1 "<< it_l1<< " name.compare(it_l1) " << name.compare(it_l1) ;
          if (name.compare(it_l1) == 0) {
            //std::cout << " it_l1 "<< it_l1 ;
            bool resultFin = (finalDecisions.at(i)).second;
            if (resultFin) {
              //std::cout << " is passed  prescale " <<(prescales.at(i)).second<< std::endl;
              hlt_psL1_map[key_hlt].push_back((prescales.at(i)).second);
            } else {
	      //std::cout << " not passed  prescale " <<(prescales.at(i)).second<< std::endl;
	      //std::cout << " not passed  prescale set to 0"<< std::endl;
	     hlt_psL1_map[key_hlt].push_back(0);
            }
          }
        }
      }
    }
  }

  //  std::cout << "removing zeros"<<std::endl;
  // for (auto&  [key_hlt , value_ps] :hlt_psHLT_map){
  //   std::cout<< "key_hlt "<< key_hlt;
  // for (auto&  ps :value_ps){
  //  std::cout<< " ps "<< ps;
  // }
  // std::cout<< std::endl;
  // }
  for (auto& [key_hlt, value_ps] : hlt_psHLT_map) {
    value_ps.erase(remove_if(value_ps.begin(), value_ps.end(), [](int i) { return i < 1; }), value_ps.end());
  }
  // std::cout << "zeros removed "<<std::endl;
  // for (auto&  [key_hlt , value_ps] :hlt_psHLT_map){
  //   std::cout<< "key_hlt "<< key_hlt;
  //   for (auto&  ps :value_ps){
  //     std::cout<< "ps "<< ps;
  //   }
  //   std::cout<< " Min ps "<<  *min_element(hlt_psHLT_map[key_hlt].begin(), hlt_psHLT_map[key_hlt].end());
  //   std::cout<< std::endl;

  // }

  for (auto& [key_hlt, value_ps] : hlt_psL1_map) {
    value_ps.erase(remove_if(value_ps.begin(), value_ps.end(), [](int i) { return i < 1; }), value_ps.end());
  }

  // std::cout << "zeros removed "<<std::endl;
  // for (auto&  [key_l1 , value_ps] :hlt_psL1_map){
  //   std::cout<< "key_l1 "<< key_l1;
  //   for (auto&  ps :value_ps){
  //     std::cout<< "ps "<< ps;
  //   }
  //   std::cout<< " Min ps "<<  *min_element(hlt_psL1_map[key_l1].begin(), hlt_psL1_map[key_l1].end());
  //   std::cout<< std::endl;

  // }

  //  if(debug) std::cout << "pass_dimuuon0 " << pass_dimuon0_path<<" pass_trk "<<pass_jpsiTrk_path<<std::endl;

  //if(pass_dimuon0_path) std::cout << "pass_dimuuon0" << std::endl;
  std::vector<bool> jpsiFromMuon_fromDimuon0_flags;
  std::vector<bool> jpsiFromMuon_fromJpsiTrk_flags;
  std::vector<bool> jpsiFromMuon_fromJpsiTrk_PsiPrime_flags;
  std::vector<bool> jpsiFromMuon_fromJpsiTrk_NonResonant_flags;
  std::vector<bool> jpsiFromMuon_fromDoubleMu_flags;
  std::vector<bool> dimuon0Flags;
  std::vector<bool> jpsiTrkFlags;
  std::vector<bool> jpsiTrk_PsiPrimeFlags;
  std::vector<bool> jpsiTrk_NonResonantFlags;

  std::vector<bool> doubleMuFlags;

  std::vector<bool> dimuon0_jpsi_flags;
  std::vector<bool> dimuon43_jpsi_displaced_flags;

  std::vector<bool> dimuon0_jpsi_4R_flags;

  std::vector<bool> muon_flags;
  std::vector<bool> iso_muon_flags;

  if (pass_dimuon01_path || pass_dimuon02_path || pass_jpsiTrk1_path || pass_jpsiTrk2_path || pass_jpsiTrk3_path ||
      pass_jpsiTrk4_path || pass_jpsiTrk5_path || pass_jpsiTrk6_path || pass_doubleMu1_path || pass_doubleMu2_path ||
      pass_dimuon0_jpsi_1_path || pass_dimuon0_jpsi_2_path || pass_dimuon43_jpsi_displaced_path ||
      pass_dimuon0_jpsi_4R_1_path || pass_dimuon0_jpsi_4R_2_path || pass_muon_1_path || pass_muon_2_path || pass_iso_muon_path) {
    //if(pass_dimuon01_path || pass_jpsiTrk1_path || pass_jpsiTrk2_path || pass_jpsiTrk3_path || pass_jpsiTrk4_path || pass_jpsiTrk5_path || pass_jpsiTrk6_path) {
    for (pat::TriggerObjectStandAlone obj :
         *triggerObjects) {  // note: not "const &" since we want to call unpackPathNames
      obj.unpackFilterLabels(iEvent, *triggerBits);
      obj.unpackPathNames(names);
      
      if ( (obj.collection()!="hltIterL3MuonCandidates::HLT")&&(obj.collection()!="hltL3MuonCandidates::HLT")) continue;
      //for (auto fl :obj.filterLabels()){
      //	std::cout << "filter label " << fl <<std::endl; 
      //}
      
      bool muonFromJpsi_fromDimuon0Path = false;
      bool muonFromJpsi_fromJpsiTrkPath = false;
      bool muonFromJpsi_fromJpsiTrk_PsiPrimePath = false;
      bool muonFromJpsi_fromJpsiTrk_NonResonantPath = false;
      bool muonFromJpsi_fromDoubleMuPath = false;
      bool dimuon0_seed = false;
      bool jpsitrk_seed = false;
      bool jpsitrk_PsiPrime_seed = false;
      bool jpsitrk_NonResonant_seed = false;

      bool doubleMu_seed = false;

      bool dimuon0_jpsi_seed = false;
      bool dimuon43_jpsi_displaced_seed = false;
      bool dimuon0_jpsi_4R_seed = false;

      bool muon_seed = false;
      bool iso_muon_seed = false;

      //if(pass_jpsiTrk5_path || pass_jpsiTrk6_path)
      //  if(debugTrg) std::cout << "pass_nonREsonant path" << std::endl;

      //Matching for lept trigger: HLT_Dimuon0_Jpsi3p5_Muon2_v4

      if (obj.hasFilterLabel("hltVertexmumuFilterJpsiMuon3p5"))
        muonFromJpsi_fromDimuon0Path = true;
      if (obj.hasFilterLabel("hltTripleMuL3PreFiltered222")) {
        dimuon0_seed = true;
      }
      // Matching for had trigger: HLT_DoubleMu4_JpsiTrk_Displaced_v15
      if (obj.hasFilterLabel("hltDisplacedmumuFilterDoubleMu4Jpsi"))
        muonFromJpsi_fromJpsiTrkPath = true;
      if (obj.hasFilterLabel("hltJpsiTkVertexFilter")) {
        jpsitrk_seed = true;  // it shouldn't be used since this matches a track
      }

      //Matching for the possible reference trigger: HLT_Dimuon0_Jpsi_NoVertexing_v
      if (obj.hasFilterLabel("hltDimuon0JpsiL3Filtered")) {
        dimuon0_jpsi_seed = true;
      }
      //Matching for the possible reference trigger: HLT_DoubleMu4_3_Jpsi_Displaced_v
      if (obj.hasFilterLabel("hltDisplacedmumuFilterDoubleMu43Jpsi")) {
        dimuon43_jpsi_displaced_seed = true;
      
	std::cout << " ---- Insidie DoubleMu4_3_Jpsi_Displaced"  <<std::endl;
	for (auto fl :obj.filterLabels()){
	  std::cout << "filter label " << fl <<std::endl;
	}
      }
      //Matching for the possible reference trigger: HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v
      if (obj.hasFilterLabel("hltDisplacedmumuFilterDimuon0JpsiL1s4R0er1p5R")) {
        dimuon0_jpsi_4R_seed = true;
      }

      if (obj.hasFilterLabel("hltDisplacedmumuFilterDoubleMu4PsiPrime")) {
        muonFromJpsi_fromJpsiTrk_PsiPrimePath = true;
      }
      if (obj.hasFilterLabel("hltPsiPrimeTkVertexFilter")) {
        jpsitrk_PsiPrime_seed = true;
      }

      if (obj.hasFilterLabel("hltDisplacedmumuFilterDoubleMu4LowMassNonResonant")) {
        muonFromJpsi_fromJpsiTrk_NonResonantPath = true;
        //if(debugTrg) std::cout << "hltDisplacedmumuFilterDoubleMu4NonResonant" << std::endl;
      }
      if (obj.hasFilterLabel("hltLowMassNonResonantTkVertexFilter")) {
        jpsitrk_NonResonant_seed = true;
      }

      // HLT_DoubleMu4_3_Jpsi  :
      if (obj.hasFilterLabel("hltmumuFilterDoubleMu43Jpsi"))
        muonFromJpsi_fromDoubleMuPath = true;

      // HLT_DoubleMu4_3_Jpsi and // HLT_DoubleMu4_3_Jpsi_Displaced filter on muons   :       
      if (obj.hasFilterLabel("hltDoubleMu43JpsiDisplacedL3Filtered")){
        doubleMu_seed = true;
	//std::cout << " ---- Inside HLT_DoubleMu4_3_Jpsi and HLT_DoubleMu4_3_Jpsi_Displaced filter on muons "  <<std::endl;
	// for (auto fl :obj.filterLabels()){
	//  std::cout << "filter label " << fl <<std::endl;
        //}
      }
      if (obj.hasFilterLabel("hltL3fL1sMu5L1f0L2f5L3Filtered8")){
       
	//	std::cout << " ---- Insidie Mu8"  <<std::endl;                                                                                                                                      
        //for (auto fl :obj.filterLabels()){
	//  std::cout << "filter label " << fl <<std::endl;
        //}
        muon_seed = true;
      }

       
       if (obj.hasFilterLabel("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07") //2017- 2018b	   
	   || obj.hasFilterLabel("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09")// 2016 b 

	   ){
	 std::cout << " ---- Insidie IsoMu24"  <<std::endl;
	 for (auto fl :obj.filterLabels()){
	   std::cout << "filter label " << fl <<std::endl;
	 }

	 iso_muon_seed = true;
       }
      //for each triggered muon I know which trigger it passes
      dimuon0Flags.push_back(dimuon0_seed);                          //unpaired muon passes the dimuon0
      jpsiTrkFlags.push_back(jpsitrk_seed);                          //unpaired muon passes the jpsitrk
      jpsiTrk_PsiPrimeFlags.push_back(jpsitrk_PsiPrime_seed);        //unpaired muon passes the psiprimetrk
      jpsiTrk_NonResonantFlags.push_back(jpsitrk_NonResonant_seed);  //unpaired muon passes the nonResonantrk
      doubleMuFlags.push_back(doubleMu_seed);

      jpsiFromMuon_fromDimuon0_flags.push_back(muonFromJpsi_fromDimuon0Path);  //the muon is from the jpsi
      jpsiFromMuon_fromJpsiTrk_flags.push_back(muonFromJpsi_fromJpsiTrkPath);  //the muon is from the jpsi
      jpsiFromMuon_fromJpsiTrk_PsiPrime_flags.push_back(
          muonFromJpsi_fromJpsiTrk_PsiPrimePath);  //the muon is from the jpsi
      jpsiFromMuon_fromJpsiTrk_NonResonant_flags.push_back(
          muonFromJpsi_fromJpsiTrk_NonResonantPath);                             //the muon is from the jpsi
      jpsiFromMuon_fromDoubleMu_flags.push_back(muonFromJpsi_fromDoubleMuPath);  //the muon is from the jpsi

      dimuon0_jpsi_flags.push_back(dimuon0_jpsi_seed);
      dimuon43_jpsi_displaced_flags.push_back(dimuon43_jpsi_displaced_seed);

      dimuon0_jpsi_4R_flags.push_back(dimuon0_jpsi_4R_seed);
      muon_flags.push_back(muon_seed);
      iso_muon_flags.push_back(iso_muon_seed);
     
      triggeringMuons.push_back(obj);
      if (debug) {
        std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi()
                  << std::endl;
        // Print trigger object collection and type
        std::cout << "\t   Collection: " << obj.collection() << std::endl;
        std::cout << "\t muonFromJpsi_fromDimuon0Path " << muonFromJpsi_fromDimuon0Path << std::endl;
        std::cout << "\t muonFromJpsi_fromJpsiTrkPath " << muonFromJpsi_fromJpsiTrkPath << std::endl;
        std::cout << "\t muonFromJpsi_fromDoubleMuPath " << muonFromJpsi_fromDoubleMuPath << std::endl;
        std::cout << "\t pass_dimuon0_path" << pass_dimuon01_path << std::endl;
        std::cout << "\t dimuon0_seed" << dimuon0_seed << std::endl;
        std::cout << "\t jpsitrk_seed" << jpsitrk_seed << std::endl;
        std::cout << "\t jpsitrk_PsiPrime_seed" << jpsitrk_PsiPrime_seed << std::endl;
        std::cout << "\t jpsitrk_NonResonant_seed" << jpsitrk_NonResonant_seed << std::endl;
        std::cout << "\t doubleMu_seed" << doubleMu_seed << std::endl;
        std::cout << "\t dimuon0_jpsi_seed" << dimuon0_jpsi_seed << std::endl;
        std::cout << "\t dimuon0_jpsi_4R_seed" << dimuon0_jpsi_4R_seed << std::endl;
        std::cout << "\t dimuon43_jpsi_displaced_seed" << dimuon43_jpsi_displaced_seed << std::endl;
        std::cout << "\t muon_seed" << muon_seed << std::endl;
	std::cout << "\t iso_muon_seed" << iso_muon_seed << std::endl;
      }

    }  //trigger objects
  }

  if (debug) {
    std::cout << "\n total n of triggering muons = " << triggeringMuons.size() << std::endl;
    for (auto ij : triggeringMuons) {
      std::cout << " >>> components (pt, eta, phi) = " << ij.pt() << " " << ij.eta() << " " << ij.phi() << std::endl;
    }
  }
  //now check for reco muons matched to triggering muons
  std::vector<int> muonIsTrigger(muons->size(), 0);
  std::vector<int> muonIsFromJpsi_dimuon0Path(muons->size(), 0);
  std::vector<int> muonIsFromJpsi_jpsiTrkPath(muons->size(), 0);
  std::vector<int> muonIsFromJpsi_jpsiTrk_PsiPrimePath(muons->size(), 0);
  std::vector<int> muonIsFromJpsi_jpsiTrk_NonResonantPath(muons->size(), 0);
  std::vector<int> muonIsFromJpsi_doubleMuPath(muons->size(), 0);
  std::vector<int> muonIsDimuon0Trg(muons->size(), 0);
  std::vector<int> muonIsJpsiTrkTrg(muons->size(), 0);
  std::vector<int> muonIsJpsiTrk_PsiPrimeTrg(muons->size(), 0);
  std::vector<int> muonIsJpsiTrk_NonResonantTrg(muons->size(), 0);

  std::vector<int> muonIsDoubleMuTrg(muons->size(), 0);
  std::vector<int> muonIsMu8_Trg(muons->size(), 0);
  std::vector<int> muonIsIsoMu24_Trg(muons->size(), 0);

  std::vector<int> muonIs_dimuon0_jpsi_Trg(muons->size(), 0);
  std::vector<int> muonIs_dimuon43_jpsi_displaced_Trg(muons->size(), 0);

  std::vector<int> muonIs_dimuon0_jpsi_4R_Trg(muons->size(), 0);

  if (debugTrg)
    std::cout << "Event -------------------- " << std::endl;

  for (const pat::Muon& muon : *muons) {
    //this is for triggering muon not really need to be configurable
    unsigned int iMuo(&muon - &(muons->at(0)));
    //if(!(muon.isLooseMuon() && muon.isSoftMuon(PV))) continue;

    bool isMuonMatchedToDimuon0Path = !(muon.triggerObjectMatchByPath("HLT_Dimuon0_Jpsi3p5_Muon2_v*") == nullptr);

    bool isMuonMatchedToJpsiTrkPath =
      !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*") == nullptr);
    bool isMuonMatchedToJpsiTrk_PsiPrimePath =
      !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v*") == nullptr);
    bool isMuonMatchedToJpsiTrk_NonResonantPath =
      !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v*") == nullptr);

    bool isMuonMatchedToDoubleMuPath = !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_3_Jpsi_v*") == nullptr);

    bool isMuonMatchedToDoubleMu43Path = !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_3_Jpsi_v*") == nullptr);
    bool isMuonMatchedTo_dimuon0_jpsi_Path =
      !(muon.triggerObjectMatchByPath("HLT_Dimuon0_Jpsi_NoVertexing_v*") == nullptr);
    bool isMuonMatchedTo_dimuon43_jpsi_displaced_Path =
      !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_3_Jpsi_Displaced_v*") == nullptr);

    bool isMuonMatchedTo_dimuon0_jpsi_4R_Path =
      !(muon.triggerObjectMatchByPath("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v*") == nullptr);

    bool isMuonMatchTo_muon_Path = !(muon.triggerObjectMatchByPath("HLT_Mu8_v*") == nullptr);
    bool isMuonMatchTo_iso_muon_Path = !(muon.triggerObjectMatchByPath("HLT_IsoMu24_v*") == nullptr);
 
    /*    
    bool isMuonMatchedToDimuon0Path = !(muon.triggerObjectMatchByPath("HLT_Dimuon0_Jpsi3p5_Muon2_v5") == nullptr) ||
                                      !(muon.triggerObjectMatchByPath("HLT_Dimuon0_Jpsi3p5_Muon2_v6") == nullptr);
    bool isMuonMatchedToJpsiTrkPath =
        !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v14") == nullptr) ||
        !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v15") == nullptr);
    bool isMuonMatchedToJpsiTrk_PsiPrimePath =
        !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v14") == nullptr) ||
        !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v15") == nullptr);
    bool isMuonMatchedToJpsiTrk_NonResonantPath =
        !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v14") == nullptr) ||
        !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15") == nullptr);

    bool isMuonMatchedToDoubleMuPath = !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_3_Jpsi_v1") == nullptr) ||
                                       !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_3_Jpsi_v2") == nullptr);

    bool isMuonMatchedTo_dimuon0_jpsi_Path =
        !(muon.triggerObjectMatchByPath("HLT_Dimuon0_Jpsi_NoVertexing_v7") == nullptr) ||
        !(muon.triggerObjectMatchByPath("HLT_Dimuon0_Jpsi_NoVertexing_v8") == nullptr);
    bool isMuonMatchedTo_dimuon43_jpsi_displaced_Path =
        !(muon.triggerObjectMatchByPath("HLT_DoubleMu4_3_Jpsi_Displaced_v7") == nullptr);

    bool isMuonMatchedTo_dimuon0_jpsi_4R_Path =
        !(muon.triggerObjectMatchByPath("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v6") == nullptr) ||
        !(muon.triggerObjectMatchByPath("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v7") == nullptr);

    bool isMuonMatchTo_muon_Path = !(muon.triggerObjectMatchByPath("HLT_Mu8_v11") == nullptr) ||
                                   !(muon.triggerObjectMatchByPath("HLT_Mu8_v12") == nullptr);
    */
    //if( isMuonMatchedToJpsiTrk_NonResonantPath ) std::cout << " isMuonMatchedToJpsiTrk_NonResonantPath-------------------- " << std::endl;

    float dRMuonMatchingDimuon0 = -1.;
    int recoMuonMatchingDimuon0_index = -1;
    int trgMuonMatchingDimuon0_index = -1;

    float dRMuonMatchingJpsiTrk = -1.;
    int recoMuonMatchingJpsiTrk_index = -1;
    int trgMuonMatchingJpsiTrk_index = -1;

    float dRMuonMatchingJpsiTrk_PsiPrime = -1.;
    int recoMuonMatchingJpsiTrk_PsiPrime_index = -1;
    int trgMuonMatchingJpsiTrk_PsiPrime_index = -1;

    float dRMuonMatchingJpsiTrk_NonResonant = -1.;
    int recoMuonMatchingJpsiTrk_NonResonant_index = -1;
    int trgMuonMatchingJpsiTrk_NonResonant_index = -1;

    float dRMuonMatchingDoubleMu = -1.;
    int recoMuonMatchingDoubleMu_index = -1;
    int trgMuonMatchingDoubleMu_index = -1;

    float dRMuonMatchingDoubleMu43 = -1.;
    int recoMuonMatchingDoubleMu43_index = -1;
    int trgMuonMatchingDoubleMu43_index = -1;

    float dRMuonMatchingDimuon0_jpsi = -1.;
    int recoMuonMatchingDimuon0_jpsi_index = -1;
    int trgMuonMatchingDimuon0_jpsi_index = -1;

    float dRMuonMatchingDimuon43_jpsi_displaced = -1.;
    int recoMuonMatchingDimuon43_jpsi_displaced_index = -1;
    int trgMuonMatchingDimuon43_jpsi_displaced_index = -1;

    float dRMuonMatchingDimuon0_jpsi_4R = -1.;
    int recoMuonMatchingDimuon0_jpsi_4R_index = -1;
    int trgMuonMatchingDimuon0_jpsi_4R_index = -1;

    float dRMuonMatchingMuon = -1.;
    int recoMuonMatchingMuon_index = -1;
    int trgMuonMatchingMuon_index = -1;

    float dRMuonMatchingMuon_iso = -1.;
    int recoMuonMatchingMuon_iso_index = -1;
    int trgMuonMatchingMuon_iso_index = -1;


    for (unsigned int iTrg = 0; iTrg < triggeringMuons.size(); ++iTrg) {
      if (!dimuon0Flags[iTrg] && !jpsiTrkFlags[iTrg] && !jpsiTrk_PsiPrimeFlags[iTrg] &&
          !jpsiTrk_NonResonantFlags[iTrg] && !doubleMuFlags[iTrg] && !dimuon0_jpsi_flags[iTrg] &&
          !dimuon43_jpsi_displaced_flags[iTrg] && !muon_flags[iTrg] && !iso_muon_flags[iTrg] && !dimuon0_jpsi_4R_flags[iTrg])
        continue;

      //it passes the dimuon0 trigger

      float dR = reco::deltaR(triggeringMuons[iTrg], muon);

      if (isMuonMatchedToDimuon0Path && dimuon0Flags[iTrg] &&
          (dR < dRMuonMatchingDimuon0 || dRMuonMatchingDimuon0 == -1) && dR < maxdR_) {
        dRMuonMatchingDimuon0 = dR;
        recoMuonMatchingDimuon0_index = iMuo;
        trgMuonMatchingDimuon0_index = iTrg;
        if (debug)
          std::cout << "trigger dimuon0 " << dimuon0Flags[iTrg] << std::endl;
        if (debug)
          std::cout << " dR = " << dR << std::endl;
        if (debug)
          std::cout << " HLT = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " "
                    << triggeringMuons[iTrg].phi() << std::endl;
        if (debug)
          std::cout << " reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << std::endl;
      }
      if (isMuonMatchTo_muon_Path && muon_flags[iTrg] && (dR < dRMuonMatchingMuon || dRMuonMatchingMuon == -1) &&
          dR < maxdR_) {
        dRMuonMatchingMuon = dR;
        recoMuonMatchingMuon_index = iMuo;
        trgMuonMatchingMuon_index = iTrg;
        if (debug){
          std::cout << "trigger Muon 8 " << muon_flags[iTrg] << std::endl;
	  if (isMuonMatchTo_muon_Path !=  muon_flags[iTrg]) std::cout << "!!!!!!! Mismatch for Mu8" << std::endl;
	}
        if (debug)
          std::cout << " dR = " << dR << std::endl;
        if (debug)
          std::cout << " HLT = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " "
                    << triggeringMuons[iTrg].phi() << std::endl;
        if (debug)
          std::cout << " reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << std::endl;
      }
      if (isMuonMatchTo_iso_muon_Path && iso_muon_flags[iTrg] && (dR < dRMuonMatchingMuon_iso || dRMuonMatchingMuon_iso == -1) &&
          dR < maxdR_) {
        dRMuonMatchingMuon_iso = dR;
        recoMuonMatchingMuon_iso_index = iMuo;
        trgMuonMatchingMuon_iso_index = iTrg;
        if (debug)
	  std::cout << "trigger Muon Iso" << iso_muon_flags[iTrg] << std::endl;
	if (isMuonMatchTo_iso_muon_Path !=  iso_muon_flags[iTrg]) std::cout << "!!!!!!! Mismatch for MuIso" << std::endl;

        if (debug)
	  std::cout << " dR = " << dR << std::endl;
        if (debug)
	  std::cout << " HLT = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " "
                    << triggeringMuons[iTrg].phi() << std::endl;
        if (debug)
	  std::cout << " reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << std::endl;
      }


      if (isMuonMatchedTo_dimuon0_jpsi_Path && dimuon0_jpsi_flags[iTrg] &&
          (dR < dRMuonMatchingDimuon0_jpsi || dRMuonMatchingDimuon0_jpsi == -1) && dR < maxdR_) {
        dRMuonMatchingDimuon0_jpsi = dR;
        recoMuonMatchingDimuon0_jpsi_index = iMuo;
        trgMuonMatchingDimuon0_jpsi_index = iTrg;
        if (debug)
          std::cout << "trigger dimuon0 jpsi " << dimuon0_jpsi_flags[iTrg] << std::endl;
        if (debug)
          std::cout << " dR = " << dR << std::endl;
        if (debug)
          std::cout << " HLT = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " "
                    << triggeringMuons[iTrg].phi() << std::endl;
        if (debug)
          std::cout << " reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << std::endl;
      }
      if (isMuonMatchedTo_dimuon43_jpsi_displaced_Path && dimuon43_jpsi_displaced_flags[iTrg] &&
          (dR < dRMuonMatchingDimuon43_jpsi_displaced || dRMuonMatchingDimuon43_jpsi_displaced == -1) && dR < maxdR_) {
        dRMuonMatchingDimuon43_jpsi_displaced = dR;
        recoMuonMatchingDimuon43_jpsi_displaced_index = iMuo;
        trgMuonMatchingDimuon43_jpsi_displaced_index = iTrg;
        if (debug){
          std::cout << "trigger dimuon43 jpsi displaced " << dimuon43_jpsi_displaced_flags[iTrg] << std::endl;
	  if (isMuonMatchedTo_dimuon43_jpsi_displaced_Path !=  dimuon43_jpsi_displaced_flags[iTrg]) std::cout << "!!!!!!! Mismatch for dimuon43_jpsi_displaced" << std::endl;
	}
        if (debug)
          std::cout << " dR = " << dR << std::endl;
        if (debug)
          std::cout << " HLT = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " "
                    << triggeringMuons[iTrg].phi() << std::endl;
        if (debug)
          std::cout << " reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << std::endl;
      }
      if (isMuonMatchedTo_dimuon0_jpsi_4R_Path && dimuon0_jpsi_4R_flags[iTrg] &&
          (dR < dRMuonMatchingDimuon0_jpsi_4R || dRMuonMatchingDimuon0_jpsi_4R == -1) && dR < maxdR_) {
        dRMuonMatchingDimuon0_jpsi_4R = dR;
        recoMuonMatchingDimuon0_jpsi_4R_index = iMuo;
        trgMuonMatchingDimuon0_jpsi_4R_index = iTrg;
        if (debug)
          std::cout << "trigger dimuon0 jpsi 4R " << dimuon0_jpsi_4R_flags[iTrg] << std::endl;
        if (debug)
          std::cout << " dR = " << dR << std::endl;
        if (debug)
          std::cout << " HLT = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " "
                    << triggeringMuons[iTrg].phi() << std::endl;
        if (debug)
          std::cout << " reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << std::endl;
      }

      if (isMuonMatchedToJpsiTrkPath && jpsiTrkFlags[iTrg] &&
          (dR < dRMuonMatchingJpsiTrk || dRMuonMatchingJpsiTrk == -1) && dR < maxdR_) {
        dRMuonMatchingJpsiTrk = dR;
        recoMuonMatchingJpsiTrk_index = iMuo;
        trgMuonMatchingJpsiTrk_index = iTrg;
        if (debug)
          std::cout << "trigger jpsitrk " << jpsiTrkFlags[iTrg] << std::endl;
        if (debug)
          std::cout << " dR = " << dR << std::endl;
        if (debug)
          std::cout << " HLT = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " "
                    << triggeringMuons[iTrg].phi() << std::endl;
        if (debug)
          std::cout << " reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << std::endl;
      }
      if (isMuonMatchedToJpsiTrk_PsiPrimePath && jpsiTrk_PsiPrimeFlags[iTrg] &&
          (dR < dRMuonMatchingJpsiTrk_PsiPrime || dRMuonMatchingJpsiTrk_PsiPrime == -1) && dR < maxdR_) {
        dRMuonMatchingJpsiTrk_PsiPrime = dR;
        recoMuonMatchingJpsiTrk_PsiPrime_index = iMuo;
        trgMuonMatchingJpsiTrk_PsiPrime_index = iTrg;
        if (debug)
          std::cout << "trigger jpsitrk " << jpsiTrk_PsiPrimeFlags[iTrg] << std::endl;
        if (debug)
          std::cout << " dR = " << dR << std::endl;
        if (debug)
          std::cout << " HLT = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " "
                    << triggeringMuons[iTrg].phi() << std::endl;
        if (debug)
          std::cout << " reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << std::endl;
      }
      if (isMuonMatchedToJpsiTrk_NonResonantPath && jpsiTrk_NonResonantFlags[iTrg] &&
          (dR < dRMuonMatchingJpsiTrk_NonResonant || dRMuonMatchingJpsiTrk_NonResonant == -1) && dR < maxdR_) {
        dRMuonMatchingJpsiTrk_NonResonant = dR;
        recoMuonMatchingJpsiTrk_NonResonant_index = iMuo;
        trgMuonMatchingJpsiTrk_NonResonant_index = iTrg;
        if (debug)
          std::cout << "trigger jpsitrk " << jpsiTrk_NonResonantFlags[iTrg] << std::endl;
        if (debug)
          std::cout << " dR = " << dR << std::endl;
        if (debug)
          std::cout << " HLT = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " "
                    << triggeringMuons[iTrg].phi() << std::endl;
        if (debug)
          std::cout << " reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << std::endl;
      }
      if (isMuonMatchedToDoubleMuPath &&   doubleMuFlags[iTrg] &&
          (dR < dRMuonMatchingDoubleMu || dRMuonMatchingDoubleMu == -1) && dR < maxdR_) {
        dRMuonMatchingDoubleMu = dR;
        recoMuonMatchingDoubleMu_index = iMuo;
        trgMuonMatchingDoubleMu_index = iTrg;
        if (debug)
          std::cout << "trigger doubleMu " << doubleMuFlags[iTrg] << std::endl;
        if (debug)
          std::cout << " dR = " << dR << std::endl;
        if (debug)
          std::cout << " HLT = " << triggeringMuons[iTrg].pt() << " " << triggeringMuons[iTrg].eta() << " "
                    << triggeringMuons[iTrg].phi() << std::endl;
        if (debug)
          std::cout << " reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << std::endl;
      }
      if (isMuonMatchedToDoubleMu43Path &&  jpsiFromMuon_fromDoubleMu_flags[iTrg] &&
	  (dR < dRMuonMatchingDoubleMu43 || dRMuonMatchingDoubleMu43 == -1) && dR < maxdR_) {
	dRMuonMatchingDoubleMu43 = dR;
	recoMuonMatchingDoubleMu43_index = iMuo;
	trgMuonMatchingDoubleMu43_index = iTrg;
	if (debug){
	  if (isMuonMatchedToDoubleMu43Path!=  jpsiFromMuon_fromDoubleMu_flags[iTrg] )
	    std::cout << "!!!!! Mismatch HLT_DoubleMu4_3_Jpsi "<< std::endl;  
	    }
      }
    }

    //save reco muon
    if (recoMuonMatchingDimuon0_index != -1 || recoMuonMatchingJpsiTrk_index != -1 ||
        recoMuonMatchingJpsiTrk_PsiPrime_index != -1 || recoMuonMatchingJpsiTrk_NonResonant_index != -1 ||
        recoMuonMatchingDoubleMu_index != -1 ||recoMuonMatchingDoubleMu43_index != -1 || recoMuonMatchingDimuon0_jpsi_index != -1 ||
        recoMuonMatchingDimuon43_jpsi_displaced_index != -1 || recoMuonMatchingDimuon0_jpsi_4R_index != -1 ||
        recoMuonMatchingMuon_index != -1 || recoMuonMatchingMuon_iso_index!=-1)

    {
      muonIsTrigger[iMuo] = 1;
      pat::Muon recoTriggerMuonCand(muon);
      recoTriggerMuonCand.addUserInt("trgMuonMuon_index", trgMuonMatchingMuon_index);
      recoTriggerMuonCand.addUserInt("trgMuonDimuon0_index", trgMuonMatchingDimuon0_index);
      recoTriggerMuonCand.addUserInt("trgMuonJpsiTrk_index", trgMuonMatchingJpsiTrk_index);
      recoTriggerMuonCand.addUserInt("trgMuonJpsiTrk_PsiPrime_index", trgMuonMatchingJpsiTrk_PsiPrime_index);
      recoTriggerMuonCand.addUserInt("trgMuonJpsiTrk_NonResonant_index", trgMuonMatchingJpsiTrk_NonResonant_index);
      recoTriggerMuonCand.addUserInt("trgMuonDoubleMu_index", trgMuonMatchingDoubleMu_index);
      recoTriggerMuonCand.addUserInt("trgMuonDimuon0_jpsi_index", trgMuonMatchingDimuon0_jpsi_index);
      recoTriggerMuonCand.addUserInt("trgMuonDimuon43_jpsi_displaced_index",
                                     trgMuonMatchingDimuon43_jpsi_displaced_index);
      recoTriggerMuonCand.addUserInt("trgMuonMu8_index", trgMuonMatchingMuon_index);
      recoTriggerMuonCand.addUserInt("trgMuonIsoMu24_index", trgMuonMatchingMuon_iso_index);
      recoTriggerMuonCand.addUserInt("trgMuonDimuon0_jpsi_4R_index", trgMuonMatchingDimuon0_jpsi_4R_index);
      trgmuons_out->emplace_back(recoTriggerMuonCand);
      //keep track of original muon index for SelectedMuons collection

      if (recoMuonMatchingMuon_index != -1) {
        muonIsMu8_Trg[iMuo] = muon_flags[trgMuonMatchingMuon_index];
      } else {
        muonIsMu8_Trg[iMuo] = 0;
      }

      if (recoMuonMatchingMuon_iso_index != -1) {
        muonIsIsoMu24_Trg[iMuo] = iso_muon_flags[trgMuonMatchingMuon_iso_index];
      } else {
        muonIsIsoMu24_Trg[iMuo] = 0;
      }

      if (recoMuonMatchingDimuon0_index != -1) {
        muonIsFromJpsi_dimuon0Path[iMuo] = jpsiFromMuon_fromDimuon0_flags[trgMuonMatchingDimuon0_index];
        muonIsDimuon0Trg[iMuo] = dimuon0Flags[trgMuonMatchingDimuon0_index];
      } else {
        muonIsFromJpsi_dimuon0Path[iMuo] = 0;  //jpsiFromMuon_fromDimuon0_flags[trgMuonMatchingDimuon0_index];
        muonIsDimuon0Trg[iMuo] = 0;            //dimuon0Flags[trgMuonMatchingDimuon0_index];
      }

      if (recoMuonMatchingDimuon0_jpsi_index != -1) {
        muonIs_dimuon0_jpsi_Trg[iMuo] = dimuon0_jpsi_flags[trgMuonMatchingDimuon0_jpsi_index];
      } else {
        muonIs_dimuon0_jpsi_Trg[iMuo] = 0;
      }

      if (recoMuonMatchingDimuon43_jpsi_displaced_index != -1) {
        muonIs_dimuon43_jpsi_displaced_Trg[iMuo] =
            dimuon43_jpsi_displaced_flags[trgMuonMatchingDimuon43_jpsi_displaced_index];
      } else {
        muonIs_dimuon43_jpsi_displaced_Trg[iMuo] = 0;
      }

      if (recoMuonMatchingDimuon0_jpsi_4R_index != -1) {
        muonIs_dimuon0_jpsi_4R_Trg[iMuo] = dimuon0_jpsi_4R_flags[trgMuonMatchingDimuon0_jpsi_4R_index];
      } else {
        muonIs_dimuon0_jpsi_4R_Trg[iMuo] = 0;
      }

      if (recoMuonMatchingJpsiTrk_index != -1) {
        if (debug)
          std::cout << "here it fills the vector for " << muon.pt() << " flag "
                    << jpsiFromMuon_fromJpsiTrk_flags[trgMuonMatchingJpsiTrk_index] << std::endl;
        muonIsFromJpsi_jpsiTrkPath[iMuo] = jpsiFromMuon_fromJpsiTrk_flags[trgMuonMatchingJpsiTrk_index];
        muonIsJpsiTrkTrg[iMuo] = jpsiTrkFlags[trgMuonMatchingJpsiTrk_index];
      } else {
        muonIsFromJpsi_jpsiTrkPath[iMuo] = 0;  //jpsiFromMuon_fromJpsiTrk_flags[trgMuonMatchingJpsiTrk_index];
        muonIsJpsiTrkTrg[iMuo] = 0;            //jpsiTrkFlags[trgMuonMatchingJpsiTrk_index];
      }
      if (recoMuonMatchingJpsiTrk_PsiPrime_index != -1) {
        muonIsFromJpsi_jpsiTrk_PsiPrimePath[iMuo] =
            jpsiFromMuon_fromJpsiTrk_PsiPrime_flags[trgMuonMatchingJpsiTrk_PsiPrime_index];
        muonIsJpsiTrk_PsiPrimeTrg[iMuo] = jpsiTrk_PsiPrimeFlags[trgMuonMatchingJpsiTrk_PsiPrime_index];
      } else {
        muonIsFromJpsi_jpsiTrk_PsiPrimePath[iMuo] = 0;  //jpsiFromMuon_fromJpsiTrk_flags[trgMuonMatchingJpsiTrk_index];
        muonIsJpsiTrk_PsiPrimeTrg[iMuo] = 0;            //jpsiTrkFlags[trgMuonMatchingJpsiTrk_index];
      }

      if (debugTrg)
        if (isMuonMatchedToJpsiTrk_NonResonantPath)
          std::cout << "recoMuonMatchingJpsiTrk_NonResonant_flag=====================  "
                    << jpsiFromMuon_fromJpsiTrk_NonResonant_flags[recoMuonMatchingJpsiTrk_NonResonant_index]
                    << std::endl;
      if (recoMuonMatchingJpsiTrk_NonResonant_index != -1) {
        muonIsFromJpsi_jpsiTrk_NonResonantPath[iMuo] =
            jpsiFromMuon_fromJpsiTrk_NonResonant_flags[trgMuonMatchingJpsiTrk_NonResonant_index];
        muonIsJpsiTrk_NonResonantTrg[iMuo] = jpsiTrk_NonResonantFlags[trgMuonMatchingJpsiTrk_NonResonant_index];
      } else {
        muonIsFromJpsi_jpsiTrk_NonResonantPath[iMuo] =
            0;                                   //jpsiFromMuon_fromJpsiTrk_flags[trgMuonMatchingJpsiTrk_index];
        muonIsJpsiTrk_NonResonantTrg[iMuo] = 0;  //jpsiTrkFlags[trgMuonMatchingJpsiTrk_index];
      }

      if (recoMuonMatchingDoubleMu43_index != -1) {
        muonIsFromJpsi_doubleMuPath[iMuo] = jpsiFromMuon_fromDoubleMu_flags[trgMuonMatchingDoubleMu43_index];
      } else {
        muonIsFromJpsi_doubleMuPath[iMuo] = 0;  //jpsiFromMuon_fromDoubleMu_flags[trgMuonMatchingDoubleMu_index];                                                                                                                                                 
      }
      

      if (recoMuonMatchingDoubleMu_index != -1) {
        muonIsDoubleMuTrg[iMuo] = doubleMuFlags[trgMuonMatchingDoubleMu_index];
      } else {
        muonIsDoubleMuTrg[iMuo] = 0;            //doubleMuFlags[trgMuonMatchingDoubleMu_index];
      }
      /*
      if(debug){      
        std::cout << "HERE" << std::endl;
        std::cout << "trgMuonMatching_index: " << trgMuonMatching_index << std::endl;
        std::cout << "jpsiFromMuon_fromDimuon0_flags.push_back(muonFromJpsi_fromDimuon0Path)"<< jpsiFromMuon_fromDimuon0_flags[trgMuonMatchingDimuon0_index] << std::endl;
        std::cout << "jpsiFromMuon_fromJpsiTrk_flags.push_back(muonFromJpsi_fromDimuon0Path)"<< jpsiFromMuon_fromJpsiTrk_flags[trgMuonMatching_index] << std::endl;
        std::cout << "dimuon0Flags.push_back(pass_dimuon0_path)"<< dimuon0Flags[trgMuonMatchingDimuon0_index] << std::endl;
        std::cout << "jpsiTrkFlags.push_back(pass_jpsiTrk_path)"<< jpsiTrkFlags[trgMuonMatching_index] << std::endl;
           std::cout  << "----- reco = " << muon.pt() << " " << muon.eta() << " " << muon.phi() << " " 
		      << " HLT Dimuon0 = " << triggeringMuons[trgMuonMatchingDimuon0_index].pt() << " " << triggeringMuons[trgMuonMatchingDimuon0_index].eta() << " " << triggeringMuons[trgMuonMatchingDimuon0_index].phi()
		      << " HLT JpsiTrk = " << triggeringMuons[trgMuonMatchingJpsiTrk_index].pt() << " " << triggeringMuons[trgMuonMatchingJpsiTrk_index].eta() << " " << triggeringMuons[trgMuonMatchingJpsiTrk_index].phi()
		      << std::endl;
      }
      */
    }
  }

  // now produce output for analysis (code simplified loop of trg inside)
  // trigger muon + all compatible in dz with any tag
  for (unsigned int muIdx = 0; muIdx < muons->size(); ++muIdx) {
    const pat::Muon& mu = (*muons)[muIdx];
    //selection cuts
    if (mu.pt() < ptMin_)
      continue;
    if (fabs(mu.eta()) > absEtaMax_)
      continue;
    //following ID is needed for trigger muons not here
    // anyway it is off in the configuration
    //G: if (softMuonsOnly_ && !mu.isSoftMuon(PV)) continue;

    /* // same PV as the tag muon, both tag and probe only dz selection
    bool SkipMuon=true;
    for (const pat::Muon & trgmu : *trgmuons_out) {
	    if( fabs(mu.vz()-trgmu.vz()) > dzTrg_cleaning_ && dzTrg_cleaning_ >0 )
	      continue;
	    SkipMuon=false;
    } 
    // needs decission: what about events without trg muon? now we SKIP them
    if (SkipMuon)  continue;
    */

    // build transient track
    const reco::TransientTrack muonTT((*(mu.bestTrack())),
                                      &(*bFieldHandle));  //sara: check, why not using inner track for muons?
    if (!muonTT.isValid())
      continue;

    muons_out->emplace_back(mu);
    muons_out->back().addUserInt("isTriggering", muonIsTrigger[muIdx]);
    muons_out->back().addUserInt("isMuonFromJpsi_dimuon0Trg", muonIsFromJpsi_dimuon0Path[muIdx]);
    muons_out->back().addUserInt("isMuonFromJpsi_jpsiTrkTrg", muonIsFromJpsi_jpsiTrkPath[muIdx]);
    muons_out->back().addUserInt("isMuonFromJpsi_jpsiTrk_PsiPrimeTrg", muonIsFromJpsi_jpsiTrk_PsiPrimePath[muIdx]);
    muons_out->back().addUserInt("isMuonFromJpsi_jpsiTrk_NonResonantTrg",
                                 muonIsFromJpsi_jpsiTrk_NonResonantPath[muIdx]);
    muons_out->back().addUserInt("isMuonFromJpsi_doubleMuTrg", muonIsFromJpsi_doubleMuPath[muIdx]);

    muons_out->back().addUserInt("isDimuon0Trg", muonIsDimuon0Trg[muIdx]);
    muons_out->back().addUserInt("isJpsiTrkTrg", muonIsJpsiTrkTrg[muIdx]);
    muons_out->back().addUserInt("isJpsiTrk_PsiPrimeTrg", muonIsJpsiTrk_PsiPrimeTrg[muIdx]);
    muons_out->back().addUserInt("isJpsiTrk_NonResonantTrg", muonIsJpsiTrk_NonResonantTrg[muIdx]);
    muons_out->back().addUserInt("isDoubleMuTrg", muonIsDoubleMuTrg[muIdx]);
    muons_out->back().addUserInt("isMu8_Trg", muonIsMu8_Trg[muIdx]);
    muons_out->back().addUserInt("isIsoMu24_Trg", muonIsIsoMu24_Trg[muIdx]);    
    muons_out->back().addUserInt("isDimuon0_jpsi_Trg", muonIs_dimuon0_jpsi_Trg[muIdx]);
    muons_out->back().addUserInt("isDimuon43_jpsi_displaced_Trg", muonIs_dimuon43_jpsi_displaced_Trg[muIdx]);
    muons_out->back().addUserInt("isDimuon0_jpsi_4R_Trg", muonIs_dimuon0_jpsi_4R_Trg[muIdx]);

    // std::cout << "removing zeros"<<std::endl;
    // for (auto&  [key_hlt , value_ps] :hlt_psHLT_map){
    //   std::cout<< "key_hlt "<< key_hlt;
    //   for (auto&  ps :value_ps){
    // 	std::cout<< " ps "<< ps;
    //   }
    //   std::cout<< std::endl;
    // }
    // for (auto&  [key_hlt , value_ps] :hlt_psHLT_map){
    //   value_ps.erase( remove_if( value_ps.begin(), value_ps.end(), []( int i ){ return i < 1; } ), value_ps.end() );
    // }
    // std::cout << "zeros removed "<<std::endl;
    // for (auto&  [key_hlt , value_ps] :hlt_psHLT_map){
    //   std::cout<< "key_hlt "<< key_hlt;
    //   for (auto&  ps :value_ps){
    // 	std::cout<< "ps "<< ps;
    //   }
    //   std::cout<< " Min ps "<<  *min_element(hlt_psHLT_map[key_hlt].begin(), hlt_psHLT_map[key_hlt].end());
    //   std::cout<< std::endl;

    // }

    // for (auto&  [key_hlt , value_ps] :hlt_psL1_map){
    //   value_ps.erase( remove_if( value_ps.begin(), value_ps.end(), []( int i ){ return i < 1; } ), value_ps.end() );
    // }

    //std::cout << "Saving PS"<<std::endl;
    muons_out->back().addUserInt("isDimuon0Trg_HLT",!(mu.triggerObjectMatchByPath("HLT_Dimuon0_Jpsi3p5_Muon2_v*") == nullptr));
				 //                                 *min_element(hlt_psHLT_map["HLT_Dimuon0_Jpsi3p5_Muon2_v"].begin(),
                                 //             hlt_psHLT_map["HLT_Dimuon0_Jpsi3p5_Muon2_v"].end()));
    muons_out->back().addUserInt("isJpsiTrkTrg_HLT",!(mu.triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*") == nullptr));
				 //                         *min_element(hlt_psHLT_map["HLT_DoubleMu4_JpsiTrk_Displaced_v"].begin(),
                                 //             hlt_psHLT_map["HLT_DoubleMu4_JpsiTrk_Displaced_v"].end()));
    muons_out->back().addUserInt("isJpsiTrk_PsiPrimeTrg_HLT",!(mu.triggerObjectMatchByPath("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v*") == nullptr));
				 //    *min_element(hlt_psHLT_map["HLT_DoubleMu4_PsiPrimeTrk_Displaced_v"].begin(),
                                 //             hlt_psHLT_map["HLT_DoubleMu4_PsiPrimeTrk_Displaced_v"].end()));
    muons_out->back().addUserInt("isJpsiTrk_NonResonantTrg_HLT", !(mu.triggerObjectMatchByPath("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v*") == nullptr));
                                 //*min_element(hlt_psHLT_map["HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v"].begin(),
                                 //             hlt_psHLT_map["HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v"].end()));
    muons_out->back().addUserInt( "isDoubleMuTrg_HLT", !(mu.triggerObjectMatchByPath("HLT_DoubleMu4_3_Jpsi_v*") == nullptr));
				 //        *min_element(hlt_psHLT_map["HLT_DoubleMu4_3_Jpsi_v"].begin(), hlt_psHLT_map["HLT_DoubleMu4_3_Jpsi_v"].end()));
    muons_out->back().addUserInt("isDimuon0_jpsi_Trg_HLT", !(mu.triggerObjectMatchByPath("HLT_Dimuon0_Jpsi_NoVertexing_v*") == nullptr));
				 //                               *min_element(hlt_psHLT_map["HLT_Dimuon0_Jpsi_NoVertexing_v"].begin(),
                                 //             hlt_psHLT_map["HLT_Dimuon0_Jpsi_NoVertexing_v"].end()));
    muons_out->back().addUserInt("isDimuon43_jpsi_displaced_Trg_HLT",!(mu.triggerObjectMatchByPath("HLT_DoubleMu4_3_Jpsi_Displaced_v*") == nullptr));
				 //   *min_element(hlt_psHLT_map["HLT_DoubleMu4_3_Jpsi_Displaced_v"].begin(),
                                 //             hlt_psHLT_map["HLT_DoubleMu4_3_Jpsi_Displaced_v"].end()));
    muons_out->back().addUserInt("isDimuon0_jpsi_4R_Trg_HLT",!(mu.triggerObjectMatchByPath("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v*") == nullptr));
				 //                                 *min_element(hlt_psHLT_map["HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v"].begin(),
                                 //             hlt_psHLT_map["HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v"].end()));
    muons_out->back().addUserInt("isMu8_Trg_HLT", !(mu.triggerObjectMatchByPath("HLT_Mu8_v*")== nullptr));

    muons_out->back().addUserInt("isIsoMu24_Trg_HLT", !(mu.triggerObjectMatchByPath("HLT_IsoMu24_v*")== nullptr));
                                
    // muons_out->back().addUserInt("isDimuon0Trg_L1ps",
    //                              *min_element(hlt_psL1_map["HLT_Dimuon0_Jpsi3p5_Muon2_v"].begin(),
    //                                           hlt_psL1_map["HLT_Dimuon0_Jpsi3p5_Muon2_v"].end()));
    // muons_out->back().addUserInt("isJpsiTrkTrg_L1ps",
    //                              *min_element(hlt_psL1_map["HLT_DoubleMu4_JpsiTrk_Displaced_v"].begin(),
    //                                           hlt_psL1_map["HLT_DoubleMu4_JpsiTrk_Displaced_v"].end()));
    // muons_out->back().addUserInt("isJpsiTrk_PsiPrimeTrg_L1ps",
    //                              *min_element(hlt_psL1_map["HLT_DoubleMu4_PsiPrimeTrk_Displaced_v"].begin(),
    //                                           hlt_psL1_map["HLT_DoubleMu4_PsiPrimeTrk_Displaced_v"].end()));
    // muons_out->back().addUserInt("isJpsiTrk_NonResonantTrg_L1ps",
    //                              *min_element(hlt_psL1_map["HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v"].begin(),
    //                                           hlt_psL1_map["HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v"].end()));
    // muons_out->back().addUserInt(
    //     "isDoubleMuTrg_L1ps",
    //     *min_element(hlt_psL1_map["HLT_DoubleMu4_3_Jpsi_v"].begin(), hlt_psL1_map["HLT_DoubleMu4_3_Jpsi_v"].end()));
    // muons_out->back().addUserInt("isDimuon0_jpsi_Trg_L1ps",
    //                              *min_element(hlt_psL1_map["HLT_Dimuon0_Jpsi_NoVertexing_v"].begin(),
    //                                           hlt_psL1_map["HLT_Dimuon0_Jpsi_NoVertexing_v"].end()));
    // muons_out->back().addUserInt("isDimuon43_jpsi_displaced_Trg_L1ps",
    //                              *min_element(hlt_psL1_map["HLT_DoubleMu4_3_Jpsi_Displaced_v"].begin(),
    //                                           hlt_psL1_map["HLT_DoubleMu4_3_Jpsi_Displaced_v"].end()));
    // muons_out->back().addUserInt("isDimuon0_jpsi_4R_Trg_L1ps",
    //                              *min_element(hlt_psL1_map["HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v"].begin(),
    //                                           hlt_psL1_map["HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v"].end()));
    // muons_out->back().addUserInt("isMu8_Trg_L1ps",
    //                              *min_element(hlt_psL1_map["HLT_Mu8_v"].begin(), hlt_psL1_map["HLT_Mu8_v"].end()));
    // //muons_out->back().addUserInt("isIsoMu24_Trg_L1ps",
    // //                             *min_element(hlt_psL1_map["HLT_IsoMu24_v"].begin(), hlt_psL1_map["HLT_IsoMu24_v"].end()));
    // muons_out->back().addUserInt("isIsoMu24_Trg_L1ps",!(mu.triggerObjectMatchByPath("HLT_IsoMu24_*")== nullptr)); 
    trans_muons_out->emplace_back(muonTT);
				 }

  iEvent.put(std::move(trgmuons_out), "trgMuons");
  iEvent.put(std::move(muons_out), "SelectedMuons");
  iEvent.put(std::move(trans_muons_out), "SelectedTransientMuons");
}

DEFINE_FWK_MODULE(MuonTriggerSelector);
