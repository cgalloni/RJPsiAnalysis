//// table to produce hlt and l1 bits that we are going to use in the analysis, to avoid saving every bit in the menu


// system include files
#include <memory>


#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/View.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TString.h"
#include <string>


#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"

#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

class TrgBitTableProducer : public edm::stream::EDProducer<> {


public:

 
 explicit TrgBitTableProducer(const edm::ParameterSet &cfg):
    hltresultsToken_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag> ("hltresults"))),
    l1resultsToken_(consumes<GlobalAlgBlkBxCollection>(cfg.getParameter<edm::InputTag> ("l1results"))),
    hltpaths_( cfg.getParameter< std::vector<std::string> >( "paths" ) ),
    l1seeds_( cfg.getParameter< std::vector<std::string> >( "seeds" ) ),
    hltPrescaleProvider_(cfg, consumesCollector(), *this) 
{
    produces<nanoaod::FlatTable>();
    fGtUtil = new l1t::L1TGlobalUtil(cfg,
				     consumesCollector(),
				     *this,
				     cfg.getParameter<edm::InputTag>("l1results"),
				     cfg.getParameter<edm::InputTag>("l1results"));
}

  ~TrgBitTableProducer() override {}

  void produce(edm::Event&, edm::EventSetup const&) override;


  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  

private:

  const edm::EDGetTokenT< edm::TriggerResults >    hltresultsToken_;
  const edm::EDGetTokenT<GlobalAlgBlkBxCollection> l1resultsToken_;
  // l1 seeds not implemented yet but can be added with litle effort
  const std::vector< std::string >                 hltpaths_;
  const std::vector< std::string >                 l1seeds_;
  TString * algoBitToName  =  new TString[512]; 
  bool loaded = false;
  l1t::L1TGlobalUtil* fGtUtil;
  HLTPrescaleProvider hltPrescaleProvider_;
  void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override;
};

void
TrgBitTableProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  bool changed(true); 
  if (!hltPrescaleProvider_.init(iRun, iSetup, "HLT", changed)) {
    edm::LogError("TriggerCandProducer") << "Error! Can't initialize HLTConfigProvider";
    throw cms::Exception("HLTConfigProvider::init() returned non 0");
  }
}
void 
TrgBitTableProducer::produce( edm::Event &evt, edm::EventSetup const &stp) 
{

  //input
  edm::Handle<GlobalAlgBlkBxCollection> l1Results;
  evt.getByToken(l1resultsToken_,l1Results);
  edm::Handle< edm::TriggerResults > hltResults;
  evt.getByToken( hltresultsToken_, hltResults);


  // returns uint8 instead of bool, because addCollumnValue<bool> is unsuported. (cmsRun error). the next "economical" class is uint8
  std::vector<uint8_t> hltbits;
  std::vector<uint8_t> l1bits, l1bits_l1;
  std::vector<float>  trigPrescl_hlt_vec,  //trigPrescl_l1_vec, 
    trigPrescl_hlt_l1_vec,trigPrescl_l1_vec_l1;
  unsigned int Npaths = hltpaths_.size();
  hltbits.reserve( Npaths );
  trigPrescl_hlt_vec.reserve( Npaths );
  trigPrescl_hlt_l1_vec.reserve( Npaths );
  unsigned int Nseeds = l1seeds_.size();
  l1bits.reserve( Nseeds );  l1bits_l1.reserve( Nseeds );
  //  trigPrescl_l1_vec.reserve( Nseeds );
  trigPrescl_l1_vec_l1.reserve( Nseeds );
  edm::TriggerNames trigName = evt.triggerNames( *hltResults );   
 
  // get L1 seeds
  if ( !l1Results.isValid() ) {

   for (unsigned int iseed=0; iseed<l1seeds_.size(); iseed++)
     l1bits.push_back( 0 );

  } else{
    if (!loaded){
      edm::ESHandle<L1TUtmTriggerMenu> menu;
      stp.get<L1TUtmTriggerMenuRcd>().get(menu);
      for (auto const & keyval: menu->getAlgorithmMap()) {
        std::string const & trigName = keyval.second.getName();
        algoBitToName[ int( keyval.second.getIndex() ) ] = TString( trigName );
      }
      loaded=true;
    }
    
    GlobalAlgBlk const &result=l1Results->at(0, 0);
    // MEthod 1 for prescales 
    L1GtUtils const& l1GtUtils = hltPrescaleProvider_.l1GtUtils();
    int iErrorCode = -1;
    L1GtUtils::TriggerCategory trigCategory = L1GtUtils::AlgorithmTrigger;
    const int pfSetIndexAlgorithmTrigger = l1GtUtils.prescaleFactorSetIndex(evt, trigCategory, iErrorCode);
    if (iErrorCode == 0) {
      // if (debug)
      //std::cout << "%Prescale set index: " << pfSetIndexAlgorithmTrigger << std::endl;
	//} else {
      //std::cout << "%Could not extract Prescale set index from event record. Error code: " << iErrorCode << std::endl;
    }  

    //Method2 for prescales:
    fGtUtil->retrieveL1(evt, stp, l1resultsToken_);
    const std::vector<std::pair<std::string, bool>> finalDecisions = fGtUtil->decisionsFinal();
    const std::vector<std::pair<std::string, int>> prescales = fGtUtil->prescales();
    //    bool l1_bit=false;
    //float l1_Prescl_L1 = 0;
    for ( auto& l1seed: l1seeds_ ){
      bool l1_bit=false;
      float l1_Prescl_L1 = 0;
      for (unsigned int i = 0; i < finalDecisions.size(); ++i) {
	std::string name = (finalDecisions.at(i)).first;

	if (name == "NULL")
	  continue;
	std::cout<<"L1 name "<<name<< std::endl;
	if (name.compare(l1seed) == 0) {
	  std::cout<<"---L1 name found in seed "<< l1seed<< std::endl; 
	  bool resultFin = (finalDecisions.at(i)).second;
	  std::cout<<"----- accepted: : "<<resultFin << std::endl;
	  if (resultFin) {
	    l1_bit=true;
	    l1_Prescl_L1 =(prescales.at(i)).second;
	  }
	}
      }
      l1bits_l1.push_back( l1_bit );
      trigPrescl_l1_vec_l1.push_back( l1_Prescl_L1 );
      std::cout<<"----- l1_Prescl_L1 : "<<l1_Prescl_L1<< std::endl;
    
    }
   

    for ( auto& l1seed: l1seeds_ ){
      bool sfire=false;
      // float l1Prescl = 0;
      for (unsigned int itrg=0; itrg<result.maxPhysicsTriggers; ++itrg){
        if (result.getAlgoDecisionFinal(itrg)!=1) continue;
        std::string l1trigName = static_cast<const char *>(algoBitToName[itrg]);
        if (l1trigName!=l1seed) continue;
        sfire=true;
	//l1Prescl = l1GtUtils.prescaleFactor(evt, l1trigName, iErrorCode);
        break;
      }
     if (sfire) {
       l1bits.push_back( 1 );
       //   trigPrescl_l1_vec.push_back( l1Prescl );
     }
     else {
       l1bits.push_back( 0 );
       //  trigPrescl_l1_vec.push_back(0);
     }
    }
  }

  // get HLT triggers
  if ( hltResults.failedToGet() ){
    for ( unsigned int ibit = 0; ibit < Npaths; ++ibit){
      hltbits.push_back( 0 );
      trigPrescl_hlt_vec.push_back( 0 );
    }
  } else {
    int Ntrg = hltResults->size();
    for ( auto& hltpath: hltpaths_ ){
      bool fire = false; 
      //float trigPrescl_hlt,trigPrescl_hlt_l1=0;
      std::pair<int,int>  trigPrescl_l1_hlt={0,0};
      for( int itrg = 0; itrg < Ntrg; ++itrg ){
        if ( !hltResults->accept( itrg ) ) continue;
        TString TrigPath = trigName.triggerName( itrg );
	if ( TrigPath.Contains( hltpath ) ){
	  trigPrescl_l1_hlt = hltPrescaleProvider_.prescaleValues(evt, stp, (string)TrigPath);
	  fire=true; 
	  //std::cout << "--- Trigger: " << TrigPath << std::endl;
	  //std::cout << "--- Trigger index: " << itrg << std::endl;
	  //std::cout <<" Prescales: L1 " << trigPrescl_l1_hlt.first <<" HLT "<<trigPrescl_l1_hlt.second<< std::endl;        
	  break; 
	} 
      }
      
      if( fire ) hltbits.push_back( 1 );
      else hltbits.push_back( 0 );
      if( fire ) trigPrescl_hlt_vec.push_back(trigPrescl_l1_hlt.second);
      else trigPrescl_hlt_vec.push_back( 0 );
      if( fire ) trigPrescl_hlt_l1_vec.push_back(trigPrescl_l1_hlt.first);
      else trigPrescl_hlt_l1_vec.push_back( 0 );
    }
  }
 
 
  auto tab  = std::make_unique<nanoaod::FlatTable>(1,"", true);
  for (unsigned int ipath = 0; ipath <Npaths; ++ipath ){
    tab->addColumnValue<uint8_t> (hltpaths_[ipath], hltbits[ipath], "hlt path", nanoaod::FlatTable::UInt8Column);
    tab->addColumnValue<float> (hltpaths_[ipath]+"_ps", trigPrescl_hlt_vec[ipath], "hlt ps", nanoaod::FlatTable::FloatColumn);
    //    tab->addColumnValue<float> (hltpaths_[ipath]+"_l1_ps", trigPrescl_hlt_l1_vec[ipath], "hlt l1 ps", nanoaod::FlatTable::FloatColumn); //NOT working well 
  }
  for (unsigned int iseed = 0; iseed <Nseeds; ++iseed ){
    tab->addColumnValue<uint8_t> (l1seeds_[iseed], l1bits[iseed], "l1 seed", nanoaod::FlatTable::UInt8Column);
    //    tab->addColumnValue<float> (l1seeds_[iseed]+"_ps", trigPrescl_l1_vec[iseed], "l1 ps", nanoaod::FlatTable::FloatColumn);    
    tab->addColumnValue<float> (l1seeds_[iseed]+"_ps_L1", trigPrescl_l1_vec_l1[iseed], "l1 ps from L1", nanoaod::FlatTable::FloatColumn);
  }
    evt.put(std::move(tab));

}


//define this as a plug-in
DEFINE_FWK_MODULE(TrgBitTableProducer);
