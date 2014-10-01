// system include files
#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>

// user include files
#include "Fastbtag/FastbAnalyzer/plugins/EventAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
// #include "DataFormats/​JetReco/​interface/​JetCollection.h"

#include <TH1.h>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// using namespace edm;
// using namespace std;
// using namespace reco; 

EventAnalyzer::EventAnalyzer(const edm::ParameterSet& iConfig)
{
  simmethod = iConfig.getParameter<std::string>("simmethod");
  outhist = iConfig.getParameter<std::string>("outhist");
  jetLabel_ = iConfig.getParameter<edm::InputTag>("jetLabel");
  genmatchingLabel_ = iConfig.getParameter<edm::InputTag>("genmatchingLabel");
  discriminatorLabel_ = iConfig.getParameter<std::vector<edm::InputTag> >("discriminatorLabel");
  discriminatorTokens_ = edm::vector_transform(discriminatorLabel_, [this](edm::InputTag const & tag){return mayConsume<reco::JetFloatAssociation::Container>(tag);});
  tagInfoLabel_ = iConfig.getParameter<std::vector<edm::InputTag> >("tagInfoLabel");
}


EventAnalyzer::~EventAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


void
EventAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


  using namespace std;


  //initialization
  bjet_pt.clear();
  bjet_eta.clear();
  bjet_phi.clear();
  bjet_disc.clear();
  bjet_flavour.clear();

  // Getting required collections to fill histograms
 
  // Get the vector of jets
  edm::Handle<edm::View<reco::Jet> > jets;
  iEvent.getByLabel(jetLabel_, jets);
 
  // for jet flavour
  edm::Handle<reco::JetFlavourInfoMatchingCollection> jetFlavInfoMatch;
  iEvent.getByLabel(genmatchingLabel_, jetFlavInfoMatch);
 
  // Get the vector of jet tags with b-tagging info
  // std::vector<edm::Handle<reco::JetFloatAssociation::Container> > jetDiscriminators;
  std::vector<edm::Handle<reco::JetTagCollection>> jetDiscriminators;
  jetDiscriminators.resize(discriminatorLabel_.size());
  for (size_t i = 0; i < discriminatorLabel_.size(); ++i) {  
      iEvent.getByLabel(discriminatorLabel_[i], jetDiscriminators[i]);
  };
  std::vector<edm::Handle<edm::View<reco::BaseTagInfo> > > jetTagInfos;
  jetTagInfos.resize(tagInfoLabel_.size());
  for (size_t i = 0; i < tagInfoLabel_.size(); ++i) {
    iEvent.getByLabel(tagInfoLabel_[i], jetTagInfos[i]);
  };

  //create map for btag and flavour
  
  FlavourMap flavourmap;
  for (reco::JetFlavourInfoMatchingCollection::const_iterator iter = jetFlavInfoMatch->begin();
       iter != jetFlavInfoMatch->end(); iter++) {
    unsigned int fl = std::abs(iter->second.getPartonFlavour());
    flavourmap.insert(FlavourMap::value_type(iter->first, fl));
  };

  std::vector< BtagMap > btagmap;
  btagmap.resize(discriminatorLabel_.size());
  for (unsigned int idisc = 0; idisc != discriminatorLabel_.size();++idisc){
    const reco::JetTagCollection & bTags = *(jetDiscriminators[idisc].product());
    for (unsigned int ibjet = 0; ibjet != bTags.size(); ++ibjet){
      float disc = bTags[ibjet].second;
      btagmap[idisc].insert(BtagMap::value_type(bTags[ibjet].first,disc));
    };
  };
 
  //Loop over jets
  bjet_pt.resize(jets->size());
  bjet_eta.resize(jets->size());
  bjet_phi.resize(jets->size());
  bjet_flavour.resize(jets->size());
  bjet_disc.resize(discriminatorLabel_.size());
  for (unsigned int i=0; i<bjet_disc.size(); ++i){
    bjet_disc[i].resize(jets->size());
  };


  for (unsigned int index = 0; index < jets->size(); ++index) {
    edm::RefToBase<reco::Jet> jetRef = jets->refAt(index);
    // edm::Ptr<reco::Jet> jetPtr = jets->ptrAt(index);
    if (flavourmap.find(jetRef) == flavourmap.end()) {
      std::cout <<" Cannot access flavour for this jet - not in the Map"<<std::endl;
      bjet_flavour.push_back(-1000);
    }else{
      unsigned int flavour = flavourmap[jetRef];
      bjet_flavour.push_back(flavour);
    };
    for (unsigned int idisc = 0; idisc != discriminatorLabel_.size();++idisc){
      if (btagmap[idisc].find(jetRef) == btagmap[idisc].end()) {
        std::cout <<" Cannot access btag discriminator for this jet - not in the Map"<<std::endl;
        bjet_disc[idisc].push_back(-1000.0);
      }else{
        float disc = btagmap[idisc][jetRef];
        bjet_disc[idisc].push_back(disc);
      };
    };

    bjet_pt.push_back((*jets)[index].pt());
    bjet_eta.push_back((*jets)[index].eta());
    bjet_phi.push_back((*jets)[index].phi());

  };

  for (unsigned int idisc = 0; idisc < bjet_disc.size(); ++idisc){
    if (discriminatorLabel_[idisc].label() == "combinedSecondaryVertexBJetTags") {
      for (unsigned int ijet=0; ijet<bjet_eta.size();++ijet){
        if (bjet_disc[idisc][ijet] > 0.898){
          CTbjet_pt_CSVT->Fill(bjet_pt[ijet]);
          CTbjet_eta_CSVT->Fill(bjet_eta[ijet]);
          CTbjet_phi_CSVT->Fill(bjet_phi[ijet]);
        };
        if (bjet_disc[idisc][ijet] > 0.679){
          CTbjet_pt_CSVM->Fill(bjet_pt[ijet]);
          CTbjet_eta_CSVM->Fill(bjet_eta[ijet]);
          CTbjet_phi_CSVM->Fill(bjet_phi[ijet]);
        };
        if (bjet_disc[idisc][ijet] > 0.244){
          CTbjet_pt_CSVL->Fill(bjet_pt[ijet]);
          CTbjet_eta_CSVL->Fill(bjet_eta[ijet]);
          CTbjet_phi_CSVL->Fill(bjet_phi[ijet]);
        };
      };
    };
      // "jetProbabilityBJetTags"
      // "combinedSecondaryVertexV1BJetTags"
      // "combinedSecondaryVertexSoftPFLeptonV1BJetTags"
      // "combinedSecondaryVertexIVFV2BJetTags"        
      // "trackCountingHighPurBJetTags"
    

  };




}

void 
EventAnalyzer::beginJob()
{
}

void 
EventAnalyzer::endJob() 
{

  TString filename = outhist;
  TFile *histfile = new TFile(filename,"RECREATE");

  CTbjet_pt_CSVT->Write();
  CTbjet_pt_CSVM->Write();
  CTbjet_pt_CSVL->Write();

  CTbjet_eta_CSVT->Write();
  CTbjet_eta_CSVM->Write();
  CTbjet_eta_CSVL->Write();

  CTbjet_phi_CSVT->Write();
  CTbjet_phi_CSVM->Write();
  CTbjet_phi_CSVL->Write();

  histfile->Close();


}

// ------------ method called when starting to processes a run  ------------
/*
void 
EventAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
EventAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
EventAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
EventAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EventAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(EventAnalyzer);
