// -*- C++ -*-
//
// Package:    Fastbtag/FastbAnalyzer
// Class:      FastbAnalyzer
// 
/**\class FastbAnalyzer FastbAnalyzer.cc Fastbtag/FastbAnalyzer/plugins/FastbAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Kin Ho Lo
//         Created:  Tue, 26 Aug 2014 12:56:25 GMT
//
//


// system include files
#include <memory>

// user include files
#include "Fastbtag/FastbAnalyzer/plugins/FastbAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace edm;
using namespace std;
using namespace reco; 

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FastbAnalyzer::FastbAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  genParticle_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticle"));

  jetLabels_ = iConfig.getParameter<std::vector<edm::InputTag> >("jetLabels");
  for ( std::vector<edm::InputTag>::const_iterator jetlabel = jetLabels_.begin(),
    jetlabelEnd = jetLabels_.end(); jetlabel != jetlabelEnd; ++jetlabel ) {
    jetTokens_.push_back( consumes<edm::View<reco::Jet>>( *jetlabel ) );
  };
  
  btaglabels_ = iConfig.getParameter<std::vector<edm::InputTag> >("btagLabels");
  for ( std::vector<edm::InputTag>::const_iterator btaglabel = btaglabels_.begin(),
    btaglabelEnd = btaglabels_.end(); btaglabel != btaglabelEnd; ++btaglabel ) {
    btagTokens_.push_back( consumes<reco::JetTagCollection>( *btaglabel ) );
  };

  matching_radius = iConfig.getParameter<double>("matchingradius");
  btagcut = iConfig.getParameter<double>("btagcut");
  taggerlabel = iConfig.getParameter<string>("taggerlabel");
  pT = iConfig.getParameter<double>("pT");
  numev = iConfig.getParameter<int>("numev");

  N_matched_bjet=0;
  N_genjet=0;

}


FastbAnalyzer::~FastbAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
FastbAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  genjet_pt.clear();
  genjet_eta.clear();
  genjet_phi.clear();
  genjet_mass.clear();
  genjet_discrim.clear();

  jet_pt.clear();
  jet_eta.clear();
  jet_phi.clear();
  jet_mass.clear();

  bjet_pt.clear();
  bjet_eta.clear();
  bjet_phi.clear();
  bjet_mass.clear();
  bjet_discrim.clear();

  double deltaR_value = matching_radius*matching_radius;

  //Get genparticles
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticle_, genParticles );
  for (reco::GenParticleCollection::const_iterator iter=genParticles->begin();iter!=genParticles->end();++iter){
    if(!(   iter->status()==3   )) continue;
    if(!(   iter->pdgId()==23   )) continue;
    if((&*iter)->pdgId() == 5 || (&*iter)->pdgId() == -5) {
      genjet_pt.push_back(iter->pt());
      genjet_eta.push_back(iter->eta());
      genjet_phi.push_back(iter->phi());
      genjet_mass.push_back(iter->mass());
      N_genjet++;
    };
  };

  //Get recojets
  for (unsigned int icoll = 0; icoll < jetLabels_.size(); ++icoll ) {
    edm::Handle<edm::View<reco::Jet>>  pfJetCollection;
    bool ValidPFJets = iEvent.getByToken(jetTokens_[icoll], pfJetCollection );
    if(!ValidPFJets) continue;
    edm::View<reco::Jet> const & pfjets = *pfJetCollection;
    for (unsigned int ijet=0; ijet<pfjets.size(); ++ijet){
      jet_pt.push_back(pfjets[ijet].pt());
      jet_eta.push_back(pfjets[ijet].eta());
      jet_phi.push_back(pfjets[ijet].phi());
      jet_mass.push_back(pfjets[ijet].mass());
    };
  };

  //Loop over each b-tagger
  for (unsigned int icoll = 0; icoll != btaglabels_.size(); ++icoll) {
    //get bjets
    edm::Handle<reco::JetTagCollection> bTagHandle;
    bool Validbtag = iEvent.getByToken(btagTokens_[icoll], bTagHandle);
    if (!Validbtag) continue;
    const reco::JetTagCollection & bTags = *(bTagHandle.product());
    for (unsigned int i = 0; i < bTags.size(); ++i) {
      if (bTags[i].second > btagcut) {
        for (unsigned int ijet=0; ijet<jet_eta.size(); ++ijet){
          if (deltaR2(jet_eta[ijet],jet_phi[ijet],bTags[i].first->eta(),bTags[i].first->phi())<deltaR_value){
            bjet_pt.push_back(bTags[i].first->pt());
            bjet_eta.push_back(bTags[i].first->eta());
            bjet_phi.push_back(bTags[i].first->phi());
            bjet_mass.push_back(bTags[i].first->mass());
            bjet_discrim.push_back(bTags[i].second);
          };
        };
      };
    };

    //matching with genjet
    for (unsigned int ibjet = 0; ibjet < bjet_eta.size(); ++ibjet) {
      for (unsigned int igenjet = 0; igenjet < genjet_eta.size(); ++igenjet) {
        if (deltaR2(bjet_eta[ibjet],bjet_phi[ibjet],genjet_eta[igenjet],genjet_phi[igenjet])<deltaR_value){
          N_matched_bjet++;
        };
      };
    };
  };

}


// ------------ method called once each job just before starting event loop  ------------
void 
FastbAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FastbAnalyzer::endJob() 
{
  cout << "The number of matched bjets: "<<N_matched_bjet<<endl;
  cout << "efficiency of btag: "<<((float)N_matched_bjet/N_genjet)<<endl;

  ofstream f_eff;

  std::string f_eff_name = taggerlabel+"_pt"+std::to_string(int(pT))+"_"+std::to_string(numev)+".txt";
  f_eff.open(f_eff_name);
  f_eff << ((float)N_matched_bjet/N_genjet)<<endl;
  f_eff.close();

}

// ------------ method called when starting to processes a run  ------------
/*
void 
FastbAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
FastbAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
FastbAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
FastbAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FastbAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FastbAnalyzer);
