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
  jetFlavourInfosToken_ = consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getParameter<edm::InputTag>("jetflavourinfo"));

  jetLabels_ = iConfig.getParameter<std::vector<edm::InputTag> >("jetLabels");
  for ( std::vector<edm::InputTag>::const_iterator jetlabel = jetLabels_.begin(),
    jetlabelEnd = jetLabels_.end(); jetlabel != jetlabelEnd; ++jetlabel ) {
    jetTokens_.push_back( consumes<edm::View<reco::Jet>>( *jetlabel ) );
  };
  
  btaglabels_ = iConfig.getParameter<std::vector<edm::InputTag> >("btagLabels");
  method = iConfig.getParameter<string>("method");
  // for ( std::vector<edm::InputTag>::const_iterator btaglabel = btaglabels_.begin(),
  //   btaglabelEnd = btaglabels_.end(); btaglabel != btaglabelEnd; ++btaglabel ) {
  //   btagTokens_.push_back( consumes<reco::JetTagCollection>( *btaglabel ) );
  // };

  pT = iConfig.getParameter<double>("pT");
  numev = iConfig.getParameter<int>("numev");

  Nbjet=0;
  Nbjet_CSVT=0;
  Nbjet_CSVM=0;
  Nbjet_CSVL=0;
  Nbjet_CSVV1T=0;
  Nbjet_CSVV1M=0;
  Nbjet_CSVV1L=0;
  Nbjet_CSVSLV1T=0;
  Nbjet_CSVSLV1M=0;
  Nbjet_CSVSLV1L=0;
  Nbjet_CSVIVFV2T=0;
  Nbjet_CSVIVFV2M=0;
  Nbjet_CSVIVFV2L=0;
  Nbjet_JPT=0;
  Nbjet_JPM=0;
  Nbjet_JPL=0;
  Nbjet_TCHPT=0;

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


  edm::Handle<reco::JetFlavourInfoMatchingCollection> theJetFlavourInfos;
  iEvent.getByToken(jetFlavourInfosToken_, theJetFlavourInfos );
  for ( reco::JetFlavourInfoMatchingCollection::const_iterator j  = theJetFlavourInfos->begin();
                                     j != theJetFlavourInfos->end();
                                     ++j ) {
    // const reco::Jet *aJet = (*j).first.get();
    reco::JetFlavourInfo aInfo = (*j).second;
    if (abs(aInfo.getPartonFlavour())==5){
      Nbjet++;
    };

  };

  //Loop over each b-tagger
  for (unsigned int icoll = 0; icoll != btaglabels_.size(); ++icoll) {
    //get bjets
    edm::Handle<reco::JetTagCollection> bTagHandle;
    bool Validbtag = iEvent.getByLabel(btaglabels_[icoll], bTagHandle);
    if (!Validbtag) continue;
    const reco::JetTagCollection & bTags = *(bTagHandle.product());
    for (unsigned int i = 0; i < bTags.size(); ++i) {

      if (btaglabels_[icoll].label() =="combinedSecondaryVertexBJetTags"){
        if (bTags[i].second > 0.898){
          Nbjet_CSVT++;
        };
        if (bTags[i].second > 0.679){
          Nbjet_CSVM++;
        };
        if (bTags[i].second > 0.244){
          Nbjet_CSVL++;
        };
      };

      if (btaglabels_[icoll].label()=="jetProbabilityBJetTags"){
        if (bTags[i].second > 0.790){
          Nbjet_JPT++;
        };
        if (bTags[i].second > 0.545){
          Nbjet_JPM++;
        };
        if (bTags[i].second > 0.275){
          Nbjet_JPL++;
        };
      };

      if (btaglabels_[icoll].label()=="combinedSecondaryVertexV1BJetTags"){
        if (bTags[i].second > 0.405){
          Nbjet_CSVV1L++;
        };
        if (bTags[i].second > 0.783){
          Nbjet_CSVV1M++;
        };
        if (bTags[i].second > 0.920){
          Nbjet_CSVV1T++;
        };
      };

      if (btaglabels_[icoll].label()=="combinedSecondaryVertexSoftPFLeptonV1BJetTags"){
        if (bTags[i].second > 0.527){
          Nbjet_CSVSLV1L++;
        };
        if (bTags[i].second > 0.756){
          Nbjet_CSVSLV1M++;
        };
        if (bTags[i].second > 0.859){
          Nbjet_CSVSLV1T++;
        };
      };

      if (btaglabels_[icoll].label()=="combinedSecondaryVertexIVFV2BJetTags"){
        if (bTags[i].second > 0.423){
          Nbjet_CSVIVFV2L++;
        };
        if (bTags[i].second > 0.814){
          Nbjet_CSVIVFV2M++;
        };
        if (bTags[i].second > 0.941){
          Nbjet_CSVIVFV2T++;
        };
      };

      if (btaglabels_[icoll].label()=="trackCountingHighPurBJetTags"){
        if (bTags[i].second > 3.41){
          Nbjet_TCHPT++;
        };
      };
    };
  };


  //Get genparticles
  // edm::Handle<reco::GenParticleCollection> genParticles;
  // iEvent.getByToken(genParticle_, genParticles );
  // for (reco::GenParticleCollection::const_iterator iter=genParticles->begin();iter!=genParticles->end();++iter){
  //   if(!(   iter->status()==3   )) continue;
  //   if(!(   iter->pdgId()==23   )) continue;
  //   if((&*iter)->pdgId() == 5 || (&*iter)->pdgId() == -5) {
  //     genjet_pt.push_back(iter->pt());
  //     genjet_eta.push_back(iter->eta());
  //     genjet_phi.push_back(iter->phi());
  //     genjet_mass.push_back(iter->mass());
  //     N_genjet++;
  //   };
  // };

  // //Filling histogram with correctly tagged b jets
  // edm::Handle< edm::View<reco::Jet> > jets;
  // edm::Handle<reco::JetTagCollection> bJetTags;
  // iEvent.getByLabel("combinedSecondaryVertexBJetTags", bJetTags);
  // iEvent.getByLabel("ak4PFJetsCHS", jets); //ak4PFJets are used to calculate CSV
  // for (edm::View<reco::Jet>::const_iterator ijet = jets->begin(); ijet != jets->end(); ++ijet ){
  //   for(unsigned int index = 0; index < jets->size(); ++index) {
  //     float disc = (*bJetTags)[jets->refAt(index)];
  //     if (disc > 0.898){
  //       CTbjet_pt_CSVT->Fill(ijet->pt());
  //       CTbjet_eta_CSVT->Fill(ijet->eta());
  //       CTbjet_phi_CSVT->Fill(ijet->phi());
  //     };
  //     if (disc > 0.679){
  //       CTbjet_pt_CSVM->Fill(ijet->pt());
  //       CTbjet_eta_CSVM->Fill(ijet->eta());
  //       CTbjet_phi_CSVM->Fill(ijet->phi());
  //     };
  //     if (disc > 0.244){
  //       CTbjet_pt_CSVL->Fill(ijet->pt());
  //       CTbjet_eta_CSVL->Fill(ijet->eta());
  //       CTbjet_phi_CSVL->Fill(ijet->phi());
  //     };
  //   };
  // };






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
  cout << "The number of CSVT bjets: "  <<Nbjet_CSVT<<endl;
  cout << "efficiency of CSVT: "        <<       ((float)Nbjet_CSVT/Nbjet)         <<endl;
  cout << "efficiency of CSVM: "        <<       ((float)Nbjet_CSVM/Nbjet)         <<endl;
  cout << "efficiency of CSVL: "        <<       ((float)Nbjet_CSVL/Nbjet)         <<endl;
  cout << "efficiency of CSVV1T: "      <<       ((float)Nbjet_CSVV1T/Nbjet)       <<endl;
  cout << "efficiency of CSVV1M: "      <<       ((float)Nbjet_CSVV1M/Nbjet)       <<endl;
  cout << "efficiency of CSVV1L: "      <<       ((float)Nbjet_CSVV1L/Nbjet)       <<endl;
  cout << "efficiency of CSVSLV1T: "    <<       ((float)Nbjet_CSVSLV1T/Nbjet)     <<endl;
  cout << "efficiency of CSVSLV1M: "    <<       ((float)Nbjet_CSVSLV1M/Nbjet)     <<endl;
  cout << "efficiency of CSVSLV1L: "    <<       ((float)Nbjet_CSVSLV1L/Nbjet)     <<endl;
  cout << "efficiency of CSVIVFV2T: "   <<       ((float)Nbjet_CSVIVFV2T/Nbjet)    <<endl;
  cout << "efficiency of CSVIVFV2M: "   <<       ((float)Nbjet_CSVIVFV2M/Nbjet)    <<endl;
  cout << "efficiency of CSVIVFV2L: "   <<       ((float)Nbjet_CSVIVFV2L/Nbjet)    <<endl;
  cout << "efficiency of JPT: "         <<       ((float)Nbjet_JPT/Nbjet)          <<endl;
  cout << "efficiency of JPM: "         <<       ((float)Nbjet_JPM/Nbjet)          <<endl;
  cout << "efficiency of JPL: "         <<       ((float)Nbjet_JPL/Nbjet)          <<endl;
  cout << "efficiency of TCHPT: "       <<       ((float)Nbjet_TCHPT/Nbjet)        <<endl;

  ofstream eff_file_CSVT;
  ofstream eff_file_CSVM;
  ofstream eff_file_CSVL;
  ofstream eff_file_CSVV1T;
  ofstream eff_file_CSVV1M;
  ofstream eff_file_CSVV1L;
  ofstream eff_file_CSVSLV1T;
  ofstream eff_file_CSVSLV1M;
  ofstream eff_file_CSVSLV1L;
  ofstream eff_file_CSVIVFV2T;
  ofstream eff_file_CSVIVFV2M;
  ofstream eff_file_CSVIVFV2L;
  ofstream eff_file_JPT;
  ofstream eff_file_JPM;
  ofstream eff_file_JPL; 
  ofstream eff_file_TCHPT;
  ofstream Nbjet_file;

  string dir = "/afs/cern.ch/work/k/klo/fastsim/validation_bjets/lsf/eff_file/";
  string filetag = std::to_string(int(pT))+"_"+std::to_string(numev);

  Nbjet_file.open(dir+"Nbjet.txt");
  eff_file_CSVT.open(dir+method+"_CSVT_"+filetag+".txt");
  eff_file_CSVM.open(dir+method+"_CSVM_"+filetag+".txt");
  eff_file_CSVL.open(dir+method+"_CSVL_"+filetag+".txt");
  eff_file_CSVV1T.open(dir+method+"_CSVV1T_"+filetag+".txt");
  eff_file_CSVV1M.open(dir+method+"_CSVV1M_"+filetag+".txt");
  eff_file_CSVV1L.open(dir+method+"_CSVV1L_"+filetag+".txt");
  eff_file_CSVSLV1T.open(dir+method+"_CSVSLV1T_"+filetag+".txt");
  eff_file_CSVSLV1M.open(dir+method+"_CSVSLV1M_"+filetag+".txt");
  eff_file_CSVSLV1L.open(dir+method+"_CSVSLV1L_"+filetag+".txt");
  eff_file_CSVIVFV2T.open(dir+method+"_CSVIVFV2T_"+filetag+".txt");
  eff_file_CSVIVFV2M.open(dir+method+"_CSVIVFV2M_"+filetag+".txt");
  eff_file_CSVIVFV2L.open(dir+method+"_CSVIVFV2L_"+filetag+".txt");
  eff_file_JPT.open(dir+method+"_JPT_"+filetag+".txt");
  eff_file_JPM.open(dir+method+"_JPM_"+filetag+".txt");
  eff_file_JPL.open(dir+method+"_JPL_"+filetag+".txt");
  eff_file_TCHPT.open(dir+method+"_TCHPT_"+filetag+".txt");

  Nbjet_file << Nbjet <<endl;
  eff_file_CSVT <<       ((float)Nbjet_CSVT/Nbjet)         <<endl;
  eff_file_CSVM <<       ((float)Nbjet_CSVM/Nbjet)         <<endl;
  eff_file_CSVL <<       ((float)Nbjet_CSVL/Nbjet)         <<endl;
  eff_file_CSVV1T <<       ((float)Nbjet_CSVV1T/Nbjet)       <<endl;
  eff_file_CSVV1M <<       ((float)Nbjet_CSVV1M/Nbjet)       <<endl;
  eff_file_CSVV1L <<       ((float)Nbjet_CSVV1L/Nbjet)       <<endl;
  eff_file_CSVSLV1T <<       ((float)Nbjet_CSVSLV1T/Nbjet)     <<endl;
  eff_file_CSVSLV1M <<       ((float)Nbjet_CSVSLV1M/Nbjet)     <<endl;
  eff_file_CSVSLV1L <<       ((float)Nbjet_CSVSLV1L/Nbjet)     <<endl;
  eff_file_CSVIVFV2T <<       ((float)Nbjet_CSVIVFV2T/Nbjet)    <<endl;
  eff_file_CSVIVFV2M <<       ((float)Nbjet_CSVIVFV2M/Nbjet)    <<endl;
  eff_file_CSVIVFV2L <<       ((float)Nbjet_CSVIVFV2L/Nbjet)    <<endl;
  eff_file_JPT <<       ((float)Nbjet_JPT/Nbjet)          <<endl;
  eff_file_JPM <<       ((float)Nbjet_JPM/Nbjet)          <<endl;
  eff_file_JPL  <<       ((float)Nbjet_JPL/Nbjet)          <<endl;
  eff_file_TCHPT <<       ((float)Nbjet_TCHPT/Nbjet)        <<endl;

  Nbjet_file.close();
  eff_file_CSVT.close();
  eff_file_CSVM.close();
  eff_file_CSVL.close();
  eff_file_CSVV1T.close();
  eff_file_CSVV1M.close();
  eff_file_CSVV1L.close();
  eff_file_CSVSLV1T.close();
  eff_file_CSVSLV1M.close();
  eff_file_CSVSLV1L.close();
  eff_file_CSVIVFV2T.close();
  eff_file_CSVIVFV2M.close();
  eff_file_CSVIVFV2L.close();
  eff_file_JPT.close();
  eff_file_JPM.close();
  eff_file_JPL.close();
  eff_file_TCHPT.close();


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
