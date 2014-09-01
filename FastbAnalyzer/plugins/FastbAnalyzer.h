#ifndef FastbAnalyzer_H
#define FastbAnalyzer_H

#include <memory>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// ParticleFlow
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"

//Jets
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

//Calculations
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/CandUtils/interface/pdgIdUtils.h"

//Genparticles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// ROOT
#include "TH1F.h"
#include "TFile.h"



class FastbAnalyzer : public edm::EDAnalyzer {

public:
  FastbAnalyzer(const edm::ParameterSet& ps);
  virtual ~FastbAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::GenParticleCollection>  genParticle_;
  std::vector<double> genjet_pt;
  std::vector<double> genjet_eta;
  std::vector<double> genjet_phi;
  std::vector<double> genjet_mass;
  std::vector<double> genjet_discrim;


  std::vector<edm::InputTag> jetLabels_;
  std::vector< edm::EDGetTokenT< edm::View<reco::Jet> > > jetTokens_;
  std::vector<double> jet_pt;
  std::vector<double> jet_eta;
  std::vector<double> jet_phi;
  std::vector<double> jet_mass;

  std::vector<edm::InputTag> btaglabels_;
  std::vector<edm::EDGetTokenT<reco::JetTagCollection>> btagTokens_;
  std::vector<double> bjet_pt;
  std::vector<double> bjet_eta;
  std::vector<double> bjet_phi;
  std::vector<double> bjet_mass;
  std::vector<double> bjet_discrim;

  unsigned int N_matched_bjet;
  unsigned int N_genjet;
  double matching_radius;
  double btagcut;
  double pT;
  int numev;
  std::string taggerlabel;

};
#endif











