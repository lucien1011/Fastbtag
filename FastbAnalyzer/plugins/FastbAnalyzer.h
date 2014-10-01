#ifndef FastbAnalyzer_H
#define FastbAnalyzer_H

#include <iostream>
#include <fstream>

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
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

//Calculations
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/CandUtils/interface/pdgIdUtils.h"

//Genparticles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// ROOT
#include "TH1F.h"
#include "TH1D.h"
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

  edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;
  
  std::vector<edm::InputTag> jetLabels_;
  std::vector< edm::EDGetTokenT< edm::View<reco::Jet> > > jetTokens_;

  std::vector<edm::InputTag> btaglabels_;
  // std::vector<edm::EDGetTokenT<reco::JetTagCollection>> btagTokens_;

  unsigned int Nbjet;
  unsigned int Nbjet_CSVT;
  unsigned int Nbjet_CSVM;
  unsigned int Nbjet_CSVL;
  unsigned int Nbjet_CSVV1T;
  unsigned int Nbjet_CSVV1M;
  unsigned int Nbjet_CSVV1L;
  unsigned int Nbjet_CSVSLV1T;
  unsigned int Nbjet_CSVSLV1M;
  unsigned int Nbjet_CSVSLV1L;
  unsigned int Nbjet_CSVIVFV2T;
  unsigned int Nbjet_CSVIVFV2M;
  unsigned int Nbjet_CSVIVFV2L;
  unsigned int Nbjet_JPT;
  unsigned int Nbjet_JPM;
  unsigned int Nbjet_JPL;
  unsigned int Nbjet_TCHPT;

  double btagcut;
  double pT;
  int numev;
  std::string method;

  //file to store histograms
  // TFile *histfile;

  // CT stands for correctly tagged
  TH1D *CTbjet_pt_CSVT = new TH1D("CTbjet_pt_CSVT","CTbjet_pt_CSVT",100,0,2000);
  TH1D *CTbjet_pt_CSVM = new TH1D("CTbjet_pt_CSVM","CTbjet_pt_CSVM",100,0,2000);
  TH1D *CTbjet_pt_CSVL = new TH1D("CTbjet_pt_CSVL","CTbjet_pt_CSVL",100,0,2000);
  TH1D CTbjet_pt_CSVV1T;
  TH1D CTbjet_pt_CSVV1M;
  TH1D CTbjet_pt_CSVV1L;
  TH1D CTbjet_pt_CSVSLV1T;
  TH1D CTbjet_pt_CSVLSLV1M;
  TH1D CTbjet_pt_CSVLSLV1L;
  TH1D CTbjet_pt_CSVLIVFV2T;
  TH1D CTbjet_pt_CSVIVFV2M;
  TH1D CTbjet_pt_CSVIVFV2L;
  TH1D CTbjet_pt_JPT;
  TH1D CTbjet_pt_JPM;
  TH1D CTbjet_pt_JPL;
  TH1D CTbjet_pt_TCHPT;

  TH1D *CTbjet_eta_CSVT = new TH1D("CTbjet_eta_CSVT","CTbjet_pt_CSVT",50,0,2);
  TH1D *CTbjet_eta_CSVM = new TH1D("CTbjet_eta_CSVM","CTbjet_pt_CSVM",50,0,2);
  TH1D *CTbjet_eta_CSVL = new TH1D("CTbjet_eta_CSVL","CTbjet_pt_CSVL",50,0,2);
  TH1D CTbjet_eta_CSVV1T;
  TH1D CTbjet_eta_CSVV1M;
  TH1D CTbjet_eta_CSVV1L;
  TH1D CTbjet_eta_CSVSLV1T;
  TH1D CTbjet_eta_CSVSLV1M;
  TH1D CTbjet_eta_CSVSLV1L;
  TH1D CTbjet_eta_CSVIVFV2T;
  TH1D CTbjet_eta_CSVIVFV2M;
  TH1D CTbjet_eta_CSVIVFV2L;
  TH1D CTbjet_eta_JPT;
  TH1D CTbjet_eta_JPM;
  TH1D CTbjet_eta_JPL;
  TH1D CTbjet_eta_TCHPT;

  TH1D *CTbjet_phi_CSVT = new TH1D("CTbjet_phi_CSVT","CTbjet_pt_CSVT",50,-3.14,3.14);
  TH1D *CTbjet_phi_CSVM = new TH1D("CTbjet_phi_CSVM","CTbjet_pt_CSVM",50,-3.14,3.14);
  TH1D *CTbjet_phi_CSVL = new TH1D("CTbjet_phi_CSVL","CTbjet_pt_CSVL",50,-3.14,3.14);
  TH1D CTbjet_phi_CSVV1T;
  TH1D CTbjet_phi_CSVV1M;
  TH1D CTbjet_phi_CSVV1L;
  TH1D CTbjet_phi_CSVSLV1T;
  TH1D CTbjet_phi_CSVSLV1M;
  TH1D CTbjet_phi_CSVSLV1L;
  TH1D CTbjet_phi_CSVIVFV2T;
  TH1D CTbjet_phi_CSVIVFV2M;
  TH1D CTbjet_phi_CSVIVFV2L;
  TH1D CTbjet_phi_JPT;
  TH1D CTbjet_phi_JPM;
  TH1D CTbjet_phi_JPL;
  TH1D CTbjet_phi_TCHPT;
  
  // MT stands for mistagged
  TH1D MTbjet_pt_CSVT;
  TH1D MTbjet_pt_CSVM;
  TH1D MTbjet_pt_CSVL;
  TH1D MTbjet_pt_CSVV1T;
  TH1D MTbjet_pt_CSVV1M;
  TH1D MTbjet_pt_CSVV1L;
  TH1D MTbjet_pt_CSVSLV1T;
  TH1D MTbjet_pt_CSVSLV1M;
  TH1D MTbjet_pt_CSVSLV1L;
  TH1D MTbjet_pt_CSVIVFV2T;
  TH1D MTbjet_pt_CSVIVFV2M;
  TH1D MTbjet_pt_CSVIVFV2L;
  TH1D MTbjet_pt_JPT;
  TH1D MTbjet_pt_JPM;
  TH1D MTbjet_pt_JPL;
  TH1D MTbjet_pt_TCHPT;

  TH1D MTbjet_eta_CSVT;
  TH1D MTbjet_eta_CSVM;
  TH1D MTbjet_eta_CSVL;
  TH1D MTbjet_eta_CSVV1T;
  TH1D MTbjet_eta_CSVV1M;
  TH1D MTbjet_eta_CSVV1L;
  TH1D MTbjet_eta_CSVSLV1T;
  TH1D MTbjet_eta_CSVSLV1M;
  TH1D MTbjet_eta_CSVSLV1L;
  TH1D MTbjet_eta_CSVIVFV2T;
  TH1D MTbjet_eta_CSVIVFV2M;
  TH1D MTbjet_eta_CSVIVFV2L;
  TH1D MTbjet_eta_JPT;
  TH1D MTbjet_eta_JPM;
  TH1D MTbjet_eta_JPL;
  TH1D MTbjet_eta_TCHPT;

  TH1D MTbjet_phi_CSVT;
  TH1D MTbjet_phi_CSVM;
  TH1D MTbjet_phi_CSVL;
  TH1D MTbjet_phi_CSVV1T;
  TH1D MTbjet_phi_CSVV1M;
  TH1D MTbjet_phi_CSVV1L;
  TH1D MTbjet_phi_CSVSLV1T;
  TH1D MTbjet_phi_CSVSLV1M;
  TH1D MTbjet_phi_CSVSLV1L;
  TH1D MTbjet_phi_CSVIVFV2T;
  TH1D MTbjet_phi_CSVIVFV2M;
  TH1D MTbjet_phi_CSVIVFV2L;
  TH1D MTbjet_phi_JPT;
  TH1D MTbjet_phi_JPM;
  TH1D MTbjet_phi_JPL;
  TH1D MTbjet_phi_TCHPT;

};
#endif











