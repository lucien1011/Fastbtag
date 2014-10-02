import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone()

from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
process.jetFlavourInfosAK4PFJets = ak4JetFlavourInfos.clone()
process.jetFlavourInfosAK4PFJets.jets = cms.InputTag("ak4PFJetsCHS")

process.demo = cms.EDAnalyzer('EventAnalyzer',
      				simmethod = cms.string("full"),
      				outhist = cms.string("/afs/cern.ch/work/k/klo/fastsim/validation_bjets/lsf/rootfile/full_hist.root"),
                              jetLabel = cms.InputTag("ak4PFJetsCHS"),
                              genmatchingLabel = cms.InputTag("jetFlavourInfosAK4PFJets"),
                              discriminatorLabel = cms.VInputTag(
                                    # cms.InputTag("jetBProbabilityBJetTags"),
                                    # cms.InputTag("jetProbabilityBJetTags"),
                                    # cms.InputTag("trackCountingHighPurBJetTags"),
                                    # cms.InputTag("trackCountingHighEffBJetTags"),
                                    # cms.InputTag("simpleSecondaryVertexHighEffBJetTags"),
                                    # cms.InputTag("simpleSecondaryVertexHighPurBJetTags"),
                                    cms.InputTag("combinedSecondaryVertexBJetTags")
                              ),
                              tagInfoLabel = cms.VInputTag(
                                    cms.InputTag("secondaryVertexTagInfos"),
                                    cms.InputTag("impactParameterTagInfos")
                              )
)

# from PhysicsTools.PatAlgos.patInputFiles_cff import filesRelValProdTTbarAODSIM
# process.source.fileNames = filesRelValProdTTbarAODSIM

process.source = cms.Source("PoolSource", 
                              fileNames = cms.untracked.vstring("file:/afs/cern.ch/work/k/klo/fastsim/validation_bjets/lsf/Full/pgun/reco_10000.root"),
                              # fileNames = cms.untracked.vstring(filesRelValProdTTbarAODSIM),
                              skipBadFiles = cms.untracked.bool(True),
                              duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                              )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.AnalyzeEvent = cms.Path(process.selectedHadronsAndPartons*process.jetFlavourInfosAK4PFJets*process.demo)

# process.outpath = cms.EndPath(process.out)






