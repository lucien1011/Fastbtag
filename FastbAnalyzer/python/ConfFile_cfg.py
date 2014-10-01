import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')
options.register ('numEv',
                  10000, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "Number of events to process")

options.register ('pt',
                  500, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "specify pt of b parton")
options.register ('jobsplit',
                  1, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "the number of files for each pt point")

options.register ('method',
                  "full", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "simulation method")

options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
jobsplit=options.jobsplit

from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone()

from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
process.jetFlavourInfosAK4PFJets = ak4JetFlavourInfos.clone()

numEv_per_job = options.numEv/options.jobsplit

templist=[]

if options.method == "full":
      # for i in range(0,jobsplit):
            # filetag = "_".join(["btag","pt"+str(options.pt),str(numEv_per_job),str(i)])
      filetag = "_".join(["btag","pt"+str(options.pt),str(options.numEv)])
      filename = "file:/afs/cern.ch/work/k/klo/fastsim/validation_bjets/lsf/Full/pgun/"+filetag+".root"
      templist.append(filename)
else:
      if options.method == "fast":
            filetag = "_".join(["btag","pt"+str(options.pt),str(options.numEv)])
            filename = "file:/afs/cern.ch/work/k/klo/fastsim/validation_bjets/lsf/Fast/pgun/"+filetag+".root"
            templist.append(filename)

filelist=cms.untracked.vstring()
filelist.extend(templist)

process.source = cms.Source("PoolSource", 
                              fileNames = filelist,
                              skipBadFiles = cms.untracked.bool(True),
                              duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
                              )

process.demo = cms.EDAnalyzer('FastbAnalyzer',
                                                GenParticle = cms.InputTag("genParticles"),
                                                jetflavourinfo = cms.InputTag("jetFlavourInfosAK4PFJets"),
                                                jetLabels = cms.VInputTag("ak4PFJets"),
                                                btagLabels = cms.VInputTag("combinedSecondaryVertexBJetTags","jetProbabilityBJetTags","combinedSecondaryVertexV1BJetTags","combinedSecondaryVertexSoftPFLeptonV1BJetTags","combinedSecondaryVertexIVFV2BJetTags","trackCountingHighPurBJetTags"),
                                                pT = cms.double(options.pt),
                                                numev = cms.int32(options.numEv),
                                                method = cms.string(options.method)
)

rootdir = "/afs/cern.ch/work/k/klo/fastsim/validation_bjets/lsf/rootfile/"
rootfile = rootdir+"Fastb_"+str(options.pt)+".root"

# process.output = cms.OutputModule("PoolOutputModule",
#     fileName = cms.untracked.string(rootfile)
# )



process.p = cms.Path(process.selectedHadronsAndPartons*process.jetFlavourInfosAK4PFJets*process.demo)

# process.out = cms.EndPath(process.output)



