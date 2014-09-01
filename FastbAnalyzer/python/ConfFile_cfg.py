import FWCore.ParameterSet.Config as cms
import sys
sys.path.insert(1,'/afs/cern.ch/work/k/klo/mypylib/')
import myutils

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
options.register ('tagger',
                  "combinedSecondaryVertexBJetTags", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "btag disciminator")
options.register ('taggerlabel',
                  "combinedSecondaryVertexTight", # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "btag disciminator label for eff file")
options.register ('btagcut',
                  0.898, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.float,          # string, int, or float
                  "btag disciminator cut")
options.register ('matchingradius',
                  0.4, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.float,          # string, int, or float
                  "matching radius for genParticles and bjets")
options.register ('jobsplit',
                  1, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.int,          # string, int, or float
                  "the number of files for each pt point")

options.parseArguments()

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.numEv) )
jobsplit=options.jobsplit


fullfilelist=myutils.getfilelistwithext("/afs/cern.ch/work/k/klo/fastsim/validation_bjets/lsf/Full/pgun/",".root")
filelist = cms.untracked.vstring()
templist = []
for i in range(0,jobsplit):
      filetag = "_".join([str(options.pt),str(options.numEv/options.jobsplit),str(i)])
      inputfile = 'file:/afs/cern.ch/work/k/klo/fastsim/validation_bjets/lsf/Full/pgun/reco_pt'+filetag+'.root'
      # inputfile = 'file:/afs/cern.ch/work/k/klo/fastsim/validation_bjets/lsf/Fast/pgun/pgun_5_pt'+filetag+'.root'
      print inputfile.replace("file:","")
      if inputfile.replace("file:","") in fullfilelist:
            templist.append(inputfile)
filelist.extend(templist)


process.source = cms.Source("PoolSource", 
                              fileNames = filelist,
                              skipBadFiles = cms.untracked.bool(True),
                              duplicateCheckMode = cms.untracked.string("noDuplicateCheck")



                              )

process.demo = cms.EDAnalyzer('FastbAnalyzer',
								GenParticle = cms.InputTag("genParticles"),
								jetLabels = cms.VInputTag("ak4PFJets"),
								btagLabels = cms.VInputTag(options.tagger),
								btagcut = cms.double(options.btagcut),
								taggerlabel = cms.string(options.taggerlabel),
								matchingradius = cms.double(options.matchingradius),
								pT = cms.double(options.pt),
								numev = cms.int32(options.numEv)
)

process.p = cms.Path(process.demo)
