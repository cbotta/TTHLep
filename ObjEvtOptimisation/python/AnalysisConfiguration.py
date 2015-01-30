import FWCore.ParameterSet.Config as cms
process = cms.Process("Demo")
#####
##   Modules for the analysis
#####
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.GlobalTag.globaltag = 'PHYS14_25_V3::All'
process.GlobalTag.globaltag = 'PLS170_V7AN1::All'
#####
##   Interaction with the analyser
#####
process.demo = cms.EDAnalyzer("ObjEvtOptimisation",
 #Collections of objects  
 pruned     = cms.InputTag("prunedGenParticles"),
 packed     = cms.InputTag("packedGenParticles"),      
 bits       = cms.InputTag("TriggerResults","","HLT"),
 objects    = cms.InputTag("selectedPatTrigger"),
 prescales  = cms.InputTag("patTrigger"),
 vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
 muons      = cms.InputTag("slimmedMuons"),
 electrons  = cms.InputTag("slimmedElectrons"),
 taus       = cms.InputTag("slimmedTaus"),
 photons    = cms.InputTag("slimmedPhotons"),
 jets       = cms.InputTag("slimmedJets"),
 fatjets    = cms.InputTag("slimmedJetsAK8"),
 mets       = cms.InputTag("slimmedMETs"),
 pfCands    = cms.InputTag("packedPFCandidates"), 
 lostTracks = cms.InputTag("lostTracks"),
 #Values for the whole analysis
 mindr_p3 = cms.untracked.double(0.3), 
 mindr_p5 = cms.untracked.double(0.5), 
)
#####
##   Input files
#####
process.source = cms.Source("PoolSource",
 fileNames = cms.untracked.vstring(
 '/store/mc/Spring14miniaod/TTbarH_HToWWTo2L2Nu_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/141029_PU40bx50_PLS170_V6AN2-v1/10000/6619F511-6065-E411-B131-0023AEFDE908.root'
 ),
 #skipEvents = cms.untracked.uint32(25) #Skip the first n evt, or comment this line if you do not want to skip evt
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) ) #Num of evt to be analysed (whatever is the starting evt)
#####
##   Output file
#####
process.TFileService = cms.Service("TFileService",
 fileName = cms.string('tthtree.root')
)
#####
##   Analysis chain
#####
process.p = cms.Path(process.demo)
