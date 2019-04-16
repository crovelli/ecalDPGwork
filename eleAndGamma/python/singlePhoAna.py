import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("singlePhoAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("SimGeneral.PileupInformation.AddPileupSummary_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )

process.source = cms.Source("PoolSource",
                            #fileNames=cms.untracked.vstring("/store/user/crovelli/Ecal/noPU/DoublePhoton_step2_fromPrivateGenSim.root"))
                            #fileNames=cms.untracked.vstring("/store/user/crovelli/Ecal/withLowPUonlyFirstBin/DoublePhoton_step2_fromPrivateGenSim.root"))
                            fileNames=cms.untracked.vstring("file:/afs/cern.ch/work/c/crovelli/dpgTests/RecoTests_10_2_5/src/noPU/DoublePhoton_step2_fromPrivateGenSim_fix6.root"))


process.TFileService = cms.Service("TFileService",fileName = cms.string("singlePhotonTree.root"))
process.singlePhoAna = cms.EDAnalyzer('SinglePhoAnalyzer',
                                      reducedBarrelRecHitCollection = cms.InputTag('reducedEcalRecHitsEB'),
                                      reducedEndcapRecHitCollection = cms.InputTag('reducedEcalRecHitsEE'),
                                      vertices = cms.InputTag("offlinePrimaryVertices"),
                                      genParticles = cms.InputTag("genParticles"),
                                      photons = cms.InputTag("gedPhotons"),
                                      PileUpTag = cms.InputTag('addPileupInfo')
                                      )

process.p = cms.Path(process.singlePhoAna)

