import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("digiAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("SimGeneral.PileupInformation.AddPileupSummary_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )

process.source = cms.Source("PoolSource",
                            #fileNames=cms.untracked.vstring("file:../../../noPU/DoublePhoton_step1__withDigis__replay2.root"))
                            #fileNames=cms.untracked.vstring("file:../../../withLowPU/DoublePhoton_step1__withDigis__replay2.root"))
                            fileNames=cms.untracked.vstring("file:../../../withLowPU/DoublePhoton_step1__withPuDigis__replay1.root"))

process.digiAna = cms.EDAnalyzer('DigiAnalyzer',
                                 # mixing?
                                 EBdigiCollection = cms.InputTag("mix",""),
                                 EEdigiCollection = cms.InputTag("mix",""),
                                 # standard digis
                                 #EBdigiCollection = cms.InputTag("simEcalDigis","ebDigis"),
                                 #EEdigiCollection = cms.InputTag("simEcalDigis","eeDigis"),
                                 genParticles = cms.InputTag("genParticles"),
                                 PileUpTag = cms.InputTag('addPileupInfo')
                                 )

process.p = cms.Path(process.digiAna)

