import FWCore.ParameterSet.Config as cms

process = cms.Process("CMSNTUPLER")

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('RecoVertex/PrimaryVertexProducer/OfflinePrimaryVertices_cfi')

process.configurationMetadata = cms.untracked.PSet(
  version    = cms.untracked.string('BACON_000'),
  annotation = cms.untracked.string('AODSIM'),
  name       = cms.untracked.string('CMSNtuplerProduction')
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(        
        '/store/mc/Summer13/DYToEE_M-20_TuneZ2star_14TeV-pythia6/GEN-SIM-RECO/UpgradePhase1Age0START_DR61SLHCx_PU140Bx25_STAR17_61_V1A-v1/10000/FE0AFAB4-A1CB-E211-AC35-002354EF3BE1.root'        
    )
)

# Global Tag
#process.GlobalTag.globaltag = 'START61_V11::All'
process.GlobalTag.globaltag = 'STAR17_61_V1A::All'
#process.GlobalTag.globaltag = 'START53_V7A::All'
#process.GlobalTag.globaltag = 'DES17_61_V5::All'

#For Transient Track Builder
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')

#For BTagging Discriminators
process.load('CMSAna.CMSNtupler.btagging_cff')

#For MC jet matching
process.load('CMSAna.CMSNtupler.jetflavorMatching_cff')

#Turn Off Alignment because it's not consistent with new geometry
#process.trackerGeometryDB.applyAlignment = cms.bool(False)

#Load Ntupler Module
process.load("CMSAna.CMSNtupler.CMSNtupler_cfi")
process.myCMSNtupler.debugLevel = cms.int32(0)
process.myCMSNtupler.FillEGRegressionVars = cms.bool(True)
process.myCMSNtupler.GenJetPtMin = cms.double(15.0)
process.myCMSNtupler.JetPtMin = cms.double(20.0)


#Define ntupler sequence
process.ntupler_sequence = cms.Sequence(
    process.myJetFlavourId*
    process.newJetBtagging*
    process.myCMSNtupler
    )

process.ntupler_step  = cms.Path(process.ntupler_sequence)

# schedule definition
process.schedule = cms.Schedule(process.ntupler_step)
