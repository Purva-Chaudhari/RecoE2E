import FWCore.ParameterSet.Config as cms

EGTagger = cms.EDProducer('EGTagger'
    , photonCollection = cms.InputTag('gedPhotons')
    , EGFrames = cms.InputTag('EGFrames','EGFrames')
    , EGModelName = cms.string('tfModels/sample.onnx')
    )
