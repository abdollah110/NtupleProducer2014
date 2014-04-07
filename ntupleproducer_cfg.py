import FWCore.ParameterSet.Config as cms

process = cms.Process("NtupleProducer")

process.out = cms.OutputModule("PoolOutputModule",
                               fileName=cms.untracked.string('PATLayer1_Output.fromAOD_full.root'),
                               SelectEvents=cms.untracked.PSet(SelectEvents=cms.vstring('p')),
                               outputCommands=cms.untracked.vstring('drop *')
                               )


process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("RecoEcal.EgammaClusterProducers.ecalClusteringSequence_cff")
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("PhysicsTools/PatAlgos/patSequences_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi")
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesHZZElectronIdentificationV06_cfi")
process.load("Analysis.NtupleProducer.New_hTozzTo4leptonsPFIsolationProducer_cff")
process.load("Analysis.NtupleProducer.JES_Uncertainty_FR_cff")
process.MessageLogger = cms.Service("MessageLogger")
process.load("RecoLocalCalo/EcalRecAlgos/EcalSeverityLevelESProducer_cfi")
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi') #new location of the code
process.load("CMGTools.External.pujetidsequence_cff") #load PU JetID sequence
process.load("EventFilter.HcalRawToDigi.hcallasereventfilter2012_cff") #laser correction
process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff") #MET correction
process.load("RecoMET.METPUSubtraction.mvaPFMET_leptons_cff") #new location of the code

#################################################   Samples and GlobalTag   ############################
#Source File
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
    'file:/afs/cern.ch/work/c/ccaillol/DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_PU_S7_START52_V9_v2.root'
    )
                            )
process.maxEvents = cms.untracked.PSet(
                                       input=cms.untracked.int32(300)
                                       )

isMC = True      #comment it for Data
#isMC = False     #comment it for MC

if isMC:
    HLTProcessName = 'HLT'
    process.GlobalTag.globaltag = cms.string('START53_V23::All') #updated JEC

else:
    HLTProcessName = 'HLT'
    process.GlobalTag.globaltag = cms.string('FT_53_V21_AN4::All') 
    process.patPFMet.addGenMET = cms.bool(False)

#################################################   EDANALYZER   ##################################
process.myanalysis = cms.EDAnalyzer("NtupleProducer",

                                    #OutPut Filles (NTuples)
                                    HistOutFile=cms.untracked.string('output_Ntuples.root'),

                                    # Include or Exclude the objects
                                    Include_HPSTau=cms.bool(True),
                                    Include_Muon=cms.bool(True),
                                    Include_Electron=cms.bool(True),
                                    Include_Jet=cms.bool(True),
                                    Include_JetCorrection=cms.bool(False),
                                    Include_MET=cms.bool( True),
                                    Include_L1=cms.bool( True),#new to have L1 trigger objects
                                    Include_MET_Uncertaity=cms.bool(True),
                                    Include_GenPartiles=cms.bool(True),
                                    Include_HLT=cms.bool(True),
                                    Include_Vertex=cms.bool(True),
                                    Is_MC=cms.bool(isMC),
                                    
                                     # storing only certain trigger strings
                                    filterTriggerResults=cms.bool(True),
                                    el_trigger_name = cms.string("Ele"),
                                    mu_trigger_name = cms.string("Mu"),

                                    #vetrex and Tracks
                                    vertices=cms.InputTag('offlinePrimaryVertices'),
                                    tracks=cms.InputTag("generalTracks"),

                                    #MET
                                    met=cms.InputTag("met"),
                                    PFmet=cms.InputTag("pfMet"),
                                    tcmet=cms.InputTag("tcMet"),
                                    Type1CorMET=cms.InputTag("patType1CorrectedPFMet"),
                                    MVAmet=cms.InputTag("pfMEtMVA"),
                                    #Jets
                                    bjets=cms.InputTag("trackCountingHighPurBJetTags"),
                                    PFAK5=cms.InputTag("selectedPatJets"),
                                    rhoJetsLabel=cms.InputTag("kt6PFJets", "rho"),
				    rhoCenChargedPU=cms.InputTag("kt6PFJetsCentralChargedPileUp", "rho"),
				    rhoCenNeutral=cms.InputTag("kt6PFJetsCentralNeutral", "rho"),
				    rhoCenNeutralTight=cms.InputTag("kt6PFJetsCentralNeutralTight", "rho"),

                                    #Leptons
                                    preselectedelectrons=cms.InputTag("selectedPatElectrons"),
                                    preselectedmuons=cms.InputTag("selectedPatMuons"),
                                    preselectedHPSTaus=cms.InputTag("selectedPatTaus"),
                                    vertexCollectionForLeptonIP=cms.InputTag("selectPrimaryVertex"),
                                    tauPtCut=cms.double(15.0),

                                    #Trigger
                                    srcTriggerResults=cms.InputTag("TriggerResults", "", HLTProcessName),

                                    # MC Information
                                    PileUpInfo=cms.InputTag("addPileupInfo"),
                                    genParticlesInfo=cms.InputTag("genParticles"),

                                    #Trigger and TriggerMatching
                                    triggerEvent=cms.InputTag("patTriggerEvent"),
                                    tauMatch_Loose=cms.string('tauTriggerMatchHLTTausLoose'),
                                    tauMatch_Medium=cms.string('tauTriggerMatchHLTTausMedium'),
                                    muonMatch_Loose=cms.string('muonTriggerMatchHLTMuonsLoose'),
				    electronMatch_Loose=cms.string('muonTriggerMatchHLTMuonsLoose'),
                                    muonMatch_Medium=cms.string('muonTriggerMatchHLTMuonsLoose'),
                                    electronMatch_Medium=cms.string('muonTriggerMatchHLTMuonsLoose'),
                                    jetMatch_Loose=cms.string('jetTriggerMatchHLTJetsLoose'),
                                    jetMatch_Medium=cms.string('jetTriggerMatchHLTJetsLoose'),

                                    puJetIdFlag=cms.InputTag("puJetMva","fullId"),
                                    #rhoProducer = cms.InputTag('kt6PFJetsForRhoComputationVoronoi','rho') #new This line was there before, I have changed it with the following one
				    rhoProducer = cms.InputTag('kt6PFJets','rho') #new
                                    )
#################################################   PFIsolation  ################################
from CommonTools.ParticleFlow.pfParticleSelection_cff import *

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso
process.eleIsoSequence = setupPFElectronIso(process, 'selectedPatElectrons')

# cone vetos as used in Htautau
process.elPFIsoValueChargedAll04NoPFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueChargedAll04PFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueChargedAll03NoPFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueChargedAll03PFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')

process.elPFIsoValueGamma04NoPFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.08)','EcalEndcaps:ConeVeto(0.08)')
process.elPFIsoValueGamma04PFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.08)','EcalEndcaps:ConeVeto(0.08)')

#################################################   scrapingVeto + hcal laser veto  ################################
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter=cms.untracked.bool(True),
                                    debugOn=cms.untracked.bool(False),
                                    numtrack=cms.untracked.uint32(10),
                                    thresh=cms.untracked.double(0.25)
                                    )
if not isMC:
    process.dataFilter = cms.Sequence(process.hcallLaserEvent2012Filter+ process.scrapingVeto )
else:
    process.dataFilter = cms.Sequence()

#################################################   Good Primary Vertex ################################
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.selectPrimaryVertex = cms.EDFilter(
                                    "PrimaryVertexObjectFilter",
                                    filterParams=pvSelector.clone(minNdof=cms.double(4.0), maxZ=cms.double(24.0)),
                                    src=cms.InputTag('offlinePrimaryVertices')
                                    )

process.GoodVertexFilter = cms.EDFilter("VertexSelector",
                                           src = cms.InputTag("offlinePrimaryVertices"),
                                           cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
                                           filter = cms.bool(True)   # otherwise it won't filter the events, just produce an empty vertex collection.
                                        )

#################################################   PAT APPLICATIONS   ##################################

#Removing MC Matching
from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['All'])

##Switch to ak5PFJets (L2 and L3 Corrections are included)
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process,
                    cms.InputTag('ak5PFJets'),
                    doJTA=True,
                    doBTagging=True,
                    jetCorrLabel=('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])),
                    doType1MET=False,
                    genJetCollection=cms.InputTag("ak5GenJets"),
                    doJetID=True,
                    jetIdLabel="ak5"
                    )
if not isMC:
    process.patJetCorrFactors.levels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual')

# Electron Id
process.mvaID = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
process.ElectronIDs = cms.Sequence(process.simpleEleIdSequence)
process.patElectrons.electronIDSources = cms.PSet(
    mvaTrigV0=cms.InputTag("mvaTrigV0"),
    mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0"),
    simpleEleId95relIso = cms.InputTag("simpleEleId95relIso")
                                                      )


#change PV source for electrons and muons
process.patElectrons.pvSrc = cms.InputTag("selectPrimaryVertex")
process.patMuons.pvSrc = cms.InputTag("selectPrimaryVertex")

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

process.pfPileUp.Enable = cms.bool(True)

process.AddPUInfo = process.myanalysis.clone()
process.AddGenInfo = process.myanalysis.clone()

if not isMC:
    process.AddPUInfo = process.myanalysis.clone(PileUpInfo=cms.InputTag(""))
    process.AddGenInfo = process.myanalysis.clone(genParticlesInfo=cms.InputTag(""))

process.elPFIsoValueChargedAll04PFId.deposits[0].vetos= cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueGamma04PFId.deposits[0].vetos= cms.vstring('EcalBarrel:ConeVeto(0.08)','EcalEndcaps:ConeVeto(0.08)')
process.elPFIsoValueChargedAll04NoPFId.deposits[0].vetos= cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueGamma04NoPFId.deposits[0].vetos= cms.vstring('EcalBarrel:ConeVeto(0.08)','EcalEndcaps:ConeVeto(0.08)')


#Trigger matching "HLT_SingleIsoTau20_Trk15_MET25_v*"
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
#new trigger matching
process.muonTriggerMatchHLTMuonsLoose = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                               src=cms.InputTag("selectedPatMuons"),
                                                               matched=cms.InputTag("patTrigger"),
                                                               matchedCuts=cms.string('path("HLT_IsoMu24_eta2p1_v*")'),
                                                               maxDPtRel=cms.double(0.5),
                                                               maxDeltaR=cms.double(0.5),
                                                               resolveAmbiguities=cms.bool(True),
                                                               resolveByMatchQuality=cms.bool(True)
                                                                                                      )

process.tauTriggerMatchHLTTausLoose = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                src=cms.InputTag("selectedPatTaus"),
                                                matched=cms.InputTag("patTrigger"),
                                                matchedCuts=cms.string('path("HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk1_eta2p1_v*") || path("HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v*")'),
                                                maxDPtRel=cms.double(0.5),
                                                maxDeltaR=cms.double(0.5),
                                                resolveAmbiguities=cms.bool(True),
                                                resolveByMatchQuality=cms.bool(True)
                                                )

process.tauTriggerMatchHLTTausMedium = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                src=cms.InputTag("selectedPatTaus"),
                                                matched=cms.InputTag("patTrigger"),
                                                matchedCuts=cms.string('path("HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v*") || path("HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30_v*") || path("HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v*") || path("HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1_v*") || path("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_v*")'),
                                                maxDPtRel=cms.double(0.5),
                                                maxDeltaR=cms.double(0.5),
                                                resolveAmbiguities=cms.bool(True),
                                                resolveByMatchQuality=cms.bool(True)
                                                )

process.jetTriggerMatchHLTJetsLoose = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                               src=cms.InputTag("selectedPatJets"),
                                                               matched=cms.InputTag("patTrigger"),
                                                               matchedCuts=cms.string('path("HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v*") || path("HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30_v*") || path("HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v*")'),
                                                               maxDPtRel=cms.double(0.5),
                                                               maxDeltaR=cms.double(0.5),
                                                               resolveAmbiguities=cms.bool(True),
                                                               resolveByMatchQuality=cms.bool(True)
                                                                                                      )

#Vertexing
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
#SoftLepton
#process.load("RecoBTag.SoftLepton.softElectronCandProducer_cfi") #new (doesn't work, is idle?)

###--------------------------------------------------------------!!!!!!!!!!!!!!!!!!
# This part is just for Met Uncertainty studies
# apply type I/type I + II PFMEt corrections to pat::MET object
# and estimate systematic uncertainties on MET
process.load("JetMETCorrections/Type1MET/pfMETsysShiftCorrections_cfi")
from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties
#from JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi import *
#runMEtUncertainties(process)
#from PhysicsTools.PatUtils.tools.metUncertaintyTools import *
#print process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_mc 
###
if isMC == False:
    runMEtUncertainties(process,
                        electronCollection = cms.InputTag('cleanPatElectrons'),
                        photonCollection = '',
                        muonCollection = 'selectedPatMuons',
                        tauCollection = 'selectedPatTaus',
                        jetCollection = cms.InputTag('selectedPatJets'),
                        jetCorrLabel = 'L2L3Residual',
                        doSmearJets = False,
                        makeType1corrPFMEt = True,
                        makeType1p2corrPFMEt = False,
                        makePFMEtByMVA = False,
                        makeNoPileUpPFMEt = False,
                        doApplyType0corr = False,
                        sysShiftCorrParameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data,#new (remove [0] at the end, runA instead of run ABC)
                        #sysShiftCorrParameter = process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_data,
                        doApplySysShiftCorr = False,
                        addToPatDefaultSequence = False,
                        )


else:
    runMEtUncertainties(process,
                        electronCollection = cms.InputTag('cleanPatElectrons'),
                        photonCollection = '',
                        muonCollection = 'selectedPatMuons',
                        tauCollection = 'selectedPatTaus',
                        jetCollection = cms.InputTag('selectedPatJets'),
                        jetCorrLabel = 'L3Absolute',
                        doSmearJets = True,
                        makeType1corrPFMEt = True,
                        makeType1p2corrPFMEt = False,
                        makePFMEtByMVA = False,
                        makeNoPileUpPFMEt = False,
                        doApplyType0corr = False,
                        sysShiftCorrParameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc,#new 
                        doApplySysShiftCorr = False,
                        addToPatDefaultSequence = False,
                        )

process.patJetsNotOverlappingWithLeptonsForMEtUncertainty.checkOverlaps.taus.preselection = cms.string('pt > 20 && abs(eta) < 2.3 && leadPFChargedHadrCand.isNonnull() && leadPFChargedHadrCand.pt() > 5 && leadPFChargedHadrCand.mva_e_pi() < 0.6 && tauID("decayModeFinding") > 0.5 && tauID("againstMuonTight") > 0.5 && tauID("againstElectronLoose") > 0.5')
process.patJetsNotOverlappingWithLeptonsForMEtUncertainty.checkOverlaps.muons.preselection = cms.string('pt > 25 && abs(eta) < 2.4')

#Making a cleanPatJet similar to smeraedJet
process.cleanPatJets.checkOverlaps.taus.src=cms.InputTag("selectedPatTaus")
process.cleanPatJets.checkOverlaps.taus.preselection= cms.string('pt > 20 && abs(eta) < 2.3 && leadPFChargedHadrCand.isNonnull() && leadPFChargedHadrCand.pt() > 5 && leadPFChargedHadrCand.mva_e_pi() < 0.6 && tauID("decayModeFinding") > 0.5 && tauID("againstMuonTight") > 0.5 && tauID("againstElectronLoose") > 0.5')
process.cleanPatJets.checkOverlaps.muons.src= cms.InputTag("selectedPatMuons")
process.cleanPatJets.checkOverlaps.muons.preselection = cms.string('pt > 25 && abs(eta) < 2.4')
process.cleanPatJets.checkOverlaps.tkIsoElectrons.requireNoOverlaps= cms.bool(True)
process.cleanPatJets.checkOverlaps.photons.requireNoOverlaps= cms.bool(True)

###--------------------------------------------------------------
## MVA MET configuration

from RecoMET.METPUSubtraction.mvaPFMET_leptons_cfi import isotaus 

process.isotaus.discriminators = cms.VPSet(
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),       selectionCut=cms.double(0.5)),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr"),           selectionCut=cms.double(0.5)),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseElectronRejection"), selectionCut=cms.double(0.5)),
    cms.PSet( discriminator=cms.InputTag("hpsPFTauDiscriminationByLooseMuonRejection"),     selectionCut=cms.double(0.5))
            )

# fix for JEC not defined above eta 5.2
process.selectedPatJets.cut = cms.string('abs(eta) < 5.0')

#################################################   PATH of CONFIG FILE ##################################
#print process.dumpPython()
process.p = cms.Path (
                      process.dataFilter *
                      process.selectPrimaryVertex *
                      process.PFTau * #prescribed by TAU POG
                      process.mvaID   *
		      process.ElectronIDs  *
                      process.inclusiveVertexing *
#                      process.softElectronCands* #new, removed
                      process.pfParticleSelectionSequence +
                      process.type0PFMEtCorrection + #new (MVA MET twiki suggested to add this these two lines if the release is CMSSW_5_3_11 or newer, should the "+" be a "*"?)
                      process.patPFMETtype0Corr + #new
                       process.patDefaultSequence +
                       process.eleIsoSequence*
                      process.pfMuonIsolationSequence                     *
                      process.metUncertaintySequence*
                      process.producePatPFMETCorrections *
                      process.puJetIdSqeuence *
                      process.pfMEtMVAsequence *
                      process.myanalysis *
                      process.AddPUInfo *
                      process.AddGenInfo
                      )

#################################################   NEEDED FOR TRIGGER MATCHING   #######################

from PhysicsTools.PatAlgos.tools.trigTools import *

#switchOnTrigger( process ) # This is optional and can be omitted.
switchOnTriggerMatching(process, ['tauTriggerMatchHLTTausLoose','tauTriggerMatchHLTTausMedium','jetTriggerMatchHLTJetsLoose','muonTriggerMatchHLTMuonsLoose'])

# Switch to selected PAT objects in the trigger matching
removeCleaningFromTriggerMatching(process)
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger(process)
#process.patTrigger.addL1Algos     = cms.bool( True ) # default: 'False'
#process.patTrigger.l1ExtraCenJet  = cms.InputTag( 'l1extraParticles', 'Central'    , 'RECO' )
process.patTrigger.processName = HLTProcessName
process.patTriggerEvent.processName = HLTProcessName
