import FWCore.ParameterSet.Config as cms

process = cms.Process("NtupleProducer")

process.out = cms.OutputModule("PoolOutputModule",
                               fileName=cms.untracked.string('PATLayer1_Output.fromAOD_full.root'),
                               # save only events passing the full path
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
#process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/EventContent/EventContent_cff')
#process.load("RecoTauTag/RecoTau/PFRecoTauProducer_cfi")
#process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
#process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.load("RecoTauTag/RecoTau/PFRecoTauProducer_cfi")
process.load("PhysicsTools/PatAlgos/patSequences_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("ElectroWeakAnalysis.WENu.simpleEleIdSequence_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi")
process.load("RecoEgamma.ElectronIdentification.cutsInCategoriesHZZElectronIdentificationV06_cfi")
#process.load("TauAnalysis.RecoTools.recoVertexSelection_cff")
process.load("Analysis.NtupleProducer.New_hTozzTo4leptonsPFIsolationProducer_cff")
process.load("Analysis.NtupleProducer.JES_Uncertainty_FR_cff")
process.MessageLogger = cms.Service("MessageLogger")
process.load("RecoLocalCalo/EcalRecAlgos/EcalSeverityLevelESProducer_cfi")
process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
process.load("CMGTools.External.pujetidsequence_cff") #load PU JetID sequence
process.load("EventFilter.HcalRawToDigi.hcallasereventfilter2012_cff") # laser correction
process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff") #MET correction
process.load('JetMETCorrections.METPUSubtraction.mvaPFMET_leptons_cff') #MVA MET

#################################################   Samples and GlobalTag   ############################
#Source File
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
                            #                            'file:/tmp/abdollah/WToTauNu_TuneZ2_AODSIM_E7TeV_FlatDist10_2011EarlyData_50ns_START311_V1G1.root')
#                                                        'file:/afs/cern.ch/work/a/abdollah/Sample/2012Data_DoubleMu.root')
                                       #                 'file:/tmp/abdollah/02B4578A-009D-E111-837D-001D0967DA3A.root')
#                                                        'file:/afs/cern.ch/work/j/jez/ntuples/2012Data_DoubleMu.root'
#    'file:/afs/cern.ch/work/j/jez/ntuples/tauID/mc/Summer12/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v2/0003/F41D7301-BE9B-E111-9C03-002481E0D6A0.root'
#    'file:/afs/cern.ch/work/j/jez/ntuples/2012C_v1/DoubleElectron/BC850526-F2C6-E111-8C96-003048CF99BA.root' #data
    'file:/afs/cern.ch/work/j/jez/ntuples/FEF4E41A-46D4-E111-9594-0025B3E06424.root' #MC
#    'file:/afs/cern.ch/work/j/jez/private/CMS/analysis/ZH/53x/PostMoriond/production/CMSSW_5_3_3_patch3/src/crab_0_130226_161836/res/pickevents_1_1_TUE.root'
#    'file:/afs/cern.ch/work/j/jez/ntuples/analysis/mc/Summer11//WH_ZH_TTH_HToTauTau_M-130_7TeV-pythia6-tauola/AODSIM/PU_S4_START42_V11-v1/0000/0EC8AE34-128E-E011-A1F3-90E6BA19A25E.root'
    )
#                            'file:/afs/cern.ch/work/a/abdollah/Sample/2012Data_DoubleMu.root')
                            #                                                                                    'file:/tmp/abdollah/skim_WtoTauNu_data_RunData42Xv4.root')
                            #                            'rfio:/castor/cern.ch/user/v/veelken/CMSSW_5_2_x/skims/data/data2012runA_doubleMu_AOD_1_1_Fzg.root')
                            )
process.maxEvents = cms.untracked.PSet(
                                       input=cms.untracked.int32(10)
                                       )

isMC = True      #comment it for Data
#isMC = False     #comment it for MC

if isMC:
    HLTProcessName = 'HLT'
    process.GlobalTag.globaltag = cms.string('START53_V20::All') #updated JEC
#    process.GlobalTag.globaltag = cms.string('START41_V0::All')
#    process.GlobalTag.globaltag = cms.string('START311_V2::All')

else:
    HLTProcessName = 'HLT'
    process.GlobalTag.globaltag = cms.string('FT53_V10A_AN3::All') #2012C_v1 re-reco Aug24
    process.patPFMet.addGenMET = cms.bool(False)
    #  process.GlobalTag.globaltag = cms.string('FT_53_V6C_AN3::All') #2012A+B re-reco Jul13 + 2012A Aug6 re-reco
    #  process.GlobalTag.globaltag = cms.string('GR_P_V42_AN3::All') #2012C_v2 +2012D Prompt Reco
    #  process.GlobalTag.globaltag = cms.string('FT_P_V42C_AN3::All') #2012C Dec11 EcalRecovery
    #  process.GlobalTag.globaltag = cms.string('FT_P_V43E_AN3::All') #2012D 16Jan2013 PixelAlignment recovery


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
                                    Include_MET_Uncertaity=cms.bool(False),
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


                                    # CIC electron Identification
                                #     eleID_VeryLooseTag=cms.InputTag("eidVeryLoose"),
#                                     eleID_LooseTag=cms.InputTag("eidLoose"),
#                                     eleID_MediumTag=cms.InputTag("eidMedium"),
#                                     eleID_TightTag=cms.InputTag("eidTight"),

                                    #Trigger
                                    srcTriggerResults=cms.InputTag("TriggerResults", "", HLTProcessName),

                                    # MC Information
                                    PileUpInfo=cms.InputTag("addPileupInfo"),
                                    genParticlesInfo=cms.InputTag("genParticles"),

                                                                      #Trigger and TriggerMatching
                                    triggerEvent=cms.InputTag("patTriggerEvent"),
                                    tauMatch_Loose=cms.string('tauTriggerMatchHLTTausLoose'),
                                    tauMatch_Medium=cms.string('tauTriggerMatchHLTTausMedium'),
                                    electronMatch_Loose=cms.string('electronTriggerMatchHLTElectronsLoose'),
                                    muonMatch_Loose=cms.string('muonTriggerMatchHLTMuonsLoose'),


                                    puJetIdFlag=cms.InputTag("puJetMva","fullId"),
                                    #rho
                                    rhoProducer = cms.InputTag('kt6PFJetsForRhoComputationVoronoi','rho')
                                    )
#################################################   PFIsolation  ################################
from CommonTools.ParticleFlow.pfParticleSelection_cff import *

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso
process.eleIsoSequence = setupPFElectronIso(process, 'selectedPatElectrons')

#from RecoParticleFlow.PFProducer.electronPFIsolationValues_cff import *
# cone vetos as used in Htautau
process.elPFIsoValueChargedAll04NoPFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueChargedAll04PFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueChargedAll03NoPFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueChargedAll03PFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')

process.elPFIsoValueGamma04NoPFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.08)','EcalEndcaps:ConeVeto(0.08)')
process.elPFIsoValueGamma04PFIdPFIso.deposits[0].vetos = cms.vstring('EcalBarrel:ConeVeto(0.08)','EcalEndcaps:ConeVeto(0.08)')

#process.elPFIsoValuePU04NoPFIdPFIso.deposits[0].vetos = cms.vstring('')
#process.elPFIsoValuePU04PFIdPFIso.deposits[0].vetos = cms.vstring('')


#process.phoIsoSequence = setupPFPhotonIso(process, 'photons')
# process.load("CommonTools.ParticleFlow.pfParticleSelection_cff")
# process.pfPileUpCandidates = cms.EDProducer(
#                                             "TPPFCandidatesOnPFCandidates",
#                                             enable=cms.bool(True),
#                                             verbose=cms.untracked.bool(False),
#                                             name=cms.untracked.string("pileUpCandidates"),
#                                             topCollection=cms.InputTag("pfNoPileUp"),
#                                             bottomCollection=cms.InputTag("particleFlow"),
#                                             )

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
#from RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_DataTuning_cfi import *
#process.CiC_Ele_Id = cms.Sequence(process.eidVeryLoose + process.eidLoose + process.eidMedium + process.eidTight + process.eidSuperTight)

# Adding different WP for electron
#process.patElectronIDs = cms.Sequence(process.simpleEleIdSequence)
#process.makePatElectrons = cms.Sequence(process.patElectronIDs * process.patElectronIsolation * process.patElectrons)
#process.patElectrons.addElectronID = cms.bool(True)
process.mvaID = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
process.patElectrons.electronIDSources = cms.PSet(
    mvaTrigV0=cms.InputTag("mvaTrigV0"),
    mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
#                                                   simpleEleId95relIso=cms.InputTag("simpleEleId95relIso"),
#                                                   simpleEleId90relIso=cms.InputTag("simpleEleId90relIso"),
#                                                   simpleEleId85relIso=cms.InputTag("simpleEleId85relIso"),
#                                                   simpleEleId80relIso=cms.InputTag("simpleEleId80relIso"),
#                                                   simpleEleId70relIso=cms.InputTag("simpleEleId70relIso"),
#                                                   simpleEleId60relIso=cms.InputTag("simpleEleId60relIso"),

#                                                   simpleEleId95cIso=cms.InputTag("simpleEleId95cIso"),
#                                                   simpleEleId90cIso=cms.InputTag("simpleEleId90cIso"),
#                                                   simpleEleId85cIso=cms.InputTag("simpleEleId85cIso"),
#                                                   simpleEleId80cIso=cms.InputTag("simpleEleId80cIso"),
#                                                   simpleEleId70cIso=cms.InputTag("simpleEleId70cIso"),
#                                                   simpleEleId60cIso=cms.InputTag("simpleEleId60cIso"),

#                                                   eidVeryLoose=cms.InputTag("eidVeryLoose"),
#                                                   eidLoose=cms.InputTag("eidLoose"),
#                                                   eidMedium=cms.InputTag("eidMedium"),
#                                                   eidTight=cms.InputTag("eidTight"),
#                                                   eidSuperTight=cms.InputTag("eidSuperTight"),
                                                      )

#add pat conversions
# process.patConversions = cms.EDProducer("PATConversionProducer",
#                                         # input collection
#                                         electronSource = cms.InputTag("cleanPatElectrons")
#                                         # this should be your last selected electron collection name since currently index is used to match with electron later. We can fix this using reference pointer. ,
#                                         )

#change PV source for electrons and muons
process.patElectrons.pvSrc = cms.InputTag("selectPrimaryVertex")
process.patMuons.pvSrc = cms.InputTag("selectPrimaryVertex")

from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

process.pfPileUp.Enable = cms.bool(True)
process.tauTriggerMatchHLTTausMedium = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
                                                                                                      src=cms.InputTag("selectedPatTaus"),
                                                                                                      matched=cms.InputTag("patTrigger"),
                                                                                                      matchedCuts=cms.string('path( "HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v*" )'),
                                                                                                      maxDPtRel=cms.double(0.5),
                                                                                                      maxDeltaR=cms.double(0.5),
                                                                                                      resolveAmbiguities=cms.bool(True),
                                                                                                      resolveByMatchQuality=cms.bool(True)
                                                                                                      )

process.AddPUInfo = process.myanalysis.clone()
process.AddGenInfo = process.myanalysis.clone()
if not isMC:
    process.AddPUInfo = process.myanalysis.clone(PileUpInfo=cms.InputTag(""))
    process.AddGenInfo = process.myanalysis.clone(genParticlesInfo=cms.InputTag(""))

#From MOnica
#process.pfPileUp.PFCandidates = 'particleFlow'
#process.pfNoPileUp.bottomCollection = 'particleFlow'
#process.pfPileUpIso.PFCandidates = 'particleFlow'
#process.pfNoPileUpIso.bottomCollection = 'particleFlow'
#process.mvaTrigV0.electronTag= cms.InputTag("selectedPatElectrons")
#process.mvaNonTrigV0.electronTag= cms.InputTag("selectedPatElectrons")
#process.pfAllPhotons.src=cms.InputTag("pfNoPileUp")
#process.pfAllChargedHadrons.src=cms.InputTag("pfNoPileUp")
#process.pfAllNeutralHadrons.src=cms.InputTag("pfNoPileUp")
#process.pfAllChargedParticles.src=cms.InputTag("pfNoPileUp")


process.elPFIsoValueChargedAll04PFId.deposits[0].vetos= cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueGamma04PFId.deposits[0].vetos= cms.vstring('EcalBarrel:ConeVeto(0.08)','EcalEndcaps:ConeVeto(0.08)')
process.elPFIsoValueChargedAll04NoPFId.deposits[0].vetos= cms.vstring('EcalBarrel:ConeVeto(0.01)','EcalEndcaps:ConeVeto(0.015)')
process.elPFIsoValueGamma04NoPFId.deposits[0].vetos= cms.vstring('EcalBarrel:ConeVeto(0.08)','EcalEndcaps:ConeVeto(0.08)')


#Trigger matching "HLT_SingleIsoTau20_Trk15_MET25_v*"
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")
process.tauTriggerMatchHLTTausLoose = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
                                                src=cms.InputTag("selectedPatTaus"),
                                                matched=cms.InputTag("patTrigger"),
                                                matchedCuts=cms.string('path( "HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v*" )'),
                                                maxDPtRel=cms.double(0.5),
                                                maxDeltaR=cms.double(0.5),
                                                resolveAmbiguities=cms.bool(True),
                                                resolveByMatchQuality=cms.bool(True)
                                                )
process.tauTriggerMatchHLTTausMedium = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
                                                src=cms.InputTag("selectedPatTaus"),
                                                matched=cms.InputTag("patTrigger"),
                                                matchedCuts=cms.string('path( "HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v*" )'),
                                                maxDPtRel=cms.double(0.5),
                                                maxDeltaR=cms.double(0.5),
                                                resolveAmbiguities=cms.bool(True),
                                                resolveByMatchQuality=cms.bool(True)
                                                )

process.electronTriggerMatchHLTElectronsLoose = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
                                                               src=cms.InputTag("selectedPatElectrons"),
                                                               matched=cms.InputTag("patTrigger"),
                                                               matchedCuts=cms.string('path( "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*" )'),
                                                               maxDPtRel=cms.double(0.5),
                                                               maxDeltaR=cms.double(0.5),
                                                               resolveAmbiguities=cms.bool(True),
                                                               resolveByMatchQuality=cms.bool(True)
                                                                                                      )

process.muonTriggerMatchHLTMuonsLoose = cms.EDProducer("PATTriggerMatcherDRDPtLessByR",
                                                     src=cms.InputTag("selectedPatMuons"),
                                                       matched=cms.InputTag("patTrigger"),
                                                       matchedCuts=cms.string('path( "HLT_Mu17_Mu8_v*" ) || path( "HLT_Mu17_TkMu8_v*" )'),
                                                       maxDPtRel=cms.double(0.5),
                                                       maxDeltaR=cms.double(0.5),
                                                       resolveAmbiguities=cms.bool(True),
                                                       resolveByMatchQuality=cms.bool(True)
                                                                                                     )
#Vertexing
process.load("RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff")
#SoftLepton
process.load("RecoBTag.SoftLepton.softElectronCandProducer_cfi")

###--------------------------------------------------------------
# This part is just for Met Uncertainty studies
# apply type I/type I + II PFMEt corrections to pat::MET object
# and estimate systematic uncertainties on MET
process.load("JetMETCorrections/Type1MET/pfMETsysShiftCorrections_cfi")
from PhysicsTools.PatUtils.tools.metUncertaintyTools import runMEtUncertainties
#from JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi import *
#runMEtUncertainties(process)
#from PhysicsTools.PatUtils.tools.metUncertaintyTools import *
print process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_mc[0]

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
                        sysShiftCorrParameter = process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_data[0],
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
                        sysShiftCorrParameter = process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_mc[0],
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

from JetMETCorrections.METPUSubtraction.mvaPFMET_leptons_cfi import isotaus

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
#                      process.scrapingVeto *
#                      process.hcallLaserEvent2012Filter+
#                      process.GoodVertexFilter *   # REJECT events without good vertex
                      process.selectPrimaryVertex *
                      process.PFTau * #prescribed by TAU POG
                      process.mvaID   *
                      process.inclusiveVertexing *
                      process.softElectronCands*
                      #                      process.kt6PFJets  *
                    
#                      process.patDefaultSequence *
                      process.pfParticleSelectionSequence +
#                      process.eleIsoSequence +
#                      process.phoIsoSequence*
                       process.patDefaultSequence +
                       process.eleIsoSequence*
#                      process.pfElectronIsolationSequence                 *
                      process.pfMuonIsolationSequence                     *
#                      process.pfPostSequence *
#                      process.produceJEC_Uncertainty_FR *
#                      process.electronPFIsoDeposits *
#                      process.electronPFIsoValues *
#                      process.muPFIsoDeposits *
#                      process.muPFIsoValues *
#                      process.patDefaultSequence*
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
switchOnTriggerMatching(process, ['tauTriggerMatchHLTTausLoose','tauTriggerMatchHLTTausMedium','electronTriggerMatchHLTElectronsLoose','muonTriggerMatchHLTMuonsLoose'])
#switchOnTriggerMatching(process, ['tauTriggerMatchHLTTausMedium'])

# Switch to selected PAT objects in the trigger matching
removeCleaningFromTriggerMatching(process)
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger(process)
process.patTrigger.processName = HLTProcessName
process.patTriggerEvent.processName = HLTProcessName
