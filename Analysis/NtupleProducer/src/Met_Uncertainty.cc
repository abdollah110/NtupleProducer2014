#include <FWCore/Framework/interface/Event.h>
//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "Analysis/NtupleProducer/interface/NtupleProducer.h"
#include "CMGTools/External/interface/PileupJetIdentifier2.h"

void NtupleProducer::DoMetUncertaityAnalysis(const edm::Event& iEvent) {
    (m->RecPFMetCorElectronEnUp).clear();
    (m->RecPFMetCorElectronEnDown).clear();
    (m->RecPFMetCorMuonEnUp).clear();
    (m->RecPFMetCorMuonEnDown).clear();
    (m->RecPFMetCorTauEnUp).clear();
    (m->RecPFMetCorTauEnDown).clear();
    (m->RecPFMetCorJetEnUp).clear();
    (m->RecPFMetCorJetEnDown).clear();
    (m->RecPFMetCorJetResUp).clear();
    (m->RecPFMetCorJetResDown).clear();
    (m->RecPFMetCorUnclusteredEnUp).clear();
    (m->RecPFMetCorUnclusteredEnDown).clear();

    (m->smearedPatJets).clear();
    (m->smearedPatJetsResUp).clear();
    (m->smearedPatJetsResDown).clear();
    (m->shiftedPatJetsEnUpForCorrMEt).clear();
    (m->shiftedPatJetsEnDownForCorrMEt).clear();
    (m->cleanPatJets).clear();

    //============== Type1CorPFMET ==================================



    typedef edm::View<reco::MET> METView;

    //    Handle<METView> Type1CorMETElectronEnUphandle;
    //    iEvent.getByLabel("patType1CorrectedPFMetElectronEnUp", Type1CorMETElectronEnUphandle);
    //    for (METView::const_iterator met_iter = Type1CorMETElectronEnUphandle->begin();
    //            met_iter != Type1CorMETElectronEnUphandle->end(); met_iter++) {
    //        myobject mymet;
    //        mymet.pt = met_iter->pt();
    //        mymet.px = met_iter->px();
    //        mymet.py = met_iter->py();
    //        mymet.phi = met_iter->phi();
    //        mymet.et = met_iter->et();
    //        (m->RecPFMetCorElectronEnUp).push_back(mymet);
    //    }
    //
    //
    //    Handle<METView> Type1CorMETElectronEnDownhandle;
    //    iEvent.getByLabel("patType1CorrectedPFMetElectronEnDown", Type1CorMETElectronEnDownhandle);
    //    for (METView::const_iterator met_iter = Type1CorMETElectronEnDownhandle->begin();
    //            met_iter != Type1CorMETElectronEnDownhandle->end(); met_iter++) {
    //        myobject mymet;
    //        mymet.pt = met_iter->pt();
    //        mymet.px = met_iter->px();
    //        mymet.py = met_iter->py();
    //        mymet.phi = met_iter->phi();
    //        mymet.et = met_iter->et();
    //        (m->RecPFMetCorElectronEnDown).push_back(mymet);
    //    }
    //
    //
    //
    //    Handle<METView> Type1CorMETMuonEnUphandle;
    //    iEvent.getByLabel("patType1CorrectedPFMetMuonEnUp", Type1CorMETMuonEnUphandle);
    //    for (METView::const_iterator met_iter = Type1CorMETMuonEnUphandle->begin();
    //            met_iter != Type1CorMETMuonEnUphandle->end(); met_iter++) {
    //        myobject mymet;
    //        mymet.pt = met_iter->pt();
    //        mymet.px = met_iter->px();
    //        mymet.py = met_iter->py();
    //        mymet.phi = met_iter->phi();
    //        mymet.et = met_iter->et();
    //        (m->RecPFMetCorMuonEnUp).push_back(mymet);
    //    }
    //
    //    Handle<METView> Type1CorMETMuonEnDownhandle;
    //    iEvent.getByLabel("patType1CorrectedPFMetMuonEnDown", Type1CorMETMuonEnDownhandle);
    //    for (METView::const_iterator met_iter = Type1CorMETMuonEnDownhandle->begin();
    //            met_iter != Type1CorMETMuonEnDownhandle->end(); met_iter++) {
    //        myobject mymet;
    //        mymet.pt = met_iter->pt();
    //        mymet.px = met_iter->px();
    //        mymet.py = met_iter->py();
    //        mymet.phi = met_iter->phi();
    //        mymet.et = met_iter->et();
    //        (m->RecPFMetCorMuonEnDown).push_back(mymet);
    //    }
    Handle<METView> Type1CorMETTauEnUphandle;
    iEvent.getByLabel("patType1CorrectedPFMetTauEnUp", Type1CorMETTauEnUphandle);
    for (METView::const_iterator met_iter = Type1CorMETTauEnUphandle->begin();
            met_iter != Type1CorMETTauEnUphandle->end(); met_iter++) {
        myobject mymet;
        mymet.pt = met_iter->pt();
        mymet.px = met_iter->px();
        mymet.py = met_iter->py();
        mymet.phi = met_iter->phi();
        mymet.et = met_iter->et();
        (m->RecPFMetCorTauEnUp).push_back(mymet);
    }

    Handle<METView> Type1CorMETTauEnDownhandle;
    iEvent.getByLabel("patType1CorrectedPFMetTauEnDown", Type1CorMETTauEnDownhandle);
    for (METView::const_iterator met_iter = Type1CorMETTauEnDownhandle->begin();
            met_iter != Type1CorMETTauEnDownhandle->end(); met_iter++) {
        myobject mymet;
        mymet.pt = met_iter->pt();
        mymet.px = met_iter->px();
        mymet.py = met_iter->py();
        mymet.phi = met_iter->phi();
        mymet.et = met_iter->et();
        (m->RecPFMetCorTauEnDown).push_back(mymet);
    }
    // JetResoultion is just done for MC
    if (IsMC) {
        Handle<METView> Type1CorMETJetResUphandle;
        iEvent.getByLabel("patType1CorrectedPFMetJetResUp", Type1CorMETJetResUphandle);
        for (METView::const_iterator met_iter = Type1CorMETJetResUphandle->begin();
                met_iter != Type1CorMETJetResUphandle->end(); met_iter++) {
            myobject mymet;
            mymet.pt = met_iter->pt();
            mymet.px = met_iter->px();
            mymet.py = met_iter->py();
            mymet.phi = met_iter->phi();
            mymet.et = met_iter->et();
            (m->RecPFMetCorJetResUp).push_back(mymet);
        }

        Handle<METView> Type1CorMETJetResDownhandle;
        iEvent.getByLabel("patType1CorrectedPFMetJetResDown", Type1CorMETJetResDownhandle);
        for (METView::const_iterator met_iter = Type1CorMETJetResDownhandle->begin();
                met_iter != Type1CorMETJetResDownhandle->end(); met_iter++) {
            myobject mymet;
            mymet.pt = met_iter->pt();
            mymet.px = met_iter->px();
            mymet.py = met_iter->py();
            mymet.phi = met_iter->phi();
            mymet.et = met_iter->et();
            (m->RecPFMetCorJetResDown).push_back(mymet);
        }
    }

    Handle<METView> Type1CorMETJetEnUphandle;
    iEvent.getByLabel("patType1CorrectedPFMetJetEnUp", Type1CorMETJetEnUphandle);
    for (METView::const_iterator met_iter = Type1CorMETJetEnUphandle->begin();
            met_iter != Type1CorMETJetEnUphandle->end(); met_iter++) {
        myobject mymet;
        mymet.pt = met_iter->pt();
        mymet.px = met_iter->px();
        mymet.py = met_iter->py();
        mymet.phi = met_iter->phi();
        mymet.et = met_iter->et();
        (m->RecPFMetCorJetEnUp).push_back(mymet);
    }

    Handle<METView> Type1CorMETJetEnDownhandle;
    iEvent.getByLabel("patType1CorrectedPFMetJetEnDown", Type1CorMETJetEnDownhandle);
    for (METView::const_iterator met_iter = Type1CorMETJetEnDownhandle->begin();
            met_iter != Type1CorMETJetEnDownhandle->end(); met_iter++) {
        myobject mymet;
        mymet.pt = met_iter->pt();
        mymet.px = met_iter->px();
        mymet.py = met_iter->py();
        mymet.phi = met_iter->phi();
        mymet.et = met_iter->et();
        (m->RecPFMetCorJetEnDown).push_back(mymet);
    }

    Handle<METView> Type1CorMETUnclusteredEnUphandle;
    iEvent.getByLabel("patType1CorrectedPFMetUnclusteredEnUp", Type1CorMETUnclusteredEnUphandle);
    for (METView::const_iterator met_iter = Type1CorMETUnclusteredEnUphandle->begin();
            met_iter != Type1CorMETUnclusteredEnUphandle->end(); met_iter++) {
        myobject mymet;
        mymet.pt = met_iter->pt();
        mymet.px = met_iter->px();
        mymet.py = met_iter->py();
        mymet.phi = met_iter->phi();
        mymet.et = met_iter->et();
        (m->RecPFMetCorUnclusteredEnUp).push_back(mymet);
    }

    Handle<METView> Type1CorMETUnclusteredEnDownhandle;
    iEvent.getByLabel("patType1CorrectedPFMetUnclusteredEnDown", Type1CorMETUnclusteredEnDownhandle);
    for (METView::const_iterator met_iter = Type1CorMETUnclusteredEnDownhandle->begin();
            met_iter != Type1CorMETUnclusteredEnDownhandle->end(); met_iter++) {
        myobject mymet;
        mymet.pt = met_iter->pt();
        mymet.px = met_iter->px();
        mymet.py = met_iter->py();
        mymet.phi = met_iter->phi();
        mymet.et = met_iter->et();
        (m->RecPFMetCorUnclusteredEnDown).push_back(mymet);
    }






    //several type of the JetCollections

    //    using reco::PFJets;
    // using pat::JetCollection;
    //edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    //get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl); 
    //JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
    //JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JetCorPar);

   Handle<edm::ValueMap<int> > puJetIdFlag;
   iEvent.getByLabel(puJetIdFlag_,puJetIdFlag);

    typedef edm::View<pat::Jet> JetView;
    unsigned int i=0;
    using pat::JetCollection;
    Handle<JetView> PFjetsAK5;
    iEvent.getByLabel("patJetsNotOverlappingWithLeptonsForMEtUncertainty", PFjetsAK5);
    for (JetView::const_iterator jet = PFjetsAK5->begin(); jet != PFjetsAK5->end(); i++,jet++) {
        myobject myjet;
        myjet.pt = jet->pt();
        myjet.eta = jet->eta();
        myjet.phi = jet->phi();
        myjet.et = jet->et();
        myjet.bDiscriminatiors_CSV = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
        myjet.bDiscriminatiors_JP = jet->bDiscriminator("jetProbabilityBJetTags");
        myjet.bDiscriminatiors_TCHPT = jet->bDiscriminator("trackCountingHighPurBJetTags");

        //int    idflag = (*puJetIdFlag)[PFjetsAK5->refAt(i)];
        //myjet.puJetIdLoose= PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kLoose );
        //myjet.puJetIdMedium = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kMedium );
        //myjet.puJetIdTight = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kTight );

//        jecUnc->setJetEta(jet->eta());
//        jecUnc->setJetPt(jet->pt()); // here you must use the CORRECTED jet pt
//        myjet.JEC = jecUnc->getUncertainty(true);

        if (jet->pt() > 10) {
            (m->cleanPatJets).push_back(myjet);
        }
    }//loop over jets


    // JetResoultion is just done for MC
    if (IsMC) {
        Handle<JetView> smearedJetHandle;
        iEvent.getByLabel("smearedPatJets", smearedJetHandle);
        i=0;
        for (JetView::const_iterator jet = smearedJetHandle->begin(); jet != smearedJetHandle->end(); i++,jet++) {
            myobject myjet;
            myjet.pt = jet->pt();
            myjet.eta = jet->eta();
            myjet.phi = jet->phi();
            myjet.et = jet->et();

            myjet.bDiscriminatiors_CSV = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
            myjet.bDiscriminatiors_JP = jet->bDiscriminator("jetProbabilityBJetTags");
            myjet.bDiscriminatiors_TCHPT = jet->bDiscriminator("trackCountingHighPurBJetTags");

            //int    idflag = (*puJetIdFlag)[smearedJetHandle->refAt(i)];
            //myjet.puJetIdLoose= PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kLoose );
            //myjet.puJetIdMedium = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kMedium );
            //myjet.puJetIdTight = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kTight );

            if (jet->pt() > 10) {
                (m->smearedPatJets).push_back(myjet);
            }
        }

        Handle<JetView> smearedJetResUpHandle;
        iEvent.getByLabel("smearedPatJetsResUp", smearedJetResUpHandle);
        i=0;
        for (JetView::const_iterator jet = smearedJetResUpHandle->begin(); jet != smearedJetResUpHandle->end(); i++,jet++) {
            myobject myjet;
            myjet.pt = jet->pt();
            myjet.eta = jet->eta();
            myjet.phi = jet->phi();
            myjet.et = jet->et();

            myjet.bDiscriminatiors_CSV = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
            myjet.bDiscriminatiors_JP = jet->bDiscriminator("jetProbabilityBJetTags");
            myjet.bDiscriminatiors_TCHPT = jet->bDiscriminator("trackCountingHighPurBJetTags");

            //int    idflag = (*puJetIdFlag)[smearedJetResUpHandle->refAt(i)];
            //myjet.puJetIdLoose= PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kLoose );
            //myjet.puJetIdMedium = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kMedium );
            //myjet.puJetIdTight = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kTight );

            if (jet->pt() > 10) {
                (m->smearedPatJetsResUp).push_back(myjet);
            }
        }


        Handle<JetView> smearedJetResDownHandle;
        iEvent.getByLabel("smearedPatJetsResDown", smearedJetResDownHandle);
        i=0;
        for (JetView::const_iterator jet = smearedJetResDownHandle->begin(); jet != smearedJetResDownHandle->end(); i++,jet++) {
            myobject myjet;
            myjet.pt = jet->pt();
            myjet.eta = jet->eta();
            myjet.phi = jet->phi();
            myjet.et = jet->et();

            myjet.bDiscriminatiors_CSV = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
            myjet.bDiscriminatiors_JP = jet->bDiscriminator("jetProbabilityBJetTags");
            myjet.bDiscriminatiors_TCHPT = jet->bDiscriminator("trackCountingHighPurBJetTags");

            //int    idflag = (*puJetIdFlag)[smearedJetResDownHandle->refAt(i)];
            //myjet.puJetIdLoose= PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kLoose );
            //myjet.puJetIdMedium = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kMedium );
            //myjet.puJetIdTight = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kTight );

            if (jet->pt() > 10) {
                (m->smearedPatJetsResDown).push_back(myjet);
            }
        }
    }
    Handle<JetView> shiftedPatJetsEnUpForCorrMEtHandle;
    iEvent.getByLabel("shiftedPatJetsEnUpForCorrMEt", shiftedPatJetsEnUpForCorrMEtHandle);
    i=0;
    for (JetView::const_iterator jet = shiftedPatJetsEnUpForCorrMEtHandle->begin(); jet != shiftedPatJetsEnUpForCorrMEtHandle->end(); ++i,jet++) {
        myobject myjet;
        myjet.pt = jet->pt();
        myjet.eta = jet->eta();
        myjet.phi = jet->phi();
        myjet.et = jet->et();

        myjet.bDiscriminatiors_CSV = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
        myjet.bDiscriminatiors_JP = jet->bDiscriminator("jetProbabilityBJetTags");
        myjet.bDiscriminatiors_TCHPT = jet->bDiscriminator("trackCountingHighPurBJetTags");

        //int    idflag = (*puJetIdFlag)[shiftedPatJetsEnUpForCorrMEtHandle->refAt(i)];
        //myjet.puJetIdLoose= PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kLoose );
        //myjet.puJetIdMedium = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kMedium );
        //myjet.puJetIdTight = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kTight );

        if (jet->pt() > 10) {
            (m->shiftedPatJetsEnUpForCorrMEt).push_back(myjet);
        }
    }

    Handle<JetView> shiftedPatJetsEnDownForCorrMEtHandle;
    iEvent.getByLabel("shiftedPatJetsEnDownForCorrMEt", shiftedPatJetsEnDownForCorrMEtHandle);
    i=0;
    for (JetView::const_iterator jet = shiftedPatJetsEnDownForCorrMEtHandle->begin(); jet != shiftedPatJetsEnDownForCorrMEtHandle->end(); i++,jet++) {
        myobject myjet;
        myjet.pt = jet->pt();
        myjet.eta = jet->eta();
        myjet.phi = jet->phi();
        myjet.et = jet->et();

        myjet.bDiscriminatiors_CSV = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
        myjet.bDiscriminatiors_JP = jet->bDiscriminator("jetProbabilityBJetTags");
        myjet.bDiscriminatiors_TCHPT = jet->bDiscriminator("trackCountingHighPurBJetTags");

        //int    idflag = (*puJetIdFlag)[shiftedPatJetsEnDownForCorrMEtHandle->refAt(i)];
        //myjet.puJetIdLoose= PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kLoose );
        //myjet.puJetIdMedium = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kMedium );
        //myjet.puJetIdTight = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kTight );

        if (jet->pt() > 10) {
            (m->shiftedPatJetsEnDownForCorrMEt).push_back(myjet);
        }
    }
}
