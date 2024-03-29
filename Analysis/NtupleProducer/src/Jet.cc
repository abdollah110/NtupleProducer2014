#include <PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h>
#include <PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h>

#include "Analysis/NtupleProducer/interface/NtupleProducer.h"
#include "CMGTools/External/interface/PileupJetIdentifier2.h"

void NtupleProducer::DoJetAnalysis(const edm::Event& iEvent) {

    (m->RecPFJetsAK5).clear();
    (m->RhoCorr) = 0;
    (m->RhoCenCharged) = 0;
    (m->RhoCenNeutral) = 0;
    (m->RhoCenNeutralTight) = 0;

    //=================== PF Jet AK5 =================================
    using pat::Jet;
    using pat::JetCollection;
    Handle<pat::JetCollection> PFjetsAK5;
    iEvent.getByLabel(JetCollection_, PFjetsAK5);

    Handle<edm::ValueMap<int> > puJetIdFlag;
    iEvent.getByLabel(puJetIdFlag_,puJetIdFlag);

    edm::Handle<edm::View<pat::Jet> > jets;
    iEvent.getByLabel("selectedPatJets",jets);
    
    edm::Handle<double> hRho;
    iEvent.getByLabel(rhoProducer_, hRho);
    double rho_ = *hRho;

    edm::Handle< pat::TriggerEvent > triggerEvent;
    iEvent.getByLabel(triggerEvent_, triggerEvent);
    // PAT trigger helper for trigger matching information
    const pat::helper::TriggerMatchHelper matchHelper;

    m->Rho = rho_;
   unsigned int i = 0;
   int qq=0;
    for (JetCollection::const_iterator jet = PFjetsAK5->begin(); jet != PFjetsAK5->end(); jet++, qq++) {
        myobject myjet;
        myjet.pt = jet->pt();
        myjet.eta = jet->eta();
        myjet.phi = jet->phi();
        myjet.px = jet->px();
        myjet.py = jet->py();
        myjet.pz = jet->pz();
        myjet.z = jet->vz();

        myjet.E = jet->p();
        myjet.Energy = jet->energy();
        myjet.mass = jet->mass();
        myjet.mt = jet->mt();

        myjet.et = jet->et();
        PFJetIDSelectionFunctor jetIDLoose(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE);
        PFJetIDSelectionFunctor jetIDTight(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT);
        pat::strbitset ret_loose = jetIDLoose.getBitTemplate();
        pat::strbitset ret_tight = jetIDLoose.getBitTemplate();
        ret_loose.set(false);
        ret_tight.set(false);

        myjet.jetId_loose = jetIDLoose(*jet, ret_loose);
        myjet.jetId_tight = jetIDTight(*jet, ret_tight);

        myjet.bDiscriminatiors_CSV = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
        myjet.bDiscriminatiors_JP = jet->bDiscriminator("jetProbabilityBJetTags");
        myjet.bDiscriminatiors_TCHPT = jet->bDiscriminator("trackCountingHighPurBJetTags");

	edm::Ref<pat::JetCollection> jetRef(PFjetsAK5, i);
	int    idflag = (*puJetIdFlag)[jetRef];
	myjet.puJetIdLoose= PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kLoose );
	myjet.puJetIdMedium = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kMedium );
	myjet.puJetIdTight = PileupJetIdentifier2::passJetId( idflag, PileupJetIdentifier2::kTight );

        const pat::TriggerObjectRef trigRef_Ditau30Jet30(matchHelper.triggerMatchObject(PFjetsAK5, qq, jetMatch_Ditau30Jet30_, iEvent, *triggerEvent));
        myjet.hasTrgObject_Ditau30Jet30 = false;
        myjet.TrgObjectEta_Ditau30Jet30 = -100;
        myjet.TrgObjectPt_Ditau30Jet30 = -100;
        myjet.TrgObjectPhi_Ditau30Jet30 = -100;
        if (trigRef_Ditau30Jet30.isAvailable()) { // check references (necessary!)
          myjet.hasTrgObject_Ditau30Jet30 = true;
          myjet.TrgObjectEta_Ditau30Jet30 = trigRef_Ditau30Jet30->eta();
          myjet.TrgObjectPt_Ditau30Jet30 = trigRef_Ditau30Jet30->pt();
          myjet.TrgObjectPhi_Ditau30Jet30 = trigRef_Ditau30Jet30->phi();
        }


    //            (m->RecPFJetsAK5).push_back(myjet);
        //        }
        if (jet->pt() > 10)
            (m->RecPFJetsAK5).push_back(myjet);
	i++;
    }//loop over jets


    edm::Handle<double> rhoHandle;
    iEvent.getByLabel(rhoCorrection_, rhoHandle);
    if (rhoHandle.isValid())
        m->RhoCorr = *rhoHandle;

    edm::Handle<double> rhoHandle1;
    iEvent.getByLabel(rhoCenChargedPU_, rhoHandle1);
    if (rhoHandle1.isValid())
        m->RhoCenCharged = *rhoHandle1;

    edm::Handle<double> rhoHandle2;
    iEvent.getByLabel(rhoCenNeutral_, rhoHandle2);
    if (rhoHandle2.isValid())
        m->RhoCenNeutral = *rhoHandle2;

    edm::Handle<double> rhoHandle3;
    iEvent.getByLabel(rhoCenNeutralTight_, rhoHandle3);
    if (rhoHandle3.isValid())
        m->RhoCenNeutral = *rhoHandle3;
}

//        //patJetCorrFactors
//        std::vector<string> jes_vec = jet->availableJECSets();
//        for (unsigned int i = 0; i < jes_vec.size(); i++) {
//            cout << jes_vec[i] << "\t";
//        }
//patJetCorrFactors
//        std::vector<string> jes_vec = jet->availableJECLevels("patJetCorrFactors");
//        std::vector<string> jes_vec = jet->availableJECLevels();
//        for (unsigned int i = 0; i < jes_vec.size(); i++) {
//            cout << jes_vec[i] << " =  ";
//            cout << jet->correctedP4(jes_vec[i]).pt()<<"\t";
//        }
//            cout<<"corrected= "<<jet->pt()<<endl;
//        cout << endl;
