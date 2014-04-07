#include "../../NtupleProducer/interface/NtupleProducer.h"
//#include "../../NtupleProducer/interface/myHelper.h"

void NtupleProducer::DoHPSTauAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    (m->PreSelectedHPSTaus).clear();
    (m->SelectedHPSTaus).clear();
    using namespace std;
    using namespace reco;
    using namespace edm;
    using namespace pat;


   // trigger event
    edm::Handle< pat::TriggerEvent > triggerEvent;
    iEvent.getByLabel(triggerEvent_, triggerEvent);
    // PAT trigger helper for trigger matching information
    const pat::helper::TriggerMatchHelper matchHelper;


    
    // Get the B-field
    ESHandle<MagneticField> B;
    iSetup.get<IdealMagneticFieldRecord > ().get(B);
    const MagneticField* magField = B.product();

    //BeamSpot
    Handle<reco::BeamSpot> theBeamSpotHandle;
    iEvent.getByLabel(std::string("offlineBeamSpot"), theBeamSpotHandle);
    const reco::BeamSpot* beamSpot = theBeamSpotHandle.product();

    // Get the geometry
    ESHandle<GlobalTrackingGeometry> geomHandle;
    iSetup.get<GlobalTrackingGeometryRecord > ().get(geomHandle);




    Handle<pat::TauCollection> tausHandle;
    iEvent.getByLabel(PreSelectedhpsCollection_, tausHandle);
    const TauCollection &tau = *(tausHandle.product());

    pat::TauCollection::const_iterator itau = tau.begin();
    pat::TauCollection::const_iterator jtau = tau.end();


    Handle< std::vector<reco::Vertex> > PrVx;
    iEvent.getByLabel(edm::InputTag("offlinePrimaryVertices"), PrVx);

    Handle< std::vector<reco::Vertex> > PrVx_match;
    iEvent.getByLabel(edm::InputTag("offlinePrimaryVertices"), PrVx_match);


    int ipftau = 0;
    for (; itau != jtau; ++itau, ipftau++) {

        myobject tauu;
        if (PrVx_match->size() > 0)
            tauu.dz_Ver_match = itau->vertex().z() - PrVx_match->front().z();
        else
            tauu.dz_Ver_match = 1000;
        tauu.pt = itau->pt();
        tauu.eta = itau->eta();
        tauu.phi = itau->phi();
        tauu.mass = itau->mass();
        tauu.mt = itau->mt();
        tauu.et = itau->et();
        tauu.px = itau->px();
        tauu.py = itau->py();
        tauu.pz = itau->pz();
        tauu.z = itau->vz();


        tauu.E = itau->p();
        tauu.Energy = itau->energy();
        tauu.emfraction = itau->emFraction();
        tauu.charge = itau->charge();

        tauu.jetPt = itau->pfJetRef().get()->pt();
        tauu.jetEta = itau->pfJetRef().get()->eta();
        tauu.jetPhi = itau->pfJetRef().get()->phi();
        tauu.jetMass = itau->pfJetRef().get()->mass();

        tauu.leadChargedParticlePt = (itau->leadPFChargedHadrCand().isNonnull() ? (itau->leadPFChargedHadrCand())->pt() : 0.);
        tauu.leadTrackD0 = (itau->leadTrack().isNonnull() ? (itau->leadTrack())->dxy(PrVx->front().position()) : 0.);

 //hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
       /*vector<reco::PFCandidatePtr> PFCons = (itau->pfJetRef())->getPFConstituents();
       vector<reco::PFCandidatePtr> PFSignal = itau->signalPFChargedHadrCands();
       int chargedMul_noPt_noDz_signalSub=0;	
       int chargedMul_noPt_Dz_signalSub=0;
       int chargedMul_Pt10_noDz_signalSub=0;
       int chargedMul_Pt10_Dz_signalSub=0;
       int chargedMul_Pt15_noDz_signalSub=0;
       int chargedMul_Pt15_Dz_signalSub=0;
       int chargedMul_noPt_noDz=0;
       int chargedMul_noPt_Dz=0;
       int chargedMul_Pt10_noDz=0;
       int chargedMul_Pt10_Dz=0;
       int chargedMul_Pt15_noDz=0;
       int chargedMul_Pt15_Dz=0;
       bool overlap_signal=false;
       for(vector<reco::PFCandidatePtr>::const_iterator iChargedHadrCand=PFCons.begin(); iChargedHadrCand!=PFCons.end();++iChargedHadrCand) {
		overlap_signal=false;
                for(vector<reco::PFCandidatePtr>::const_iterator iSignalHadrCand=PFSignal.begin(); iSignalHadrCand!=PFSignal.end();++iSignalHadrCand) {
                        if (deltaR((**iChargedHadrCand).p4(),(**iSignalHadrCand).p4()) <0.01){
                                overlap_signal=true;
                        }
                }
		if (fabs((**iChargedHadrCand).vertex().z() - PrVx_match->front().z())<0.2){
			chargedMul_noPt_Dz++;
			if ((**iChargedHadrCand).pt()>1.0)
				chargedMul_Pt10_Dz++;
			if ((**iChargedHadrCand).pt()>1.5)
                                chargedMul_Pt15_Dz++;
			if (overlap_signal==false){
                               	chargedMul_noPt_Dz_signalSub++;
				if ((**iChargedHadrCand).pt()>1.0)
                                	chargedMul_Pt10_Dz_signalSub++;
                        	if ((**iChargedHadrCand).pt()>1.5)
                                	chargedMul_Pt15_Dz_signalSub++;
			}
                }
		chargedMul_noPt_noDz++;
                if ((**iChargedHadrCand).pt()>1.0)
                	chargedMul_Pt10_noDz++;
                if ((**iChargedHadrCand).pt()>1.5)
                        chargedMul_Pt15_noDz++;
		if (overlap_signal==false){
                	chargedMul_noPt_noDz_signalSub++;
                        if ((**iChargedHadrCand).pt()>1.0)
                        	chargedMul_Pt10_noDz_signalSub++;
                        if ((**iChargedHadrCand).pt()>1.5)
                                chargedMul_Pt15_noDz_signalSub++;
               }
       } 
       tauu.chargedMul_noPt_noDz_signalSub=chargedMul_noPt_noDz_signalSub;
       tauu.chargedMul_noPt_noDz=chargedMul_noPt_noDz;
       tauu.chargedMul_Pt10_noDz_signalSub=chargedMul_Pt10_noDz_signalSub;
       tauu.chargedMul_Pt10_noDz=chargedMul_Pt10_noDz;
       tauu.chargedMul_Pt15_noDz_signalSub=chargedMul_Pt15_noDz_signalSub;
       tauu.chargedMul_Pt15_noDz=chargedMul_Pt15_noDz;
       tauu.chargedMul_noPt_Dz_signalSub=chargedMul_noPt_Dz_signalSub;
       tauu.chargedMul_noPt_Dz=chargedMul_noPt_Dz;
       tauu.chargedMul_Pt10_Dz_signalSub=chargedMul_Pt10_Dz_signalSub;
       tauu.chargedMul_Pt10_Dz=chargedMul_Pt10_Dz;
       tauu.chargedMul_Pt15_Dz_signalSub=chargedMul_Pt15_Dz_signalSub;
       tauu.chargedMul_Pt15_Dz=chargedMul_Pt15_Dz;
*/
       //vector<reco::PFCandidatePtr> IsolationHadr = itau->isolationPFChargedHadrCands();
       //for(vector<reco::PFCandidatePtr>::const_iterator iIsolationHadrCand=IsolationHadr.begin(); iIsolationHadrCand!=IsolationHadr.end();++iIsolationHadrCand) {
       //        (tauu.isolationpt).push_back((**iIsolationHadrCand).pt());
       //}
  

        //tauu.chargedMultiplicity = taille;



        tauu.numChargedParticlesSignalCone = itau->signalPFChargedHadrCands().size();
        tauu.numNeutralHadronsSignalCone = itau->signalPFNeutrHadrCands().size();
        tauu.numPhotonsSignalCone = itau->signalPFGammaCands().size();
        tauu.numParticlesSignalCone = itau->signalPFCands().size();
	tauu.signalPiZeroCandidates = itau->signalPiZeroCandidates().size();

        tauu.numChargedParticlesIsoCone = itau->isolationPFChargedHadrCands().size();
        tauu.numNeutralHadronsIsoCone = itau->isolationPFNeutrHadrCands().size();
        tauu.numPhotonsIsoCone = itau->isolationPFGammaCands().size();
        tauu.numParticlesIsoCone = itau->isolationPFCands().size();

        tauu.ptSumChargedParticlesIsoCone = itau->isolationPFChargedHadrCandsPtSum();
        tauu.ptSumPhotonsIsoCone = itau->isolationPFGammaCandsEtSum();

        tauu.mva_e_pi = (itau->leadPFChargedHadrCand().isNonnull() ? itau->leadPFChargedHadrCand()->mva_e_pi() : 0.);
        tauu.mva_pi_mu = (itau->leadPFChargedHadrCand().isNonnull() ? itau->leadPFChargedHadrCand()->mva_pi_mu() : 0.);
        tauu.mva_e_mu = (itau->leadPFChargedHadrCand().isNonnull() ? itau->leadPFChargedHadrCand()->mva_e_mu() : 0.);
        tauu.hcalEnergy = (itau->leadPFChargedHadrCand().isNonnull() ? itau->leadPFChargedHadrCand()->hcalEnergy() : 0.);
        tauu.ecalEnergy = (itau->leadPFChargedHadrCand().isNonnull() ? itau->leadPFChargedHadrCand()->ecalEnergy() : 0.);
        tauu.trackRefPt = (itau->leadPFChargedHadrCand().isNonnull() ? itau->leadPFChargedHadrCand()->pt() : 0.);

        tauu.discriminationByDecayModeFinding = itau->tauID("decayModeFinding") > 0.5 ? true : false;

        tauu.byVLooseCombinedIsolationDeltaBetaCorr = itau->tauID("byVLooseCombinedIsolationDeltaBetaCorr") > 0.5 ? true : false;
        tauu.byLooseCombinedIsolationDeltaBetaCorr = itau->tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 ? true : false;
        tauu.byMediumCombinedIsolationDeltaBetaCorr = itau->tauID("byMediumCombinedIsolationDeltaBetaCorr") > 0.5 ? true : false;
        tauu.byTightCombinedIsolationDeltaBetaCorr = itau->tauID("byTightCombinedIsolationDeltaBetaCorr") > 0.5 ? true : false;

	// 3Hits isolation
        tauu.byLooseCombinedIsolationDeltaBetaCorr3Hits = itau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 ? true : false;
        tauu.byMediumCombinedIsolationDeltaBetaCorr3Hits = itau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") > 0.5 ? true : false;
        tauu.byTightCombinedIsolationDeltaBetaCorr3Hits = itau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits") > 0.5 ? true : false;
	tauu.byRawCombinedIsolationDeltaBetaCorr3Hits = itau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");

	// electron rejection
        tauu.discriminationByElectronLoose = itau->tauID("againstElectronLoose") > 0.5 ? true : false;
        tauu.discriminationByElectronMedium = itau->tauID("againstElectronMedium") > 0.5 ? true : false;
        tauu.discriminationByElectronTight = itau->tauID("againstElectronTight") > 0.5 ? true : false;
        //tauu.discriminationByElectronMVA = itau->tauID("againstElectronMVA") > 0.5 ? true : false;

        tauu.discriminationByElectronMVA5Loose = itau->tauID("againstElectronLooseMVA5") > 0.5 ? true : false;
        tauu.discriminationByElectronMVA5Medium = itau->tauID("againstElectronMediumMVA5") > 0.5 ? true : false;
        tauu.discriminationByElectronMVA5Tight = itau->tauID("againstElectronTightMVA5") > 0.5 ? true : false;
	tauu.discriminationByElectronMVA5VTight = itau->tauID("againstElectronVTightMVA5") > 0.5 ? true : false;
	tauu.discriminationByElectronMVA5VLoose = itau->tauID("againstElectronVLooseMVA5") > 0.5 ? true : false;

	// muon rejection
        tauu.discriminationByMuonLoose = itau->tauID("againstMuonLoose") > 0.5 ? true : false;
        tauu.discriminationByMuonMedium = itau->tauID("againstMuonMedium") > 0.5 ? true : false;
        tauu.discriminationByMuonTight = itau->tauID("againstMuonTight") > 0.5 ? true : false;

	tauu.discriminationByMuonLoose2 = itau->tauID("againstMuonLoose2") > 0.5 ? true : false;
        tauu.discriminationByMuonMedium2 = itau->tauID("againstMuonMedium2") > 0.5 ? true : false;
	tauu.discriminationByMuonTight2 = itau->tauID("againstMuonTight2") > 0.5 ? true : false;

	tauu.discriminationByMuonLoose3 = itau->tauID("againstMuonLoose3") > 0.5 ? true : false;
	tauu.discriminationByMuonTight3 = itau->tauID("againstMuonTight3") > 0.5 ? true : false;

	tauu.discriminationByMuonMVALoose = itau->tauID("againstMuonLooseMVA") > 0.5 ? true : false;
        tauu.discriminationByMuonMVAMedium = itau->tauID("againstMuonMediumMVA") > 0.5 ? true : false;
        tauu.discriminationByMuonMVATight = itau->tauID("againstMuonTightMVA") > 0.5 ? true : false;

        tauu.discriminationByMuonMVAraw = itau->tauID("againstMuonMVAraw");


        //MVA isolation
        //tauu.byIsolationMVAraw = itau->tauID("byIsolationMVAraw");
        //tauu.byLooseIsolationMVA = itau->tauID("byLooseIsolationMVA") > 0.5 ? true : false;
        //tauu.byMediumIsolationMVA = itau->tauID("byMediumIsolationMVA") > 0.5 ? true : false;
        //tauu.byTightIsolationMVA = itau->tauID("byTightIsolationMVA") > 0.5 ? true : false;
	
 	//tauu.byIsolationMVA2raw = itau->tauID("byIsolationMVA2raw");
        //tauu.byLooseIsolationMVA2 = itau->tauID("byLooseIsolationMVA2") > 0.5 ? true : false;
        //tauu.byMediumIsolationMVA2 = itau->tauID("byMediumIsolationMVA2") > 0.5 ? true : false;
        //tauu.byTightIsolationMVA2 = itau->tauID("byTightIsolationMVA2") > 0.5 ? true : false;

        tauu.byIsolationMVA3oldDMwLTraw = itau->tauID("byIsolationMVA3oldDMwLTraw");
        tauu.byVLooseIsolationMVA3oldDMwLT = itau->tauID("byVLooseIsolationMVA3oldDMwLT") > 0.5 ? true : false;
        tauu.byLooseIsolationMVA3oldDMwLT = itau->tauID("byLooseIsolationMVA3oldDMwLT") > 0.5 ? true : false;
        tauu.byMediumIsolationMVA3oldDMwLT = itau->tauID("byMediumIsolationMVA3oldDMwLT") > 0.5 ? true : false;
        tauu.byTightIsolationMVA3oldDMwLT = itau->tauID("byTightIsolationMVA3oldDMwLT") > 0.5 ? true : false;
        tauu.byVTightIsolationMVA3oldDMwLT = itau->tauID("byVTightIsolationMVA3oldDMwLT") > 0.5 ? true : false;
        tauu.byVVTightIsolationMVA3oldDMwLT = itau->tauID("byVVTightIsolationMVA3oldDMwLT") > 0.5 ? true : false;

        tauu.byIsolationMVA3oldDMwoLTraw = itau->tauID("byIsolationMVA3oldDMwoLTraw");
        tauu.byVLooseIsolationMVA3oldDMwoLT = itau->tauID("byVLooseIsolationMVA3oldDMwoLT") > 0.5 ? true : false;
        tauu.byLooseIsolationMVA3oldDMwoLT = itau->tauID("byLooseIsolationMVA3oldDMwoLT") > 0.5 ? true : false;
        tauu.byMediumIsolationMVA3oldDMwoLT = itau->tauID("byMediumIsolationMVA3oldDMwoLT") > 0.5 ? true : false;
        tauu.byTightIsolationMVA3oldDMwoLT = itau->tauID("byTightIsolationMVA3oldDMwoLT") > 0.5 ? true : false;
        tauu.byVTightIsolationMVA3oldDMwoLT = itau->tauID("byVTightIsolationMVA3oldDMwoLT") > 0.5 ? true : false;
        tauu.byVVTightIsolationMVA3oldDMwoLT = itau->tauID("byVVTightIsolationMVA3oldDMwoLT") > 0.5 ? true : false;

        tauu.byIsolationMVA3newDMwLTraw = itau->tauID("byIsolationMVA3newDMwLTraw");
        tauu.byVLooseIsolationMVA3newDMwLT = itau->tauID("byVLooseIsolationMVA3newDMwLT") > 0.5 ? true : false;
        tauu.byLooseIsolationMVA3newDMwLT = itau->tauID("byLooseIsolationMVA3newDMwLT") > 0.5 ? true : false;
        tauu.byMediumIsolationMVA3newDMwLT = itau->tauID("byMediumIsolationMVA3newDMwLT") > 0.5 ? true : false;
        tauu.byTightIsolationMVA3newDMwLT = itau->tauID("byTightIsolationMVA3newDMwLT") > 0.5 ? true : false;
        tauu.byVTightIsolationMVA3newDMwLT = itau->tauID("byVTightIsolationMVA3newDMwLT") > 0.5 ? true : false;
        tauu.byVVTightIsolationMVA3newDMwLT = itau->tauID("byVVTightIsolationMVA3newDMwLT") > 0.5 ? true : false;

        tauu.byIsolationMVA3newDMwoLTraw = itau->tauID("byIsolationMVA3newDMwoLTraw");
        tauu.byVLooseIsolationMVA3newDMwoLT = itau->tauID("byVLooseIsolationMVA3newDMwoLT") > 0.5 ? true : false;
        tauu.byLooseIsolationMVA3newDMwoLT = itau->tauID("byLooseIsolationMVA3newDMwoLT") > 0.5 ? true : false;
        tauu.byMediumIsolationMVA3newDMwoLT = itau->tauID("byMediumIsolationMVA3newDMwoLT") > 0.5 ? true : false;
        tauu.byTightIsolationMVA3newDMwoLT = itau->tauID("byTightIsolationMVA3newDMwoLT") > 0.5 ? true : false;
        tauu.byVTightIsolationMVA3newDMwoLT = itau->tauID("byVTightIsolationMVA3newDMwoLT") > 0.5 ? true : false;
        tauu.byVVTightIsolationMVA3newDMwoLT = itau->tauID("byVVTightIsolationMVA3newDMwoLT") > 0.5 ? true : false;


        tauu.z_expo = 0;
        if (itau->leadPFChargedHadrCand().isNonnull()) {

            if (itau->leadPFChargedHadrCand()->trackRef().isNonnull()) {
                reco::TransientTrack track(itau->leadPFChargedHadrCand()->trackRef(), magField, geomHandle);
                //FreeTrajectoryState state = track.impactPointTSCP().theState();
                //	states.push_back(state);
                //calculate Z at beamspot
                TransverseImpactPointExtrapolator extrapolator(magField);
                TrajectoryStateOnSurface closestOnTransversePlaneState = extrapolator.extrapolate(track.impactPointState(), GlobalPoint(beamSpot->position().x(), beamSpot->position().y(), 0.0));
                tauu.z_expo = closestOnTransversePlaneState.globalPosition().z();
            }

        }



                //TriggerObjectmatching
        const pat::TriggerObjectRef trigRef_loose(matchHelper.triggerMatchObject(tausHandle, ipftau, tauMatch_Loose_, iEvent, *triggerEvent));
        tauu.hasTrgObject_loose = false;
        tauu.TrgObjectEta_loose = -100;
        tauu.TrgObjectPt_loose = -100;
        tauu.TrgObjectPhi_loose = -100;
        // finally we can fill the histograms
        if (trigRef_loose.isAvailable()) { // check references (necessary!)

            tauu.hasTrgObject_loose = true;
            tauu.TrgObjectEta_loose = trigRef_loose->eta();
            tauu.TrgObjectPt_loose = trigRef_loose->pt();
            tauu.TrgObjectPhi_loose = trigRef_loose->phi();
        }

        const pat::TriggerObjectRef trigRef_medium(matchHelper.triggerMatchObject(tausHandle, ipftau, tauMatch_Medium_, iEvent, *triggerEvent));
        tauu.hasTrgObject_medium = false;
        tauu.TrgObjectEta_medium = -100;
        tauu.TrgObjectPt_medium = -100;
        tauu.TrgObjectPhi_medium = -100;
        // finally we can fill the histograms
        if (trigRef_medium.isAvailable()) { // check references (necessary!)
             tauu.hasTrgObject_medium = true;
             tauu.TrgObjectEta_medium = trigRef_medium->eta();
             tauu.TrgObjectPt_medium = trigRef_medium->pt();
             tauu.TrgObjectPhi_medium = trigRef_medium->phi();
        }


        //        if (itau->pt() > 10)
        (m->PreSelectedHPSTaus).push_back(tauu);

	//filling selected collection
	if(itau->pt() > tauPtcut_ &&  itau->tauID("decayModeFinding") > 0.5 &&  itau->tauID("againstMuonLoose") > 0.5 &&  itau->tauID("againstElectronLoose") > 0.5)
	  {
	    myobject seltau = tauu;
	    seltau.gen_index=ipftau;
	    for(unsigned int iCharged=0; iCharged < itau->signalPFChargedHadrCands().size(); iCharged++)
	      {
		const reco::PFCandidatePtr& cand = itau->signalPFChargedHadrCands().at(iCharged);
		math::XYZTLorentzVector candP4 = cand->p4();
		switch(iCharged)
		  {
		  case 0:
		    seltau.sig_track1_pt = candP4.pt();
		    seltau.sig_track1_eta = candP4.eta();
		    seltau.sig_track1_phi = candP4.phi();
		    seltau.sig_track1_m = candP4.M();
		    break;
		  case 1:
		    seltau.sig_track2_pt = candP4.pt();
                    seltau.sig_track2_eta = candP4.eta();
                    seltau.sig_track2_phi = candP4.phi();
                    seltau.sig_track2_m = candP4.M();
                    break;
		  case 2:
		    seltau.sig_track3_pt = candP4.pt();
                    seltau.sig_track3_eta = candP4.eta();
                    seltau.sig_track3_phi = candP4.phi();
                    seltau.sig_track3_m = candP4.M();
                    break;
		  default:
		    break;
		  }
	      }
	    for(unsigned int iPi0=0; iPi0 < itau->signalPiZeroCandidates().size(); iPi0++)
	      {
		const reco::RecoTauPiZero& cand = itau->signalPiZeroCandidates().at(iPi0);
		math::XYZTLorentzVector candP4 = cand.p4();
		switch(iPi0)
                  {
                  case 0:
                    seltau.sig_pi0_1_pt = candP4.pt();
                    seltau.sig_pi0_1_eta = candP4.eta();
                    seltau.sig_pi0_1_phi = candP4.phi();
                    seltau.sig_pi0_1_m = candP4.M();
                    break;
		  case 1:
		    seltau.sig_pi0_2_pt = candP4.pt();
                    seltau.sig_pi0_2_eta = candP4.eta();
                    seltau.sig_pi0_2_phi = candP4.phi();
                    seltau.sig_pi0_2_m = candP4.M();
		    break;
		  default:
		    break;
		  }

	      }
	    (m->SelectedHPSTaus).push_back(seltau);
	  }

    }

}
