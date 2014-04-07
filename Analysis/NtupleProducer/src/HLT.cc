#include <FWCore/Framework/interface/Event.h>
#include "Analysis/NtupleProducer/interface/NtupleProducer.h"

void NtupleProducer::DoHLTAnalysis(const edm::Event& iEvent) {

    (m->L2Particles).clear();
    (m->HLT).clear();
    (m->HLT_DiElectron) = 0;
    (m->HLT_DiMuon) = 0;

    edm::Handle<TriggerResults> triggerResults;
    iEvent.getByLabel(srcTriggerResults_, triggerResults);


//ooooooooooooooooooooooooooooooo
    edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT"); //make sure have correct process on MC
    edm::Handle<trigger::TriggerEvent> trigEvent; 
    iEvent.getByLabel(trigEventTag,trigEvent);

    std::string filterName("hltL2Tau25eta2p1"); 

    trigger::size_type filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())); 
    if(filterIndex<trigEvent->sizeFilters()){ 
    const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
    const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
    for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
      const trigger::TriggerObject& obj = trigObjColl[*keyIt];
      myobject triggo;
      triggo.pt = obj.pt();
      triggo.px = obj.px();
      triggo.py = obj.py();
      triggo.pz = obj.pz();
      triggo.eta = obj.eta();
      triggo.phi = obj.phi();

      (m->L2Particles).push_back(triggo);
    }
    
}//end filter size check

//oooooooooooooooooooooooooooooooooooooooooooooooo
    
    
    string eleTrigger (el_trigger_name);
    string muTrigger (mu_trigger_name);


    if (triggerResults.isValid()) {
        int ntrigs = triggerResults->size();
        TriggerNames const &triggerNames = iEvent.triggerNames(*triggerResults);

        for (int itrig = 0; itrig < ntrigs; itrig++) {
            string name = triggerNames.triggerName(itrig);
            bool result = triggerResults->accept(itrig);
	    size_t foundEl=name.find(eleTrigger);
	    size_t foundMu=name.find(muTrigger);
	    if(filterTriggerResults)
             (m->HLT)[name] = result;
	    else if(!filterTriggerResults)
	      (m->HLT)[name] = result;


            if (name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v1" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v2" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v3" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v4" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v5" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v6" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v7" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v8" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v9" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v10" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v11" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v12" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v13" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v14" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v15" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v16" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v17" or name == "HLT_IsoMu18_eta2p1_MediumIsoPFTau25_Trk5_eta2p1_v18") {
                m->HLT_DiElectron = result;
                //                cout << "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2" << endl;
            }
            if (name == "HLT_DoubleMu7_v1") {
                m->HLT_DiMuon = result;
                //                cout << "HLT_DoubleMu7_v1" << endl;
            }


        }//for itrig
    }//if triggerResults valid



}
