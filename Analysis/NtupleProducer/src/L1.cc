#include "Analysis/NtupleProducer/interface/NtupleProducer.h"

void NtupleProducer::DoL1Analysis(const edm::Event& iEvent) {

    (m->L1Tau).clear();
    (m->L1Jet).clear();

    edm::Handle<edm::View<l1extra::L1JetParticle> > l1Jet;
    iEvent.getByLabel("l1extraParticles","Central",l1Jet);

    for (edm::View<l1extra::L1JetParticle>::const_iterator jetl1 = l1Jet->begin(); jetl1 != l1Jet->end(); jetl1++) {
        myobject myl1;
        myl1.pt = jetl1->pt();
        myl1.eta = jetl1->eta();
        myl1.phi = jetl1->phi();
        myl1.et = jetl1->et();
	myl1.px = jetl1->px();
	myl1.py = jetl1->py();
	myl1.pz = jetl1->pz();
        (m->L1Jet).push_back(myl1);
    }//loop over L1 central jets

    edm::Handle<edm::View<l1extra::L1JetParticle> > l1Tau;
    iEvent.getByLabel("l1extraParticles","Tau",l1Tau);

    for (edm::View<l1extra::L1JetParticle>::const_iterator taul1 = l1Tau->begin(); taul1 != l1Tau->end(); taul1++) {
        myobject myl1;
        myl1.pt = taul1->pt();
        myl1.eta = taul1->eta();
        myl1.phi = taul1->phi();
        myl1.et = taul1->et();
        myl1.px = taul1->px();
        myl1.py = taul1->py();
        myl1.pz = taul1->pz();
        (m->L1Tau).push_back(myl1);
    }//loop over L1 taus

}



