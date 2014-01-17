#!/bin/csh

#eval `scramv1 runtime -csh`
rootcint -f eventdict.cc -c -I${PWD}/../../.. \
            Analysis/NtupleProducer/interface/myobject.h \
            Analysis/NtupleProducer/interface/myevent.h \
            Analysis/NtupleProducer/interface/LinkDef.h

