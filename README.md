NtupleProducer2014
==================
#Directory with all packages
mkdir NtupleProducer_help
cd NtupleProducer_help/
cmsrel CMSSW_5_3_14_patch2
cd CMSSW_5_3_14_patch2/src/
cmsenv

#MVA MET
git cms-addpkg PhysicsTools/PatAlgos
git cms-merge-topic cms-analysis-tools:5_3_14-updateSelectorUtils
git cms-merge-topic cms-analysis-tools:5_3_13_patch2-testNewTau
git cms-merge-topic -u TaiSakuma:53X-met-131120-01
git-cms-merge-topic -u cms-met:53X-MVaNoPuMET-20131217-01

#EGamma tools
cvs co -r V00-00-09 EgammaAnalysis/ElectronTools
#cvs co -r jakob19April2013_2012ID EgammaAnalysis/ElectronTools
cvs co -r V09-00-01 RecoEgamma/EgammaTools
cd EgammaAnalysis/ElectronTools/data
cat download.url | xargs wget
cd ../../..

#CMG tools
wget --no-check-certificate https://jez.web.cern.ch/jez/CMGTools.tgz
tar xzvf CMGTools.tgz
cd ../../..

#Remove some conflicting packages
rm -rf DataFormats/TauReco
rm -rf RecoTauTag/RecoTau
rm -rf RecoTauTag/Configuration
rm -rf RecoTauTag/ImpactParameter
rm -rf RecoTauTag/TauTagTools
rm -rf PhysicsTools/PatAlgos
rm -rf DataFormats/PatCandidates

#Directory with tau ID packages
cd ../../..
mkdir NtupleProducer_main
cd NtupleProducer_main
cmsrel CMSSW_5_3_14_patch2
cd CMSSW_5_3_14_patch2/src/
cmsenv

#New tau ID
git cms-merge-topic -u cms-tau-pog:CMSSW_5_3_X_boostedTaus
cp -rf ../../NtupleProducer_help/CMSSW_5_3_14_patch2/src/CMGTools
cp -rf ../../NtupleProducer_help/CMSSW_5_3_14_patch2/src/CMGTools .
cp -rf ../../NtupleProducer_help/CMSSW_5_3_14_patch2/src/EgammaAnalysis/ .
cp -rf ../../NtupleProducer_help/CMSSW_5_3_14_patch2/src/GeneratorInterface/ .
cp -rf ../../NtupleProducer_help/CMSSW_5_3_14_patch2/src/RecoBTag/ .
cp -rf ../../NtupleProducer_help/CMSSW_5_3_14_patch2/src/RecoEgamma/ .
cp -rf ../../NtupleProducer_help/CMSSW_5_3_14_patch2/src/RecoMET/ .
cp -rf ../../NtupleProducer_help/CMSSW_5_3_14_patch2/src/DataFormats/ .
cp -rf ../../NtupleProducer_help/CMSSW_5_3_14_patch2/src/RecoJets/ .
cp -rf ../../NtupleProducer_help/CMSSW_5_3_14_patch2/src/PhysicsTools/ .
cp -rf ../../NtupleProducer_help/CMSSW_5_3_14_patch2/src/JetMETCorrections/ .

#Correct TauMETAlgo.cc
#vi JetMETCorrections/Type1MET/src/TauMETAlgo.cc 
#Type :%s/PFCandidateRefVector/vector<reco::PFCandidatePtr>

#Correct Ntuple producer
#vi Analysis/NtupleProducer/src/HPSTauAnalysis.cc 
