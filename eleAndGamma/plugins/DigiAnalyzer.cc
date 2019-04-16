#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TSystem.h>
#include "TTree.h"
#include "TLorentzVector.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDigi/interface/EBDataFrame.h"  
#include "DataFormats/EcalDigi/interface/EEDataFrame.h"  
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h" 
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

using namespace std;
using namespace edm;
using namespace reco;

class DigiAnalyzer : public edm::EDAnalyzer {
public:
  explicit DigiAnalyzer(const edm::ParameterSet&);
  ~DigiAnalyzer();
  
private:

  edm::Service<TFileService> fs_;
  
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // collections
  edm::EDGetTokenT<EBDigiCollection> ecalDigiEBToken_;
  edm::EDGetTokenT<EEDigiCollection> ecalDigiEEToken_;
  edm::EDGetTokenT<GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<vector< PileupSummaryInfo > > PileUpToken_; 
};

DigiAnalyzer::DigiAnalyzer(const edm::ParameterSet& iConfig) {

  ecalDigiEBToken_   = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("EBdigiCollection"));
  ecalDigiEEToken_   = consumes<EEDigiCollection>(iConfig.getParameter<edm::InputTag>("EEdigiCollection"));
  genParticlesToken_ = consumes<GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  PileUpToken_       = consumes<vector< PileupSummaryInfo > >(iConfig.getParameter<edm::InputTag>("PileUpTag"));
}

DigiAnalyzer::~DigiAnalyzer() { } 

void DigiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // access edm objects
  Handle< GenParticleCollection > genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);

  Handle< EBDigiCollection > EcalBarrelDigis;
  iEvent.getByToken(ecalDigiEBToken_, EcalBarrelDigis); 

  Handle< EEDigiCollection > EcalEndcapDigis;
  iEvent.getByToken(ecalDigiEEToken_, EcalEndcapDigis); 

  Handle<vector< PileupSummaryInfo> > PileupInfos;
  iEvent.getByToken(PileUpToken_,PileupInfos);

  cout << endl;

  // PU-info
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  float pu_n = 0.;
  int pu_obs = 0.;
  for(PVI = PileupInfos->begin(); PVI != PileupInfos->end(); ++PVI) {
    int pu_bunchcrossing = PVI->getBunchCrossing();
    if( pu_bunchcrossing == 0 ) {
      pu_n   = PVI->getTrueNumInteractions();
      pu_obs = PVI->getPU_NumInteractions(); 
      continue;
    }
  }

  cout << "Run = " << iEvent.id().run() << ", event = " << iEvent.id().event() << ", pu_true = " << pu_n << ", pu_obs = " << pu_obs << endl;
  cout << "All BXs: " << endl;
  for(PVI = PileupInfos->begin(); PVI != PileupInfos->end(); ++PVI) {
    cout << PVI->getBunchCrossing() << ", true = " << PVI->getTrueNumInteractions() << ", obs = " << PVI->getPU_NumInteractions() << endl;
  }  

  // EB
  for (unsigned int digis=0; digis<EcalBarrelDigis->size(); ++digis) {
    EBDataFrame ebdf = (*EcalBarrelDigis)[digis]; 
    int nrSamples = ebdf.size();  
    EBDetId ebid  = ebdf.id () ; 
    for (int sample = 0 ; sample < nrSamples; ++sample) {
      EcalMGPASample thisSample = ebdf[sample]; 

      bool toPlot = 
	(iEvent.id().event()==1 && ebid.ieta()==42 && ebid.iphi()==278) || 
	(iEvent.id().event()==1 && ebid.ieta()==-40 && ebid.iphi()==98) || 
	(iEvent.id().event()==2 && ebid.ieta()==-73 && ebid.iphi()==103) || 
	(iEvent.id().event()==2 && ebid.ieta()==-74 && ebid.iphi()==282) || 
	(iEvent.id().event()==3 && ebid.ieta()==-49 && ebid.iphi()==293) || 
	(iEvent.id().event()==3 && ebid.ieta()==46 && ebid.iphi()==114) || 
	(iEvent.id().event()==4 && ebid.ieta()==53 && ebid.iphi()==216) || 
	(iEvent.id().event()==4 && ebid.ieta()==-53 && ebid.iphi()==37) || 
	(iEvent.id().event()==6 && ebid.ieta()==59 && ebid.iphi()==185) || 
	(iEvent.id().event()==6 && ebid.ieta()==-61 && ebid.iphi()==6) || 
	(iEvent.id().event()==8 && ebid.ieta()==-21 && ebid.iphi()==300) || 
	(iEvent.id().event()==8 && ebid.ieta()==18 && ebid.iphi()==119) || 
	(iEvent.id().event()==11 && ebid.ieta()==56 && ebid.iphi()==55) || 
	(iEvent.id().event()==11 && ebid.ieta()==-57 && ebid.iphi()==233);

      if (toPlot) 
	cout << "Sample " << thisSample << ", adc = " << thisSample.adc() 
	     << ", gain = " << thisSample.gainId() 
	     << ", ieta = " << ebid.ieta() << ", iphi = " << ebid.iphi() << endl;
    }
  }

  // EE
  /*
  for (unsigned int digis=0; digis<EcalEndcapDigis->size(); ++digis) {
    EEDataFrame eedf = (*EcalEndcapDigis)[digis]; 
    int nrSamples = eedf.size();  
    EEDetId eeid  = eedf.id () ; 
    for (int sample = 0 ; sample < nrSamples; ++sample) {
      EcalMGPASample thisSample = eedf[sample]; 
      //cout << "Sample " << thisSample << ", adc = " << thisSample.adc() << ", gain = " << thisSample.gainId() << endl;
    }
  }
  */

}

void DigiAnalyzer::beginJob() { }

void DigiAnalyzer::endJob() { }

DEFINE_FWK_MODULE(DigiAnalyzer);
