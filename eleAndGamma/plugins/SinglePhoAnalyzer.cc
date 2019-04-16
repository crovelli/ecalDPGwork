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
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

using namespace std;
using namespace edm;
using namespace reco;

// per event tree
struct eventTree_struc_ {  

  int nPho;
  int nMatchedPho;
};

// per photon tree
struct phoTree_struc_ {    

  int run;
  int event;
  float nputrue;  
  int npuobs;  

  float pt;
  float eta;
  float phi;
  float scEta;
  float scPhi;
  float scRawEnergy;

  float eMax;  
  float e5x5;
  float r9;  

  float trueEnergy;
  float truePt;
  float trueEta;
  float truePhi;

  float amplit[25];
  int ieta[25];
  int iphi[25];
  int ix[25];
  int iy[25];
  int iz[25];
};

class SinglePhoAnalyzer : public edm::EDAnalyzer {
public:
  explicit SinglePhoAnalyzer(const edm::ParameterSet&);
  ~SinglePhoAnalyzer();
  
private:

  edm::Service<TFileService> fs_;
  
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  void initEvtTreeStructure();
  void initPhoTreeStructure();

  // declarations: trees per event and photon
  TTree *eventTree;
  TTree *phoTree;  
  eventTree_struc_ treeEv_;
  phoTree_struc_ tree_;  

  // collections
  edm::EDGetTokenT<PhotonCollection> photonToken_;   
  edm::EDGetTokenT<EcalRecHitCollection> ecalHitEBToken_;
  edm::EDGetTokenT<EcalRecHitCollection> ecalHitEEToken_;
  edm::EDGetTokenT<GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<vector< PileupSummaryInfo > > PileUpToken_; 
};

SinglePhoAnalyzer::SinglePhoAnalyzer(const edm::ParameterSet& iConfig) {

  photonToken_       = consumes<PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"));  
  ecalHitEBToken_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedBarrelRecHitCollection"));
  ecalHitEEToken_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEndcapRecHitCollection"));
  genParticlesToken_ = consumes<GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  PileUpToken_       = consumes<vector< PileupSummaryInfo > >(iConfig.getParameter<edm::InputTag>("PileUpTag"));
}

SinglePhoAnalyzer::~SinglePhoAnalyzer() { } 

void SinglePhoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // access edm objects
  Handle< PhotonCollection > photons;   
  iEvent.getByToken(photonToken_, photons);
  
  Handle< GenParticleCollection > genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);

  Handle< EcalRecHitCollection > EcalBarrelRecHits;
  iEvent.getByToken(ecalHitEBToken_, EcalBarrelRecHits); 

  Handle< EcalRecHitCollection > EcalEndcapRecHits;
  iEvent.getByToken(ecalHitEEToken_, EcalEndcapRecHits); 

  Handle<vector< PileupSummaryInfo> > PileupInfos;
  iEvent.getByToken(PileUpToken_,PileupInfos);

  const CaloSubdetectorTopology* theSubdetTopologyEB_;
  const CaloSubdetectorTopology* theSubdetTopologyEE_;
  edm::ESHandle<CaloTopology> theCaloTopology;
  iSetup.get<CaloTopologyRecord>().get(theCaloTopology);
  theSubdetTopologyEB_ = theCaloTopology->getSubdetectorTopology(DetId::Ecal,EcalBarrel);
  theSubdetTopologyEE_ = theCaloTopology->getSubdetectorTopology(DetId::Ecal,EcalEndcap);
  


  // real analysis
  int allGamma = 0;
  int matchedGamma = 0;

  // event tree re-initialization
  initEvtTreeStructure();

  cout << endl;
  cout << endl;

  // loop over photons
  PhotonCollection::const_iterator g1;
  for (g1 = photons->begin(); g1 != photons->end(); ++g1) {

    // minimal selection to speed up: photon eta and pT cuts
    float gammaPt = g1->pt();
    if (gammaPt<10) continue;
    
    // photon tree re-initialization
    initPhoTreeStructure();

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

    // filling the photon tree with basic infos
    tree_.run     = iEvent.id().run(); 
    tree_.event   = iEvent.id().event(); 
    tree_.nputrue = pu_n; 
    tree_.npuobs  = pu_obs; 

    tree_.pt  = g1->pt();
    tree_.eta = g1->eta();
    tree_.phi = g1->phi();
    tree_.scEta = (g1->superCluster())->eta();
    tree_.scPhi = (g1->superCluster())->phi();
    tree_.scRawEnergy = (g1->superCluster())->rawEnergy();

    tree_.eMax = g1->maxEnergyXtal();  
    tree_.e5x5 = g1->e5x5();
    tree_.r9   = g1->r9();

    cout << "Run = " << iEvent.id().run() << ", event = " << iEvent.id().event() << ", pu_true = " << pu_n << ", pu_obs = " << pu_obs << endl;

    // we save all photons, and we decide using these flags if they are true or not
    float matchedEne = -999.;
    float matchedPt  = -999.;
    float matchedEta = -999.;
    float matchedPhi = -999.;
    double dR = 999;

    //cout << endl;
    GenParticleCollection::const_iterator genP;
    for( genP = genParticles->begin(); genP<genParticles->end(); ++genP ){
      
      if( abs(genP->pdgId()) != 22 || genP->status() != 1 ) continue;
      
      TLorentzVector theGen(0,0,0,0);
      TLorentzVector theRec(0,0,0,0);
      theGen.SetPtEtaPhiM(genP->pt(), genP->eta(), genP->phi(), 0.);
      theRec.SetPtEtaPhiM(g1->pt(), g1->eta(), g1->phi(), 0.);
      float dRtmp = theGen.DeltaR( theRec );
      if( dRtmp < dR && dRtmp<0.1 ){
	dR = dRtmp;
	matchedEne = genP->energy();
	matchedPt  = genP->pt();
	matchedEta = genP->eta();
	matchedPhi = genP->phi();
      }
    }
    tree_.trueEnergy = matchedEne;
    tree_.truePt = matchedPt;
    tree_.trueEta = matchedEta;
    tree_.truePhi = matchedPhi;

    cout << "Etrue = " << matchedEne << ", eMax = " << g1->maxEnergyXtal() << endl;
    
    // extra info on rechits for xtals in the 5x5 matrix around the seed
    DetId seedDetId = ( (g1->superCluster())->seed() )->seed();
    
    if(seedDetId.subdetId()==EcalEndcap) {
      
      int iNeigh=0; 
      
      CaloNavigator<DetId> cursorE = CaloNavigator<DetId>(seedDetId, theSubdetTopologyEE_ );

      for(int ix=-2; ix<3; ++ix) {
	for(int iy=-2; iy<3; ++iy) {
	  cursorE.home();
	  cursorE.offsetBy( ix, iy );
	  DetId cryId = cursorE.pos();

	  if(cryId.subdetId()!=EcalEndcap) { 
	    tree_.amplit[iNeigh] = -5000.;
	    tree_.ieta[iNeigh] = -5000; 
	    tree_.iphi[iNeigh] = -5000;
	    tree_.ix[iNeigh] = -5000; 
	    tree_.iy[iNeigh] = -5000; 
	    tree_.iz[iNeigh] = -5000; 
	    iNeigh++;
	    continue;  
	  }

	  EcalRecHitCollection::const_iterator itneigh = EcalEndcapRecHits->find( cryId );

	  if( itneigh != EcalEndcapRecHits->end() ) {
	    tree_.amplit[iNeigh] = itneigh->energy();
	    tree_.ieta[iNeigh] = -999; 
	    tree_.iphi[iNeigh] = -999; 
	    tree_.ix[iNeigh] = ((EEDetId)itneigh->detid()).ix();
	    tree_.iy[iNeigh] = ((EEDetId)itneigh->detid()).iy();
	    tree_.iz[iNeigh] = ((EEDetId)itneigh->detid()).zside();

	    if (iNeigh==12)
	      cout << "In EE: ix = "    << ((EEDetId)itneigh->detid()).ix() 
		   << ", iy = "         << ((EEDetId)itneigh->detid()).iy() 
		   << ", iz = "         << ((EEDetId)itneigh->detid()).zside()
		   << ", amplit[12] = " << itneigh->energy() << endl;

	  } else {
	    tree_.amplit[iNeigh] = -2000.;
	    tree_.ieta[iNeigh]   = -2000; 
	    tree_.iphi[iNeigh]   = -2000;
	    tree_.ix[iNeigh] = -2000; 
	    tree_.iy[iNeigh] = -2000; 
	    tree_.iz[iNeigh] = -2000; 
	  }
	  
	  iNeigh++;
	}
      }
      if (iNeigh!=25) cout << "problem: not 25 crystals!  ==> " << iNeigh << endl;
      
    } else if (seedDetId.subdetId()==EcalBarrel) {

      int iNeigh=0; 

      CaloNavigator<DetId> cursorE = CaloNavigator<DetId>(seedDetId, theSubdetTopologyEB_ );

      for(int ix=-2; ix<3; ++ix) {
	for(int iy=-2; iy<3; ++iy) {
	  cursorE.home();
	  cursorE.offsetBy( ix, iy );
	  DetId cryId = cursorE.pos();

	  if(cryId.subdetId()!=EcalBarrel) { 
	    tree_.amplit[iNeigh] = -5000.;
	    tree_.ieta[iNeigh] = -5000; 
	    tree_.iphi[iNeigh] = -5000;
	    tree_.ix[iNeigh] = -5000; 
	    tree_.iy[iNeigh] = -5000; 
	    tree_.iz[iNeigh] = -5000; 
	    iNeigh++;
	    continue;  
	  }
	  
	  EcalRecHitCollection::const_iterator itneigh = EcalBarrelRecHits->find( cryId );

	  if( itneigh != EcalBarrelRecHits->end() ) { 
	    tree_.amplit[iNeigh] = itneigh->energy();
	    tree_.ieta[iNeigh] = ((EBDetId)itneigh->detid()).ieta();
	    tree_.iphi[iNeigh] = ((EBDetId)itneigh->detid()).iphi();
	    tree_.ix[iNeigh] = -999;
	    tree_.iy[iNeigh] = -999;
	    tree_.iz[iNeigh] = -999;

	    if (iNeigh==12)
	      cout << "In EB: ieta = "  << ((EBDetId)itneigh->detid()).ieta() 
		   << ", iphi = "       << ((EBDetId)itneigh->detid()).iphi() 
		   << ", amplit[12] = " << itneigh->energy() << endl;

	  } else {
	    tree_.amplit[iNeigh] = -2000.;
	    tree_.ieta[iNeigh] = -2000; 
	    tree_.iphi[iNeigh] = -2000;
	    tree_.ix[iNeigh] = -2000; 
	    tree_.iy[iNeigh] = -2000; 
	    tree_.iz[iNeigh] = -2000; 
	  }

	  iNeigh++;
	}
      }
      if (iNeigh!=25) cout << "problem: not 25 crystals!  ==> " << iNeigh << endl;
    }

    // to count the number of reco photons and reco matching gen ones
    allGamma++;
    if (matchedEne>-1) matchedGamma++;

    phoTree->Fill();

  } // loop over reco photons

  // filling the tree with per-event infos
  treeEv_.nPho = allGamma;
  treeEv_.nMatchedPho = matchedGamma;
  eventTree->Fill();
}

void SinglePhoAnalyzer::beginJob() {

  eventTree = fs_->make<TTree>("eventTree","per-event tree");
  phoTree   = fs_->make<TTree>("singlePhotons","single photon tree");  
  
  // tree per event
  TString treeEvent = "nPho/I:nMatchedPho/I";
  eventTree->Branch("event",&(treeEv_.nPho),treeEvent); 

  // tree per photon
  TString treeKine = "pt/F:eta/F:phi/F";
  TString treeSc   = "scEta/F:scPhi/F:scRawEnergy/F";
  TString treeEne  = "eMax/F:e5x5/F:r9/F";
  TString treeTrue = "trueEnergy/F:truePt/F:trueEta/F:truePhi/F";
  TString tree5x5  = "amplit[25]/F:ieta[25]/I:iphi[25]/I:ix[25]/I:iy[25]/I:iz[25]/I";
  TString treeEve  = "run/I:event/I:nputrue/F:npuobs/I";
  phoTree->Branch("kinematics",&(tree_.pt),treeKine);
  phoTree->Branch("supercluster",&(tree_.scEta),treeSc);
  phoTree->Branch("energy",&(tree_.eMax),treeEne);
  phoTree->Branch("mctruth",&(tree_.trueEnergy),treeTrue);
  phoTree->Branch("tree5x5",&(tree_.amplit),tree5x5);  
  phoTree->Branch("eventInfo",&(tree_.run),treeEve);
}

void SinglePhoAnalyzer::endJob() { }

void SinglePhoAnalyzer::initPhoTreeStructure() {

  tree_.run=-500; 
  tree_.event=-500; 
  tree_.nputrue=-500.;
  tree_.npuobs=-500;
  tree_.pt=-500.;
  tree_.eta=-500.;
  tree_.phi=-500.;
  tree_.scEta=-500.;
  tree_.scPhi=-500.;
  tree_.scRawEnergy=-500.;
  tree_.eMax=-500.;
  tree_.e5x5=-500.;
  tree_.r9=-500.; 
  tree_.trueEnergy=-500.;
  tree_.truePt=-500.;
  tree_.trueEta=-500.;
  tree_.truePhi=-500.;
  for (uint iNeigh=0; iNeigh<25; iNeigh++) {
    tree_.amplit[iNeigh]=-500.;
    tree_.ieta[iNeigh]=-500;
    tree_.iphi[iNeigh]=-500;
    tree_.ix[iNeigh]=-500;
    tree_.iy[iNeigh]=-500;
    tree_.iz[iNeigh]=-500;
  }
}

void SinglePhoAnalyzer::initEvtTreeStructure() {

  treeEv_.nPho=-500;
  treeEv_.nMatchedPho=-500;
}

DEFINE_FWK_MODULE(SinglePhoAnalyzer);
