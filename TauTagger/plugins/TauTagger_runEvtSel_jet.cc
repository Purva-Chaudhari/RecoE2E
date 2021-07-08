#include "RecoE2E/TauTagger/interface/TauTagger.h"
#include "RecoE2E/FrameProducers/interface/DetFrameProducer.h"

const int search_window = 7;
const int image_padding = 12;
vector<int>   vFailedJetIdx_;
extern unsigned int jet_runId_;
extern unsigned int jet_lumiId_;
extern unsigned long long jet_eventId_;


void TauTagger::runEvtSel_jet ( const edm::Event& iEvent, const edm::EventSetup& iSetup, e2e::Frame2D& vJetSeeds ) {

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles );
 
  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();
	
  edm::Handle<HBHERecHitCollection> HBHERecHitsH_;
  iEvent.getByToken( HBHERecHitCollectionT_, HBHERecHitsH_ );
	
  if ( debug ) std::cout << " >> PFJetCol.size: " << jets->size() << std::endl;

  //edm::Handle<reco::PFTauCollection> taus;
  //iEvent.getByToken(tauCollectionT_, taus);

  float seedE;
  int iphi_, ieta_, ietaAbs_;
  int nJet = 0;
  vFailedJetIdx_.clear();

  vJetIdxs.clear();
  v_tau_jetPdgIds_.clear();
  v_jetIsTau.clear();
  v_jetdR.clear();
  v_TaudR.clear();
  v_TaupT1.clear();
  v_TaupT2.clear();

  /*
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJetCollectionT_, genJets);
  std::vector<float> v_tau_jetFakePhoIdxs;
  */

  unsigned int nMatchedJets = 0;
  unsigned int PdgId        = 0;
  float jetdR               = -99.;
  float taudR               = -99.;
  float taupT1              = -99.;
  float taupT2              = -99.;
  bool JetIsTau             = false;

  edm::LogInfo("TopTagger") << " >> Reading and selecting Jets from " << jets->size() << " jet seeds: ";

 if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
 if ( debug ) std::cout << " JETS IN THE EVENT = " << jets->size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;
  // Loop over jets
  for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {

    bool keepJet = true;
    	iphi_ = -1;
    	ieta_ = -1;

    reco::PFJetRef iJet( jets, iJ );
    if ( std::abs(iJet->pt())  < minJetPt_ ) {keepJet = false; 
		//std::cout<<"     * Selection failed at Jet index: "<<iJ<<" because abs(pt) < minJetPt_ --> pt: "<<std::abs(iJet->pt())<<" minJetPt_: "<<minJetPt_<<". Adding "<<iphi_<<" to JetSeediphi and "<<ieta_<<" JetSeedieta vectors."<<std::endl;
		edm::LogInfo("TauTagger") << "     * Selection failed at Jet index: "<<iJ<<" because abs(pt) < minJetPt_ --> pt: "<<std::abs(iJet->pt())<<" minJetPt_: "<<minJetPt_<<". Adding "<<iphi_<<" to JetSeediphi and "<<ieta_<<" JetSeedieta vectors.";
    	nJet++;	
        continue;
    }
    if ( std::abs(iJet->eta()) > maxJetEta_ ) {keepJet = false; 
		//std::cout<<"     * Selection failed at Jet index: "<<iJ<<" because abs(eta) > maxJetEta_ --> eta: "<<std::abs(iJet->eta())<<" maxJetEta_: "<<maxJetEta_<<". Adding "<<iphi_<<" to JetSeediphi and "<<ieta_<<" JetSeedieta vectors."<<std::endl;
		edm::LogInfo("TauTagger") << "     * Selection failed at Jet index: "<<iJ<<" because abs(eta) > maxJetEta_ --> eta: "<<std::abs(iJet->eta())<<" maxJetEta_: "<<maxJetEta_<<". Adding "<<iphi_<<" to JetSeediphi and "<<ieta_<<" JetSeedieta vectors.";
    	nJet++;	
        continue;
    }
    
    if (debug ) std::cout << "  >>>>>> Jet [" << iJ << "] ->  Pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi() << std::endl;
    bool passedGenSel = false;
    unsigned int iGenParticle = 0;
    for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
      float dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );
      if ( dR > 0.4 ) continue;
      
      //if ( !(  (std::abs(iGen->pdgId()) == 15 && iGen->status() == 2 && iGen->numberOfMothers() == 0 ) || iGen->status() == 23 || iGen->status() == 43 || iGen->status() == 43) ) continue;
      
      if ( std::abs(iGen->pdgId()) != 25 ) continue;
      if ( iGen->numberOfDaughters() != 2 ) continue;
      passedGenSel = true;
      taudR = reco::deltaR( iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi() );
        
      ++iGenParticle;

      if ( debug ) std::cout << "   GEN particle " << iGenParticle << " -> status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << " nMoms: " <<iGen->numberOfMothers() << " | pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << " | dR: " << dR << std::endl;

      PdgId = std::abs(iGen->pdgId());
      jetdR = dR;
      if ( iGen->daughter(0)->pt() > iGen->daughter(1)->pt() ){
        taupT1 = iGen->daughter(0)->pt(); 
        taupT2 = iGen->daughter(1)->pt();
      } else {
        taupT1 = iGen->daughter(1)->pt(); 
        taupT2 = iGen->daughter(0)->pt();
      }

      if (PdgId == 15)  JetIsTau = true;
       
    } // primary gen particles

    if (passedGenSel) { 
      ++nMatchedJets;
      vJetIdxs.push_back( iJ );
      v_tau_jetPdgIds_.push_back( PdgId );
      v_jetdR.push_back( jetdR );
      v_TaudR.push_back( taudR );
      v_TaupT1.push_back( taupT1 );
      v_TaupT2.push_back( taupT2 );
      v_jetIsTau.push_back( JetIsTau );

    }

    if (keepJet){
	if ( debug ) std::cout << " >> jet[" << iJ << "]Pt:" << iJet->pt()  << " Eta:" << iJet->eta()  << " Phi:" << iJet->phi() 
			   << " jetE:" << iJet->energy() << " jetM:" << iJet->mass() << std::endl;
	HcalDetId hId( spr::findDetIdHCAL( caloGeom, iJet->eta(), iJet->phi(), false ) );
    	if ( hId.subdet() != HcalBarrel && hId.subdet() != HcalEndcap ){
      		vFailedJetIdx_.push_back(iJ);
      		continue;
    	}
	   
	HBHERecHitCollection::const_iterator iRHit( HBHERecHitsH_->find(hId) );
    	seedE = ( iRHit == HBHERecHitsH_->end() ) ? 0. : iRHit->energy() ;
    	HcalDetId seedId = hId;
    	if ( debug ) std::cout << " >> hId.ieta:" << hId.ieta() << " hId.iphi:" << hId.iphi() << " E:" << seedE << std::endl;
	
	// Look for the most energetic HBHE tower deposit within a search window
	for ( int ieta = 0; ieta < search_window; ieta++ ) {

      		ieta_ = hId.ieta() - (search_window/2)+ieta;

      		if ( std::abs(ieta_) > HBHE_IETA_MAX_HE-1 ) continue;
      		if ( std::abs(ieta_) < HBHE_IETA_MIN_HB ) continue;

      		HcalSubdetector subdet_ = std::abs(ieta_) > HBHE_IETA_MAX_HB ? HcalEndcap : HcalBarrel;

      		for ( int iphi = 0; iphi < search_window; iphi++ ) {

        		iphi_ = hId.iphi() - (search_window/2)+iphi;

        		// iphi should wrap around
        		if ( iphi_ > HBHE_IPHI_MAX ) {
          		iphi_ = iphi_-HBHE_IPHI_MAX;
        		} else if ( iphi_ < HBHE_IPHI_MIN ) {
          		iphi_ = HBHE_IPHI_MAX-abs(iphi_); 
        		}
        		// Skip non-existent and lower energy towers 
        		HcalDetId hId_( subdet_, ieta_, iphi_, 1 );
        		HBHERecHitCollection::const_iterator iRHit( HBHERecHitsH_->find(hId_) );
        		if ( iRHit == HBHERecHitsH_->end() ) continue;
        		if ( iRHit->energy() <= seedE ) continue;
        		if ( debug ) std::cout << " !! hId.ieta:" << hId_.ieta() << " hId.iphi:" << hId_.iphi() << " E:" << iRHit->energy() << std::endl;

        		seedE = iRHit->energy();
        		seedId = hId_;

      		} // iphi 
    	} // ieta
	
    	// NOTE: HBHE iphi = 1 does not correspond to EB iphi = 1!
    	// => Need to shift by 2 HBHE towers: HBHE::iphi: [1,...,71,72]->[3,4,...,71,72,1,2]
    	iphi_  = seedId.iphi() + 2; // shift
    	iphi_  = iphi_ > HBHE_IPHI_MAX ? iphi_-HBHE_IPHI_MAX : iphi_; // wrap-around
    	iphi_  = iphi_ - 1; // make histogram-friendly
    	ietaAbs_  = seedId.ietaAbs() == HBHE_IETA_MAX_HE ? HBHE_IETA_MAX_HE-1 : seedId.ietaAbs(); // last HBHE ieta embedded
    	ieta_  = seedId.zside() > 0 ? ietaAbs_-1 : -ietaAbs_;
    	ieta_  = ieta_+HBHE_IETA_MAX_HE-1;

    	// If the seed is too close to the edge of HE, discard event
    	// Required to keep the seed at the image center
    	if ( HBHE_IETA_MAX_HE-1 - ietaAbs_ < image_padding ) { 
      	if ( debug ) std::cout << " Fail HE edge cut " << std::endl;
      	vFailedJetIdx_.push_back(iJ);
	//std::cout<<"     * Failed Jet seed at index: "<<iJ<<" seed is too close to the edge of HE. Adding -1 to JetSeediphi and JetSeedieta vectors."<<std::endl;
	edm::LogInfo("TopTagger") << "     * Failed Jet seed at index: " << iJ << " seed is too close to the edge of HE. Adding -1 to JetSeediphi and JetSeedieta vectors.";	
	ieta_=-1;
	iphi_=-1;
    	nJet++;
      	continue;
    	}
    	// Save position of most energetic HBHE tower
    	// in EB-aligned coordinates
    	if ( debug ) std::cout << " !! ieta_:" << ieta_ << " iphi_:" << iphi_ << " ietaAbs_:" << ietaAbs_ << " E:" << seedE << std::endl;
	
    	vJetSeeds[iJ][0] = ieta_;
	vJetSeeds[iJ][1] = iphi_;
    	nJet++;
      }

  } // reco jets

  if ( debug ) std::cout << " Matched jets " << nMatchedJets << std::endl;

  // Check jet multiplicity (return false)
  if ( nMatchedJets < 1 ) std::cout << " Matched jets Less than 1 " << nMatchedJets << std::endl;

  if ( debug ) std::cout << " >> Event contains a tau jet" << std::endl;
  return ;

} // runEvtSel_jet()
