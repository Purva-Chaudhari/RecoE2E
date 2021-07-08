#include "RecoE2E/FrameProducers/interface/DetFrameProducer.h"
#include "RecoE2E/FrameProducers/interface/JetFrameProducer.h"


using std::vector;

const unsigned nJets = 10; //TODO: use cfg level nJets_

TH1D *h_tau_jet_nJet;

vector<float> v_tau_jet_m0_;
vector<float> v_tau_jet_pt_;
vector<float> v_tau_jetPdgIds_;
vector<float> v_tau_jetIsTau_;
vector<float> v_tau_jetdR_;
vector<float> v_tau_TaudR_;
vector<float> v_tau_TaupT1_;
vector<float> v_tau_TaupT2_;

vector<float> v_tau_subJetE_[nJets];
vector<float> v_tau_subJetPx_[nJets];
vector<float> v_tau_subJetPy_[nJets];
vector<float> v_tau_subJetPz_[nJets];


void DetFrameProducer::fillEvtSeljetDijetTau ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  h_tau_jet_nJet->Fill( vJetIdxs.size() );

  v_tau_jet_pt_.clear();
  v_tau_jet_m0_.clear();
  v_tau_jetIsTau_.clear();
  v_tau_jetdR_.clear();
  v_tau_TaudR_.clear();
  v_tau_TaupT1_.clear();
  v_tau_TaupT2_.clear();
 
  for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ ) {

    reco::PFJetRef iJet( jets, vJetIdxs[iJ] );

  

    // Fill branches 
    v_tau_jet_pt_.push_back( iJet->pt() );
    v_tau_jet_m0_.push_back( iJet->mass() );
    v_tau_jetIsTau_.push_back( v_jetIsTau[iJ] );
    v_tau_jetdR_.push_back( v_jetdR[iJ] );

    // Gen jet constituents
    v_tau_subJetE_[iJ].clear();
    v_tau_subJetPx_[iJ].clear();
    v_tau_subJetPy_[iJ].clear();
    v_tau_subJetPz_[iJ].clear();
    //std::vector<reco::PFCandidatePtr> jetConstituents = iJet->getPFConstituents();
    unsigned int nConstituents = iJet->getPFConstituents().size();
    //if ( debug ) std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
    for ( unsigned int j = 0; j < nConstituents; j++ ) {
      const reco::PFCandidatePtr subJet = iJet->getPFConstituent( j );
      //if ( debug ) std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << subJet->energy() << " px:" << subJet->px() << " py:" << subJet->py() << " pz:" << subJet->pz() << std::endl;
      v_tau_subJetE_[iJ].push_back( subJet->energy() );
      v_tau_subJetPx_[iJ].push_back( subJet->px() );
      v_tau_subJetPy_[iJ].push_back( subJet->py() );
      v_tau_subJetPz_[iJ].push_back( subJet->pz() );
    }
  }

} // fillEvtSeljetDijetTau()
