// -*- C++ -*-
//
// Package:    PFAnalysis/PFAnalyzers
// Class:      PFTrackHFAnalyzer
//
/**\class PFTrackHFAnalyzer PFTrackHFAnalyzer.cc PFAnalysis/PFAnalyzers/plugins/PFTrackHFAnalyzer.cc

 Description: Analyzer of PFTracks/clusters in HF region
              The input step3 files will need: --outputCommand 'keep recoPFRecHits_particleFlow*_*_*','keep *_*pfTrack*_*_*'

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Kenichi Hatakeyama
//         Created:  Wed, 17 Jul 2019 15:24:34 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "boost/format.hpp"

#include "TH1.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class PFTrackHFAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit PFTrackHFAnalyzer(const edm::ParameterSet&);
      ~PFTrackHFAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::GenParticleCollection> genparToken_; 
      edm::EDGetTokenT<CaloParticleCollection> caloparToken_; 
      edm::EDGetTokenT<reco::VertexCollection> vertexToken_; 
      edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfcandToken_; 
      edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfcandPFAlgoToken_; 
      edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfcandPFTICLToken_; 
      edm::EDGetTokenT<std::vector<reco::PFCluster>> pfclusterHFToken_; 
      edm::EDGetTokenT<std::vector<reco::PFRecHit>> pfrechitHFToken_; 
      edm::EDGetTokenT<std::vector<reco::PFRecTrack>> pftrackToken_; 
      edm::EDGetTokenT<reco::TrackCollection> trackToken_;  //used to select what tracks to read from configuration file

      bool debug_;
  
      TH1I * histo;
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
PFTrackHFAnalyzer::PFTrackHFAnalyzer(const edm::ParameterSet& iConfig)
 :
  genparToken_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("source_genpars"))),
  caloparToken_(consumes<CaloParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("source_calopars"))),
  vertexToken_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("source_vertices"))),
  pfcandToken_(consumes<std::vector<reco::PFCandidate>>(iConfig.getUntrackedParameter<edm::InputTag>("source_pfcands"))),
  pfcandPFAlgoToken_(consumes<std::vector<reco::PFCandidate>>(iConfig.getUntrackedParameter<edm::InputTag>("source_pfcands_pfalgo"))),
  pfcandPFTICLToken_(consumes<std::vector<reco::PFCandidate>>(iConfig.getUntrackedParameter<edm::InputTag>("source_pfcands_pfticl"))),
  pfclusterHFToken_(consumes<std::vector<reco::PFCluster>>(iConfig.getUntrackedParameter<edm::InputTag>("source_pfclustersHF"))),
  pfrechitHFToken_(consumes<std::vector<reco::PFRecHit>>(iConfig.getUntrackedParameter<edm::InputTag>("source_pfrechitsHF"))),
  pftrackToken_(consumes<std::vector<reco::PFRecTrack>>(iConfig.getUntrackedParameter<edm::InputTag>("source_pftracks"))),
  trackToken_(consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("source_tracks"))),
  debug_(iConfig.getUntrackedParameter<bool>("debug"))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   histo = fs->make<TH1I>("charge" , "Charges" , 3 , -1 , 2 );

}


PFTrackHFAnalyzer::~PFTrackHFAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
PFTrackHFAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   edm::Handle<reco::GenParticleCollection> genpars; iEvent.getByToken(genparToken_, genpars);
   edm::Handle<CaloParticleCollection> calopars; iEvent.getByToken(caloparToken_, calopars);
   edm::Handle<reco::VertexCollection> vertices; iEvent.getByToken(vertexToken_, vertices);

   edm::Handle<std::vector<reco::PFCandidate>> pfcands; iEvent.getByToken(pfcandToken_, pfcands);
   edm::Handle<std::vector<reco::PFCandidate>> pfcands_pfalgo; iEvent.getByToken(pfcandPFAlgoToken_, pfcands_pfalgo);
   edm::Handle<std::vector<reco::PFCandidate>> pfcands_pfticl; iEvent.getByToken(pfcandPFTICLToken_, pfcands_pfticl);

   edm::Handle<std::vector<reco::PFCluster>> pfclustersHF; iEvent.getByToken(pfclusterHFToken_, pfclustersHF);
   edm::Handle<std::vector<reco::PFRecHit>> pfrechitsHF; iEvent.getByToken(pfrechitHFToken_, pfrechitsHF);
   edm::Handle<std::vector<reco::PFRecTrack>> pftracks; iEvent.getByToken(pftrackToken_, pftracks);

   //------------------------------ ------------------------------ ------------------------------
   if (debug_) {
   
     LogPrint("PFTrackHFAnalyzer") << "genpars: "      << genpars->size();
     LogPrint("PFTrackHFAnalyzer") << "calopars: "     << calopars->size();
     LogPrint("PFTrackHFAnalyzer") << "vertices: "     << vertices->size();
     LogPrint("PFTrackHFAnalyzer") << "pfcands: "      << pfcands->size();
     LogPrint("PFTrackHFAnalyzer") << "pfcands_pfalgo: "      << pfcands_pfalgo->size();
     LogPrint("PFTrackHFAnalyzer") << "pfcands_pfticl: "      << pfcands_pfticl->size();
     LogPrint("PFTrackHFAnalyzer") << "pfclustersHF: " << pfclustersHF->size();
     LogPrint("PFTrackHFAnalyzer") << "pfrechitsHF: "  << pfrechitsHF->size();
     LogPrint("PFTrackHFAnalyzer") << "pftracks: "     << pftracks->size();

     LogPrint("PFTrackHFAnalyzer") << "genpar: ";
     for(const auto& genpar : *(genpars.product()) )
       LogPrint("PFTrackHFAnalyzer") << boost::format("genpar (pt,eta,phi): (%6.1f, %6.2f, %6.2f)") % genpar.pt() % genpar.eta() % genpar.phi();
    
     LogPrint("PFTrackHFAnalyzer") << "calopar: ";
     for(const auto& calopar : *(calopars.product()) )
       LogPrint("PFTrackHFAnalyzer") << boost::format("calopar (pt,eta,phi): (%6.1f, %6.2f, %6.2f)") % calopar.pt() % calopar.eta() % calopar.phi();

     std::vector<unsigned> vtrackKey;
     std::map<unsigned, unsigned> duplicate;
     LogPrint("PFTrackHFAnalyzer") << "pfcand: ";
     for(const auto& pfcand : *(pfcands.product()) ){
       if (fabs(pfcand.charge()) && pfcand.trackRef().isNonnull())
	 LogPrint("PFTrackHFAnalyzer") << boost::format("pfcand (pt,eta,phi): (%6.1f, %6.2f, %6.2f) %5d %5d") % pfcand.pt() % pfcand.eta() % pfcand.phi() % pfcand.pdgId() % pfcand.trackRef().key();
       if (pfcand.charge() && pfcand.trackRef().isNonnull()) {  // otherwise PFCandidate throws
	 const auto trackKey = pfcand.trackRef().key();
	 vtrackKey.push_back(trackKey);
       }
     }
     std::sort(vtrackKey.begin(), vtrackKey.end());

     auto beg = begin(vtrackKey) + 1;
     for (;beg != end(vtrackKey); ++beg) {
       if (*beg == *(beg - 1)) {
	 duplicate[*beg]++;
       }
     }

     for (const auto& i : duplicate){
       std::cout << i.first << " appear " << i.second+1 << " times" << std::endl;
       LogWarning("PFTrackHFAnalyzer") << boost::format("duplicates");

       for(const auto& pfcand : *(pfcands.product()) ){
	 if (pfcand.charge()) {  // otherwise PFCandidate throws
	   const auto trackKey = pfcand.trackRef().key();
	   if (i.first==trackKey && pfcand.trackRef().isNonnull())
	     std::cout << " a " 
		       << pfcand.trackRef()->pt() << " " 
		       << pfcand.trackRef()->eta() << " " 
		       << pfcand.trackRef()->phi() << " " 
		       << pfcand.trackRef().key() << " " 
		       << pfcand.trackRef().id() << std::endl;
	 }
       }
       for(const auto& pfcand : *(pfcands_pfticl.product()) ){
	 if (pfcand.charge()) {  // otherwise PFCandidate throws
	   const auto trackKey = pfcand.trackRef().key();
	   if (i.first==trackKey && pfcand.trackRef().isNonnull())
	     std::cout << " pfticl " 
		       << pfcand.pt() << " " 
		       << pfcand.eta() << " " 
		       << pfcand.phi() << " " 
		       << pfcand.pdgId() << " " 
		       << pfcand.trackRef()->pt() << " " 
		       << pfcand.trackRef()->eta() << " " 
		       << pfcand.trackRef()->phi() << " " 
		       << pfcand.trackRef()->charge() << " " 
		       << pfcand.trackRef().key() << " " 
		       << pfcand.trackRef().id() << std::endl;
	 }
       }
       for(const auto& pfcand : *(pfcands_pfalgo.product()) ){
	 if (pfcand.charge()) {  // otherwise PFCandidate throws
	   const auto trackKey = pfcand.trackRef().key();
	   if (i.first==trackKey && pfcand.trackRef().isNonnull())
	     std::cout << " pfalgo " 
		       << pfcand.pt() << " " 
		       << pfcand.eta() << " " 
		       << pfcand.phi() << " " 
		       << pfcand.pdgId() << " " 
		       << pfcand.trackRef()->pt() << " " 
		       << pfcand.trackRef()->eta() << " " 
		       << pfcand.trackRef()->phi() << " " 
		       << pfcand.trackRef()->charge() << " " 
		       << pfcand.trackRef().key() << " " 
		       << pfcand.trackRef().id() << std::endl;
	 }
       }
     }

     vtrackKey.clear();
     LogPrint("PFTrackHFAnalyzer") << "pfticl: ";
     for(const auto& pfcand : *(pfcands_pfticl.product()) ){
       if (fabs(pfcand.charge()))
	 LogPrint("PFTrackHFAnalyzer") << boost::format("pfcand (pt,eta,phi): (%6.1f, %6.2f, %6.2f) %5d %5d") % pfcand.pt() % pfcand.eta() % pfcand.phi() % pfcand.pdgId() % pfcand.trackRef().key();
       if (pfcand.charge()) {  // otherwise PFCandidate throws
	 const auto trackKey = pfcand.trackRef().key();
	 vtrackKey.push_back(trackKey);
       }
     }

     vtrackKey.clear();
     LogPrint("PFTrackHFAnalyzer") << "pfalgo: ";
     for(const auto& pfcand : *(pfcands_pfalgo.product()) ){
       if (fabs(pfcand.charge()))
	 LogPrint("PFTrackHFAnalyzer") << boost::format("pfcand (pt,eta,phi): (%6.1f, %6.2f, %6.2f) %5d %5d") % pfcand.pt() % pfcand.eta() % pfcand.phi() % pfcand.pdgId() % pfcand.trackRef().key();
       if (pfcand.charge()) {  // otherwise PFCandidate throws
	 const auto trackKey = pfcand.trackRef().key();
	 vtrackKey.push_back(trackKey);
       }
     }
   
     LogPrint("PFTrackHFAnalyzer") << "pfclusterHF: ";
     for(const auto& pfclus : *(pfclustersHF.product()) ){
       if (pfclus.pt()>1.) {
	 LogPrint("PFTrackHFAnalyzer") << boost::format("pfclus (pt,eta,phi,E): (%6.1f, %6.2f, %6.2f, %6.1f)")
	   % pfclus.pt() % pfclus.eta() % pfclus.phi() % pfclus.energy();
	 const std::vector<reco::PFRecHitFraction> &fracs = pfclus.recHitFractions();
	 const std::vector<std::pair<DetId, float>> &hfracs = pfclus.hitsAndFractions();
	 for(unsigned i=0; i<fracs.size(); i++) {
	   const auto& id = hfracs[i].first.rawId();
	   const reco::PFRecHitRef& pfRecHits = fracs[i].recHitRef();
	   double rawenergy = pfRecHits->energy();
	   LogPrint("PFTrackHFAnalyzer") << boost::format(" rechit (ieta,iphi,depth,frac,E): (%3d, %3d, %2d, %6.1f, %6.1f)")
	     % HcalDetId(id).ieta() % HcalDetId(id).iphi() % HcalDetId(id).depth() % hfracs[i].second % rawenergy;
	 }
       }
     }
     
     LogPrint("PFTrackHFAnalyzer") << "pftrack: ";
     for(const auto& pftrack : *(pftracks.product()) ){

       const reco::TrackRef trackref = pftrack.trackRef();
       std::vector<reco::PFTrajectoryPoint> trajectoryPoints = pftrack.trajectoryPoints();
       
       constexpr reco::PFTrajectoryPoint::LayerType VFcalEntrance =
       	 reco::PFTrajectoryPoint::VFcalEntrance;
       
       const reco::PFTrajectoryPoint& tkAtHF =
	 pftrack.extrapolatedPoint( VFcalEntrance );

       const double tracketa = tkAtHF.positionREP().Eta();
       const double trackphi = tkAtHF.positionREP().Phi();

       LogPrint("PFTrackHFAnalyzer") << boost::format("pftrack (pt,eta,phi)@origin (eta,phi)@HF: (%6.1f +- %4.1f, %6.2f, %6.2f) (%6.2f, %6.2f)")
	 % trackref->pt() % trackref->ptError() % trackref->eta() % trackref->phi() % tracketa % trackphi;

     }
     
   } // debug_ printout ends here
   //------------------------------ ------------------------------ ------------------------------
   
   for(const auto& track : iEvent.get(trackToken_) ) {
     // do something with track parameters, e.g, plot the charge.
     // int charge = track.charge();
     //std::cout << track.charge() << std::endl;
     histo->Fill( track.charge() );
   }
   
}


// ------------ method called once each job just before starting event loop  ------------
void
PFTrackHFAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
PFTrackHFAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PFTrackHFAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PFTrackHFAnalyzer);
