// -*- C++ -*-
//
// Package:    TTHLep/ObjEvtOptimisation
// Class:      ObjEvtOptimisation
// 
/**\class ObjEvtOptimisation ObjEvtOptimisation.cc TTHLep/ObjEvtOptimisation/plugins/ObjEvtOptimisation.cc

 Description:
 This class aims to map miniAOD infos into flat tree for stand-alone analyses (using only root)

 User defined classes:
 TreeVariables.h:      It save the variables in the tree
 HelpingFunctions.h:   It define some methods that are sometimes used in the main workflow (e.g. to order object by pT, or cout object pT,eta)
 ObjEvtFunctions.h:    It define some methods for object id and event selection
 LeptonJetVariables.h: It pick up the variables related to the jet associated to the lepton (for discriminating TTHMultilep from ttbar)

 Implementation:
 The implementation follows some rules. 
 If you modify part of the code, you are kindly invited to be coherent with the style it is written.
 For instance, pay attention to  
  the way you write comments
  the indentation
  the use of methods and functions in the main code 

 Usefull reference:
 https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
*/
//
// Original Author:  Francesco Romeo
//         Created:  Mon, 19 Jan 2015 08:54:10 GMT
//
//
/////
//   Headers
/////
//System and event
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TStopwatch.h"
//Gen info
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//Trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
//Vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//Muon
#include "DataFormats/PatCandidates/interface/Muon.h"
//Electron
#include "DataFormats/PatCandidates/interface/Electron.h"
//Tau
#include "DataFormats/PatCandidates/interface/Tau.h"
//Lepton
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
//Photon
#include "DataFormats/PatCandidates/interface/Photon.h"
//Jet and Met
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
//Packed
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
//Track builder infos
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
//Math
#include "DataFormats/Math/interface/deltaR.h"
//Store info
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandidateWithRef.h"
//User defined classes
#include "TTHLep/ObjEvtOptimisation/interface/TreeVariables.h"
#include "TTHLep/ObjEvtOptimisation/interface/ObjEvtFunctions.h"
#include "TTHLep/ObjEvtOptimisation/interface/LeptonJetVariables.h"
/////
//   Namespaces
/////
using namespace reco;
using namespace edm;
using namespace std;
/////
//   Class declaration
/////
class ObjEvtOptimisation : public edm::EDAnalyzer {
 public:
 explicit ObjEvtOptimisation(const edm::ParameterSet&);
 ~ObjEvtOptimisation();
 static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 private:
 //Default methods
 virtual void beginJob() override;
 virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
 virtual void endJob() override;
 //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
 //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
 //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
 //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
 //Collections of objects
 edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
 edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
 edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
 edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
 edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
 edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
 edm::EDGetTokenT<pat::MuonCollection> muonToken_;
 edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
 edm::EDGetTokenT<pat::TauCollection> tauToken_;
 edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
 edm::EDGetTokenT<pat::JetCollection> jetToken_;
 edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
 edm::EDGetTokenT<pat::METCollection> metToken_;
 edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
 edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
 //Values for begin job
 int evt_totnum;
 //Values for the whole analysis
 const double mindr_p3;
 const double mindr_p5;
 //Watch time and cpu for the analysis
 TStopwatch* stopwatch;
 //Tree
 CTree *tree;
 const edm::Service<TFileService> fs;
};
/////
//   Constructors and destructor
/////
ObjEvtOptimisation::ObjEvtOptimisation(const edm::ParameterSet& iConfig):
 //Collections of objects
 prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
 packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
 triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
 triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
 triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
 vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
 muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
 electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
 tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
 photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
 jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
 fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
 metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
 pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
 lostTracksToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("lostTracks"))),
 //Values for the whole analysis
 mindr_p3(iConfig.getUntrackedParameter<double>("mindr_p3")),
 mindr_p5(iConfig.getUntrackedParameter<double>("mindr_p5")),
 //TTree 
 tree(new CTree(fs->make<TTree>("tree", "tree")))   
{
 //Now do what ever initialization is needed
 tree->make_branches();
 stopwatch = new TStopwatch();
}

ObjEvtOptimisation::~ObjEvtOptimisation()
{
 // do anything here that needs to be done at desctruction time
 // (e.g. close files, deallocate resources etc.)
 delete stopwatch;
}
/////
//   Member functions
/////
// ------------ method called for each event  ------------
void ObjEvtOptimisation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
 using namespace edm;
 evt_totnum++; 
 tree->loop_initialize();
 /////
 //   Handle iEvent.getByToken
 /////
 Handle<edm::View<reco::GenParticle> > pruned;
 iEvent.getByToken(prunedGenToken_,pruned);
 Handle<edm::View<pat::PackedGenParticle> > packed;
 iEvent.getByToken(packedGenToken_,packed);
 edm::Handle<edm::TriggerResults> triggerBits;
 iEvent.getByToken(triggerBits_, triggerBits);
 edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
 iEvent.getByToken(triggerObjects_, triggerObjects);
 edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
 iEvent.getByToken(triggerPrescales_, triggerPrescales);
 edm::Handle<reco::VertexCollection> vertices;
 iEvent.getByToken(vtxToken_, vertices);
 edm::Handle<pat::MuonCollection> muons;
 iEvent.getByToken(muonToken_, muons);
 edm::Handle<pat::ElectronCollection> electrons;
 iEvent.getByToken(electronToken_, electrons);
 edm::Handle<pat::TauCollection> taus;
 iEvent.getByToken(tauToken_, taus);
 edm::Handle<pat::PhotonCollection> photons;
 iEvent.getByToken(photonToken_, photons);
 edm::Handle<pat::JetCollection> jets;
 iEvent.getByToken(jetToken_, jets);
 PFJetIDSelectionFunctor pfLooseJetID(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE),
                         pfTightJetID(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT);
 pat::strbitset passLooseCuts(pfLooseJetID.getBitTemplate()),
                passTightCuts(pfTightJetID.getBitTemplate());
 edm::Handle<pat::JetCollection> fatjets;
 iEvent.getByToken(fatjetToken_, fatjets);
 edm::Handle<pat::METCollection> mets;
 iEvent.getByToken(metToken_, mets);
 edm::Handle<pat::PackedCandidateCollection> pfs;
 iEvent.getByToken(pfToken_, pfs);
 edm::Handle<pat::PackedCandidateCollection> lostrks;
 iEvent.getByToken(lostTracksToken_, lostrks);
 edm::ESHandle<TransientTrackBuilder> ttrkbuilder;
 iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrkbuilder);
 /////
 //   Trigger 
 /////
 //bool istriggeredevent = false;
 //const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
 //for(unsigned int i = 0; i<triggerBits->size(); ++i){
 // if(names.triggerName(i)=="mytriggername" && triggerBits->accept(i)) istriggeredevent = true;
 //} 
 //if(!istriggeredevent) return;
 /////
 //   Primary vertex
 /////
 if(vertices->empty()) return; // skip the event if no PV found
 const reco::Vertex &PV = vertices->front();
 /////
 //   Initial skim
 /////
 /*
 int nmuevt = 0;
 for(const pat::Muon &mu : *muons) if(is_tth_loosemu(mu,PV)) nmuevt++; 
 int neleevt = 0;
 for(const pat::Electron &ele : *electrons) if(is_tth_looseele(ele,PV)) neleevt++;
 int njetevt   = 0;
 int nlbjetevt = 0;
 int nmbjetevt = 0;
 for(const pat::Jet &j : *jets){
  if(is_tth_jet(j)){
   njetevt++;
   if(j.bDiscriminator("combinedSecondaryVertexBJetTags")>0.244) nlbjetevt++;
   if(j.bDiscriminator("combinedSecondaryVertexBJetTags")>0.679) nmbjetevt++;
  }
 }
 if((nmuevt+neleevt)<2 || !(njetevt>=4 && (nlbjetevt>=2 || nmbjetevt>=1))) return;
 */
 /////
 //   Start the analysis 
 /////
 //Take the good lepton candidates
 vector<const reco::Candidate*> goodleps;
 //vector<const pat::PackedCandidate*> packedcandgoodleps;
 vector<int> goodcandpfcindices; //Save the index of the PackedCandidateCollection element associated to a good lep 
 /////
 //   Muons
 /////
 int mu_num = 0;
 vector<double> goodmu_pt;
 for(const pat::Muon &mu : *muons){
  //Select good mu
  if(!is_tth_loosemu(mu,PV)) continue;
  if(!is_cand_inpv(mindr_p3, mu, *pfs, goodcandpfcindices)) continue;
  //Save info
  goodleps.push_back((const reco::Candidate*)&mu);
  goodmu_pt.push_back(mu.pt());
  //packedcandgoodleps.push_back( &(*pfs)[goodcandpfcindices[mu_num]] );
  tree->mu_pt[mu_num]  = mu.pt();
  tree->mu_eta[mu_num] = mu.eta();
  tree->mu_phi[mu_num] = mu.phi();
  tree->mu_en[mu_num]  = mu.energy();
  tree->mu_charge[mu_num] = mu.charge();
  tree->mu_relisodbc[mu_num]    = rel_iso_dbc(mu);
  tree->mu_neurelisodbc[mu_num] = rel_iso_dbc(mu) - mu.chargedHadronIso()/mu.pt();
  tree->mu_chrelisodbc[mu_num]  = mu.chargedHadronIso()/mu.pt();
  tree->mu_chhadiso[mu_num]     = mu.chargedHadronIso();
  tree->mu_neuhadiso[mu_num]    = mu.neutralHadronIso();
  tree->mu_photoniso[mu_num]    = mu.photonIso();
  tree->mu_puchhadiso[mu_num]   = mu.puChargedHadronIso();
  const pat::PackedCandidate &pf = (*pfs)[goodcandpfcindices[mu_num]];
  Track trk = Track(pf.pseudoTrack());
  TransientTrack ttrk = ttrkbuilder->build(&trk);
  Measurement1D mu_3dpv = IPTools::absoluteImpactParameter3D(ttrk,PV).second;
  tree->mu_3dip[mu_num] = mu_3dpv.value()/mu_3dpv.error();
  mu_num++;
 }
 tree->mu_num = mu_num;
 /////
 //   Electrons 
 /////
 int ele_num = 0;
 vector<double> goodele_pt;
 for(const pat::Electron &ele : *electrons){
  //Select good ele
  if(!is_tth_looseele(ele,PV)) continue;
  bool matchelemu = false;
  for(uint gl=0; gl<goodleps.size(); gl++) if(deltaR(goodleps[gl]->p4(),ele.p4())<mindr_p3) matchelemu = true;  
  if(matchelemu) continue;
  if(!is_cand_inpv(mindr_p3, ele, *pfs, goodcandpfcindices)) continue;
  //Save info
  goodleps.push_back((const reco::Candidate*)&ele);
  goodele_pt.push_back(ele.pt());
  //packedcandgoodleps.push_back( &(*pfs)[goodcandpfcindices[mu_num+ele_num]] );
  tree->ele_pt[ele_num]  = ele.pt();
  tree->ele_eta[ele_num] = ele.eta();
  tree->ele_phi[ele_num] = ele.phi();
  tree->ele_en[ele_num]  = ele.energy();
  tree->ele_charge[ele_num] = ele.charge();
  tree->ele_relisodbc[ele_num]    = rel_iso_dbc(ele);
  tree->ele_neurelisodbc[ele_num] = rel_iso_dbc(ele) - ele.chargedHadronIso()/ele.pt();
  tree->ele_chrelisodbc[ele_num]  = ele.chargedHadronIso()/ele.pt();
  tree->ele_chhadiso[ele_num]     = ele.chargedHadronIso();
  tree->ele_neuhadiso[ele_num]    = ele.neutralHadronIso();
  tree->ele_photoniso[ele_num]    = ele.photonIso();
  tree->ele_puchhadiso[ele_num]   = ele.puChargedHadronIso();
  ele_num++;
 }
 tree->ele_num = ele_num;
 /////
 //   Jets
 /////
 //Now jet selection for the event
 int njet     = 0;
 int nlbjet   = 0;
 int nmbjet   = 0;
 for(const pat::Jet &j : *jets){
  if(!is_tth_jet(j)) continue;
  bool jetmatchedlepts = false;
  for(uint gl=0; gl<goodleps.size(); gl++) if(deltaR(goodleps[gl]->p4(),j.p4())<mindr_p5) jetmatchedlepts = true;
  if(jetmatchedlepts) continue;
  njet++;
  if(j.bDiscriminator("combinedSecondaryVertexBJetTags")>0.244) nlbjet++;
  if(j.bDiscriminator("combinedSecondaryVertexBJetTags")>0.679) nmbjet++;
 }//Loop on jets
 /////
 //   Event info
 /////
 if(!(goodleps.size()==2 && (goodleps[0]->charge()==goodleps[1]->charge()) && njet>=4 && (nlbjet>=2 || nmbjet>=1))) return;
 //if(!(goodleps.size()>1 && njet>=4 && (nlbjet>=2 || nmbjet>=1))) return;
 /////
 //   Lepton
 /////
 int lepnum = goodleps.size();
 tree->lep_num = lepnum; 
 sort(goodleps.begin(), goodleps.end(), [](const reco::Candidate* p1, const reco::Candidate* p2) {return p1->pt() > p2->pt();});
 //sort(packedcandgoodleps.begin(), packedcandgoodleps.end(), [](const pat::PackedCandidate* p1, const pat::PackedCandidate* p2) {return p1->pt() > p2->pt();});
 //Lepton
 for(int gl=0; gl<lepnum; gl++){
  //Flavor
  int lepismu = 0;
  for(uint gm=0; gm<goodmu_pt.size(); gm++) if(goodleps[gl]->pt()==goodmu_pt[gm]) lepismu = 1;   
  tree->lep_mu[gl] = lepismu; 
  int lepisele = 0;
  for(uint ge=0; ge<goodele_pt.size(); ge++) if(goodleps[gl]->pt()==goodele_pt[ge]) lepisele = 1;   
  tree->lep_ele[gl] = lepisele; 
  //Kinematics
  tree->lep_pt[gl]  = goodleps[gl]->pt();
  tree->lep_eta[gl] = goodleps[gl]->eta();
  tree->lep_phi[gl] = goodleps[gl]->phi();
  tree->lep_en[gl]  = goodleps[gl]->energy(); 
  tree->lep_charge[gl] = goodleps[gl]->charge();
  //Variables related to the jet associated to the lepton
  int lepjetindex = lepjet_index(mindr_p5, goodleps[gl], *jets); //Take the index of the jet associated to the lepton
  if(lepjetindex!=-1){//-1 is the default value that means not jet associated to leptons (see lepjetindex in interface/LeptonJetVariables.h) 
   const pat::Jet & lepjet = (*jets)[lepjetindex];
   GlobalVector lepjetgv(lepjet.px(),lepjet.py(),lepjet.pz());
   double lepmaxjettrkvtx3Daipsig = -1; //M(|IP3D/Delta_IP3D|(jettrks,PV)): the maximum significance of the 3D IP between a tracks of the jet and the PV
   double lepmaxjettrkvtx3Dsipsig = -1; //M(|sIP3D/Delta_sIP3D|(jettrks,PV)): as M(|IP3D/Delta_IP3D|(jettrks,PV)), but signed 
   double lepmaxjettrkvtx2Daipsig = -1; //M(|IP2D/Delta_IP2D|(jettrks,PV)): as M(|IP3D/Delta_IP3D|(jettrks,PV)), but 2D (x,y plane)
   double lepmaxjettrkvtx2Dsipsig = -1; //M(|sIP2D/Delta_sIP2D|(jettrks,PV)): as M(|IP2D/Delta_IP2D|(jettrks,PV)), but signed 
   int lepjetndaus    = 0; //# of jet daughters
   int lepjetnchtrks  = 0; //# of jet charged tracks
   int lepjetntrkspv  = 0; //# of jet tracks from PV 
   int lepjetntrksnpv = 0; //# of jet tracks not from PV 
   lepjet_variables(lepjet,PV,lepjetgv,*ttrkbuilder,
    lepmaxjettrkvtx3Daipsig, lepmaxjettrkvtx3Dsipsig, lepmaxjettrkvtx2Daipsig, lepmaxjettrkvtx2Dsipsig,                       
    lepjetndaus, lepjetnchtrks, lepjetntrkspv, lepjetntrksnpv);
   tree->lep_maxjettrkvtx3Daipsig[gl] = lepmaxjettrkvtx3Daipsig;
   tree->lep_maxjettrkvtx3Dsipsig[gl] = lepmaxjettrkvtx3Dsipsig;
   tree->lep_maxjettrkvtx2Daipsig[gl] = lepmaxjettrkvtx2Daipsig;
   tree->lep_maxjettrkvtx2Dsipsig[gl] = lepmaxjettrkvtx2Dsipsig;
   tree->lep_jetndaus[gl]    = lepjetndaus;
   tree->lep_jetnchtrks[gl]  = lepjetnchtrks;
   tree->lep_jetntrkspv[gl]  = lepjetntrkspv;
   tree->lep_jetntrksnpv[gl] = lepjetntrksnpv;
  }
 }
 tree->evt_id   = (unsigned int)iEvent.id().event(); 
 //tree->evt_json = //Fill when json will be used
 tree->evt_lumi = (unsigned int)iEvent.id().luminosityBlock();
 tree->evt_run  = (unsigned int)iEvent.id().run();
 /////
 //   Fill tree
 /////
 tree->tree->Fill();
 /////
 //   End of event
 /////
#ifdef THIS_IS_AN_EVENT_EXAMPLE
 Handle<ExampleData> pIn;
 iEvent.getByLabel("example",pIn);
#endif
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
 ESHandle<SetupData> pSetup;
 iSetup.get<SetupRecord>().get(pSetup);
#endif
}
// ------------ method called once each job just before starting event loop  ------------
void 
ObjEvtOptimisation::beginJob()
{
 stopwatch->Start();
 evt_totnum = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ObjEvtOptimisation::endJob() 
{
 stopwatch->Stop(); 
 cout<<endl;
 cout<<"Rapid job summary "<<endl;
 cout<<evt_totnum<<" events analysed in "<<stopwatch->RealTime()<<" seconds"<<endl;
 cout<<endl;
}
// ------------ method called when starting to processes a run  ------------
/*
void 
ObjEvtOptimisation::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
ObjEvtOptimisation::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
ObjEvtOptimisation::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
ObjEvtOptimisation::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ObjEvtOptimisation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(ObjEvtOptimisation);
