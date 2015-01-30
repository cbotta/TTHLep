//Jet and Met
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
//Math
#include <algorithm>
//Namespaces
using namespace reco;
using namespace edm;
using namespace std;
/////
//   Lepton ID 
/////
//Ask candidate to come from PV
bool is_cand_inpv(const double mindR, const reco::Candidate& cand, const pat::PackedCandidateCollection& pcc, vector<int>& goodcandpfcindices){
 bool iscandinpv = false;
 double mindr = mindR;
 double goodcandpfcindex = -1;
 for(uint i=0; i<pcc.size(); i++){
  const pat::PackedCandidate cpf = pcc[i];
  if(deltaR(cand.p4(),cpf.p4())<mindr && //dR is the standard geometrical way to associate 
    (fabs(cand.pt()-cpf.pt())/cand.pt())<0.05 && //Check in pT, because ele,tau are usually faked by jets (many trks) and dR may not be enough
    cpf.charge()!=0 && cpf.numberOfHits()>0 && //Leptons are charged and built from tracks, also to be consistent with PV tracks  
    cpf.fromPV()==pat::PackedCandidate::PVUsedInFit //Coming from PV
    ){
   bool wasnotagoodcandindex = false; //Check that cpf was not already used for association
   for(uint vgi=0; vgi<goodcandpfcindices.size(); vgi++) if(goodcandpfcindices[vgi]==int(i)) wasnotagoodcandindex = true;
   if(!wasnotagoodcandindex){
    iscandinpv = true;
    mindr = deltaR(cand.p4(),cpf.p4());
    goodcandpfcindex = int(i);
   } 
  }
 }
 if(goodcandpfcindex!=-1) goodcandpfcindices.push_back(goodcandpfcindex);
 return iscandinpv;
}
//Relative isolation with delta-beta correction
double rel_iso_dbc(const pat::Muon& lepton){
 return ( (lepton.chargedHadronIso() + 
           std::max(0.0, lepton.neutralHadronIso() + lepton.photonIso() - 0.5*lepton.puChargedHadronIso()))/lepton.pt() );
}
double rel_iso_dbc(const pat::Electron& lepton){
 return ( (lepton.chargedHadronIso() + 
           std::max(0.0, lepton.neutralHadronIso() + lepton.photonIso() - 0.5*lepton.puChargedHadronIso()))/lepton.pt() );
}
/////
//   Muon ID
/////
bool is_tth_loosemu(const pat::Muon& mu, const reco::Vertex& vtx){
 return( mu.isLooseMuon() && mu.pt()>5 && fabs(mu.eta())<2.5 && mu.muonBestTrack().isNonnull() &&
  fabs(mu.muonBestTrack()->dxy(vtx.position()))<0.05 && fabs(mu.muonBestTrack()->dz(vtx.position()))<0.2 &&
  rel_iso_dbc(mu)<0.4 );
}
/////
//   Ele ID
/////
bool is_tth_looseele(const pat::Electron& ele, const reco::Vertex& vtx){
 return( ele.electronID("eidLoose")>0.5 && ele.pt()>7 && fabs(ele.eta())<2.5 && ele.gsfTrack().isNonnull() &&
  fabs(ele.gsfTrack()->dxy(vtx.position()))<0.05 && fabs(ele.gsfTrack()->dz(vtx.position()))<0.2 &&
  ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=1 &&
  rel_iso_dbc(ele)<0.4 );
}
/////
//   Jet id
/////
//Namespace to identify the pileup nature of a jet 
//working points from RecoJets/JetProducers/python/PileupJetIDCutParams_cfi.py
//Eta Categories  0-2.5 2.5-2.75 2.75-3.0 3.0-5.0
namespace pu_mva{
 float full_chs_loose[4][4] = {
  //full_5x_chs_wp Loose id
  {-0.98,-0.95,-0.94,-0.94}, //pt 0-10
  {-0.98,-0.95,-0.94,-0.94}, //pt 10-20
  {-0.89,-0.77,-0.69,-0.75}, //pt 20-30
  {-0.89,-0.77,-0.69,-0.57}, //pt 30-50
 };
 int eta_idx(float eta){//Note that in principle this is not needed as we always require jet to be within eta<2.4
  const float ae = TMath::Abs(eta);
  if(ae < 2.5){
   return 0;
  }else if(ae < 2.75){
   return 1;
  }else if(ae < 3.0){
   return 2;
  }else if(ae < 5.0) {
   return 3;
  }else{
   return 3;
  }
 }
 int pt_idx(float pt){
  if(pt < 10){
   return 0;
  }else if(pt < 20){
   return 1;
  }else if(pt < 30){
   return 2;
  }else if(pt < 50){
   return 3;
  }else{
   //edm::LogWarning("jet_pu_id") << "pt outside range " << pt;
   return -1;
  }
 }
 bool pass_id(const pat::Jet& x, float mva) {
  int vpt_idx  = pt_idx(x.pt());
  int veta_idx = eta_idx(x.eta());
  //hard jet
  if(vpt_idx==-1) return true;
  if(mva>full_chs_loose[vpt_idx][veta_idx]) return true;
  return false;
 }
}
bool is_tth_jet(const pat::Jet &j){
 PFJetIDSelectionFunctor pfLooseJetID(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE),
                         pfTightJetID(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT);
 pat::strbitset passLooseCuts(pfLooseJetID.getBitTemplate()),
                passTightCuts(pfTightJetID.getBitTemplate());
 return( j.pt()>25 && fabs(j.eta())<2.5 && pfLooseJetID(j,passLooseCuts) && pu_mva::pass_id(j, j.userFloat("pileupJetId:fullDiscriminant")) );
}
