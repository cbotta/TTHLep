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
KalmanVertexFitter vtxFitter(true);
//Namespaces
using namespace reco;
using namespace edm;
using namespace std;
//Return the index of the jet associated to the lepton
int lepjet_index(double mindr_objetmatch, const reco::Candidate* cand, const pat::JetCollection& jets){
 int index = -1;
 double mindr = mindr_objetmatch;
 int currjetpos = 0;
 for(const pat::Jet &j : jets){
  if(deltaR(j.p4(),cand->p4())<mindr){
   mindr = deltaR(j.p4(),cand->p4());
   index = currjetpos;
  }
  currjetpos++;
 }
 return index; 
}
//Access lepton jet information related to jet trks IP and trk multiplicity
void lepjet_variables(const pat::Jet& lepjet, const reco::Vertex& vtx, GlobalVector lepjetgv, const TransientTrackBuilder& ttrkbuilder,
 double& lepmaxjettrkvtx3Daipsig, double& lepmaxjettrkvtx3Dsipsig, double& lepmaxjettrkvtx2Daipsig, double& lepmaxjettrkvtx2Dsipsig,
 int& lepjetndaus, int& lepjetnchtrks, int& lepjetntrkspv, int& lepjetntrksnpv){
 //Declare following values to get max IP values 
 Measurement1D aip3Dtljpv; bool is_aip3Dtljpv = false; double maxaip3Dtljpv = -9999;
 Measurement1D sip3Dtljpv; bool is_sip3Dtljpv = false; double maxsip3Dtljpv = -9999;
 Measurement1D aip2Dtljpv; bool is_aip2Dtljpv = false; double maxaip2Dtljpv = -9999; 
 Measurement1D sip2Dtljpv; bool is_sip2Dtljpv = false; double maxsip2Dtljpv = -9999; 
 //Access jet daughters
 vector<CandidatePtr> jdaus(lepjet.daughterPtrVector());
 sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});
 lepjetndaus = jdaus.size();
 for(int jd=0; jd<lepjetndaus; jd++){
  const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
  //Minimal conditions for a track 
  if(jcand.charge()!=0 && jcand.numberOfHits()>0){
   lepjetnchtrks++; 
   //Take maximum IP values 
   Track trk = Track(jcand.pseudoTrack());
   TransientTrack ttrk = ttrkbuilder.build(&trk);   
   pair<bool, Measurement1D> currIP = IPTools::absoluteImpactParameter3D(ttrk,vtx);
   if(currIP.second.value()>maxaip3Dtljpv){
    maxaip3Dtljpv = currIP.second.value();
    aip3Dtljpv    = currIP.second;
    is_aip3Dtljpv = currIP.first; 
   }
   currIP = IPTools::signedImpactParameter3D(ttrk,lepjetgv,vtx);
   if(currIP.second.value()>maxsip3Dtljpv){
    maxsip3Dtljpv = currIP.second.value();
    sip3Dtljpv    = currIP.second;
    is_sip3Dtljpv = currIP.first; 
   }
   currIP = IPTools::absoluteTransverseImpactParameter(ttrk,vtx);
   if(currIP.second.value()>maxaip2Dtljpv){
    maxaip2Dtljpv = currIP.second.value();
    aip2Dtljpv    = currIP.second;
    is_aip2Dtljpv = currIP.first; 
   }
   currIP = IPTools::signedTransverseImpactParameter(ttrk,lepjetgv,vtx);
   if(currIP.second.value()>maxsip2Dtljpv){
    maxsip2Dtljpv = currIP.second.value();
    sip2Dtljpv    = currIP.second;
    is_sip2Dtljpv = currIP.first; 
   }
   //Other conditions on jet daughters
   if(jcand.fromPV()==pat::PackedCandidate::PVUsedInFit){
    lepjetntrkspv++;
   }else{
    lepjetntrksnpv++;
   }
  } 
 }
 //Save max IP values 
 if(is_aip3Dtljpv && aip3Dtljpv.value()>=0 && aip3Dtljpv.error()>=0) lepmaxjettrkvtx3Daipsig = aip3Dtljpv.value()/aip3Dtljpv.error();
 if(is_sip3Dtljpv && aip3Dtljpv.value()>=0 && aip3Dtljpv.error()>=0) lepmaxjettrkvtx3Dsipsig = sip3Dtljpv.value()/sip3Dtljpv.error();
 if(is_aip2Dtljpv && aip2Dtljpv.value()>=0 && aip2Dtljpv.error()>=0) lepmaxjettrkvtx2Daipsig = aip2Dtljpv.value()/aip2Dtljpv.error(); 
 if(is_sip2Dtljpv && aip2Dtljpv.value()>=0 && aip2Dtljpv.error()>=0) lepmaxjettrkvtx2Dsipsig = sip2Dtljpv.value()/sip2Dtljpv.error(); 
}
