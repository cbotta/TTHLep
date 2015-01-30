//New class
#ifndef TREEVARIABLES
#define TREEVARIABLES
/////
//   Headers
/////
#include <TTree.h>
#include <string>
#include <map>
#include <cmath>
#include <iostream>
/////
//   Constants
/////
#define DEF_SIZE1D 25
#define DEF_VAL_INT -9999
#define DEF_VAL_FLOAT -9999.0f
#define DEF_VAL_DOUBLE -9999.0d
#define FLOAT_EPS 0.0000001f
#define DOUBLE_EPS 0.0000001d
/////
//   Functions
/////
#define INIT_1DARRAY(x,n,y) for(int i=0;i<n;i++) {x[i]=y;}
#define INIT_2DARRAY(x,n,m,y) for(int i=0;i<n;i++) { for(int j=0;j<m;j++) { x[i][j]=y; } }
inline bool is_undef(int x) { return x==DEF_VAL_INT; };
inline bool is_undef(float x) { return fabs(x-DEF_VAL_FLOAT) < FLOAT_EPS; };
inline bool is_undef(double x) { return fabs(x-DEF_VAL_DOUBLE) < DOUBLE_EPS; }
/////
//   Class declaration
/////
class CTree{
public:
 CTree(TTree* _tree) { tree = _tree; };
 TTree* tree;
 /////
 //   Helper functions for accessing branches
 /////
 template <typename T>
 T get_address(const std::string name){
  auto* br = tree->GetBranch(name.c_str());
  if(br==0){
   std::cerr << "ERROR: get_address CTree " << "branch " << name << " does not exist" << std::endl;
   throw std::exception();
  }
  auto* p = br->GetAddress();
  return reinterpret_cast<T>(p);
 }
 /////
 //   Declare variables
 /////
 //Event
 int evt_id, evt_json, evt_lumi, evt_run;
 //Muon
 int mu_num;
 int mu_charge[DEF_SIZE1D];
 double mu_pt[DEF_SIZE1D], mu_eta[DEF_SIZE1D], mu_phi[DEF_SIZE1D], mu_en[DEF_SIZE1D];
 double mu_relisodbc[DEF_SIZE1D], mu_neurelisodbc[DEF_SIZE1D], mu_chrelisodbc[DEF_SIZE1D], mu_chhadiso[DEF_SIZE1D], mu_neuhadiso[DEF_SIZE1D], mu_photoniso[DEF_SIZE1D], mu_puchhadiso[DEF_SIZE1D];
 double mu_3dip[DEF_SIZE1D];
 //Electron
 int ele_num;
 int ele_charge[DEF_SIZE1D];
 double ele_pt[DEF_SIZE1D], ele_eta[DEF_SIZE1D], ele_phi[DEF_SIZE1D], ele_en[DEF_SIZE1D];
 double ele_relisodbc[DEF_SIZE1D], ele_neurelisodbc[DEF_SIZE1D], ele_chrelisodbc[DEF_SIZE1D], ele_chhadiso[DEF_SIZE1D], ele_neuhadiso[DEF_SIZE1D], ele_photoniso[DEF_SIZE1D], ele_puchhadiso[DEF_SIZE1D];
 //Lepton
 int lep_num;
 int lep_charge[DEF_SIZE1D];
 int lep_mu[DEF_SIZE1D], lep_ele[DEF_SIZE1D];
 double lep_pt[DEF_SIZE1D], lep_eta[DEF_SIZE1D], lep_phi[DEF_SIZE1D], lep_en[DEF_SIZE1D];
 double lep_maxjettrkvtx3Daipsig[DEF_SIZE1D], lep_maxjettrkvtx3Dsipsig[DEF_SIZE1D], lep_maxjettrkvtx2Daipsig[DEF_SIZE1D], lep_maxjettrkvtx2Dsipsig[DEF_SIZE1D];
 int lep_jetndaus[DEF_SIZE1D], lep_jetnchtrks[DEF_SIZE1D], lep_jetntrkspv[DEF_SIZE1D], lep_jetntrksnpv[DEF_SIZE1D];
 /////
 //   Initialise
 /////
 void loop_initialize(void){
  //Event
  evt_id   = DEF_VAL_INT;
  evt_json = DEF_VAL_INT;
  evt_lumi = DEF_VAL_INT;
  evt_run  = DEF_VAL_INT;
  //Muon
  mu_num = DEF_VAL_INT;
  INIT_1DARRAY(mu_charge,DEF_SIZE1D,DEF_VAL_INT);
  INIT_1DARRAY(mu_pt,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(mu_eta,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(mu_phi,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(mu_en,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(mu_neurelisodbc,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(mu_chrelisodbc,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(mu_chhadiso,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(mu_neuhadiso,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(mu_photoniso,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(mu_puchhadiso,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(mu_3dip,DEF_SIZE1D,DEF_VAL_DOUBLE);
  //Electron
  ele_num = DEF_VAL_INT;
  INIT_1DARRAY(ele_charge,DEF_SIZE1D,DEF_VAL_INT);
  INIT_1DARRAY(ele_pt,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(ele_eta,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(ele_phi,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(ele_en,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(ele_relisodbc,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(ele_neurelisodbc,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(ele_chrelisodbc,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(ele_chhadiso,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(ele_neuhadiso,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(ele_photoniso,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(ele_puchhadiso,DEF_SIZE1D,DEF_VAL_DOUBLE);
  //Lepton
  lep_num  = DEF_VAL_INT;
  INIT_1DARRAY(lep_charge,DEF_SIZE1D,DEF_VAL_INT); 
  INIT_1DARRAY(lep_mu,DEF_SIZE1D,DEF_VAL_INT); 
  INIT_1DARRAY(lep_ele,DEF_SIZE1D,DEF_VAL_INT); 
  INIT_1DARRAY(lep_pt,DEF_SIZE1D,DEF_VAL_DOUBLE); 
  INIT_1DARRAY(lep_eta,DEF_SIZE1D,DEF_VAL_DOUBLE); 
  INIT_1DARRAY(lep_phi,DEF_SIZE1D,DEF_VAL_DOUBLE); 
  INIT_1DARRAY(lep_en,DEF_SIZE1D,DEF_VAL_DOUBLE); 
  INIT_1DARRAY(lep_maxjettrkvtx3Daipsig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(lep_maxjettrkvtx3Dsipsig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(lep_maxjettrkvtx2Daipsig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(lep_maxjettrkvtx2Dsipsig,DEF_SIZE1D,DEF_VAL_DOUBLE);
  INIT_1DARRAY(lep_jetndaus,DEF_SIZE1D,DEF_VAL_INT);
  INIT_1DARRAY(lep_jetnchtrks,DEF_SIZE1D,DEF_VAL_INT);
  INIT_1DARRAY(lep_jetntrkspv,DEF_SIZE1D,DEF_VAL_INT);
  INIT_1DARRAY(lep_jetntrksnpv,DEF_SIZE1D,DEF_VAL_INT);
 } 
 /////
 //   Set branches
 /////
 void make_branches(void){
  //Event
  tree->Branch("evt_id", &evt_id, "evt_id/I");
  tree->Branch("evt_json", &evt_json, "evt_json/I");
  tree->Branch("evt_lumi", &evt_lumi, "evt_lumi/I");
  tree->Branch("evt_run", &evt_run, "evt_run/I");
  //Muon
  tree->Branch("mu_num", &mu_num, "mu_num/I");
  tree->Branch("mu_charge", &mu_charge, "mu_charge[mu_num]/I");
  tree->Branch("mu_pt", &mu_pt, "mu_pt[mu_num]/D");
  tree->Branch("mu_eta", &mu_eta, "mu_eta[mu_num]/D");
  tree->Branch("mu_phi", &mu_phi, "mu_phi[mu_num]/D");
  tree->Branch("mu_en", &mu_en, "mu_en[mu_num]/D");
  tree->Branch("mu_relisodbc", &mu_relisodbc, "mu_relisodbc[mu_num]/D");
  tree->Branch("mu_neurelisodbc", &mu_neurelisodbc, "mu_neurelisodbc[mu_num]/D");
  tree->Branch("mu_chrelisodbc", &mu_chrelisodbc, "mu_chrelisodbc[mu_num]/D");
  tree->Branch("mu_chhadiso", &mu_chhadiso, "mu_chhadiso[mu_num]/D");
  tree->Branch("mu_neuhadiso", &mu_neuhadiso, "mu_neuhadiso[mu_num]/D");
  tree->Branch("mu_photoniso", &mu_photoniso, "mu_photoniso[mu_num]/D");
  tree->Branch("mu_puchhadiso", &mu_puchhadiso, "mu_puchhadiso[mu_num]/D");
  tree->Branch("mu_3dip", &mu_3dip, "mu_3dip[mu_num]/D");
  //Electron
  tree->Branch("ele_num", &ele_num, "ele_num/I");
  tree->Branch("ele_charge", &ele_charge, "ele_charge[ele_num]/I");
  tree->Branch("ele_pt", &ele_pt, "ele_pt[ele_num]/D");
  tree->Branch("ele_eta", &ele_eta, "ele_eta[ele_num]/D");
  tree->Branch("ele_phi", &ele_phi, "ele_phi[ele_num]/D");
  tree->Branch("ele_en", &ele_en, "ele_en[ele_num]/D");
  tree->Branch("ele_relisodbc", &ele_relisodbc, "ele_relisodbc[ele_num]/D");
  tree->Branch("ele_neurelisodbc", &ele_neurelisodbc, "ele_neurelisodbc[ele_num]/D");
  tree->Branch("ele_chrelisodbc", &ele_chrelisodbc, "ele_chrelisodbc[ele_num]/D");
  tree->Branch("ele_chhadiso", &ele_chhadiso, "ele_chhadiso[ele_num]/D");
  tree->Branch("ele_neuhadiso", &ele_neuhadiso, "ele_neuhadiso[ele_num]/D");
  tree->Branch("ele_photoniso", &ele_photoniso, "ele_photoniso[ele_num]/D");
  tree->Branch("ele_puchhadiso", &ele_puchhadiso, "ele_puchhadiso[ele_num]/D");
  //Lepton
  tree->Branch("lep_num", &lep_num, "lep_num/I");
  tree->Branch("lep_charge", &lep_charge, "lep_charge[lep_num]/I");
  tree->Branch("lep_mu", &lep_mu, "lep_mu[lep_num]/I");
  tree->Branch("lep_ele", &lep_ele, "lep_ele[lep_num]/I");
  tree->Branch("lep_pt", &lep_pt, "lep_pt[lep_num]/D");
  tree->Branch("lep_eta", &lep_eta, "lep_eta[lep_num]/D");
  tree->Branch("lep_phi", &lep_phi, "lep_phi[lep_num]/D");
  tree->Branch("lep_en", &lep_en, "lep_en[lep_num]/D");
  tree->Branch("lep_maxjettrkvtx3Daipsig", &lep_maxjettrkvtx3Daipsig, "lep_maxjettrkvtx3Daipsig[lep_num]/D");
  tree->Branch("lep_maxjettrkvtx3Dsipsig", &lep_maxjettrkvtx3Dsipsig, "lep_maxjettrkvtx3Dsipsig[lep_num]/D");
  tree->Branch("lep_maxjettrkvtx2Daipsig", &lep_maxjettrkvtx2Daipsig, "lep_maxjettrkvtx2Daipsig[lep_num]/D");
  tree->Branch("lep_maxjettrkvtx2Dsipsig", &lep_maxjettrkvtx2Dsipsig, "lep_maxjettrkvtx2Dsipsig[lep_num]/D");
  tree->Branch("lep_jetndaus", &lep_jetndaus, "lep_jetndaus[lep_num]/I");
  tree->Branch("lep_jetnchtrks", &lep_jetnchtrks, "lep_jetnchtrks[lep_num]/I");
  tree->Branch("lep_jetntrkspv", &lep_jetntrkspv, "lep_jetntrkspv[lep_num]/I");
  tree->Branch("lep_jetntrksnpv", &lep_jetntrksnpv, "lep_jetntrksnpv[lep_num]/I");
 }
 /////
 //   Set branch address
 /////
 //Connects the branches of an existing TTree to variables used when loading the file
 void set_branch_addresses(void){
  //Event
  tree->SetBranchAddress("evt_id", &evt_id);
  tree->SetBranchAddress("evt_json", &evt_json);
  tree->SetBranchAddress("evt_lumi", &evt_lumi);
  tree->SetBranchAddress("evt_run", &evt_run);
  //Muon
  tree->SetBranchAddress("mu_num", &mu_num);
  tree->SetBranchAddress("mu_charge", &mu_charge);
  tree->SetBranchAddress("mu_pt", &mu_pt);
  tree->SetBranchAddress("mu_eta", &mu_eta);
  tree->SetBranchAddress("mu_phi", &mu_phi);
  tree->SetBranchAddress("mu_en", &mu_en);
  tree->SetBranchAddress("mu_relisodbc", &mu_relisodbc);
  tree->SetBranchAddress("mu_neurelisodbc", &mu_neurelisodbc);
  tree->SetBranchAddress("mu_chrelisodbc", &mu_chrelisodbc);
  tree->SetBranchAddress("mu_chhadiso", &mu_chhadiso);
  tree->SetBranchAddress("mu_neuhadiso", &mu_neuhadiso);
  tree->SetBranchAddress("mu_photoniso", &mu_photoniso);
  tree->SetBranchAddress("mu_puchhadiso", &mu_puchhadiso);
  tree->SetBranchAddress("mu_3dip", &mu_3dip);
  //Electron
  tree->SetBranchAddress("ele_num", &ele_num);
  tree->SetBranchAddress("ele_charge", &ele_charge);
  tree->SetBranchAddress("ele_pt", &ele_pt);
  tree->SetBranchAddress("ele_eta", &ele_eta);
  tree->SetBranchAddress("ele_phi", &ele_phi);
  tree->SetBranchAddress("ele_en", &ele_en);
  tree->SetBranchAddress("ele_relisodbc", &ele_relisodbc);
  tree->SetBranchAddress("ele_neurelisodbc", &ele_neurelisodbc);
  tree->SetBranchAddress("ele_chrelisodbc", &ele_chrelisodbc);
  tree->SetBranchAddress("ele_chhadiso", &ele_chhadiso);
  tree->SetBranchAddress("ele_neuhadiso", &ele_neuhadiso);
  tree->SetBranchAddress("ele_photoniso", &ele_photoniso);
  tree->SetBranchAddress("ele_puchhadiso", &ele_puchhadiso);
  //Lepton
  tree->SetBranchAddress("lep_num", &lep_num);
  tree->SetBranchAddress("lep_charge", &lep_charge);
  tree->SetBranchAddress("lep_mu", &lep_mu);
  tree->SetBranchAddress("lep_ele", &lep_ele);
  tree->SetBranchAddress("lep_pt", &lep_pt);
  tree->SetBranchAddress("lep_eta", &lep_eta);
  tree->SetBranchAddress("lep_phi", &lep_phi);
  tree->SetBranchAddress("lep_en", &lep_en);
  tree->SetBranchAddress("lep_maxjettrkvtx3Daipsig", &lep_maxjettrkvtx3Daipsig);
  tree->SetBranchAddress("lep_maxjettrkvtx3Dsipsig", &lep_maxjettrkvtx3Dsipsig);
  tree->SetBranchAddress("lep_maxjettrkvtx2Daipsig", &lep_maxjettrkvtx2Daipsig);
  tree->SetBranchAddress("lep_maxjettrkvtx2Dsipsig", &lep_maxjettrkvtx2Dsipsig);
  tree->SetBranchAddress("lep_jetndaus", &lep_jetndaus);
  tree->SetBranchAddress("lep_jetnchtrks", &lep_jetnchtrks);
  tree->SetBranchAddress("lep_jetntrkspv", &lep_jetntrkspv);
  tree->SetBranchAddress("lep_jetntrksnpv", &lep_jetntrksnpv);
 }
};
#endif
