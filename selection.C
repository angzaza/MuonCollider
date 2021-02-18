//macro to select good muons and events 
//Pt resolution of selected muons is also evaluated

//   good muons:
//   eta<2.5
//   Pt>5 GeV
//   D0<2 mm
//   Z0 <10 mm


#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "Math/Vector4D.h"
#include "TRotation.h"

void selection(){

// declaration of file paths
	TChain* c_sig = new TChain("c_sig");
  c_sig->Add("path/ntuple_mue_Zmumuvv_10000_newrel_hits.root/MyLCTuple");
  
  ofstream myfile;
  myfile.open ("selection_output_bkg1_1_5_tev.txt");
  
  
  
  
  
// cross section and luminosity
  //double xs_sig=9.29e-3;
  double xs=0.0797;
  double cm_energy=1.5;
  double lumi=500 ;
  
// VARIABLES DECLARATION
  int n_events_sig=0; // number of generated events 
  int n_events_bkg=0;
  
  int n=0;
  int n_rc_muon=0; //number of reco muons for each event
  int counter_eta_mu=0, counter_pt_mu=0, counter_d0_mu=0, counter_z0_mu=0; //counter for good muons
  int count_eta_mu_evt=0, count_pt_mu_evt=0, count_d0_mu_evt=0, count_z0_mu_evt=0, count_mup_evt=0, count_mum_evt=0; //counter for the number of muons in each events passing the selections
  int counter_eta_evt=0, counter_pt_evt=0, counter_d0_evt=0, counter_z0_evt=0, counter_OS_evt=0; //counter for events with at least 4 muons passing the corresponing selections
  int counter_0mu_evt=0, counter_1mu_evt=0, counter_2mu_evt=0, counter_3mu_evt=0, counter_4mu_evt=0, counter_5mu_evt=0, counter_6mu_evt=0, counter_6mup_evt=0; //counter for events with 0,1,2,3,4,5,6,>6 reco muons
  
  
  
  int counter_dR=0, counter_pt10=0, counter_pt20=0;
  int condition_4mu_sel=0;
  int condition_dR=0, condition_pt=0;
  int condition_Z=0, condition_Z1=0, condition_ZZ=0, condition_mass_4mu=0;
  int counter_h=0;
  
  int counter_dR_evt=0, counter_2pt_evt=0, counter_ZZ_evt=0, counter_Z1_evt=0, counter_mass_evt=0;
  
  int za_index=0,zb_index=0,z1_index=0,z2_index=0;
  
  TLorentzVector rcmu,mcmu,z1,z2,higgs;
  TLorentzVector mup[4];
  TLorentzVector mum[4];
  TLorentzVector mup_sel[2];
  TLorentzVector mum_sel[2];
  
  double mass_z=91.1876; 
  double deltaR;
  double delta_m1=0,delta_m2=0,delta_m=0;
  double delta_min=500, delta_min_2=500, pt_max=-1;
  double inv_mass_z1, inv_mass_z2, inv_mass_H, pt_H;
  double inv_mass[4]={0};
  double inv_mass_4mu=0;
  double pt_sum_z2=0;
// VARIABLES FROM BRANCHES
   Int_t n_mcp=0;
   Float_t *mc_px;
  Float_t *mc_py;
  Float_t *mc_pz;
  Float_t *mc_e;
  Float_t *mc_vtx;
  Float_t *mc_vty;
  Float_t *mc_vtz;
  Int_t *mc_pdg;
  Int_t *mc_pa0;
  Int_t *mc_gst;
  Int_t *match_trk;
  
  int num=800000;
  mc_px = (float *) malloc(sizeof(float) * num);
  mc_py = (float *) malloc(sizeof(float) * num);
  mc_pz = (float *) malloc(sizeof(float) * num);
  mc_e = (float *) malloc(sizeof(float) * num);
  mc_vtx = (float *) malloc(sizeof(float) * num);
  mc_vty = (float *) malloc(sizeof(float) * num);
  mc_vtz = (float *) malloc(sizeof(float) * num);
  mc_pdg = (int *) malloc(sizeof(int) * num);
  mc_pa0 = (int *) malloc(sizeof(int) * num);
  mc_gst = (int *) malloc(sizeof(int) * num);
  match_trk = (int *) malloc(sizeof(int) * num);
  
  Float_t trk_d0[30000]={0};
  Float_t trk_curv[30000]={0.};
  Float_t trk_phi[30000]={0};
  Float_t trk_z0[30000]={0};
  Float_t trk_tanl[30000];
  Int_t trk_atIP[30000]={0};
  Int_t tr_fts[30000]={0};
  Int_t trk_n=0;
  
  Int_t n_rec;
  Int_t rc_typ[30000];
  Int_t rc_ntr[30000];
  Int_t rc_ftr[30000];
  Float_t rc_mox[30000];
  Float_t rc_moy[30000];
  Float_t rc_moz[30000];
  Float_t rc_ene[30000];
  Float_t rc_cha[1000];
  
    Int_t r2m_nrel;
  Int_t r2m_t[1000];
  Int_t r2m_f[1000];
  Float_t r2m_w[1000];
  
  ////MC particles branches/////////////////////////////////
  c_sig->SetBranchAddress("mcmox", mc_px);
  c_sig->SetBranchAddress("mcmoy", mc_py);
  c_sig->SetBranchAddress("mcmoz", mc_pz);
  c_sig->SetBranchAddress("mcpdg", mc_pdg);
  c_sig->SetBranchAddress("mcpa0", mc_pa0);
  c_sig->SetBranchAddress("mcvtx", mc_vtx);
  c_sig->SetBranchAddress("mcvty", mc_vty);
  c_sig->SetBranchAddress("mcvtz", mc_vtz);
  c_sig->SetBranchAddress("mcene", mc_e);
  c_sig->SetBranchAddress("nmcp", &n_mcp);
  
  /// reco tracks branches/////////////////////////////
  c_sig->SetBranchAddress("trsip", trk_atIP);
  c_sig->SetBranchAddress("tsdze", trk_d0);
  c_sig->SetBranchAddress("tsphi", trk_phi);
  c_sig->SetBranchAddress("tsome", trk_curv);
  c_sig->SetBranchAddress("tszze", trk_z0);
  c_sig->SetBranchAddress("tstnl", trk_tanl);
  c_sig->SetBranchAddress("trfts", tr_fts);
  c_sig->SetBranchAddress("ntrk", &trk_n);
  
  //reco particles branches///////////////
  c_sig->SetBranchAddress("nrec", &n_rec);
  c_sig->SetBranchAddress("rctyp", rc_typ);
  c_sig->SetBranchAddress("rcntr", rc_ntr);
  c_sig->SetBranchAddress("rcmox", rc_mox);
  c_sig->SetBranchAddress("rcmoy", rc_moy);
  c_sig->SetBranchAddress("rcmoz", rc_moz);
  c_sig->SetBranchAddress("rcene", rc_ene);
  c_sig->SetBranchAddress("rcftr", rc_ftr);
  c_sig->SetBranchAddress("rccha", rc_cha);
  
   //relation branches///////////////////
  c_sig->SetBranchAddress("r2mnrel", &r2m_nrel);
  c_sig->SetBranchAddress("r2mf", r2m_f);
  c_sig->SetBranchAddress("r2mt", r2m_t);
  c_sig->SetBranchAddress("r2mw", r2m_w);
  
 
  
  // VARIABLES DECLARATION FOR PT RESOLUTION
  	double pT_bins[7]={5,20,30,40,60,100,200};
	double eta_bins[51]={0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5};
	double res_bins[101];
	for(int i=0;i<101;i++){
	res_bins[i]=-0.11+i*0.0022;
	}

   TH1F* deltaPt_histo = new TH1F("deltaPt","deltaPt",100,-3,3);
    TH3F* Pt_res_pt = new TH3F("res pt eta=0","res pt eta=0",6,pT_bins,50,eta_bins,100,res_bins);
   TH2F* Pt_res_eta= new TH2F("res #eta","res #eta",50,0,2.5,80,-0.14,0.14);
   
   TH2F* Pt_res_pt_1;
   TH2F* Pt_res_pt_2;
   TH2F* Pt_res_pt_3;

   double ptBin_array_1[6]={12.5,25,35,50,80,150};
   double ptErr_array_1[6]={7.5,5,5,10,20,50};
   double sigma_array_1[6]={0};
   double sigmaErr_array_1[6]={0};
   double ptBin_array_2[6]={12.5,25,35,50,80,150};
   double ptErr_array_2[6]={7.5,5,5,10,20,50};
   double sigma_array_2[6]={0};
   double sigmaErr_array_2[6]={0};
   double ptBin_array_3[6]={12.5,25,35,50,80,150};
   double ptErr_array_3[6]={7.5,5,5,10,20,50};
   double sigma_array_3[6]={0};
   double sigmaErr_array_3[6]={0};
    double ptBin_array_4[6]={12.5,25,35,50,80,150};
   double ptErr_array_4[6]={7.5,5,5,10,20,50};
   double sigma_array_4[6]={0};
   double sigmaErr_array_4[6]={0};
   
   
   double etaBin_array[54]={0};
   double etaErr_array[54]={0};
   double sigma_eta_array[54]={0};
   double sigmaErr_eta_array[54]={0};
   
   double stdBin_array1[20]={0};
   double stdErr_array1[20]={0};
   
   TLorentzVector mc_mu,rc_mu;
  
   double deltaPhi=0,deltaPt=0;
  // int i=0,k=0,n=0;
   
  
   
   TF1 *g[20];
   TF1 *g_eta[20];
   TF1 *gaus_fit=new TF1("gaus","gaus",-0.05,0.05);
   gaus_fit->SetParameter(1,0.5);
   gaus_fit->SetParameter(2,0.5);
   gaus_fit->SetParLimits(1,-1,1);
   gaus_fit->SetParLimits(2,0.001,1);
   char *title; // or val[3]
   
    TF1 *gaus_fit2=new TF1("gaus","gaus",-0.09,0.09);
		gaus_fit2->SetParameter(1,0.5);
   gaus_fit2->SetParameter(2,0.5);
   gaus_fit2->SetParLimits(1,-1,1);
   gaus_fit2->SetParLimits(2,0.001,1);
   
   ///////////////////////////////////////////////////////////////////////////////////////////////////////
   
   //////////// HISTOGRAMS //////////////////////////////////////////
  TH1F* dR_4mu_histo=new TH1F("#DeltaR betwwen each of the 4 muons","#DeltaR between each of the 4 muons",30,0,7.5);
  TH1F* invM_Z1=new TH1F("Z1 mass","Z1 mass",75,0,150);
  TH1F* invM_Z2=new TH1F("Z2 mass","Z2 mass",75,0,150);
  TH1F* invM_H=new TH1F("H mass","H mass",200,50,2000);
  TH1F* pt_H_histo=new TH1F("H Pt","H Pt",100,0,500);
  

  
  /////// CUTS //////////////////////////////////////////////
  // for muons//////
  double eta_cut=2.5;
  double pt_cut=5.;
  double d0_cut=2.;
  double z0_cut=10.;
  
  //for events ////
  double dR_cut=0.02;
  double Pt_1=10.;
  double Pt_2=20.;
  double inv_mass_Z_low=12.;
  double inv_mass_Z_high=120.;
  double inv_mass_Z1_cut=40.;
  double inv_mass_4mu_cut=70.;
  ////////////////////////////////////////////////////////////////////////
  

  //***************** SIGNAL SELECTION ********************************************
  myfile << "mu+e-->Zmumuvv"<<endl;
  myfile << "SELECTION"<<endl;
  myfile<<"cm energy(TeV): "<< cm_energy << "   lumi(fb): "<<lumi<<endl;
  myfile<<" "<<endl;
  myfile<< "xs: "<<xs<<endl;
  myfile<<"events: "<<c_sig->GetEntries()<<endl;
  myfile<<"MUONS SELECTION"<<endl;
  
  n_events_sig=c_sig->GetEntries(); 
   
  //cycle on events
  for(unsigned int ientry=0; ientry<c_sig->GetEntries(); ++ientry){
  	c_sig->GetEntry(ientry);
   	n_rc_muon=0;
   	condition_4mu_sel=0; //this variable is switched to one if the event contains at least 4 muons passing the selection
   	
   	// selection of muons /////////////////////////////////
   	
   	count_eta_mu_evt=0; count_pt_mu_evt=0; count_d0_mu_evt=0; count_z0_mu_evt=0; //counters for the number of muons that pass the selections for each event -> they are set to 0 for each event
   	count_mup_evt=0; count_mum_evt=0;
   	
   	//cycle on reco particles
   	for(int i_rc=0;i_rc<n_rec;i_rc++){
   		//condition on muons with an associated track
   		if(abs(rc_typ[i_rc])==13 && rc_ntr[i_rc]==1){
   			n_rc_muon++;
   			rcmu.SetPxPyPzE(rc_mox[i_rc],rc_moy[i_rc],rc_moz[i_rc],rc_ene[i_rc]);
   			
   			if(abs(rcmu.Eta())<eta_cut){
   				counter_eta_mu++;
   				count_eta_mu_evt++;
   				
   				if(rcmu.Pt()>pt_cut){
   					counter_pt_mu++;
   					count_pt_mu_evt++;
   					
   					if(abs(trk_d0[tr_fts[rc_ftr[i_rc]]])<d0_cut){
   						counter_d0_mu++;
   						count_d0_mu_evt++;
   						
   						if(abs(trk_z0[tr_fts[rc_ftr[i_rc]]])<z0_cut){
   							counter_z0_mu++;
   							count_z0_mu_evt++;
   							
   							
   							//  cycle on relations between RECO and GEN particles in order to evaluate Pt resolution
   							for (unsigned int u=0; u<r2m_nrel; ++u){
     							if(r2m_f[u]==i_rc){
     								mc_mu.SetPxPyPzE(mc_px[r2m_t[u]],mc_py[r2m_t[u]],mc_pz[r2m_t[u]],mc_e[r2m_t[u]]);
										deltaPt=mc_mu.Pt()-rc_mu.Pt();
										deltaPt_histo->Fill(deltaPt);
										//cout<<"deltaPT: "<<deltaPt<<endl;
										//cout<<"eta: "<<mc_mu.Eta()<<endl;
										Pt_res_eta->Fill(abs(rc_mu.Eta()),deltaPt/rc_mu.Pt());
										Pt_res_pt->Fill(rc_mu.Pt(),abs(rc_mu.Eta()),deltaPt/rc_mu.Pt());
     							}
     						}
     						//////////////////////////////////////////////////////////////////////////////////////
     						
     						// muons TLorentz vectors will be stored in two different arrays according to the charge sign
     						if(rc_cha[i_rc]>0){
     							mup[count_mup_evt].SetPxPyPzE(rc_mox[i_rc],rc_moy[i_rc],rc_moz[i_rc],rc_ene[i_rc]);
     							count_mup_evt++;
     						}
   							else if(rc_cha[i_rc]<0){
   								mum[count_mum_evt].SetPxPyPzE(rc_mox[i_rc],rc_moy[i_rc],rc_moz[i_rc],rc_ene[i_rc]);
   								count_mum_evt++;
   							}
   							
   						} // end condition on z0
   					} // end condition on d0
   			  } //end condition on pt
   	  	} //end condition on eta
   	} //end selection of muons
   } // end cycle on reco particles
   
  
   
   // selection of events //////////////////////////////////////
  
  if(n_rc_muon>6){
   		counter_6mup_evt++;
   	}
  else if(n_rc_muon==6){
   		counter_6mu_evt++;
   	}
  else if(n_rc_muon==5){
   		counter_5mu_evt++;
   	}
  else if(n_rc_muon==4){
   		counter_4mu_evt++;
   	}
  else if(n_rc_muon==3){
   		counter_3mu_evt++;
   	}
  else if(n_rc_muon==2){
   		counter_2mu_evt++;
   	}
  else if(n_rc_muon==1){
   		counter_1mu_evt++;
   	}
  else if(n_rc_muon==0){
   		counter_0mu_evt++;
   	}
  
  
  if(count_eta_mu_evt>=4){
  	counter_eta_evt++;
  	if(count_pt_mu_evt>=4){
  		counter_pt_evt++;
  		if(count_d0_mu_evt>=4){
  			counter_d0_evt++;
  			if(count_z0_mu_evt>=4){
  				counter_z0_evt++;
  				if(count_mup_evt>=2 && count_mum_evt>=2){
   					counter_OS_evt++;	
 			 			condition_4mu_sel=1;
 					}
  			}
  		}
  	}
  }
  
 
  
  
  if(condition_4mu_sel==1){
  
  counter_dR=0; counter_pt10=0; counter_pt20=0;
  condition_dR=0; condition_pt=0; condition_Z=0; condition_Z1=0; condition_mass_4mu=0; condition_ZZ=0;
  
  delta_min_2=500.;
  pt_max=-1.;
  
  //in general we have count_mup_evt positive muons and count_mum_evt negative muons,  with count_mup_evt>=4 and count_mum_evt>=4
  //in order to reconstruct the ZZ candidates 4 muons combinations are constructed with the following cycles in i,j,k,l
  for(int i=0;i<count_mup_evt;i++){
  	for (int j=0;j<count_mum_evt;j++){
    	for(int k=i;k<count_mup_evt;k++){
     		for (int l=0;l<count_mum_evt;l++){	
     			if(k!=j && l!=i && k!=i && l!=j){	
     				counter_dR=0; counter_pt10=0; counter_pt20=0;
     				
     				mup_sel[0]=mup[i];
     				mup_sel[1]=mup[k];
     				mum_sel[0]=mum[j];
     				mum_sel[1]=mum[l];	
     				
     				//ZZ candidates are required to have deltaR between each of the four leptons > 0.02
     				for(int a=0;a<2;a++){
     					//deltaR between same sign muons
     					for(int b=a+1;b<2;b++){
     						//deltaR=sqrt((mum_sel[a].Eta()-mum_sel[b].Eta())*(mum_sel[a].Eta()-mum_sel[b].Eta())+(mum_sel[a].Phi()-mum_sel[b].Phi())*(mum_sel[a].Phi()-mum_sel[b].Phi()));
     						deltaR=mum_sel[a].DeltaR(mum_sel[b]);
     						dR_4mu_histo->Fill(deltaR);
     						if(deltaR>dR_cut){ //condition on deltaR between negative sign muons
   								counter_dR++;
   							}
   							//deltaR=sqrt((mup_sel[a].Eta()-mup_sel[b].Eta())*(mup_sel[a].Eta()-mup_sel[b].Eta())+(mup_sel[a].Phi()-mup_sel[b].Phi())*(mup_sel[a].Phi()-mup_sel[b].Phi()));
   							deltaR=mup_sel[a].DeltaR(mup_sel[b]);
   							dR_4mu_histo->Fill(deltaR);
   							if(deltaR>dR_cut){ //condition on deltaR between positive sign muons
   								counter_dR++;
   							}
     					}
     					//deltaR between opposite sign muons
     					for(int b=0;b<2;b++){
     						//deltaR=sqrt((mum_sel[a].Eta()-mup_sel[b].Eta())*(mum_sel[a].Eta()-mup_sel[b].Eta())+(mum_sel[a].Phi()-mup_sel[b].Phi())*(mum_sel[a].Phi()-mup_sel[b].Phi()));
     						deltaR=mum_sel[a].DeltaR(mup_sel[b]);
     						dR_4mu_histo->Fill(deltaR);
   							if(deltaR>dR_cut){ //condition on deltaR between opposite sign muons
   								counter_dR++;
   							}
     					}
     				}
     				
     				if(counter_dR==6){
     					condition_dR++;
     					
     					for(int a=0;a<2;a++){
   							if((mup_sel[a].Pt())>Pt_2){counter_pt20++;} //20
   							if((mup_sel[a].Pt())>Pt_1){counter_pt10++;}
   							if((mum_sel[a].Pt())>Pt_1){counter_pt10++;}
   							if((mum_sel[a].Pt())>Pt_2){counter_pt20++;} //20
   						}
   							
   							//selected events should have at least 2 muons with Pt,i>10 GeV and Pt,j>20 GeV
   							if(counter_pt10>=2 && counter_pt20>=1){
   								condition_pt++;
   								
   								for(int a=0;a<2;a++){
   									for(int b=0;b<2;b++){	
   										inv_mass[b+2*a]=(mum_sel[a]+mup_sel[b]).M(); //invariant mass of Z candidates
   									}
   								}
   								//non-overlapping ZZ invariant mass candidates are inv_mass[0]&inv_mass[3]  and inv_mass[1]&inv_mass[2
   								
   								delta_min=500.;
   								condition_mass_4mu=0;
   								//cycle on non-overlapping ZZ candidates
   								for(n=0;n<2;n++){
   									//Z candidates are required to have 12<invMass<120 GeV
   									if((inv_mass[n]>inv_mass_Z_low && inv_mass[n]<inv_mass_Z_high) && (inv_mass[3-n]>inv_mass_Z_low && inv_mass[3-n]<inv_mass_Z_high)){  //selection ZZ candidates 
   										condition_Z++; //condition_Z is increased by 1 if both Z candidates have mass between 12 and 120 GeV
   										delta_m1=abs(inv_mass[n]-mass_z);
   										delta_m2=abs(inv_mass[3-n]-mass_z);
   										za_index=n;
   										zb_index=3-n;
   										delta_m=delta_m1;
   										if(delta_m2<delta_m1){
   											za_index=3-n;
   											zb_index=n;
   											delta_m=delta_m2;
   										}
   										//Z1 candidates are required to have an invariant mass larger than 40 GeV
   										if(inv_mass[za_index]>inv_mass_Z1_cut){
   											condition_Z1++; //condition_Z1 is increased by 1 if the Z1 candidate has a mass greater than 40 GeV
   											if(delta_m<delta_min){
   												delta_min=delta_m;
   												z1_index=za_index;
   												z2_index=zb_index;
   												pt_sum_z2=mum[(int)(z2_index/2)].Pt()+mup[(int)(z2_index%2)].Pt();
   												//ZZcandidates are required to have 4mu invariant mass > 70 GeV
   												if(((mum_sel[0]+mum_sel[1]+mup_sel[0]+mup_sel[1]).M())>70){
   													condition_mass_4mu=1; //condition_mass_4mu is switched to 1 if the invariant mass of the selected 4 muons is larger than 70 GeV
   												}
   											}
   										}
   									}
   								} // end cycle on the ZZ combinations 
   								
   						
   								if(condition_mass_4mu==1){
   									condition_ZZ++; //condition_ZZ is increased by 1 if there is at least one ZZ candidate passing the selection
   									if(fabs(inv_mass[z1_index]-mass_z)<delta_min_2){
   										delta_min_2=fabs(inv_mass[z1_index]-mass_z);
   										inv_mass_z1=inv_mass[z1_index];
   										inv_mass_z2=inv_mass[z2_index];
   										inv_mass_H=(mum_sel[0]+mum_sel[1]+mup_sel[0]+mup_sel[1]).M();
   										pt_H=(mum_sel[0]+mum_sel[1]+mup_sel[0]+mup_sel[1]).Pt();
   												
   									}
   									else if(fabs(inv_mass[z1_index]-mass_z)==delta_min_2 && pt_sum_z2>pt_max){
   										delta_min_2=fabs(inv_mass[z1_index]-mass_z);
   										pt_max=pt_sum_z2;
   										inv_mass_z1=inv_mass[z1_index];
   										inv_mass_z2=inv_mass[z2_index];
   										inv_mass_H=(mum_sel[0]+mum_sel[1]+mup_sel[0]+mup_sel[1]).M();
   									}
   								}
   						
   								
   								
   							}// end condition on pt
     					
     				}//end condition on dR
     				
     				
     				
     				
     			}
     		}
     	}
    }
  }
  
  if(condition_ZZ>=1){
  	invM_Z1->Fill(inv_mass_z1);
   	invM_Z2->Fill(inv_mass_z2);
   	invM_H->Fill(inv_mass_H);
    pt_H_histo->Fill(pt_H);
   			
		if(inv_mass_H>100 && inv_mass_H<140){
			counter_h++;
		}
  }
  
  
  if(condition_dR>=1){
  	counter_dR_evt++;
  	if(condition_pt>=1){
  		counter_2pt_evt++;
  		if(condition_Z>=1){
  			counter_ZZ_evt++;
  			if(condition_Z1>=1){
  				counter_Z1_evt++;
  				if(condition_ZZ>=1){
  					counter_mass_evt++;
  				}
  			}
  	}
  }
  
  }
  
  
  	
  }
  
 
  

}//end cycle on events

 myfile<<"eta < "<<eta_cut<<" :  "<<counter_eta_mu<<endl;
 myfile<<"pt > "<<pt_cut<<" :  "<<counter_pt_mu<<endl;
 myfile<<"d0 < "<<d0_cut<<" :  "<<counter_d0_mu<<endl;
 myfile<<"z0 < "<<z0_cut<<" :  "<<counter_z0_mu<<endl<<endl;
 myfile<<" EVENTS SELECTION"<<endl;
 
 myfile<<"evts with 0 rc muons: "<<counter_0mu_evt<<endl;
  myfile<<"evts with 1 rc muons: "<<counter_1mu_evt<<endl;
  myfile<<"evts with 2 rc muons: "<<counter_2mu_evt<<endl;
  myfile<<"evts with 3 rc muons: "<<counter_3mu_evt<<endl;
  myfile<<"evts with 4 rc muons: "<<counter_4mu_evt<<endl;
  myfile<<"evts with 5 rc muons: "<<counter_5mu_evt<<endl;
  myfile<<"evts with 6 rc muons: "<<counter_6mu_evt<<endl;
  myfile<<"evts with rc muons > 6: "<<counter_6mup_evt<<endl<<endl;
  
   myfile<<"evts with at least 4mu passing eta sel: "<<counter_eta_evt<<endl;
  myfile<<"evts with at least 4mu passing pt sel: "<<counter_pt_evt<<endl;
  myfile<<"evts with at least 4mu passing d0 sel: "<<counter_d0_evt<<endl;
  myfile<<"evts with at least 4mu passing z0 sel: "<<counter_z0_evt<<endl;
  myfile<<"evts with at least 2mu+ and 2mu-: "<<counter_OS_evt<<endl;
  
   

myfile<<"number of events with : "<<endl;
myfile<<"dR>"<<dR_cut<<" beteen each of the 4 muons: "<<counter_dR_evt<<endl;
myfile<<"at least 2 muons with Pt,i>20 and Pt,j>10: "<<counter_2pt_evt<<endl;
myfile<<"at least 1 ZZ candidate with 12<Mass(Z)>120 GeV: "<<counter_ZZ_evt<<endl;
myfile<<"at least 1 ZZ candidate with Mass(Z1)>40 GeV: "<<counter_Z1_evt<<endl;
myfile<<"the invariant mass of the 4 muons>70 GeV: "<<counter_mass_evt<<endl;
myfile<<"Higgs mass between 100 and 140 GeV: "<<counter_h<<endl;

myfile.close();

TFile *rootFile = new TFile("Selection_bkg1_1_5tev_plots.root","RECREATE");
invM_Z1->Write();
invM_Z2->Write();
invM_H->Write();
pt_H_histo->Write();




}
