#include "ostream"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TMath.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TChain.h"
#include "Fit/FitResult.h"

void NanoReader1(){
 gROOT->Reset();
    //TFile *hist_sv = new TFile("hist.root", "Recreate"); 
    TFile *file = new TFile("/eos/user/m/mjamshid/MY_ANALYSIS/RunIISummer16NanoAODv6_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_NANOAODSIM_PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1.root "); 
   
   //TFile *file = new TFile("/eos/user/m/mjamshid/MY_ANALYSIS/F9A90ED4-23E8-664E-97FD-205F8072210A.root");
    TTreeReader reader("Events", file);
    TTreeReaderArray<Float_t>   Muon_pt(reader, "Muon_pt"    );
    TTreeReaderArray<Float_t>  Muon_eta(reader, "Muon_eta"   );
    TTreeReaderArray<Float_t>  Muon_phi(reader, "Muon_phi"   );
    TTreeReaderArray<Float_t> Muon_mass(reader, "Muon_mass"  );
    TTreeReaderArray<UInt_t>      nMuon(reader, "nMuon"      );
    TTreeReaderArray<Int_t> Muon_charge(reader, "Muon_charge");
    TTreeReaderArray<Bool_t> Muon_softId(reader, "Muon_softId");
    TTreeReaderValue<Bool_t> HLT_IsoMu24(reader, "HLT_IsoMu24");
//    TTreeReaderValue<Bool_t> HLT_IsoMu27(reader, "HLT_IsoMu27");
    Int_t nEvents = reader.GetEntries(0);
    Double_t ptCut = 5;
    Int_t nMuonCut = 2;
    Int_t count = 0;
// Defining Histograms
    
    TH1F *h1 = new TH1F("h1", "Muon0_pt", 50, 0,150);
    TH1F *h2 = new TH1F("h2", "Muon0_eta", 50, -10, 10);
    TH1F *h3 = new TH1F("h3", "Muon0_phi", 51, -6,6);
    TH1F *h4 = new TH1F("h4", "Muon0_charge", 50, -10,10);
    TH1F *h5 = new TH1F("h5", "Muon1_pt", 50, 0,150);
    TH1F *h6 = new TH1F("h6", "Muon1_eta", 50, -10, 10);
    TH1F *h7 = new TH1F("h7", "Muon1_phi", 51, -6,6);
    TH1F *h8 = new TH1F("h8", "Muon1_charge", 50, -10,10);
    TH1F *h9 = new TH1F("h9", "Muon0_mass", 20, 0,1);
    TH1F *h10 = new TH1F("h10", "Muon1_mass", 20, 0,1);
    TH1F *h11 = new TH1F("h11", "di-#mu Invariant Mass", 50, 60,120);
    while(reader.Next()){

        int n = nMuon[0];
  //   Double_t mu_id = Muon_softId[0];
        Float_t ptMax = Muon_pt[0];
        if(ptMax <  ptCut || n< nMuonCut ){continue;}
        count++;
        std::cout << "==============================================" << std::endl;
        std::cout << "Event Number: " << reader.GetCurrentEntry() << std::endl;
        std::cout << "Number of Muons: " << n << std::endl;
        TLorentzVector pMuon[n];
        for (Int_t i=0; i<2; i++) {
            pMuon[i].SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],Muon_mass[i]);

          /*std::cout << Muon_charge[i]     << "       ";
            std::cout << Muon_pt[i]         << "       ";
            std::cout << pMuon[i].E()       << "       ";
            std::cout << pMuon[i].Px()      << "       ";
            std::cout << pMuon[i].Py()      << "       ";
            std::cout << pMuon[i].Pz()      << std::endl; */
        //  Muon_softId[i];
         if(Muon_softId[i] == false) {continue;}
         count++;
	
	 std::cout << " Muons of Soft ID = " << Muon_softId[i] <<std::endl;

        }
// Accessing to pMuon[0] and pMuon[1] since they are the highest pt Muons in the event
      std::cout<<" Muon 0 PT = "<<Muon_pt[0]<<" Muon 0 eta =  "<<Muon_eta[0]<<" Muon 0 Phi = "<<Muon_phi[0]<<" Muon 0 Charge = "
      <<Muon_charge[0]<<std::endl;
      std::cout<<" Muon 1 PT = "<<Muon_pt[1]<<" Muon 1 eta = "<<Muon_eta[1]<<" Muon 1 Phi = "<<Muon_phi[1]<<" Muon 1 Charge = "
      <<Muon_charge[1]<<std::endl;
// Invariant Masses of Muon[0] & Muon[1]
  /*    std::cout<<" Muon 0 Px = "<< pMuon[0].Px()<<" Muon 0 Py = "<<pMuon[0].Py()<<" Muon 0 Pz = "<<pMuon[0].Pz()<<" Muon 0 E = "          <<pMuon[0].E()<<std::endl;
      std::cout<<" Muon 0 invariant Mass = "<<pMuon[0].Mag()<<std::endl;
      std::cout<<" Muon 1 Px = "<< pMuon[1].Px()<<" Muon 1 Py = "<<pMuon[1].Py()<<" Muon 1 Pz = "<<pMuon[1].Pz()<<" Muon 1 E = "
      <<pMuon[1].E()<<std::endl;
      std::cout<<" Muon 1 invariant Mass = "<<pMuon[1].Mag()<<std::endl;
 */
// Invariant mass of Di-Muon
    //  std::cout<<" Di-Muon Invariant Mass = "<<(pMuon[0]+pMuon[1]).Mag()<<std::endl;
    

  Double_t dimu = (pMuon[0]+pMuon[1]).Mag();
    Double_t cut0 = 60;
    Double_t cut1 = 120;
  Double_t zchar = Muon_charge[0]+Muon_charge[1] ;
    Double_t cut2 = 0;
  Double_t high_pt0 = Muon_pt[0];
  Double_t high_pt1 = Muon_pt[1];
 Double_t cut4 = 26;
 Double_t cut5 = 3;
 Bool_t softId = Muon_softId[0] && Muon_softId[1];
  if(dimu < cut0 || dimu > cut1 || zchar != cut2 || high_pt0 < cut4 ||  high_pt1 < cut5 || *HLT_IsoMu24 == false || softId == false) {continue;}
  count++;

std::cout<< " Largest pT of Muon = "<< high_pt0 <<std::endl;



  std::cout<< " Di_Muon Invariant Mass from Z decay = "<< dimu <<std::endl;
  std::cout << " Electric charge of Z = "<< zchar<<std::endl;
  std::cout << " Pass HLT_IsoMu24 ? = " << *HLT_IsoMu24 <<std::endl;
  std::cout << " Muons of Soft ID = " << Muon_softId[i] <<std::endl;
//  std::cout<< " Largest pT of Muon = "<< high_pt <<std::endl;

// Fiilling Histograms
      h1->Fill(Muon_pt[0]);
      h2->Fill(Muon_eta[0]);
      h3->Fill(Muon_phi[0]);
      h4->Fill(Muon_charge[0]);
      h5->Fill(Muon_pt[1]);
      h6->Fill(Muon_eta[1]);
      h7->Fill(Muon_phi[1]);
      h8->Fill(Muon_charge[1]);
      h9->Fill(pMuon[0].Mag());
      h10->Fill(pMuon[1].Mag());
    //  h11->Fill((pMuon[0]+pMuon[1]).Mag());
        h11 ->Fill(dimu);
    }

    Float_t eff = float(count)/float(nEvents);
    std::cout << "==============================================" << std::endl;
    std::cout << "Number of events on the dataset: " << nEvents << std::endl;
    std::cout << "Number of events that passed the cut Pt > " << ptCut
              << " GeV and number of muons = "<< nMuonCut <<": "<< count << std::endl;
    std::cout << "Cut efficiency: " << eff << std::endl;
    std::cout << "==============================================" << std::endl;
   
 TFile *hist_sv = new TFile("hist.root", "RECREATE");
  // Defining Canvases
    TCanvas *c1 = new TCanvas("c1", "c1");
    c1->cd();
    h1->Draw();
    h1->Write();
c1->SaveAs("mu_pt0mc.root");
    TCanvas *c2 = new TCanvas("c2", "c2");
    c2->cd();
    h2->Draw();
    h2->Write();
c2->SaveAs("mu_eta0mc.root");

   TCanvas *c3 = new TCanvas("c3", "c3");
    c3->cd();
    h3->Draw();
    h3->Write();
c3->SaveAs("mu_phi0mc.root");

    TCanvas *c4 = new TCanvas("c4", "c4");
    c4->cd();
    h4->Draw();
     h4->Write();

c4->SaveAs("mu_charge0mc.root");
 

  TCanvas *c5 = new TCanvas("c5", "c5");
    c5->cd();
    h5->Draw();
    h5->Write();
c5->SaveAs("mu_pt1mc.root");

    TCanvas *c6 = new TCanvas("c6", "c6");
    c6->cd();
    h6->Draw();
     h6->Write();
c6->SaveAs("mu_eta1mc.root");


   TCanvas *c7 = new TCanvas("c7", "c7");
    c7->cd();
    h7->Draw();
     h7->Write();
c7->SaveAs("mu_phi1mc.root");


    TCanvas *c8 = new TCanvas("c8", "c8");
    c8->cd();
    h8->Draw();
    h8->Write();
c8->SaveAs("mu_charge1mc.root");


TCanvas *c9 = new TCanvas("c9", "c9");
    c9->cd();
    h9->Draw();
    h9->Write();
c9->SaveAs("mu_0_massmc.root");


    TCanvas *c10 = new TCanvas("c10", "c10");
    c10->cd();
    h10->Draw();
    h10->Write();
c10->SaveAs("mu_1_massmc.root");



 TCanvas *c11 = new TCanvas("c11", "di-muon-mass");
    c11->cd();
    h11->Draw();
    h11->Write();
c11->SaveAs("di_mu-massmc.root");


hist_sv->Write();
hist_sv->Close();	


}
