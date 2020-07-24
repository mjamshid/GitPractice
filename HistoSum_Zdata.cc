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
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TMath.h"
#include "BParticle.h"
#include "BEvent.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooChebychev.h"
#include "RooPlot.h"
#include "RooDataHist.h"

using namespace RooFit;
void HistoSum_Zdata(){
 gROOT->Reset();

//  gStyle->SetOptTitle(kFALSE);
//  gStyle->SetOptStat(000000000);

 // TFile *f0 = TFile::Open("hist.root");
 TFile *f1 = TFile::Open("hist_data.root");

/*
  f0->cd();
 gDirectory->ls();
 TH1F* h_z = (TH1F*)f0->Get("h11");
*/
// TH1F* h_pt1 = (TH1F*)f0->Get("h5");
 f1->cd();
 gDirectory->ls();
 TH1F* d_z = (TH1F*)f1->Get("k11");
 // TH1F* d_pt1 = (TH1F*)f1->Get("k5");
 
// Assigning Histogram dataset to di_mu object

RooDataHist data("data", "dataset", di_mu ,d_z);

// Assigning new variable 
// RooRealVar x("x", "x", 60,120);

// Signal Model and PArameters

RooRealVar ms("ms", "ms", 91, 80,100);
RooRealVar w("w","w", 5);
RooGaussian G("G", "G", di_mu, ms, w);

// Background Model and parameters

RooRealVar mb("mb", "mb", 91, 80, 100);
RooRealVar k("k", "k", 5);
RooArgusBG A("A", "A", di_mu, mb, k);

// Composite Model and Parameters
/*
RooRealVar f("f", "signal fraction", 0,5);
RooAddPdf M("M", "G+A", RooArgList(G,A), f);
*/

// Fitting to data

G.fitTo(data);
A.fitTo(data);

 
//  TLegend *leg =  new TLegend(.70,.70,.80,.80);
 TCanvas *c5 = new TCanvas("c5", "c5");
 c5->cd();
 c5->SetTicks(1, 1);
/*
f0->cd(); 
 h_z->Draw("HIST");
Double_t norm = 1;
Double_t scale = norm/(h_z->Integral());
h_z->Scale(scale);
 h_z->SetLineColor(kBlue);
 h_z->SetLineWidth(3);
*/


/*  h_pt1->Draw("SAME");
 h_pt1->SetLineColor(kRed);
 h_pt1->SetLineWidth(3); */


//  leg->AddEntry(h_z, "Z_mass_MC" , "L");
//  leg->AddEntry(h_pt1, "mu1pT_MC", "L");


 d_z->Draw("HIST");
Double_t norm0 = 1;
Double_t scale0 = norm0/(d_z->Integral());
d_z->Scale(scale0);
 d_z->SetLineColor(kMagenta);
 d_z->SetLineWidth(4);
d_z->GetXaxis()->SetTitle("#mu#mu Mass (GeV)");
d_z->GetYaxis()->SetTitle("Events");


// Plotting

 RooPlot* frame = di_mu.frame();
 data.plotOn(frame);
 G.plotOn(frame, LineColor(kBlue));
 A.plotOn(frame, LineColor(kSpring), LineStyle(kDashed));
frame->Draw();

 d_z->SetTitle("Z-Boson Mass");
 
 /* d_pt1->Draw("SAME");
 d_pt1->SetLineColor(kMagenta);
 d_pt1->SetLineWidth(4); */
 

//  leg->AddEntry(d_z, "Z_mass_Data", "L");

//  leg->AddEntry(d_pt1, "mu1pT_Data", "L");

 /*
 leg->SetBorderSize(0);
 leg->SetTextSize(0.035);
 leg->SetTextFont(40);
 leg->Draw();
*/
c5->SaveAs("Z1_mass.pdf");
}  
