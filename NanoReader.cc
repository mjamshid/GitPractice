#include "ostream"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"

void NanoReader(){
    TFile *file = new TFile("/afs/cern.ch/user/m/mjamshid/RunIISummer16NanoAODv6_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_NANOAODSIM_PUMoriond17_Nano25Oct2019_102X_mcRun2_asymptotic_v7_ext1-v1.root");
    TTreeReader reader("Events", file);
    TTreeReaderArray<Float_t>   Muon_pt(reader, "Muon_pt"    );
    TTreeReaderArray<Float_t>  Muon_eta(reader, "Muon_eta"   );
    TTreeReaderArray<Float_t>  Muon_phi(reader, "Muon_phi"   );
    TTreeReaderArray<Float_t> Muon_mass(reader, "Muon_mass"  );
    TTreeReaderArray<UInt_t>      nMuon(reader, "nMuon"      );
    TTreeReaderArray<Int_t> Muon_charge(reader, "Muon_charge");

    Int_t nEvents = reader.GetEntries(0);
    Double_t ptCut = 5;
    Int_t nMuonCut = 2;
    Int_t count = 0;
    while(reader.Next()){
        int n = nMuon[0];
        Float_t ptMax = Muon_pt[0];
        if(ptMax <  ptCut || n < nMuonCut ){continue;}
        count++;
        std::cout << "==============================================" << std::endl;
        std::cout << "Event Number: " << reader.GetCurrentEntry() << std::endl;
        std::cout << "Number of Muons: " << n << std::endl;
        TLorentzVector pMuon[n];
        for (Int_t i=0; i<n; i++) {
           pMuon[i].SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],Muon_mass[i]);
            std::cout << Muon_charge[i]     << "       ";
            std::cout << Muon_pt[i]         << "       ";
            std::cout << pMuon[i].E()       << "       ";
            std::cout << pMuon[i].Px()      << "       ";
            std::cout << pMuon[i].Py()      << "       ";
            std::cout << pMuon[i].Pz()      << std::endl;
        }
    }
    Float_t eff = float(count)/float(nEvents);
    std::cout << "==============================================" << std::endl;
    std::cout << "Number of events on the dataset: " << nEvents << std::endl;
    std::cout << "Number of events that passed the cut Pt > " << ptCut
              << " GeV and number of muons = "<< nMuonCut <<": "<< count << std::endl;
    std::cout << "Cut efficiency: " << eff << std::endl;
    std::cout << "==============================================" << std::endl;
}
