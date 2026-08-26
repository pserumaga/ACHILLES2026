#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TKey.h"
#include "TCollection.h"
#include <iostream>
#include <string>

TH1D* getHist(TFile* f, const std::vector<std::string>& candidates) {
    for (const auto& name : candidates) {
        TH1D* h = (TH1D*)f->Get(name.c_str());
        if (h) return h;
    }
    std::cerr << "  [warn] none of the candidate names matched directly, "
                 "scanning file contents...\n";
    return nullptr;
}

void printAvailableKeys(TFile* f) {
    std::cerr << "  Available histogram keys in " << f->GetName() << ":\n";
    TIter next(f->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        std::string cname = key->GetClassName();
        if (cname.find("TH1") != std::string::npos)
            std::cerr << "    " << key->GetName() << "\n";
    }
}

TH1D* getRatio(TFile* f, const char* tag) {
    TH1D* hQE = getHist(f, {"hIn_nu_E_qe"});
    TH1D* hIntf = getHist(f, {"hIn_nu_E_intf"});

    if (!hQE || !hIntf) {
        std::cerr << "[error] Could not find QE and/or Interference "
                     "hIn_nu_E histograms in " << f->GetName() << "\n";
        printAvailableKeys(f);
        return nullptr;
    }

    TH1D* ratio = (TH1D*)hIntf->Clone(Form("hRatio_IntfQE_%s", tag));
    ratio->SetTitle(Form("Interference/QE vs E^{in}_{#nu} (%s); "
                          "E^{in}_{#nu} [MeV];Interference / QE", tag));
    ratio->Sumw2();
    ratio->Divide(hQE);
    return ratio;
}

void plot_ratio_intf_qe(const char* sbndFile = "achilles_Ar_SBNDFLUX_output.root",
                         const char* duneFile = "achilles_Ar_DUNEFLUX_output.root",
                         const char* outFile  = "ratio_intf_qe.root") {

    TFile* fSBND = TFile::Open(sbndFile, "READ");
    TFile* fDUNE = TFile::Open(duneFile, "READ");
    if (!fSBND || fSBND->IsZombie()) { std::cerr << "Cannot open " << sbndFile << "\n"; return; }
    if (!fDUNE || fDUNE->IsZombie()) { std::cerr << "Cannot open " << duneFile << "\n"; return; }

    TH1D* ratioSBND = getRatio(fSBND, "SBND");
    TH1D* ratioDUNE = getRatio(fDUNE, "DUNE");
    if (!ratioSBND || !ratioDUNE) {
        std::cerr << "[error] Aborting - fix histogram names above and rerun.\n";
        return;
    }

    // Build a comparison canvas (kept as a TCanvas object inside the file so it can be reopened/edited directly in ROOT
    ratioSBND->SetStats(0);
    ratioDUNE->SetStats(0);  // if it also overlays this one
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("cRatioComparison", "Interference/QE ratio", 800, 600);
    ratioSBND->SetLineColor(kOrange+1);
    ratioSBND->SetLineWidth(2);
    ratioDUNE->SetLineColor(kBlue+1);
    ratioDUNE->SetLineWidth(2);
    ratioDUNE->SetLineStyle(2);

    ratioSBND->Draw("HIST E");
    ratioDUNE->Draw("HIST E SAME");

    TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->AddEntry(ratioSBND, "SBND flux", "l");
    leg->AddEntry(ratioDUNE, "DUNE flux", "l");
    leg->SetBorderSize(0);
    leg->Draw();

    // Write everything to the output ROOT file: the two ratio histograms
    // individually (fully editable/re-drawable) plus the canvas.
    TFile* fOut = TFile::Open(outFile, "RECREATE");
    ratioSBND->Write();
    ratioDUNE->Write();
    c->Write();
    fOut->Close();

    std::cout << "Wrote ratio histograms + canvas to " << outFile << "\n";

    fSBND->Close();
    fDUNE->Close();
}
