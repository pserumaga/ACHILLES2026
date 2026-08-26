// plot_prism_onoff_ratio.C
// Ratio of the on-axis PRISM bin (0.0-0.2 deg) to the far off-axis bin (1.4-1.6 deg),
// for QE-only and for QE+Interference combined, as a function of outgoing muon energy.
//
// CAVEAT: TH1::Divide()'s default error propagation assumes the numerator and
// denominator are statistically independent. Here they're not - both bins come from
// the same underlying flat-flux event sample, just reweighted differently - so the
// drawn error bars are an upper-bound estimate, not a rigorous uncertainty. Fine for
// exploratory/poster use; flag to your supervisor before using these errors in the
// paper, since a proper treatment needs the covariance between the two reweightings.

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include <string>
#include <iostream>
#include <algorithm>

namespace {

TH1D* GetHistFlex(TFile* f, const std::string& base, const std::string& cat) {
  if (!f) return nullptr;
  const std::string name_us = base + "_" + cat;
  const std::string name_sp = base + " " + cat;
  if (auto h = dynamic_cast<TH1D*>(f->Get(name_us.c_str()))) return h;
  if (auto h = dynamic_cast<TH1D*>(f->Get(name_sp.c_str()))) return h;
  return nullptr;
}

} // anon

void plot_prism_onoff_ratio(const char* onaxisFile     = "prism_outputs/output_flux_oaa_numu_0_0_0_2.root",
                             const char* faroffaxisFile = "prism_outputs/output_flux_oaa_numu_1_4_1_6.root",
                             const char* outdir         = "PRISMplots",
                             const char* histbase       = "hEmu")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gSystem->mkdir(outdir, kTRUE);

  TFile* fOn  = TFile::Open(onaxisFile);
  TFile* fOff = TFile::Open(faroffaxisFile);
  if (!fOn || fOn->IsZombie())   { std::cerr << "ERROR: cannot open " << onaxisFile << "\n"; return; }
  if (!fOff || fOff->IsZombie()) { std::cerr << "ERROR: cannot open " << faroffaxisFile << "\n"; return; }

  TH1D* onQE    = GetHistFlex(fOn,  histbase, "qe");
  TH1D* onIntf  = GetHistFlex(fOn,  histbase, "intf");
  TH1D* offQE   = GetHistFlex(fOff, histbase, "qe");
  TH1D* offIntf = GetHistFlex(fOff, histbase, "intf");

  if (!onQE || !offQE) {
    std::cerr << "ERROR: missing " << histbase << "_qe in one of the input files.\n";
    fOn->Close(); fOff->Close();
    return;
  }

  // QE-only ratio: on-axis / far-off-axis
  TH1D* ratioQE = (TH1D*)onQE->Clone("ratioQE");
  ratioQE->SetDirectory(nullptr);
  ratioQE->Divide(offQE);

  // QE+Intf ratio: build combined histograms first, then divide
  TH1D* ratioQEIntf = nullptr;
  if (onIntf && offIntf) {
    TH1D* onCombo = (TH1D*)onQE->Clone("onCombo");
    onCombo->SetDirectory(nullptr);
    onCombo->Add(onIntf);

    TH1D* offCombo = (TH1D*)offQE->Clone("offCombo");
    offCombo->SetDirectory(nullptr);
    offCombo->Add(offIntf);

    ratioQEIntf = (TH1D*)onCombo->Clone("ratioQEIntf");
    ratioQEIntf->SetDirectory(nullptr);
    ratioQEIntf->Divide(offCombo);
  } else {
    std::cerr << "WARNING: " << histbase << "_intf not found in one of the files; QE+Intf ratio skipped.\n";
  }

  ratioQE->SetLineColor(kAzure + 2);
  ratioQE->SetLineWidth(3);
  ratioQE->SetMarkerColor(kAzure + 2);
  ratioQE->SetMarkerStyle(20);

  if (ratioQEIntf) {
    ratioQEIntf->SetLineColor(kOrange + 7);
    ratioQEIntf->SetLineStyle(2);
    ratioQEIntf->SetLineWidth(3);
    ratioQEIntf->SetMarkerColor(kOrange + 7);
    ratioQEIntf->SetMarkerStyle(21);
  }

  double ymax = 1.25 * std::max(ratioQE->GetMaximum(), ratioQEIntf ? ratioQEIntf->GetMaximum() : 0.0);
  if (ymax <= 0) ymax = 2.0;

  TCanvas* c = new TCanvas("c_prism_onoff_ratio", "c_prism_onoff_ratio", 1000, 750);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.13);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.07);

  ratioQE->GetXaxis()->SetTitle("E_{#mu}^{out} [MeV]");
  ratioQE->GetXaxis()->SetRangeUser(0, 4000);
  ratioQE->GetYaxis()->SetTitle("On-axis / Far off-axis");
  ratioQE->SetMinimum(0.0);
  ratioQE->SetMaximum(ymax);
  ratioQE->Draw("E1");
  if (ratioQEIntf) ratioQEIntf->Draw("E1 same");

  TLine* unity = new TLine(ratioQE->GetXaxis()->GetXmin(), 1.0, ratioQE->GetXaxis()->GetXmax(), 1.0);
  unity->SetLineColor(kGray + 1);
  unity->SetLineStyle(3);
  unity->Draw();

  TLegend* leg = new TLegend(0.63, 0.75, 0.92, 0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.032);
  leg->AddEntry(ratioQE, "QE only", "lpe");
  if (ratioQEIntf) leg->AddEntry(ratioQEIntf, "QE + Interference", "lpe");
  leg->Draw();

  TString outBase = TString::Format("%s/prism_onoff_ratio_%s", outdir, histbase);
  c->SaveAs(outBase + ".png");
  c->SaveAs(outBase + ".pdf");

  fOn->Close(); fOff->Close();
  delete fOn; delete fOff;
  delete leg; delete unity; delete c;
}
