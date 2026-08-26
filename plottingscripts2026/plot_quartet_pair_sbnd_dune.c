#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TAxis.h"
#include "TMath.h"
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

namespace {

TH1D* GetHistFlex(TFile* f, const std::string& base, const std::string& cat) {
  if (!f) return nullptr;
  const std::string name_us = base + "_" + cat; // e.g. hEmu_res
  const std::string name_sp = base + " " + cat; // e.g. hEmu res
  if (auto h = dynamic_cast<TH1D*>(f->Get(name_us.c_str()))) return h;
  if (auto h = dynamic_cast<TH1D*>(f->Get(name_sp.c_str()))) return h;
  return nullptr;
}

void StyleCurve(TH1* h, Color_t c, int lstyle = 1) {
  if (!h) return;
  h->SetLineColor(c);
  h->SetLineStyle(lstyle);
  h->SetMarkerColor(c);
  h->SetLineWidth(3);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.8);
}

double MaxY(const std::vector<TH1D*>& hs) {
  double m = 0.0;
  for (auto* h : hs) if (h) m = std::max(m, h->GetMaximum());
  return m;
}

// Last bin (with padding) where any histogram in hs still has real content,
// as an x-axis upper edge. frac_of_max is relative to that group's own max.
double AutoXMax(const std::vector<TH1D*>& hs, double frac_of_max = 0.003) {
  double gmax = MaxY(hs);
  if (gmax <= 0) return 0;
  double threshold = frac_of_max * gmax;

  TH1D* ref = nullptr;
  for (auto* h : hs) if (h) { ref = h; break; }
  if (!ref) return 0;
  int nbins = ref->GetNbinsX();

  int lastBin = 1;
  for (auto* h : hs) {
    if (!h) continue;
    for (int b = h->GetNbinsX(); b >= 1; --b) {
      if (h->GetBinContent(b) > threshold) {
        lastBin = std::max(lastBin, b);
        break;
      }
    }
  }

  int pad = std::max(2, nbins / 40); // ~2.5% padding past the last real content
  int paddedBin = std::min(nbins, lastBin + pad);
  return ref->GetXaxis()->GetBinUpEdge(paddedBin);
}

bool RightSideIsHigh(const std::vector<TH1D*>& hs, double xmaxVisible, double ymax) {
  TH1D* ref = nullptr;
  for (auto* h : hs) if (h) { ref = h; break; }
  if (!ref || ymax <= 0) return false;

  int lastVisibleBin = ref->GetXaxis()->FindBin(xmaxVisible);
  lastVisibleBin = std::min(lastVisibleBin, ref->GetNbinsX());
  int startCheck = std::max(1, (int)(lastVisibleBin * 0.72));

  for (auto* h : hs) {
    if (!h) continue;
    for (int b = startCheck; b <= lastVisibleBin; ++b) {
      if (h->GetBinContent(b) > 0.45 * ymax) return true;
    }
  }
  return false;
}

void SetGlobalStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelSize(0.038, "XYZ");
  gStyle->SetTitleSize(0.048, "XYZ");
  gStyle->SetTitleOffset(1.15, "X");
  gStyle->SetTitleOffset(1.45, "Y");
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(42);
}

enum class AxisMode { AutoShared, ManualShared, Angle, Cosine };

struct Fam {
  std::string base;
  std::string xtitle;
  std::string ytitle;
  AxisMode mode;
  double manualCap = 0.0; // used when mode == ManualShared
};

// One quartet (Total/QE/Interference/Resonance) drawn on the current pad.
// xmin/xmax define the visible x-range (already resolved by the caller).
void DrawQuartet(TFile* f, const Fam& fam, const std::vector<std::string>& cats,
                  const Color_t cols[4], const int lstyles[4], const char* legLabels[4],
                  double xmin, double xmax, double ymax,
                  const char* fileTag, const char* displayLabel, const char* outdir) {
  std::vector<TH1D*> h(4, nullptr);
  for (int i = 0; i < 4; i++) {
    h[i] = GetHistFlex(f, fam.base, cats[i]);
    if (!h[i]) std::cerr << "WARNING: " << fam.base << "_" << cats[i]
                          << " not found in " << fileTag << ".\n";
  }
  if (!h[0] && !h[1] && !h[2] && !h[3]) {
    std::cerr << "WARNING: none of " << fam.base << " {total,qe,intf,res} found for "
               << fileTag << "; skipping.\n";
    return;
  }

  TH1D* frame = nullptr;
  for (auto* hist : h) { if (hist) { frame = hist; break; } }
  for (int i = 0; i < 4; i++) StyleCurve(h[i], cols[i], lstyles[i]);

  TString cname = TString::Format("c_%s_%s", fam.base.c_str(), fileTag);
  TCanvas* c = new TCanvas(cname, cname, 1000, 750);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.13);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.07);

  frame->GetXaxis()->SetTitle(fam.xtitle.c_str());
  frame->GetYaxis()->SetTitle(fam.ytitle.c_str());
  frame->SetMinimum(0.0);
  frame->SetMaximum(ymax);
  frame->SetLineWidth(2);
  frame->GetXaxis()->SetRangeUser(xmin, xmax);

  frame->Draw("hist");
  for (int i = 0; i < 4; i++) {
    if (h[i] && h[i] != frame) h[i]->Draw("hist same");
  }

  bool useLeft = RightSideIsHigh(h, xmax, ymax);
  double lx1 = useLeft ? 0.18 : 0.63;
  double lx2 = useLeft ? 0.45 : 0.90;
  TLegend* leg = new TLegend(lx1, 0.68, lx2, 0.90);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  for (int i = 0; i < 4; i++) if (h[i]) leg->AddEntry(h[i], legLabels[i], "l");
  leg->Draw();

  // On-plot corner label (e.g. detector name) -- only drawn if requested.
  if (displayLabel && std::string(displayLabel).length() > 0) {
    TLatex lt;
    lt.SetNDC();
    lt.SetTextFont(42);
    lt.SetTextSize(0.035);
    lt.SetTextAlign(11);
    lt.DrawLatex(0.15, 0.925, displayLabel);
  }

  TString outBase = TString::Format("%s/%s_%s", outdir, fam.base.c_str(), fileTag);
  c->SaveAs(outBase + ".png");
  c->SaveAs(outBase + ".pdf");

  delete leg;
  delete c;
}

} // anon

// infileSBND / infileDUNE: the two ROOT files to compare.
// outdir: single directory both SBND_*.png and DUNE_*.png get written to.
void plot_quartet_pair(const char* infileSBND = "sbnd_output.root",
                        const char* infileDUNE = "dune_output.root",
                        const char* outdir = "paperplots")
{
  SetGlobalStyle();
  gSystem->mkdir(outdir, kTRUE);

  TFile* fS = TFile::Open(infileSBND);
  TFile* fD = TFile::Open(infileDUNE);
  if (!fS || fS->IsZombie()) { std::cerr << "ERROR: cannot open " << infileSBND << "\n"; return; }
  if (!fD || fD->IsZombie()) { std::cerr << "ERROR: cannot open " << infileDUNE << "\n"; return; }

  std::vector<Fam> fams = {
    {"hIn_nu_E",         "E_{#nu}^{in} [MeV]",             "Cross Section [a.u.]", AxisMode::ManualShared, 5000.0},
    {"hOut_mu_cosTheta", "cos#theta(#mu)",                 "Cross Section [a.u.]", AxisMode::Cosine},
    {"hOut_p_cosTheta",  "cos#theta(proton)",              "Cross Section [a.u.]", AxisMode::Cosine},
    {"homega",           "#omega = E_{#nu}-E_{#mu} [MeV]", "Cross Section [a.u.]", AxisMode::AutoShared},
    {"hQ2",              "Q^{2} [MeV^{2}]",                "Cross Section [a.u.]", AxisMode::AutoShared},
    {"htheta_mu_p",      "#theta(#mu,p) [rad]",            "Cross Section [a.u.]", AxisMode::Angle},
    {"hdeltapt",         "#delta p_{T} [MeV]",             "Cross Section [a.u.]", AxisMode::AutoShared},
    {"halphat",          "#delta #alpha_{T} [rad]",        "Cross Section [a.u.]", AxisMode::Angle},
    {"hphit",            "#delta #phi_{T} [rad]",          "Cross Section [a.u.]", AxisMode::Angle},
    {"hEmu",             "E_{#mu}^{out} [MeV]",            "Cross Section [a.u.]", AxisMode::ManualShared, 5000.0},
  };

  const std::vector<std::string> cats = {"total", "qe", "intf", "res"};
  const Color_t colTotal = kBlack;
  const Color_t colQE    = kBlue + 2;
  const Color_t colIntf  = kOrange + 7;
  const Color_t colRes   = kGreen + 3;
  const Color_t cols[4] = { colTotal, colQE, colIntf, colRes };
  const int lstyles[4]  = { 1, 1, 2, 3 };
  const char* legLabels[4] = { "Total", "QE", "Interference", "Resonance" };

  for (const auto& fam : fams) {
    std::vector<TH1D*> hS(4, nullptr), hD(4, nullptr);
    for (int i = 0; i < 4; i++) {
      hS[i] = GetHistFlex(fS, fam.base, cats[i]);
      hD[i] = GetHistFlex(fD, fam.base, cats[i]);
    }
    bool anyS = hS[0] || hS[1] || hS[2] || hS[3];
    bool anyD = hD[0] || hD[1] || hD[2] || hD[3];
    if (!anyS && !anyD) {
      std::cerr << "WARNING: " << fam.base << " not found in either file; skipping.\n";
      continue;
    }

    // Resolve the shared x-range for this family.
    double xmin = 0.0, xmax = 0.0;
    if (fam.mode == AxisMode::Angle) {
      xmin = 0.0;
      xmax = TMath::Pi();
    } else if (fam.mode == AxisMode::Cosine) {
      xmin = -1.0;
      xmax = 1.0;
    } else if (fam.mode == AxisMode::ManualShared) {
      xmin = 0.0;
      xmax = fam.manualCap;
    } else { // AutoShared
      double xS = anyS ? AutoXMax(hS) : 0.0;
      double xD = anyD ? AutoXMax(hD) : 0.0;
      xmax = std::max(xS, xD);
      if (xmax <= 0) xmax = 1.0;
      xmin = 0.0;
    }

    // Each detector still gets its own y-axis scale (SBND and DUNE cross
    // sections differ by orders of magnitude -- that's expected, not a bug).
    double ymaxS = anyS ? 1.25 * MaxY(hS) : 0.0;
    double ymaxD = anyD ? 1.25 * MaxY(hD) : 0.0;
    if (ymaxS <= 0) ymaxS = 1.0;
    if (ymaxD <= 0) ymaxD = 1.0;

    if (anyS) DrawQuartet(fS, fam, cats, cols, lstyles, legLabels, xmin, xmax, ymaxS, "SBND", "", outdir);
    if (anyD) DrawQuartet(fD, fam, cats, cols, lstyles, legLabels, xmin, xmax, ymaxD, "DUNE", "", outdir);
  }

  fS->Close();
  fD->Close();
  delete fS;
  delete fD;
}
