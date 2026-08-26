/* plot_quartet.C
Draws total, QE, intf, res for every histogram family
booked in analyze_hepmcCC.cpp.
Changes to Nathan's version:
- X-axis is now auto-trimmed per histogram so trailing empty bins get cut off automatically.
- Legend position is now chosen based on where the curves are
- If the content is high on the right side of the visible range, the legend moves to the top-left instead of overlapping the lines which was an issue.*/

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TAxis.h"
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

// Find upper x-edge that trims near-empty bins.
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

  int pad = std::max(2, nbins / 40); // ~2.5% padding past the last real content to avoid cutting last bin edge.
  int paddedBin = std::min(nbins, lastBin + pad);
  return ref->GetXaxis()->GetBinUpEdge(paddedBin);
}

// Decide whether curves are tall on the right-hand side of the visible and determines legend placement.
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

}

void plot_quartet(const char* infile = "my_output.root", const char* outdir = "plots",
                   const char* label = "")
{
  SetGlobalStyle();
  gSystem->mkdir(outdir, kTRUE);

  TFile* f = TFile::Open(infile);
  if (!f || f->IsZombie()) { std::cerr << "ERROR: cannot open " << infile << "\n"; return; }

  /* fam.xmax left at 0 for every family now -- x-range is auto-computed per histogram instead. 
  Set a nonzero value here only if you want to force a manual cut for a specific variable. Some of them needed manual cuts because the auto 
  trim didn't function properly and there was still a lot of empty space.*/
  struct Fam { std::string base; std::string xtitle; std::string ytitle; double xmax; };
  std::vector<Fam> fams = {
    {"hIn_nu_E",         "E_{#nu}^{in} [MeV]",             "Cross Section [a.u.]", 0},
    {"hOut_mu_pz",       "p_{z}^{#mu} [MeV]",              "Cross Section [a.u.]", 0},
    {"hOut_mu_cosTheta", "cos#theta(#mu)",                 "Cross Section [a.u.]", 0},
    {"hOut_p_cosTheta",  "cos#theta(hadron)",              "Cross Section [a.u.]", 0},
    {"homega",           "#omega = E_{#nu}-E_{#mu} [MeV]", "Cross Section [a.u.]", 0},
    {"hmu_pmag",         "|#vec{p}_{#mu}| [MeV]",          "Cross Section [a.u.]", 0},
    {"hp_pmag",          "|#vec{p}_{h}| [MeV]",            "Cross Section [a.u.]", 0},
    {"hQ2",              "Q^{2} [MeV^{2}]",                "Cross Section [a.u.]", 0},
    {"htheta_mu_p",      "#theta(#mu,h) [rad]",            "Cross Section [a.u.]", 0},
    {"hdeltapt",         "#delta p_{T} [MeV]",             "Cross Section [a.u.]", 0},
    {"halphat",          "#delta #alpha_{T} [rad]",        "Cross Section [a.u.]", 0},
    {"hphit",            "#delta #phi_{T} [rad]",          "Cross Section [a.u.]", 0},
    {"hEmu",             "E_{#mu}^{out} [MeV]",            "Cross Section [a.u.]", 0}
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
    std::vector<TH1D*> h(4, nullptr);
    for (int i = 0; i < 4; i++) {
      h[i] = GetHistFlex(f, fam.base, cats[i]);
      if (!h[i]) {
        std::cerr << "WARNING: " << fam.base << "_" << cats[i] << " not found.\n";
      }
    }
    if (!h[0] && !h[1] && !h[2] && !h[3]) {
      std::cerr << "WARNING: none of " << fam.base << " {total,qe,intf,res} found; skipping.\n";
      continue;
    }

    TH1D* frame = nullptr;
    for (auto* hist : h) { if (hist) { frame = hist; break; } }

    for (int i = 0; i < 4; i++) StyleCurve(h[i], cols[i], lstyles[i]);

    // Manual x-axis overrides for specific detector combos, because the auto-trim doesn't always work well for some variables.
    std::string labelStr = label ? label : "";
    bool isDUNE = labelStr.find("DUNE") != std::string::npos;
    bool isSBND = labelStr.find("SBND") != std::string::npos;
    double manualXMax = 0.0;
    if (fam.base == "hIn_nu_E") {
      if (isSBND) manualXMax = 3000.0;
      else if (isDUNE) manualXMax = 6000.0;
    } else if (fam.base == "hEmu") {
      if (isSBND) manualXMax = 3000.0;
      else if (isDUNE) manualXMax = 4000.0;
    }

    TString cname = TString::Format("c_%s", fam.base.c_str());
    TCanvas* c = new TCanvas(cname, cname, 1000, 750);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.13);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.07);

    double ymax = 1.25 * MaxY(h);
    if (ymax <= 0) ymax = 1.0;

    // Auto-trim the x-axis unless a manual cap (above) or fam.xmax applies.
    double xmaxUse = manualXMax > 0 ? manualXMax : fam.xmax;
    if (xmaxUse <= 0) xmaxUse = AutoXMax(h);

    frame->GetXaxis()->SetTitle(fam.xtitle.c_str());
    frame->GetYaxis()->SetTitle(fam.ytitle.c_str());
    frame->SetMinimum(0.0);
    frame->SetMaximum(ymax);
    frame->SetLineWidth(2);
    if (xmaxUse > 0) frame->GetXaxis()->SetRangeUser(0, xmaxUse);

    frame->Draw("hist");
    for (int i = 0; i < 4; i++) {
      if (h[i] && h[i] != frame) h[i]->Draw("hist same");
    }

    // Pick a legend corner that avoids the curves. Default is top-right;
    // switch to top-left if the content is tall near the right edge.
    bool useLeft = RightSideIsHigh(h, xmaxUse > 0 ? xmaxUse : frame->GetXaxis()->GetXmax(), ymax);
    double lx1 = useLeft ? 0.18 : 0.63;
    double lx2 = useLeft ? 0.45 : 0.90;
    TLegend* leg = new TLegend(lx1, 0.68, lx2, 0.90);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    for (int i = 0; i < 4; i++) {
      if (h[i]) leg->AddEntry(h[i], legLabels[i], "l");
    }
    leg->Draw();

    if (label && std::string(label).length() > 0) {
      TLatex lt;
      lt.SetNDC();
      lt.SetTextFont(42);
      lt.SetTextSize(0.035);
      lt.SetTextAlign(11);
      lt.DrawLatex(0.15, 0.925, label);
    }

    TString outBase = TString::Format("%s/%s", outdir, fam.base.c_str());
    c->SaveAs(outBase + ".png");
    c->SaveAs(outBase + ".pdf");

    delete leg;
    delete c;
  }

  f->Close();
  delete f;
}
