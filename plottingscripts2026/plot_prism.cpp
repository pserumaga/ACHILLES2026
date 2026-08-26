/*plot_prism_bins -- this is for Emu right now but can be changed to any variable by changing the histbase and cat arguments.  
The histograms are expected to be in files named like output_flux_oaa_numu_0_0_0_2.root,
where the four numbers are the low and high angles of the bin in degrees (e.g., 0.0-0.2 degrees). 
The macro will plot all bins found in the specified input directory and save the plots to the specified output directory.*/

#include "TFile.h"
#include "TROOT.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TColor.h"
#include "TAxis.h"
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <regex>
namespace {

struct BinFile {
  std::string path;
  double angleLow;
  double angleHigh;
};

// Pulls the four underscore-separated numbers out of a filename like
// "output_flux_oaa_numu_0_0_0_2.root" -> angleLow=0.0, angleHigh=0.2
bool ParseAngleBin(const std::string& fname, double& lo, double& hi) {
  static const std::regex re(R"((\d+)_(\d+)_(\d+)_(\d+)\.root$)");
  std::smatch m;
  if (!std::regex_search(fname, m, re)) return false;
  lo = std::stod(m[1].str() + "." + m[2].str());
  hi = std::stod(m[3].str() + "." + m[4].str());
  return true;
}

TH1D* GetHistFlex(TFile* f, const std::string& base, const std::string& cat) {
  if (!f) return nullptr;
  const std::string name_us = base + "_" + cat;
  const std::string name_sp = base + " " + cat;
  if (auto h = dynamic_cast<TH1D*>(f->Get(name_us.c_str()))) return h;
  if (auto h = dynamic_cast<TH1D*>(f->Get(name_sp.c_str()))) return h;
  return nullptr;
}

// Maps each PRISM bin to a distinct rainbow color,
// `slot` is a fixed, caller-assigned index (0, 1, 2, ...)
int GradientColor(double t, int slot) { // t in [0,1]
  t = std::max(0.0, std::min(1.0, t));
  Float_t hue = 270.0 * t;
  Float_t r, g, b;
  TColor::HSV2RGB(hue, 1.0, 0.85, r, g, b); // value slightly dimmed for legibility on white
  const Int_t idx = 9000 + slot; // fixed custom-color range, unlikely to collide with anything else
  if (TColor* existing = gROOT->GetColor(idx)) {
    existing->SetRGB(r, g, b); // reuse the slot if this macro already ran once this session
  } else {
    new TColor(idx, r, g, b);
  }
  return idx;
}

}

void plot_prism_emu(const char* indir   = "prism_outputs",
                     const char* outdir = "PRISMplots",
                     const char* histbase = "hEmu",
                     const char* cat      = "total")
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gSystem->mkdir(outdir, kTRUE);

  void* dirp = gSystem->OpenDirectory(indir);
  if (!dirp) { std::cerr << "ERROR: cannot open directory " << indir << "\n"; return; }

  std::vector<BinFile> bins;
  const char* entry;
  while ((entry = gSystem->GetDirEntry(dirp))) {
    std::string fname = entry;
    if (fname.size() < 5 || fname.substr(fname.size() - 5) != ".root") continue;
    double lo, hi;
    if (!ParseAngleBin(fname, lo, hi)) {
      std::cerr << "WARNING: could not parse angle bin from " << fname << "; skipping.\n";
      continue;
    }
    bins.push_back({std::string(indir) + "/" + fname, lo, hi});
  }
  gSystem->FreeDirectory(dirp);

  if (bins.empty()) {
    std::cerr << "ERROR: no parseable OAA output files found in " << indir << "\n";
    return;
  }

  std::sort(bins.begin(), bins.end(),
            [](const BinFile& a, const BinFile& b) { return a.angleLow < b.angleLow; });

  const double globalLo = bins.front().angleLow;
  const double globalHi = bins.back().angleHigh;
  const double span = (globalHi - globalLo > 0) ? (globalHi - globalLo) : 1.0;

  TCanvas* c = new TCanvas("c_prism_emu", "c_prism_emu", 1000, 750);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.13);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.07);

  TLegend* leg = new TLegend(0.63, 0.50, 0.92, 0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.028);

  std::vector<TFile*> openFiles; 
  std::vector<TH1D*> keepHists;  
  TH1D* frame = nullptr;
  double ymax = 0.0;
  int slot = 0; 

  for (const auto& b : bins) {
    TFile* f = TFile::Open(b.path.c_str());
    if (!f || f->IsZombie()) {
      std::cerr << "WARNING: cannot open " << b.path << "; skipping.\n";
      continue;
    }
    openFiles.push_back(f);

    TH1D* h = GetHistFlex(f, histbase, cat);
    if (!h) {
      std::cerr << "WARNING: " << histbase << "_" << cat << " not found in " << b.path << "; skipping.\n";
      continue;
    }
    h->SetDirectory(nullptr); 
    keepHists.push_back(h);

    const double t = (b.angleLow - globalLo) / span;
    const int color = GradientColor(t, slot++);
    h->SetLineColor(color);
    h->SetLineWidth(3);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.6);

    ymax = std::max(ymax, h->GetMaximum());
    if (!frame) frame = h;

    TString label = TString::Format("%.1f#circ - %.1f#circ", b.angleLow, b.angleHigh);
    leg->AddEntry(h, label, "l");
  }

  if (!frame) { std::cerr << "ERROR: no histograms drawn.\n"; delete c; delete leg; return; }

  frame->GetXaxis()->SetTitle("E_{#mu}^{out} [MeV]");
  frame->GetXaxis()->SetRangeUser(0, 2000); // adjust accordingly
  frame->GetYaxis()->SetTitle("Cross Section [a.u.]");
  frame->SetMinimum(0.0);
  frame->SetMaximum(1.25 * ymax);
  frame->SetLineWidth(3);
  frame->Draw("hist");

  for (auto* h : keepHists) if (h != frame) h->Draw("hist same");

  leg->Draw();

  TString outBase = TString::Format("%s/prism_%s_%s", outdir, histbase, cat);
  c->SaveAs(outBase + ".png");
  c->SaveAs(outBase + ".pdf");

  for (auto* f : openFiles) { f->Close(); delete f; }
  delete leg;
  delete c;
}
