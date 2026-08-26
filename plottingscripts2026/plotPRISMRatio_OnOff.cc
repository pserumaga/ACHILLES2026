void plotPRISMRatio_OnOff(std::string onAxisFile, std::string offAxisFile, std::string outDir,
                           std::string detectorLabel = "SBND"){

  /* PRISM on-axis/off-axis ratio, muon energy (E_mu), argon target.
  // Two ratio curves on one plot: QE only, QE + Interference (summed)
  Usage:
  root -l 'plotPRISMRatio_OnOff.cc("/path/output_flux_oaa_numu_0_0_0_2.root",
                                    "/path/output_flux_oaa_numu_1_4_1_6.root",
                                    "/path/output/dir/",
                                    "SBND")'*/

  if (!outDir.empty() && outDir.back() != '/') outDir += "/";

  TFile* fOn  = new TFile(onAxisFile.c_str());
  TFile* fOff = new TFile(offAxisFile.c_str());

  if (!fOn || fOn->IsZombie() || !fOff || fOff->IsZombie()) {
    std::cerr << "Error: could not open one or both input files." << std::endl;
    return;
  }

  // QE only
  TH1D* hOn_qe  = (TH1D*)fOn->Get("hEmu_qe");
  TH1D* hOff_qe = (TH1D*)fOff->Get("hEmu_qe");

  //QE + Interference (summed)
  TH1D* hOn_intf  = (TH1D*)fOn->Get("hEmu_intf");
  TH1D* hOff_intf = (TH1D*)fOff->Get("hEmu_intf");

  if (!hOn_qe || !hOff_qe || !hOn_intf || !hOff_intf){
    std::cerr << "Error: hEmu_qe / hEmu_intf not found in one or both files." << std::endl;
    return;
  }

  TH1D* hOn_qeIntf  = (TH1D*)hOn_qe->Clone("hOn_qeIntf");
  hOn_qeIntf->Add(hOn_intf);
  TH1D* hOff_qeIntf = (TH1D*)hOff_qe->Clone("hOff_qeIntf");
  hOff_qeIntf->Add(hOff_intf);

  // Ratios: off-axis / on-axis
  // (swap Clone base and Divide argument here to flip to on/off)
  TH1D* hRatio_qe = (TH1D*)hOff_qe->Clone("hRatio_qe");
  hRatio_qe->Divide(hOn_qe);

  TH1D* hRatio_qeIntf = (TH1D*)hOff_qeIntf->Clone("hRatio_qeIntf");
  hRatio_qeIntf->Divide(hOn_qeIntf);

  Int_t colorQE     = TColor::GetColor("#E69F00"); // QE only
  Int_t colorQEIntf = TColor::GetColor("#56B4E9"); // QE + Interference

  hRatio_qe->SetLineColor(colorQE);         hRatio_qe->SetLineWidth(2);     hRatio_qe->SetLineStyle(1);
  hRatio_qeIntf->SetLineColor(colorQEIntf); hRatio_qeIntf->SetLineWidth(2); hRatio_qeIntf->SetLineStyle(2);

  hRatio_qe->SetStats(0);
  hRatio_qe->SetTitle("");
  hRatio_qe->GetXaxis()->SetTitle("E_{#mu}^{out} [MeV]");
  hRatio_qe->GetYaxis()->SetTitle("Off-Axis / On-Axis Ratio");
  hRatio_qe->GetXaxis()->SetTitleSize(0.05);
  hRatio_qe->GetYaxis()->SetTitleSize(0.05);

  // Emu has manual cap consistent with the Ar/C to avoid empty space
  hRatio_qe->GetXaxis()->SetRangeUser(0, 3000);

  double yMax = std::max(hRatio_qe->GetMaximum(), hRatio_qeIntf->GetMaximum()) * 1.2;
  hRatio_qe->SetMinimum(0.0);
  hRatio_qe->SetMaximum(yMax);

  TCanvas* c1 = new TCanvas("c1", (detectorLabel + " PRISM on/off-axis ratio, E_mu").c_str(), 800, 600);
  hRatio_qe->Draw("HIST");
  hRatio_qeIntf->Draw("HIST same");

  TLegend* l1 = new TLegend(0.55, 0.72, 0.88, 0.88);
  l1->AddEntry(hRatio_qe, "QE", "l");
  l1->AddEntry(hRatio_qeIntf, "QE + Interference", "l");
  l1->SetBorderSize(0);
  l1->Draw();

  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.045);
  lat.DrawLatex(0.15, 0.83, detectorLabel.c_str());

  string outName = "PRISM_OnOffRatio_" + detectorLabel + "_Emu";
  c1->SaveAs((outDir + outName + ".png").c_str());
  c1->SaveAs((outDir + outName + ".pdf").c_str());
}
