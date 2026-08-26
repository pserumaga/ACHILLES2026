void plotArCRatioEmu(std::string argonFile, std::string carbonFile, std::string outDir){

  //---------------------------------------------------------------
  // Argon/Carbon cross-section ratio: outgoing muon energy (E_mu)
  // SBND flux, cascades ON
  //
  // Usage:
  //   root -l 'plotArCRatioEmu.cc("/path/to/achilles_Ar_..._output.root",
  //                                "/path/to/achilles_C_..._output.root",
  //                                "/path/to/output/dir/")'
  //---------------------------------------------------------------

  // Make sure outDir ends with a trailing slash
  if (!outDir.empty() && outDir.back() != '/') outDir += "/";

  // Read in files
  TFile* fArgon  = new TFile(argonFile.c_str());
  TFile* fCarbon = new TFile(carbonFile.c_str());

  if (!fArgon || fArgon->IsZombie() || !fCarbon || fCarbon->IsZombie()) {
    std::cerr << "Error: could not open one or both input files." << std::endl;
    return;
  }

  // Variable to plot: outgoing muon energy.
  // NOTE: the real histogram base name in these files is "hEmu"
  // (last year's script assumed "hDelta_phi_T", which doesn't exist
  // in this output -- the actual TKI names are hdeltapt/halphat/hphit).
  string variable = "hEmu";

  // Build the names of the 3 component histograms
  string totalName = variable + "_total";
  string qeName    = variable + "_qe";
  string intfName  = variable + "_intf";

  TH1D* hArgon_total  = (TH1D*)fArgon->Get(totalName.c_str());
  TH1D* hCarbon_total = (TH1D*)fCarbon->Get(totalName.c_str());
  TH1D* hRatio_total  = (TH1D*)hArgon_total->Clone("hRatio_total");
  hRatio_total->Divide(hCarbon_total);

  TH1D* hArgon_qe  = (TH1D*)fArgon->Get(qeName.c_str());
  TH1D* hCarbon_qe = (TH1D*)fCarbon->Get(qeName.c_str());
  TH1D* hRatio_qe  = (TH1D*)hArgon_qe->Clone("hRatio_qe");
  hRatio_qe->Divide(hCarbon_qe);

  TH1D* hArgon_intf  = (TH1D*)fArgon->Get(intfName.c_str());
  TH1D* hCarbon_intf = (TH1D*)fCarbon->Get(intfName.c_str());
  TH1D* hRatio_intf  = (TH1D*)hArgon_intf->Clone("hRatio_intf");
  hRatio_intf->Divide(hCarbon_intf);

  //-------------------------------------------------
  // Colorblind-friendly palette (Okabe-Ito), pinned to fixed
  // hex-derived indices rather than TColor::GetFreeColorIndex(),
  // so colors stay stable across repeated macro runs.
  //-------------------------------------------------
  Int_t colorBlack   = TColor::GetColor("#000000"); // Total
  Int_t colorOrange  = TColor::GetColor("#E69F00"); // QE
  Int_t colorSkyBlue = TColor::GetColor("#56B4E9"); // Interference

  hRatio_total->SetLineColor(colorBlack);
  hRatio_qe->SetLineColor(colorOrange);
  hRatio_intf->SetLineColor(colorSkyBlue);

  hRatio_total->SetLineWidth(2);
  hRatio_qe->SetLineWidth(2);
  hRatio_intf->SetLineWidth(2);

  // Different line styles too, so the plot still reads in
  // grayscale printouts on top of the color coding.
  hRatio_total->SetLineStyle(1); // solid
  hRatio_qe->SetLineStyle(2);    // dashed
  hRatio_intf->SetLineStyle(3);  // dotted

  // Cosmetics
  hRatio_total->SetStats(0);
  hRatio_total->SetTitle("");
  hRatio_total->GetYaxis()->SetTitle("Ratio of Cross-Sections (#sigma(Ar)/#sigma(C))");
  hRatio_total->GetXaxis()->SetTitle("E_{#mu}^{out} [MeV]");
  // E_mu is an energy variable -> manual cap rather than auto-trim.
  // 99% of total-mode content sits below ~2700 MeV for this SBND sample,
  // so 3000 MeV keeps the physically populated region without dead space.
  hRatio_total->GetXaxis()->SetRangeUser(0, 3000);
  hRatio_total->SetMinimum(0.0);
  hRatio_total->SetMaximum(5.0);
  hRatio_total->GetXaxis()->SetTitleSize(0.05);
  hRatio_total->GetYaxis()->SetTitleSize(0.05);

  // Canvas
  TCanvas* c1 = new TCanvas("c1", "Ar/C ratio - outgoing muon energy", 800, 600);
  hRatio_total->Draw("HIST");
  hRatio_qe->Draw("HIST same");
  hRatio_intf->Draw("HIST same");

  TLegend* l1 = new TLegend(0.6, 0.7, 0.9, 0.9);
  l1->AddEntry(hRatio_total, "Total", "l");
  l1->AddEntry(hRatio_qe, "QE", "l");
  l1->AddEntry(hRatio_intf, "Interference", "l");
  l1->SetBorderSize(0);
  l1->Draw();

  c1->SaveAs((outDir + "ArC_ratio_Emu.png").c_str());
  c1->SaveAs((outDir + "ArC_ratio_Emu.pdf").c_str());
}
