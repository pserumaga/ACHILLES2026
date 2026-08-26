void plotQECosTheta_SBND_DUNE(std::string sbndFile, std::string duneFile, std::string outDir, bool normalizeShape = true){

  /*QE-only comparison: outgoing muon cos(theta) and outgoing proton
  cos(theta), SBND flux overlaid on DUNE ND flux (argon target,
  cascades ON). Produces two plots, each with SBND/DUNE on the
  same axes and a shared legend.
  Usage:
  root -l 'plotQECosTheta_SBND_DUNE.cc("/path/achilles_Ar_..._SBNDFLUX_cascadesON_output.root",
                                       "/path/achilles_Ar_..._DUNEFLUX_cascadesON_output.root",
                                       "/path/output/dir/")' */

  if (!outDir.empty() && outDir.back() != '/') outDir += "/";

  TFile* fSBND = new TFile(sbndFile.c_str());
  TFile* fDUNE = new TFile(duneFile.c_str());

  if (!fSBND || fSBND->IsZombie() || !fDUNE || fDUNE->IsZombie()) {
    std::cerr << "Error: could not open one or both input files." << std::endl;
    return;
  }

  // QE-only histograms
  TH1D* hSBND_mu = (TH1D*)fSBND->Get("hOut_mu_cosTheta_qe");
  TH1D* hDUNE_mu = (TH1D*)fDUNE->Get("hOut_mu_cosTheta_qe");
  TH1D* hSBND_p  = (TH1D*)fSBND->Get("hOut_p_cosTheta_qe");
  TH1D* hDUNE_p  = (TH1D*)fDUNE->Get("hOut_p_cosTheta_qe");

  if (!hSBND_mu || !hDUNE_mu || !hSBND_p || !hDUNE_p){
    std::cerr << "Error: one or more QE histograms (hOut_mu_cosTheta_qe / hOut_p_cosTheta_qe) not found." << std::endl;
    return;
  }

  hSBND_mu = (TH1D*)hSBND_mu->Clone("hSBND_mu_qe");
  hDUNE_mu = (TH1D*)hDUNE_mu->Clone("hDUNE_mu_qe");
  hSBND_p  = (TH1D*)hSBND_p->Clone("hSBND_p_qe");
  hDUNE_p  = (TH1D*)hDUNE_p->Clone("hDUNE_p_qe");

  if (normalizeShape){
    if (hSBND_mu->Integral() > 0) hSBND_mu->Scale(1.0 / hSBND_mu->Integral());
    if (hDUNE_mu->Integral() > 0) hDUNE_mu->Scale(1.0 / hDUNE_mu->Integral());
    if (hSBND_p->Integral()  > 0) hSBND_p->Scale(1.0 / hSBND_p->Integral());
    if (hDUNE_p->Integral()  > 0) hDUNE_p->Scale(1.0 / hDUNE_p->Integral());
  }

  Int_t colorSBND = TColor::GetColor("#0072B2");
  Int_t colorDUNE = TColor::GetColor("#D55E00");

  hSBND_mu->SetLineColor(colorSBND); hSBND_mu->SetLineWidth(2); hSBND_mu->SetLineStyle(1);
  hDUNE_mu->SetLineColor(colorDUNE); hDUNE_mu->SetLineWidth(2); hDUNE_mu->SetLineStyle(2);
  hSBND_p->SetLineColor(colorSBND);  hSBND_p->SetLineWidth(2);  hSBND_p->SetLineStyle(1);
  hDUNE_p->SetLineColor(colorDUNE);  hDUNE_p->SetLineWidth(2);  hDUNE_p->SetLineStyle(2);

  for (TH1D* h : {hSBND_mu, hDUNE_mu, hSBND_p, hDUNE_p}){
    h->SetStats(0);
    h->SetTitle("");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
  }

  // Cosine variables: use each histogram's native fixed binning
  double muXmin = std::min(hSBND_mu->GetXaxis()->GetXmin(), hDUNE_mu->GetXaxis()->GetXmin());
  double muXmax = std::max(hSBND_mu->GetXaxis()->GetXmax(), hDUNE_mu->GetXaxis()->GetXmax());
  double pXmin  = std::min(hSBND_p->GetXaxis()->GetXmin(),  hDUNE_p->GetXaxis()->GetXmin());
  double pXmax  = std::max(hSBND_p->GetXaxis()->GetXmax(),  hDUNE_p->GetXaxis()->GetXmax());

  hSBND_mu->GetXaxis()->SetRangeUser(muXmin, muXmax);
  hDUNE_mu->GetXaxis()->SetRangeUser(muXmin, muXmax);
  hSBND_p->GetXaxis()->SetRangeUser(pXmin, pXmax);
  hDUNE_p->GetXaxis()->SetRangeUser(pXmin, pXmax);

  hSBND_mu->GetXaxis()->SetTitle("cos(#theta_{#mu})");
  hSBND_p->GetXaxis()->SetTitle("cos(#theta_{p})");

  string yTitle = normalizeShape ? "Normalized Events (QE)" : "Events (QE)";
  hSBND_mu->GetYaxis()->SetTitle(yTitle.c_str());
  hSBND_p->GetYaxis()->SetTitle(yTitle.c_str());

  hSBND_mu->SetMaximum(std::max(hSBND_mu->GetMaximum(), hDUNE_mu->GetMaximum()) * 1.3);
  hSBND_p->SetMaximum(std::max(hSBND_p->GetMaximum(), hDUNE_p->GetMaximum()) * 1.3);
  hSBND_mu->SetMinimum(0.0);
  hSBND_p->SetMinimum(0.0);

  //Muon cos(theta) plot
  TCanvas* cMu = new TCanvas("cMu", "QE only: muon cos(theta), SBND vs DUNE", 800, 600);
  hSBND_mu->Draw("HIST");
  hDUNE_mu->Draw("HIST same");

  TLegend* lMu = new TLegend(0.15, 0.7, 0.45, 0.88);
  lMu->AddEntry(hSBND_mu, "SBND, QE", "l");
  lMu->AddEntry(hDUNE_mu, "DUNE, QE", "l");
  lMu->SetBorderSize(0);
  lMu->Draw();

  cMu->SaveAs((outDir + "SBND_DUNE_QE_MuonCosTheta.png").c_str());
  cMu->SaveAs((outDir + "SBND_DUNE_QE_MuonCosTheta.pdf").c_str());

  //Proton cos(theta) plot
  TCanvas* cP = new TCanvas("cP", "QE only: proton cos(theta), SBND vs DUNE", 800, 600);
  hSBND_p->Draw("HIST");
  hDUNE_p->Draw("HIST same");

  TLegend* lP = new TLegend(0.15, 0.7, 0.45, 0.88);
  lP->AddEntry(hSBND_p, "SBND, QE", "l");
  lP->AddEntry(hDUNE_p, "DUNE, QE", "l");
  lP->SetBorderSize(0);
  lP->Draw();

  cP->SaveAs((outDir + "SBND_DUNE_QE_ProtonCosTheta.png").c_str());
  cP->SaveAs((outDir + "SBND_DUNE_QE_ProtonCosTheta.pdf").c_str());
}
