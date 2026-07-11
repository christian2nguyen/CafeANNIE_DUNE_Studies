// Make oscillated predictions
// cafe demo1.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Cuts/AnaCuts.h"
#include "CAFAna/Analysis/Plots.h"
#include "StandardRecord/StandardRecord.h"
#include "TCanvas.h"
#include "TH1.h"

#include "TPad.h"
#include <string>
#include "TAxis.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"

// New includes for this macro
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Analysis/Calcs.h"
#include "OscLib/func/OscCalculatorPMNSOpt.h"
#include "OscLib/func/IOscCalculator.h"
//#include "OscLib/func/IOscCalculatorAdjustable.h"
//#include "OscLib/OscCalcPMNS.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TLatex.h"
#include <vector>
#include "TROOT.h"

/**
 * Draws a DUNE nu_e Appearance parameter label box on a histogram.
 *
 * @param hist         The TH1D to draw the label on (must already be drawn)
 * @param ordering     "Normal" or "Inverted" 
 * @param sin2_2th13   Value of sin^2(2theta_13)
 * @param sin2_th23    Value of sin^2(theta_23)
 * @param years        Exposure in years (e.g. 3.5)
 * @param staged       Whether the run is staged (default true)
 * @param x1, y1       NDC lower-left corner of the box (default top-left)
 * @param x2, y2       NDC upper-right corner of the box
 */
 TPaveText* DrawDUNELabel(
    bool isDisapp,
    const char* ordering    = "Normal",
    double      sin2_2th13  = 0.088,
    double      sin2_th23   = 0.580,
    double       dmsq32     = .002451,
    double      years       = 3.5,
    bool        staged      = true,
    double      x1          = 0.5,
    double      y1          = 0.62,
    double      x2          = 0.8,
    double      y2          = 0.88
) ;


 TPaveText* DrawDUNELabel(
    bool isDisapp,
    const char* ordering    ,
    double      s22th13     ,      // renamed
    double      s2th23      ,      // renamed
    double       dmsq32 ,
    double      years       ,
    bool        staged      ,
    double      x1          ,
    double      y1          ,
    double      x2          ,
    double      y2          
) {
    TPaveText* pt = new TPaveText(x1, y1, x2, y2, "NDC");
    pt->SetFillColor(0);
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetTextFont(62);
    pt->SetTextSize(0.035);
   if(isDisapp) pt->AddText("DUNE #nu_{#mu} Disapperance");
   else pt->AddText("DUNE #nu_{e} Appearance");
    
    pt->AddText(Form("%s Ordering", ordering));
    pt->AddText(Form("sin^{2}2#theta_{13} = %.3f", s22th13));   // fixed
    pt->AddText(Form("sin^{2}#theta_{23} = %.3f",  s2th23));    // fixed
    pt->AddText(Form("#Delta m^{2}_{32} = %.3f #times 10^{-3} eV^{2}",  dmsq32*1000));    // fixed
    pt->AddText(staged ? Form("%.1f years (staged)", years)
                       : Form("%.1f years", years));

    pt->Draw("same");
    return pt;
}

using namespace ana;




void demo1()
{ 
  gROOT->SetBatch(kTRUE);
  TLegend* lg1 = new TLegend( 0.5, 0.3, 0.8, 0.6 );
    lg1->SetNColumns(1);
    lg1->SetBorderSize(0);
    lg1->SetTextSize(.03); //  
    
    

    
      TLegend* lg2 = new TLegend( 0.35, 0.7, 0.8, 0.88 );
    lg2->SetNColumns(1);
    lg2->SetBorderSize(0);
    lg2->SetTextSize(.03); //  
    
    
          TLegend* lg3 = new TLegend( 0.55, 0.28, 0.8, 0.6);
    lg3->SetNColumns(1);
    lg3->SetBorderSize(0);
    lg3->SetTextSize(.028); //  
    
    
    
    
      TLegend* lg4 = new TLegend( 0.2, 0.7, 0.8, 0.88 );
    lg4->SetNColumns(1);
    lg4->SetBorderSize(0);
    lg4->SetTextSize(.03); //
    
 char pdf_title[1024];
 std::string Pdf_name = "DUNEexample_demo1_RHC_v3";
  TCanvas *c1 = new TCanvas("c1");
 sprintf(pdf_title, "%s.pdf(", Pdf_name.c_str());
 c1 -> Print(pdf_title);
 sprintf(pdf_title, "%s.pdf", Pdf_name.c_str());

 auto HistRate = [](TH1* h){
   return h ? h->Integral(1, h->GetNbinsX()) : 0.0;
 };

 auto HistRateRange = [](TH1* h, double xmin, double xmax){
   if(!h) return 0.0;
   const int bmin = h->GetXaxis()->FindBin(xmin);
   const int bmax = h->GetXaxis()->FindBin(xmax);
   return h->Integral(bmin, bmax);
 };

 auto MakeRateLabel = [&](const char* base, TH1* h, double total){
   const double rate = HistRate(h);
   const double frac = (total > 0) ? (100.0 * rate / total) : 0.0;
   return std::string(Form("%s: %.1f (%.1f%%)", base, rate, frac));
 };

 auto MakeRateLabelRange = [&](const char* base, TH1* h, double total, double xmin, double xmax){
   const double rate = HistRateRange(h, xmin, xmax);
   const double frac = (total > 0) ? (100.0 * rate / total) : 0.0;
   return std::string(Form("%s: %.1f (%.1f%%)", base, rate, frac));
 };

 auto PrintFigureSummary = [&](const char* figName,
                               const std::vector<std::pair<std::string, TH1*>>& comps,
                               double total){
   std::cout << "\n=== " << figName << " ===\n";
   std::cout << "Total events = " << total << "\n";
   for(const auto& c: comps){
     const double rate = HistRate(c.second);
     const double frac = (total > 0) ? (100.0 * rate / total) : 0.0;
     std::cout << "  " << c.first << ": " << rate << " (" << frac << "%)\n";
   }
 };

  // See demo0.C for explanation of these repeated parts
  //const std::string fnameNonSwap = "/pnfs/dune/persistent/users/cnguyen/old_statefiles_cp_uboonearea/FD_FHC_nonswap.root";
  //const std::string fnameNueSwap = "/pnfs/dune/persistent/users/cnguyen/old_statefiles_cp_uboonearea/FD_FHC_nueswap.root";
  //const std::string fnameTauSwap = "/pnfs/dune/persistent/users/cnguyen/old_statefiles_cp_uboonearea/FD_FHC_tauswap.root";
  
    const std::string fnameNonSwap = "/pnfs/dune/persistent/users/LBL_TDR/CAFs/v3/FD_FHC_nonswap.root";
    const std::string fnameNueSwap = "/pnfs/dune/persistent/users/LBL_TDR/CAFs/v3/FD_FHC_nueswap.root";
    const std::string fnameTauSwap = "/pnfs/dune/persistent/users/LBL_TDR/CAFs/v3/FD_FHC_tauswap.root";
  
  
  
  SpectrumLoader loaderNonSwap(fnameNonSwap);
  SpectrumLoader loaderNueSwap(fnameNueSwap);
  SpectrumLoader loaderTauSwap(fnameTauSwap);
  const Var kRecoEnergy = SIMPLEVAR(dune.Ev_reco_numu);
  const Var kCVNNumu = SIMPLEVAR(dune.cvnnumu);
  const Var kCVNNue = SIMPLEVAR(dune.cvnnue);
  const Binning binsEnergy = Binning::Simple(40, 0, 10);
  const Binning binsPro = Binning::Simple(60, 0, 1.0);
  const HistAxis axEnergy("Reco #nu_{Energy} [GeV]", binsEnergy, kRecoEnergy);
  const HistAxis axCVN("CVN  #nu_{#mu} Probability", binsPro, kCVNNumu);
  const HistAxis axCVN_nue("CVN  #nu_{e} Probability", binsPro, kCVNNue);
 // const double pot = 3.5 * 1.47e21 * 40/1.13;
 //const double pot = 3.5 * 1.1e21 * 40/1.13;
 const double pot = (3.5) * 1.1e21 * 40/1.13 ;
  // A cut is structured like a Var, but returning bool
  const Cut kPassesCVN({},
                       [](const caf::StandardRecord* sr)
                       {
                         return sr->dune.cvnnumu > 0.5;
                       });

 
  // In many cases it's easier to form them from existing Vars like this
  //  const Cut kPassesCVN = kCVNNumu > 0;

  // A Prediction is an objects holding a variety of "OscillatableSpectrum"
  // objects, one for each original and final flavour combination.
  PredictionNoExtrap pred(loaderNonSwap, loaderNueSwap, loaderTauSwap,
                          axEnergy, kPassFD_CVN_NUMU&&kIsTrueFV);

  PredictionNoExtrap pred_nue(loaderNonSwap, loaderNueSwap, loaderTauSwap,
                          axEnergy, kPassFD_CVN_NUE&&kIsTrueFV);


  PredictionNoExtrap pred_nusignal(loaderNonSwap, loaderNueSwap, loaderTauSwap,
                          axCVN, kIsTrueFV);

    PredictionNoExtrap pred_CVN_nue(loaderNonSwap, loaderNueSwap, loaderTauSwap,
                          axCVN_nue, kIsTrueFV);

    PredictionNoExtrap pred_tao(loaderNonSwap, loaderNueSwap, loaderTauSwap,
                          axCVN, kIsTauFromMu&&kIsTrueFV||kIsTauFromE&&kIsTrueFV);
  
  // These calls will fill all of the constituent parts of the prediction
  loaderNonSwap.Go();
  loaderNueSwap.Go();
  loaderTauSwap.Go();

  // We can extract a total prediction unoscillated
  const Spectrum sUnosc = pred.PredictUnoscillated();
  // Or oscillated, in this case using reasonable parameters from
  // Analysis/Calcs.h
  //osc::IOscCalculator* calc = DefaultOscCalc();
  osc::IOscCalculatorAdjustable* calc = DefaultOscCalc();
  
calc->SetDmsq32(2.451e-3);
calc->SetDmsq21(7.39e-5);
double s2th23 = 0.580; // target sin^2(theta23)
calc->SetTh23(asin(sqrt(s2th23)));
double ss2th12 = 0.088;
calc->SetTh12(0.5 * asin(sqrt(ss2th12)));
calc->SetTh13(.151);
calc->SetdCP(0.0); // radians


  osc::IOscCalculatorAdjustable* calc_pi_2 = DefaultOscCalc();
  double half_pi = TMath::Pi() / 2.0;
calc_pi_2->SetDmsq32(2.451e-3);
calc_pi_2->SetDmsq21(7.39e-5);
calc_pi_2->SetTh23(asin(sqrt(s2th23)));
calc_pi_2->SetTh12(0.5 * asin(sqrt(ss2th12)));
calc_pi_2->SetTh13(.151);
calc_pi_2->SetdCP(half_pi); // radians

  osc::IOscCalculatorAdjustable* calc_pi_3 = DefaultOscCalc();
calc_pi_3->SetDmsq32(2.451e-3);
calc_pi_3->SetDmsq21(7.39e-5);
calc_pi_3->SetTh23(asin(sqrt(s2th23)));
calc_pi_3->SetTh12(0.5 * asin(sqrt(ss2th12)));
calc_pi_3->SetTh13(.151);
calc_pi_3->SetdCP(-half_pi); // radians

  //osc::OscCalcPMNS* pmns = dynamic_cast<osc::OscCalcPMNS*>(calc);
double dmsq32 = calc->GetDmsq32();
double th23 = calc->GetTh23();
double th13 = calc->GetTh13();
 //double th23 = pmns->GetTh23();
//
 double sin2_2theta13 = pow(sin(2.0 * th13), 2);
  double sin2_theta23  = pow(sin(th23), 2);

 // const Spectrum sOsc = pred.Predict(calc); // old 
   const Spectrum sOsc = pred.PredictComponent(calc,
                                                  Flavors::kAllNuMu,
                                                  Current::kCC,
                                                  Sign::kNu);
  //const Spectrum sOsc_nue = pred_nue.Predict(calc);



  const Spectrum sOsc_nusignal = pred_nusignal.Predict(calc);
  
  //const Spectrum sOsc_NC = pred_NC.Predict(calc);
  const Spectrum sOsc_ta = pred_tao.Predict(calc);


  // And we can break things down by flavour
  const Spectrum sOsc_NC = pred.PredictComponent(calc,
                                                  Flavors::kAll,
                                                  Current::kNC,
                                                  Sign::kBoth);
                                                  
   const Spectrum sOsc_ANTInu = pred.PredictComponent(calc,
                                                  Flavors::kAllNuMu,
                                                  Current::kCC,
                                                  Sign::kAntiNu);
                                                  
    const Spectrum sOsc_ALLNue = pred.PredictComponent(calc,
                                                  Flavors::kAllNuE,
                                                  Current::kCC,
                                                  Sign::kBoth);
    const Spectrum sOsc_ALLtao = pred.PredictComponent(calc,
                                                  Flavors::kAllNuTau,
                                                  Current::kCC,
                                                  Sign::kBoth);
  const Spectrum sOsc_Data = pred.PredictComponent(calc,
                                                  Flavors::kAll,
                                                  Current::kBoth,
                                                  Sign::kBoth);                                                  

                                                  
                                 //Flavors::kNuEToNuEx
                                 //Flavors::kNuMuToNuEx                 
                                 //Flavors::kNuEToNuE
                                 //Sign::kNu
/*
      kNuEToNuE    = 1<<0, ///< \f$\nu_e\to\nu_e\f$ ('beam \f$\nu_e \f$')
      kNuEToNuMu   = 1<<1, ///< \f$\nu_e\to\nu_\mu\f$ ('\f$\nu_\mu\f$ appearance')
      kNuEToNuTau  = 1<<2, ///< \f$\nu_e\to\nu_\tau\f$
      kNuMuToNuE   = 1<<3, ///< \f$\nu_\mu\to\nu_e\f$ ('\f$\nu_e\f$ appearance')
      kNuMuToNuMu  = 1<<4, ///< \f$\nu_\mu\to\nu_\mu\f$ ('\f$\nu_\mu\f$ survival')
      kNuMuToNuTau = 1<<5, ///< \f$\nu_\mu\to\nu_\tau\f$

      kAllNuE   = kNuEToNuE   | kNuMuToNuE,   ///< All \f$\nu_e\f$
      kAllNuMu  = kNuEToNuMu  | kNuMuToNuMu,  ///< All \f$\nu_\mu\f$
      kAllNuTau = kNuEToNuTau | kNuMuToNuTau, ///< All \f$\nu_\tau\f$

      kAll = kAllNuE | kAllNuMu | kAllNuTau   ///< All neutrinos, any flavor


*/                                 
                                 
                                 
  // Plot what we have so far
  TH1D* h_sUnosc = sUnosc.ToTH1(pot);
  //sUnosc.ToTH1(pot)->Draw("hist");
  TH1D* h_sUnoscNC = sOsc_NC.ToTH1(pot);
  TH1D* h_sUnosc_ANTInumu = sOsc_ANTInu.ToTH1(pot);
  TH1D* h_sUnosc_allnue = sOsc_ALLNue.ToTH1(pot);
  TH1D* h_sUnosc_alltao = sOsc_ALLtao.ToTH1(pot);
  
  
  THStack* hs = new THStack("hs", "");

// Fill colors to match your example
  h_sUnosc_alltao   ->SetFillColor(kCyan+1);   // ντ CC
  h_sUnoscNC        ->SetFillColor(kGreen+2);     // NC
  h_sUnosc_ANTInumu ->SetFillColor(kMagenta+1);      // anti-numu CC
  h_sUnosc_allnue   ->SetFillColor(kBlue+1);      // (νe + anti-νe) CC
  
  h_sUnosc_alltao   ->SetLineWidth(0);
  h_sUnoscNC        ->SetLineWidth(0);
  h_sUnosc_ANTInumu ->SetLineWidth(0);
  h_sUnosc_allnue  ->SetLineWidth(0);
  
  
  hs->Add(h_sUnosc_allnue);
  hs->Add(h_sUnoscNC);
  hs->Add(h_sUnosc_ANTInumu);
  hs->Add(h_sUnosc_alltao);
  //sUnoscNC.ToTH1(pot, kBlue)->Draw("hist same");
  TH1D* h_sOsc = sOsc.ToTH1(pot, kRed);
  //sOsc.ToTH1(pot, kRed)->Draw("hist same");
  TH1D* h_sOsc_fakedata = sOsc_Data.FakeData(pot).ToTH1(pot, kBlack);
  
  hs->Add(h_sOsc);
  TH1D* h_sOsc_mockdata = sOsc.MockData(pot).ToTH1(pot);
  h_sUnosc->SetTitle("applied cut: cvnnumu > 0.5 && cvnnue < 0.55"); 
  h_sOsc->SetMaximum(h_sOsc_mockdata->GetMaximum()* 1.45); 
  //h_sUnosc->Draw("hist");
  //h_sOsc->Draw("hist");
  //h_sUnoscNC->Draw("same hist");

  
  //h_sOsc_fakedata->Draw("Axis");
  //h_sOsc_mockdata->Draw("same pe");
  h_sOsc_fakedata->GetXaxis()->SetTitle("Reco E_{#nu} [GeV]");
  h_sOsc_fakedata->GetYaxis()->SetTitle("Events per 0.25 GeV");
  h_sOsc_fakedata->SetMaximum(h_sOsc_fakedata->GetMaximum()* 1.45); 
  h_sOsc_fakedata->Draw("pe");
  hs->Draw("hist same");
  //h_sOsc_fakedata->SetMarkerStyle(kFullCircle);
  h_sOsc_fakedata->SetMarkerColor(kBlack);
  h_sOsc_fakedata->SetMarkerSize(1.0);
  h_sOsc_fakedata->Draw("same pe");
 //lg1->AddEntry(h_sUnosc, "Unosc", "l" );
 //lg1->AddEntry(h_sUnoscNC, "Unosc NC", "l" );
  const double total_numu = HistRate(h_sOsc) + HistRate(h_sUnosc_ANTInumu) +
                            HistRate(h_sUnoscNC) + HistRate(h_sUnosc_alltao) +
                            HistRate(h_sUnosc_allnue);
  PrintFigureSummary("Figure 1: FD numu selection",
                     {{"Signal #nu_{#mu} CC", h_sOsc},
                      {"BG: #bar{#nu}_{#mu} CC", h_sUnosc_ANTInumu},
                      {"BG: NC", h_sUnoscNC},
                      {"BG: (#nu_{e} + #bar{#nu}_{e}) CC", h_sUnosc_allnue},
                      {"BG: (#nu_{#tau} + #bar{#nu}_{#tau}) CC", h_sUnosc_alltao}},
                     total_numu);

  lg1->AddEntry(h_sOsc_fakedata, Form("FakeData: %.1f (100.0%%)", HistRate(h_sOsc_fakedata)), "pe" );
  lg1->AddEntry(h_sOsc, MakeRateLabel("Signal #nu_{#mu} CC", h_sOsc, total_numu).c_str(), "l" );

 
  lg1->AddEntry(h_sUnosc_ANTInumu, MakeRateLabel("BG: #bar{#nu}_{#mu} CC", h_sUnosc_ANTInumu, total_numu).c_str(), "f");
  lg1->AddEntry(h_sUnoscNC,        MakeRateLabel("BG: NC", h_sUnoscNC, total_numu).c_str(), "f");
  lg1->AddEntry(h_sUnosc_allnue,   MakeRateLabel("BG: (#nu_{e} + #bar{#nu}_{e}) CC", h_sUnosc_allnue, total_numu).c_str(), "f");
  lg1->AddEntry(h_sUnosc_alltao,   MakeRateLabel("BG: (#nu_{#tau} + #bar{#nu}_{#tau}) CC", h_sUnosc_alltao, total_numu).c_str(), "f");
 
  

 
 
 
 
 //lg1->AddEntry(h_sOsc_mockdata, "mock data", "pe" );
 lg1->Draw("same");
 
     
   TPaveText* label  = DrawDUNELabel(
   true,
     "Normal",
    sin2_2theta13,
    sin2_theta23,dmsq32 ) ;
    //label->Draw("same");
  gPad->Update();
  c1 -> Print(pdf_title);

  // Extra page: same first plot, zoomed x-range
  const double total_numu_zoom = HistRateRange(h_sOsc, 0.5, 8.0) +
                                 HistRateRange(h_sUnosc_ANTInumu, 0.5, 8.0) +
                                 HistRateRange(h_sUnoscNC, 0.5, 8.0) +
                                 HistRateRange(h_sUnosc_alltao, 0.5, 8.0) +
                                 HistRateRange(h_sUnosc_allnue, 0.5, 8.0);
  PrintFigureSummary("Figure 1 (zoom 0.5-8.0): FD numu selection",
                     {{"Signal #nu_{#mu} CC", h_sOsc},
                      {"BG: #bar{#nu}_{#mu} CC", h_sUnosc_ANTInumu},
                      {"BG: NC", h_sUnoscNC},
                      {"BG: (#nu_{#tau} + #bar{#nu}_{#tau}) CC", h_sUnosc_alltao},
                      {"BG: (#nu_{e} + #bar{#nu}_{e}) CC", h_sUnosc_allnue}},
                     total_numu_zoom);

  TLegend* lg1_zoom = new TLegend(0.4, 0.5, 0.8, 0.8);
  lg1_zoom->SetNColumns(1);
  lg1_zoom->SetBorderSize(0);
  lg1_zoom->SetTextSize(.04);
  lg1_zoom->AddEntry(h_sOsc_fakedata, Form("FakeData %.1f (100.0%%)", HistRateRange(h_sOsc_fakedata, 0.5, 8.0)), "pe" );
  lg1_zoom->AddEntry(h_sOsc, MakeRateLabelRange("Signal #nu_{#mu} CC", h_sOsc, total_numu_zoom, 0.5, 8.0).c_str(), "l" );
  lg1_zoom->AddEntry(h_sUnosc_ANTInumu, MakeRateLabelRange("BG: #bar{#nu}_{#mu} CC", h_sUnosc_ANTInumu, total_numu_zoom, 0.5, 8.0).c_str(), "f");
  lg1_zoom->AddEntry(h_sUnoscNC, MakeRateLabelRange("BG: NC", h_sUnoscNC, total_numu_zoom, 0.5, 8.0).c_str(), "f");
  lg1_zoom->AddEntry(h_sUnosc_allnue, MakeRateLabelRange("BG: (#nu_{e} + #bar{#nu}_{e}) CC", h_sUnosc_allnue, total_numu_zoom, 0.5, 8.0).c_str(), "f");
  lg1_zoom->AddEntry(h_sUnosc_alltao, MakeRateLabelRange("BG: (#nu_{#tau} + #bar{#nu}_{#tau}) CC", h_sUnosc_alltao, total_numu_zoom, 0.5, 8.0).c_str(), "f");
  h_sOsc->GetXaxis()->SetRangeUser(0.5, 8.0);
  h_sOsc->GetXaxis()->SetTitle("Reco E_{#nu} [GeV]");
  h_sOsc->GetYaxis()->SetTitle("Events per 0.25 GeV");
  h_sOsc->Draw("Axis");
  h_sOsc_fakedata->Draw("same pe");
  hs->Draw("hist same");
  lg1_zoom->Draw("same");
  c1->Print(pdf_title);
  h_sOsc->GetXaxis()->SetRangeUser(0, 10.0);
  
  
///////////////////////////////////////////////////////////
   const Spectrum sOsc_nue = pred_nue.PredictComponent(calc,
                                                  Flavors::kNuMuToNuE,
                                                  Current::kCC,
                                                  Sign::kBoth);
                                                  
   const Spectrum sOsc_nue_nuonly = pred_nue.PredictComponent(calc,
                                                  Flavors::kNuMuToNuE,
                                                  Current::kCC,
                                                  Sign::kNu);                                                  
                                                  
     const Spectrum sOsc_nue_antinuonly = pred_nue.PredictComponent(calc,
                                                        Flavors::kNuMuToNuE,
                                                        Current::kCC,
                                                        Sign::kAntiNu);

  TH1D* h_sOsc_nue = sOsc_nue.ToTH1(pot, kRed);
  TH1D* h_sOsc_nue_nu_only =sOsc_nue_nuonly.ToTH1(pot, kRed+1);
  TH1D* h_sOsc_nue_antinu_only =sOsc_nue_antinuonly.ToTH1(pot, kOrange+7);

  
  const Spectrum noscillated_nue_BG_Beam = pred_nue.PredictComponent(
    calc,
    Flavors::kNuEToNuE,
    Current::kCC,
    Sign::kBoth);
  
  
  const Spectrum sOsc_NC_nue = pred_nue.PredictComponent(calc,
                                                  Flavors::kAll,
                                                  Current::kNC,
                                                  Sign::kBoth);
                                                  
  const Spectrum sOsc_ANTInu_nue = pred_nue.PredictComponent(calc,
                                                  Flavors::kAllNuMu,
                                                  Current::kCC,
                                                  Sign::kBoth);
                                                  
    //const Spectrum sOsc_ALLNue_nue = pred_nue.PredictComponent(calc,
    //                                              Flavors::kAllNuE,
    //                                              Current::kCC,
    //                                              Sign::kBoth);
                                                  
    const Spectrum sOsc_ALLtao_nue = pred_nue.PredictComponent(calc,
                                                  Flavors::kAllNuTau,
                                                  Current::kCC,
                                                  Sign::kBoth);
                                                  
  const Spectrum fakedata_nue = pred_nue.PredictComponent(calc,
        Flavors::kAll,
        Current::kBoth,
        Sign::kBoth);
        
  const Spectrum sOsc_Data_pi_2 = pred_nue.PredictComponent(calc_pi_2,
                                                  Flavors::kAll,
                                                  Current::kBoth,
                                                  Sign::kBoth);  
  const Spectrum sOsc_Data_pi_2_neg = pred_nue.PredictComponent(calc_pi_3,
                                                  Flavors::kAll,
                                                  Current::kBoth,
                                                  Sign::kBoth);                                                    
                                                      
        
        
  TH1D* h_sOsc_fakedata_nue = fakedata_nue.FakeData(pot).ToTH1(pot, kBlack);
  TH1D* h_sUnoscNC_nue = sOsc_NC_nue.ToTH1(pot);
  TH1D* h_sUnosc_numu_nue = sOsc_ANTInu_nue.ToTH1(pot);
  TH1D* h_sUnosc_allnue_nue = noscillated_nue_BG_Beam.ToTH1(pot);
  TH1D* h_sUnosc_alltao_nue = sOsc_ALLtao_nue.ToTH1(pot);
  
  h_sUnosc_allnue_nue->SetMaximum(h_sOsc_fakedata_nue->GetMaximum()* 1.6); 
  THStack* hs_nue = new THStack("hs_nue", "");

// Fill colors to match your example
  h_sUnosc_alltao_nue   ->SetFillColor(kCyan+1);   // ντ CC
  h_sUnoscNC_nue        ->SetFillColor(kGreen+2);     // NC
  h_sUnosc_numu_nue ->SetFillColor(kMagenta+1);      // anti-numu CC
  h_sUnosc_allnue_nue   ->SetFillColor(kBlue+1);      // (νe + anti-νe) CC
  
  h_sUnosc_alltao_nue   ->SetLineWidth(0);
  h_sUnoscNC_nue        ->SetLineWidth(0);
  h_sUnosc_numu_nue ->SetLineWidth(0);
  h_sUnosc_allnue_nue  ->SetLineWidth(0);
  
  
  hs_nue->Add(h_sUnosc_allnue_nue);
  hs_nue->Add(h_sUnosc_numu_nue);
  hs_nue->Add(h_sUnoscNC_nue);
  hs_nue->Add(h_sUnosc_alltao_nue);
  hs_nue->Add(h_sOsc_nue);
  
    //h_sOsc_nue->Draw("hist");
  
   //h_sOsc_nue_mockdata->Draw("same pe");
   h_sOsc_fakedata_nue->GetXaxis()->SetTitle("Reco E_{#nu} [GeV]");
   h_sOsc_fakedata_nue->GetYaxis()->SetTitle("Events per 0.25 GeV");
   h_sOsc_fakedata_nue->SetMaximum(h_sOsc_fakedata_nue->GetMaximum()* 1.45); 
   h_sOsc_fakedata_nue ->Draw("pe");
   hs_nue->Draw("hist same");
   gPad->Update();
   h_sOsc_fakedata_nue ->SetMarkerColor(kBlack);
   h_sOsc_fakedata_nue ->SetMarkerSize(1.0);
   h_sOsc_fakedata_nue ->Draw("same pe");
   
  TH1D* h_sOsc_fakedata_pi_2 = sOsc_Data_pi_2.FakeData(pot).ToTH1(pot, kBlack);
  TH1D* h_sOsc_fakedata_pi_2_neg = sOsc_Data_pi_2_neg.FakeData(pot).ToTH1(pot, kBlack);
  
  
  h_sOsc_fakedata_pi_2->SetLineStyle(4);
  h_sOsc_fakedata_pi_2_neg->SetLineStyle(5);
  
     h_sOsc_fakedata_pi_2->Draw("hist same");
     h_sOsc_fakedata_pi_2_neg->Draw("hist same");
   
  const double total_nue = HistRate(h_sOsc_nue) +
                           HistRate(h_sUnoscNC_nue) + HistRate(h_sUnosc_numu_nue) +
                           HistRate(h_sUnosc_alltao_nue)+ HistRate(h_sUnosc_allnue_nue);
  PrintFigureSummary("Figure 2: FD nue selection",
                     {{"Signal (#nu_{e} + #bar{#nu}_{e}) CC", h_sOsc_nue},
                      {"Signal #nu_{e} only CC", h_sOsc_nue_nu_only},
                      {"Signal #bar{#nu}_{e} only CC", h_sOsc_nue_antinu_only},
                      {"BG: Beam (#nu_{e} + #bar{#nu}_{e}) CC", h_sUnosc_allnue_nue},
                      {"BG: NC", h_sUnoscNC_nue},
                      {"BG: (#nu_{#mu} + #bar{#nu}_{#mu}) CC", h_sUnosc_numu_nue},
                      {"BG: (#nu_{#tau} + #bar{#nu}_{#tau}) CC", h_sUnosc_alltao_nue}},
                     total_nue);

  lg3->AddEntry(h_sOsc_fakedata_nue, Form("Fake Data (total): %.1f (100.0%%)", HistRate(h_sOsc_fakedata_nue )), "pe" );
  lg3->AddEntry(h_sOsc_nue , MakeRateLabel("Signal(#nu_{e} + #bar{#nu}_{e})CC", h_sOsc_nue, total_nue).c_str(), "l" );
  lg3->AddEntry(h_sOsc_nue_nu_only, MakeRateLabel("Signal #nu_{e} only CC", h_sOsc_nue_nu_only, total_nue).c_str(), "l" );
  lg3->AddEntry(h_sOsc_nue_antinu_only, MakeRateLabel("Signal #bar{#nu}_{e} only CC", h_sOsc_nue_antinu_only, total_nue).c_str(), "l" );
  lg3->AddEntry(h_sUnosc_allnue_nue,   MakeRateLabel("BG: Beam (#nu_{e} + #bar{#nu}_{e}) CC", h_sUnosc_allnue_nue, total_nue).c_str(), "f");
  lg3->AddEntry(h_sUnoscNC_nue,        MakeRateLabel("BG: NC", h_sUnoscNC_nue, total_nue).c_str(), "f");
  lg3->AddEntry(h_sUnosc_numu_nue, MakeRateLabel("BG:(#nu_{#mu} + #bar{#nu}_{#mu}) CC", h_sUnosc_numu_nue, total_nue).c_str(), "f");
  lg3->AddEntry(h_sUnosc_alltao_nue,   MakeRateLabel("BG:(#nu_{#tau} + #bar{#nu}_{#tau}) CC", h_sUnosc_alltao_nue, total_nue).c_str(), "f");
  lg3->AddEntry(h_sOsc_fakedata_pi_2,"cp = #pi/2"   , "l");
   lg3->AddEntry(h_sOsc_fakedata_pi_2_neg,"cp = -#pi/2"   , "l");
    TPaveText* label_nue  = DrawDUNELabel(
    false,
     "Normal",
    sin2_2theta13,
    sin2_theta23,dmsq32 ) ;
 
 lg3->Draw("same");
  gPad->Update();
  c1 -> Print(pdf_title);

  // Extra page: same second plot, zoomed x-range
  const double total_nue_zoom = HistRateRange(h_sOsc_nue, 0.5, 8.0) +
                                HistRateRange(h_sUnosc_allnue_nue, 0.5, 8.0) +
                                HistRateRange(h_sUnoscNC_nue, 0.5, 8.0) +
                                HistRateRange(h_sUnosc_numu_nue, 0.5, 8.0) +
                                HistRateRange(h_sUnosc_alltao_nue, 0.5, 8.0);
  PrintFigureSummary("Figure 2 (zoom 0.5-8.0): FD nue selection",
                     {{"Signal (#nu_{e} + #bar{#nu}_{e}) CC", h_sOsc_nue},
                      {"Signal #nu_{e} only CC", h_sOsc_nue_nu_only},
                      {"Signal #bar{#nu}_{e} only CC", h_sOsc_nue_antinu_only},
                      {"BG: Beam (#nu_{e} + #bar{#nu}_{e}) CC", h_sUnosc_allnue_nue},
                      {"BG: NC", h_sUnoscNC_nue},
                      {"BG: (#nu_{#mu} + #bar{#nu}_{#mu}) CC", h_sUnosc_numu_nue},
                      {"BG: (#nu_{#tau} + #bar{#nu}_{#tau}) CC", h_sUnosc_alltao_nue}},
                     total_nue_zoom);

  TLegend* lg3_zoom = new TLegend(0.42, 0.6, 0.8, 0.84);
  lg3_zoom->SetNColumns(1);
  lg3_zoom->SetBorderSize(0);
  lg3_zoom->SetTextSize(.03);
  lg3_zoom->AddEntry(h_sOsc_fakedata_nue, Form("Fake Data (total): %.1f (100.0%%)", HistRateRange(h_sOsc_fakedata_nue, 0.5, 8.0)), "pe" );
  lg3_zoom->AddEntry(h_sOsc_nue, MakeRateLabelRange("Signal(#nu_{e} + #bar{#nu}_{e})CC", h_sOsc_nue, total_nue_zoom, 0.5, 8.0).c_str(), "l");
  lg3_zoom->AddEntry(h_sOsc_nue_nu_only, MakeRateLabelRange("Signal #nu_{e} only CC", h_sOsc_nue_nu_only, total_nue_zoom, 0.5, 8.0).c_str(), "l");
  lg3_zoom->AddEntry(h_sOsc_nue_antinu_only, MakeRateLabelRange("Signal #bar{#nu}_{e} only CC", h_sOsc_nue_antinu_only, total_nue_zoom, 0.5, 8.0).c_str(), "l");
  lg3_zoom->AddEntry(h_sUnosc_allnue_nue, MakeRateLabelRange("BG: Beam (#nu_{e} + #bar{#nu}_{e}) CC", h_sUnosc_allnue_nue, total_nue_zoom, 0.5, 8.0).c_str(), "f");
  lg3_zoom->AddEntry(h_sUnoscNC_nue, MakeRateLabelRange("BG:NC", h_sUnoscNC_nue, total_nue_zoom, 0.5, 8.0).c_str(), "f");
  lg3_zoom->AddEntry(h_sUnosc_numu_nue, MakeRateLabelRange("BG:(#nu_{#mu} + #bar{#nu}_{#mu}) CC", h_sUnosc_numu_nue, total_nue_zoom, 0.5, 8.0).c_str(), "f");
  lg3_zoom->AddEntry(h_sUnosc_alltao_nue, MakeRateLabelRange("BG:(#nu_{#tau} + #bar{#nu}_{#tau}) CC", h_sUnosc_alltao_nue, total_nue_zoom, 0.5, 8.0).c_str(), "f");
  lg3_zoom->AddEntry(h_sOsc_fakedata_pi_2,"cp = #pi/2"   , "l");
  lg3_zoom->AddEntry(h_sOsc_fakedata_pi_2_neg,"cp = -#pi/2"   , "l");
  h_sOsc_nue->GetXaxis()->SetRangeUser(0.5, 8.0);
  h_sOsc_nue->SetMaximum(h_sOsc_nue->GetMaximum()* 1.55); 
  h_sOsc_nue->GetXaxis()->SetTitle("Reco E_{#nu} [GeV]");
  h_sOsc_nue->GetYaxis()->SetTitle("Events per 0.25 GeV");
  h_sOsc_nue->Draw("Axis");
  h_sOsc_fakedata_nue->Draw("same pe");
   hs_nue->Draw("hist same");
       h_sOsc_fakedata_pi_2->Draw("hist same");
     h_sOsc_fakedata_pi_2_neg->Draw("hist same");
  lg3_zoom->Draw("same");


  //hs_nue->GetXaxis()->SetRangeUser(0.5, 8.0);
    //h_sOsc_nue->SetMaximum(h_sOsc_nue_fakedata->GetMaximum()* 1.8); 
  //h_sOsc_nue->Draw("hist");
  //
  //hs_nue->Draw("hist");
  //h_sOsc_fakedata_nue->Draw("same pe");
  //lg3_zoom->Draw("same");
  //label_nue->Draw("same");
  c1->Print(pdf_title);
  //h_sOsc_nue->GetXaxis()->SetRangeUser(0, 10.0);



gPad->SetLogy();
///////////////////////////////////////////
const Spectrum sOsc_num_signal_cvn = pred_nusignal.PredictComponent(calc,
                                                  Flavors::kAllNuMu,
                                                  Current::kCC,
                                                  Sign::kBoth);

const Spectrum sOsc_num_NC_cvn = pred_nusignal.PredictComponent(calc,
                                                   Flavors::kAll,
                                                  Current::kNC,
                                                  Sign::kBoth);

const Spectrum sOsc_num_tao_cvn = pred_nusignal.PredictComponent(calc,
                                                   Flavors::kAllNuTau,
                                                  Current::kCC,
                                                  Sign::kBoth);   



TH1D* h_sOsca_CVN_nusignal = sOsc_num_signal_cvn.ToTH1(pot, kBlue);
TH1D* h_sUnosc_CVN_NC = sOsc_num_NC_cvn.ToTH1(pot);
TH1D* h_sUnosc_CVN_taC = sOsc_num_tao_cvn.ToTH1(pot,kGreen);


//TH1D* h_sOsca_CVN_nusignal = sOsc_nusignal.ToTH1(pot, kBlue);
//TH1D* h_sUnosc_CVN_NC = sOsc_NC.ToTH1(pot);
//TH1D* h_sUnosc_CVN_taC = sOsc_ta.ToTH1(pot,kGreen);
h_sOsca_CVN_nusignal->SetMinimum(2);
h_sOsca_CVN_nusignal->Draw("hist");
  h_sUnosc_CVN_NC->Draw("same hist");
  h_sUnosc_CVN_taC->Draw("same hist");
 const double total_cvn_numu = HistRate(h_sOsca_CVN_nusignal) + HistRate(h_sUnosc_CVN_NC) + HistRate(h_sUnosc_CVN_taC);
 PrintFigureSummary("Figure 3: CVN numu composition",
                    {{"(#nu_{#mu} + #bar{#nu}_{#mu}) CC", h_sOsca_CVN_nusignal},
                     {"NC", h_sUnosc_CVN_NC},
                     {"(#nu_{#tau} + #bar{#nu}_{#tau}) CC", h_sUnosc_CVN_taC}},
                    total_cvn_numu);
 lg2->AddEntry(h_sOsca_CVN_nusignal, MakeRateLabel("(#nu_{#mu} + #bar{#nu}_{#mu}) CC", h_sOsca_CVN_nusignal, total_cvn_numu).c_str(), "l" );
 lg2->AddEntry(h_sUnosc_CVN_NC, MakeRateLabel("NC", h_sUnosc_CVN_NC, total_cvn_numu).c_str(), "l" );
 lg2->AddEntry(h_sUnosc_CVN_taC, MakeRateLabel("(#nu_{#tau} + #bar{#nu}_{#tau}) CC", h_sUnosc_CVN_taC, total_cvn_numu).c_str(), "l" );
 lg2->Draw("same");
  c1 -> Print(pdf_title);
  
  
  ///////////////////////////////////////////
const Spectrum sOsc_num_numu_cvn_nue =   pred_CVN_nue.PredictComponent(calc,
                                                  Flavors::kAllNuMu,
                                                  Current::kCC,
                                                  Sign::kBoth);

const Spectrum sOsc_num_NC_cvn_nue =   pred_CVN_nue.PredictComponent(calc,
                                                   Flavors::kAll,
                                                  Current::kNC,
                                                  Sign::kBoth);

const Spectrum sOsc_num_tao_cvn_nue =   pred_CVN_nue.PredictComponent(calc,
                                                   Flavors::kAllNuTau,
                                                  Current::kCC,
                                                  Sign::kBoth);   

const Spectrum sOsc_num_nue_cvn_nue =   pred_CVN_nue.PredictComponent(calc,
                                                  Flavors::kAllNuE,
                                                  Current::kCC,
                                                  Sign::kBoth);


TH1D* h_sOsca_CVN_nusignal_nue = sOsc_num_numu_cvn_nue.ToTH1(pot, kBlue);
TH1D* h_sUnosc_CVN_NC_nue = sOsc_num_NC_cvn_nue.ToTH1(pot);
TH1D* h_sUnosc_CVN_taC_nue = sOsc_num_tao_cvn_nue.ToTH1(pot,kGreen);
TH1D* h_sUnosc_CVN_nue_nue = sOsc_num_nue_cvn_nue.ToTH1(pot,kRed);

//TH1D* h_sOsca_CVN_nusignal = sOsc_nusignal.ToTH1(pot, kBlue);
//TH1D* h_sUnosc_CVN_NC = sOsc_NC.ToTH1(pot);
//TH1D* h_sUnosc_CVN_taC = sOsc_ta.ToTH1(pot,kGreen);
h_sOsca_CVN_nusignal_nue->SetMinimum(2);
h_sOsca_CVN_nusignal_nue->Draw("hist");
  h_sUnosc_CVN_NC_nue->Draw("same hist");
  h_sUnosc_CVN_taC_nue->Draw("same hist");
  h_sUnosc_CVN_nue_nue->Draw("same hist");
 const double total_cvn_nue = HistRate(h_sUnosc_CVN_nue_nue) + HistRate(h_sOsca_CVN_nusignal_nue) +
                              HistRate(h_sUnosc_CVN_NC_nue) + HistRate(h_sUnosc_CVN_taC_nue);
 PrintFigureSummary("Figure 4: CVN nue composition",
                    {{"Signal CC (#nu_{e} + #bar{#nu}_{e})", h_sUnosc_CVN_nue_nue},
                     {"(#nu_{#mu} + #bar{#nu}_{#mu}) CC", h_sOsca_CVN_nusignal_nue},
                     {"NC", h_sUnosc_CVN_NC_nue},
                     {"(#nu_{#tau} + #bar{#nu}_{#tau}) CC", h_sUnosc_CVN_taC_nue}},
                    total_cvn_nue);
 lg4->AddEntry(h_sUnosc_CVN_nue_nue, MakeRateLabel("Signal CC (#nu_{e} + #bar{#nu}_{e})", h_sUnosc_CVN_nue_nue, total_cvn_nue).c_str(), "l" );
 lg4->AddEntry(h_sOsca_CVN_nusignal_nue, MakeRateLabel("(#nu_{#mu} + #bar{#nu}_{#mu}) CC", h_sOsca_CVN_nusignal_nue, total_cvn_nue).c_str(), "l" );
 lg4->AddEntry(h_sUnosc_CVN_NC_nue, MakeRateLabel("NC", h_sUnosc_CVN_NC_nue, total_cvn_nue).c_str(), "l" );
 lg4->AddEntry(h_sUnosc_CVN_taC_nue, MakeRateLabel("(#nu_{#tau} + #bar{#nu}_{#tau}) CC", h_sUnosc_CVN_taC_nue, total_cvn_nue).c_str(), "l" );
 lg4->Draw("same");
  c1 -> Print(pdf_title);
  
  


  
  //DataMCComparisonComponents(sUnoscNC,
  //                                &pred,
  //                                calc);
 //
 
  // "Fake" data is synonymous with the Asimov data sample
  //new TCanvas;
  //sOsc.ToTH1(pot, kRed)->Draw("hist");
  //sUnoscNC.ToTH1(pot, kBlue)->Draw("hist same");
  //sOsc.FakeData(pot).ToTH1(pot)->Draw("ep same");

  // While "mock" data has statistical fluctuations in
  //new TCanvas;
  //sOsc.ToTH1(pot, kRed)->Draw("hist");
  //sUnoscNC.ToTH1(pot, kBlue)->Draw("hist same");
  //sOsc.MockData(pot).ToTH1(pot)->Draw("ep same");
  
  
  
  sprintf(pdf_title, "%s.pdf)", Pdf_name.c_str());
  c1 -> Print(pdf_title);
  
}
