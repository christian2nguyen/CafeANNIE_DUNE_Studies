// Make oscillated predictions
// cafe demo1.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
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

using namespace ana;

void demo1()
{ 
  TLegend* lg1 = new TLegend( 0.35, 0.7, 0.8, 0.88 );
    lg1->SetNColumns(2);
    lg1->SetBorderSize(0);
    lg1->SetTextSize(.04); //  
    
 char pdf_title[1024];
 std::string Pdf_name = "DUNEexample_demo1";
  TCanvas *c1 = new TCanvas("c1");
 sprintf(pdf_title, "%s.pdf(", Pdf_name.c_str());
 c1 -> Print(pdf_title);
 sprintf(pdf_title, "%s.pdf", Pdf_name.c_str());

  // See demo0.C for explanation of these repeated parts
  const std::string fnameNonSwap = "/pnfs/dune/persistent/users/cnguyen/cafana_example/FD_FHC_nonswap.root";
  const std::string fnameNueSwap = "/pnfs/dune/persistent/users/cnguyen/cafana_example/FD_FHC_nueswap.root";
  const std::string fnameTauSwap = "/pnfs/dune/persistent/users/cnguyen/cafana_example/FD_FHC_tauswap.root";
  SpectrumLoader loaderNonSwap(fnameNonSwap);
  SpectrumLoader loaderNueSwap(fnameNueSwap);
  SpectrumLoader loaderTauSwap(fnameTauSwap);
  const Var kRecoEnergy = SIMPLEVAR(dune.Ev_reco_numu);
  const Var kCVNNumu = SIMPLEVAR(dune.cvnnumu);
  const Binning binsEnergy = Binning::Simple(40, 0, 10);
  const HistAxis axEnergy("Reco #nu_{Energy} [GeV]", binsEnergy, kRecoEnergy);
  const double pot = 3.5 * 1.47e21 * 40/1.13;

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
                          axEnergy, kPassesCVN);

  // These calls will fill all of the constituent parts of the prediction
  loaderNonSwap.Go();
  loaderNueSwap.Go();
  loaderTauSwap.Go();

  // We can extract a total prediction unoscillated
  const Spectrum sUnosc = pred.PredictUnoscillated();
  // Or oscillated, in this case using reasonable parameters from
  // Analysis/Calcs.h
  osc::IOscCalculator* calc = DefaultOscCalc();
  const Spectrum sOsc = pred.Predict(calc);

  // And we can break things down by flavour
  const Spectrum sUnoscNC = pred.PredictComponent(calc,
                                                  Flavors::kAll,
                                                  Current::kNC,
                                                  Sign::kBoth);

  // Plot what we have so far
  TH1D* h_sUnosc = sUnosc.ToTH1(pot);
  //sUnosc.ToTH1(pot)->Draw("hist");
  TH1D* h_sUnoscNC = sUnoscNC.ToTH1(pot, kBlue);
  //sUnoscNC.ToTH1(pot, kBlue)->Draw("hist same");
  TH1D* h_sOsc = sOsc.ToTH1(pot, kRed);
  //sOsc.ToTH1(pot, kRed)->Draw("hist same");
  TH1D* h_sOsc_fakedata = sOsc.FakeData(pot).ToTH1(pot);
  
  TH1D* h_sOsc_mockdata = sOsc.MockData(pot).ToTH1(pot);
  h_sUnosc->SetTitle("applied cut: cvnnumu > 0.5"); 
  h_sUnosc->SetMaximum(h_sUnosc->GetMaximum()* 1.55); 
  h_sUnosc->Draw("hist");
  h_sUnoscNC->Draw("same hist");
  h_sOsc->Draw("same hist");
  
  h_sOsc_fakedata->Draw("same pe");
  h_sOsc_mockdata->Draw("same pe");

 lg1->AddEntry(h_sUnosc, "Unosc", "l" );
 lg1->AddEntry(h_sUnoscNC, "Unosc NC", "l" );
 lg1->AddEntry(h_sOsc, "Osc", "l" );
 lg1->AddEntry(h_sOsc_fakedata, "fake data", "pe" );
 lg1->AddEntry(h_sOsc_mockdata, "mock data", "pe" );
 lg1->Draw("same");
  c1 -> Print(pdf_title);
  gPad->SetLogy();

  DataMCComparisonComponents(sUnoscNC,
                                  &pred,
                                  calc);
 
 
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
