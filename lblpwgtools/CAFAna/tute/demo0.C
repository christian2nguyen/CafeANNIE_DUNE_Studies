// Make a simple spectrum plot
// cafe demo0.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"

#include "CAFAna/Cuts/TruthCuts.h"

#include "StandardRecord/StandardRecord.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TPad.h"
#include <string>
#include "TAxis.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"
using namespace ana;

void demo0()
{

    TLegend* lg1 = new TLegend( 0.35, 0.7, 0.8, 0.88 );
    lg1->SetNColumns(2);
    lg1->SetBorderSize(0);
    lg1->SetTextSize(.04); //
    
    
    
    TLegend* lg2 = new TLegend( 0.35, 0.7, 0.8, 0.88 );
    lg2->SetNColumns(2);
    lg2->SetBorderSize(0);
    lg2->SetTextSize(.04); //
 char pdf_title[1024];
 std::string Pdf_name = "DUNEexample_demo0";
  TCanvas *c1 = new TCanvas("c1");
 sprintf(pdf_title, "%s.pdf(", Pdf_name.c_str());
 c1 -> Print(pdf_title);
 sprintf(pdf_title, "%s.pdf", Pdf_name.c_str());
  // Environment variables and wildcards work. As do SAM datasets.
  const std::string fname = "/pnfs/dune/persistent/users/cnguyen/cafana_example/FD_FHC_nonswap.root";

  // Source of events
  SpectrumLoader loader(fname);

  // A Var is a little snippet of code that takes a record representing the
  // event record and returns a single number to plot.
  const Var kRecoEnergy({},
                        [](const caf::StandardRecord* sr)
                        {
                          return sr->dune.Ev_reco_numu;
                        });

  // For such a simple variable you can use a shortcut like this
  const Var kCVNNumu = SIMPLEVAR(dune.cvnnumu);

  // Define a spectrum, ie a histogram with associated POT information
  const Binning binsEnergy = Binning::Simple(40, 0, 10);
  const HistAxis axEnergy("Reco energy (GeV)", binsEnergy, kRecoEnergy);
  // kIsNumuCC here is a "Cut". Same as a Var but returning a boolean. In this
  // case, we're only keeping events that are truly numu CC interactions.
  Spectrum sEnergy(loader, axEnergy, kIsNumuCC);

  Spectrum sEnergyNC(loader, axEnergy, kIsNC);

  // And another
  const Binning binsCVN = Binning::Simple(50, 0, 1);
  const HistAxis axCVN("CVN_{#mu}", binsCVN, kCVNNumu);
  Spectrum sCVN(loader, axCVN, kIsNumuCC);
  Spectrum sCVNNC(loader, axCVN, kIsNC);

  // This is the call that actually fills in those spectra
  loader.Go();

  // POT/yr * 3.5yrs * mass correction for the workspace geometry
  const double pot = 3.5 * 1.47e21 * 40/1.13;

  // For plotting purposes we can convert to TH1s
  TH1D* h_ccnumu = sEnergy.ToTH1(pot); //sEnergy.ToTH1(pot)->Draw("hist");
  TH1D* h_NC = sEnergyNC.ToTH1(pot, kBlue); //sEnergyNC.ToTH1(pot, kBlue)->Draw("hist same");
  
 lg1->AddEntry(h_ccnumu, "CC#nu_{#mu}", "l" );
 lg1->AddEntry(h_NC, "NC", "l" );
  
  h_ccnumu->SetMaximum(h_ccnumu->GetMaximum()* 1.45); 
  h_ccnumu->Draw("hist");
  h_NC->Draw("hist same");
  lg1->Draw("same");
  c1 -> Print(pdf_title);
  
  //sCVN.ToTH1(pot)->Draw("hist");
  //sCVNNC.ToTH1(pot, kBlue)->Draw("hist same");
  
  
  TH1D* h_sCVN = sCVN.ToTH1(pot);
  TH1D* h_sCVNNC = sCVNNC.ToTH1(pot, kBlue);
  h_sCVN->SetMaximum(h_sCVN->GetMaximum()* 1.45); 
  h_sCVN->Draw("hist");
  h_sCVNNC->Draw("hist same");
  lg2->AddEntry(h_sCVN, "CVN", "l" );
  lg2->AddEntry(h_sCVNNC, "CVN NC", "l" );
  lg2->Draw("same");
  gPad->SetLogy();
  c1 -> Print(pdf_title);
  
 sprintf(pdf_title, "%s.pdf)", Pdf_name.c_str());
 c1 -> Print(pdf_title);
  
}
