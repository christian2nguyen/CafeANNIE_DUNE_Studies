// Make a simple spectrum plot
// cafe demo0.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"

#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Cuts/AnaCuts.h"

#include "StandardRecord/StandardRecord.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TPad.h"
#include <string>
#include "TAxis.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TColor.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"

using namespace ana;

void Input_distributions_forOscinputs()
{

    TLegend* lg1 = new TLegend( 0.35, 0.7, 0.8, 0.88 );
    lg1->SetNColumns(2);
    lg1->SetBorderSize(0);
    lg1->SetTextSize(.04); //
    
    
    
    TLegend* lg2 = new TLegend( 0.35, 0.7, 0.8, 0.88 );
    lg2->SetNColumns(2);
    lg2->SetBorderSize(0);
    lg2->SetTextSize(.04); //
    
    
        TLegend* lg3 = new TLegend( 0.35, 0.7, 0.8, 0.88 );
    lg3->SetNColumns(2);
    lg3->SetBorderSize(0);
    lg3->SetTextSize(.04); //
    
            TLegend* lg4 = new TLegend( 0.35, 0.7, 0.8, 0.88 );
    lg4->SetNColumns(2);
    lg4->SetBorderSize(0);
    lg4->SetTextSize(.04); //
    
                TLegend* lg5 = new TLegend( 0.35, 0.7, 0.8, 0.88 );
    lg5->SetNColumns(2);
    lg5->SetBorderSize(0);
    lg5->SetTextSize(.04); //
    
    
 char pdf_title[1024];
 std::string Pdf_name = "Inputs_distribution_Theia_OscAnaylsis_Orginalfiles";
  TCanvas *c1 = new TCanvas("c1");
 sprintf(pdf_title, "%s.pdf(", Pdf_name.c_str());
 c1 -> Print(pdf_title);
 sprintf(pdf_title, "%s.pdf", Pdf_name.c_str());
  // Environment variables and wildcards work. As do SAM datasets.
  const std::string fname = "/pnfs/dune/persistent/users/cnguyen/old_statefiles_cp_uboonearea/FD_FHC_nonswap.root";
const std::string fname5 =  "/pnfs/dune/persistent/users/cnguyen/old_statefiles_cp_uboonearea/FD_FHC_tauswap.root";
const std::string fname3 =  "/pnfs/dune/persistent/users/cnguyen/old_statefiles_cp_uboonearea/FD_FHC_nueswap.root";

 const std::string fname1 =  "/pnfs/dune/persistent/users/cnguyen/old_statefiles_cp_uboonearea/FD_FHC_nonswap.root";
 const std::string fname2 =  "/pnfs/dune/persistent/users/cnguyen/old_statefiles_cp_uboonearea/FD_RHC_nonswap.root";
 
 const std::string fname4 =  "/pnfs/dune/persistent/users/cnguyen/old_statefiles_cp_uboonearea/FD_RHC_nueswap.root";
 
 const std::string fname6 =  "/pnfs/dune/persistent/users/cnguyen/old_statefiles_cp_uboonearea/FD_RHC_tauswap.root";


  // Source of events
  SpectrumLoader loader(fname);
  SpectrumLoader loadertao(fname5);
  SpectrumLoader loadernue(fname3);
  
  

  // A Var is a little snippet of code that takes a record representing the
  // event record and returns a single number to plot.
  const Var kRecoEnergy({},
                        [](const caf::StandardRecord* sr)
                        {
                          return sr->dune.Ev_reco_numu;
                        });

 const Var kRecoEnergy_nue({},
                        [](const caf::StandardRecord* sr)
                        {
                          return sr->dune.Ev_reco_nue;
                        });

 const Var kEnu({},
                        [](const caf::StandardRecord* sr)
                        {
                          return sr->dune.Ev;
                        });

 const Var klepE_numu({},
                        [](const caf::StandardRecord* sr)
                        {
                          return sr->dune.Elep_reco;
                        });
                    
                    
 const Var klepAng_numu({},
                        [](const caf::StandardRecord* sr)
                        {
                          return sr->dune.theta_reco;
                        });
                      
 const Var klepEhad_numu({},
                        [](const caf::StandardRecord* sr)
                        {
                          return sr->dune.RecoHadEnNumu;
                        });
/*

 const Var klepE_nue({},
                        [](const caf::StandardRecord* sr)
                        {
                          return sr->dune.RecoLepEnNue;
                        });
 const Var klepAng_nue({},
                        [](const caf::StandardRecord* sr)
                        {
                          return sr->dune.RecoLepAngNue;
                        });
  
   const Var klepEhad_nue({},
                        [](const caf::StandardRecord* sr)
                        {
                          return sr->dune.RecoHadEnNue;
                        });
*/

  // For such a simple variable you can use a shortcut like this
  const Var kCVNNumu = SIMPLEVAR(dune.cvnnumu);
  const Var kCVNNue = SIMPLEVAR(dune.cvnnue);
  // Define a spectrum, ie a histogram with associated POT information
  const Binning binsEnergy = Binning::Simple(32, 0, 8);
  const Binning binsAngle = Binning::Simple(40, 0, 3.14);
  const HistAxis axEnergy("Reco energy (E_{#nu})(GeV)", binsEnergy, kRecoEnergy);
  const HistAxis axEnergy_nue("Reco energy (E_{#nu}) (GeV)", binsEnergy, kRecoEnergy_nue);
  
  const HistAxis axEnergy_LepE("Reco Muon Energy (GeV) [PassFD_CVN_NUMU]", binsEnergy, klepE_numu);
  const HistAxis axEnergy_LepAngle("Reco Muon Angle )(rad)[PassFD_CVN_NUMU]", binsAngle, klepAng_numu);
  const HistAxis axEnergy_LepEHad("Reco Muon E Had )(GeV)", binsEnergy, klepEhad_numu);
  // kIsNumuCC here is a "Cut". Same as a Var but returning a boolean. In this
  // case, we're only keeping events that are truly numu CC interactions.
  Spectrum sEnergy(loader, axEnergy, kIsNumuCC);
  Spectrum sEnergyNC(loader, axEnergy, kIsNC);
  Spectrum sEnergy_nue(loader, axEnergy_nue, kIsBeamNue);
  
  Spectrum sEnergy_NDFHC(loader, axEnergy, kPassND_FHC_NUMU);
  Spectrum sEnergy_NDRHC(loader, axEnergy, kPassND_RHC_NUMU);
  
  Spectrum sEnergy_FDNUMU(loader, axEnergy, kPassFD_CVN_NUMU);
  Spectrum sEnergy_FDNUE(loader, axEnergy, kPassFD_CVN_NUE);
  
  Spectrum sEnergy_FDNUMU_LepE_FD(loader, axEnergy_LepE, kPassFD_CVN_NUMU);
  Spectrum sEnergy_FDNUMU_LepAngle_FD(loader, axEnergy_LepAngle, kPassFD_CVN_NUMU);
  Spectrum sEnergy_FDNUMU_LepEhad_FD(loader, axEnergy_LepEHad, kPassFD_CVN_NUMU);
  

  
  //Spectrum sEnergy_FDNUMU_LepE_FD_FHC(loader, axEnergy_LepE, kPassND_FHC_NUMU);
  //Spectrum sEnergy_FDNUMU_LepAngle_FD_FHC(loader, axEnergy_LepE, kPassND_FHC_NUMU);
  //Spectrum sEnergy_FDNUMU_LepEhad_FD_FHC(loader, axEnergy_LepEHad, kPassND_FHC_NUMU);
  //
  //Spectrum sEnergy_FDNUMU_LepE_FD_RHC(loader, axEnergy_LepE, kPassND_RHC_NUMU);
  //Spectrum sEnergy_FDNUMU_LepAngle_FD_RHC(loader, axEnergy_LepE, kPassND_RHC_NUMU);
  //Spectrum sEnergy_FDNUMU_LepEhad_FD_RHC(loader, axEnergy_LepEHad, kPassND_RHC_NUMU);
  
  
  
  
  // And anothe
  const Binning binsCVN = Binning::Simple(50, 0, 1);
  const HistAxis axCVN("CVN_{#mu}", binsCVN, kCVNNumu);
  const HistAxis axCVN_nue("CVN_{e}", binsCVN, kCVNNue);
  
  PredictionNoExtrap pred(loader, loadertao, loadernue,
                          axCVN, kIsTauFromMu);
  
  const Spectrum sCVNNC_tau = pred.PredictUnoscillated();
  
  Spectrum sCVN(loader, axCVN, kIsNumuCC);
  Spectrum sCVNNC(loader, axCVN, kIsNC);
  //Spectrum sCVNNC_tau(loader, axCVN, kIsTauFromMu);
  
  Spectrum sCVNnue_nue(loader, axCVN_nue, kIsBeamNue);
  Spectrum sCVNnue_numuCC(loader, axCVN_nue, kIsNumuCC);
  Spectrum sCVNnue_NC(loader, axCVN_nue, kIsNC);
  // This is the call that actually fills in those spectra
  loader.Go();

  // POT/yr * 3.5yrs * mass correction for the workspace geometry
  const double pot = 3.5 * 1.47e21 * 40/1.13;
  gPad->SetLogy();
  // For plotting purposes we can convert to TH1s
  TH1D* h_ccnumu = sEnergy.ToTH1(pot); //sEnergy.ToTH1(pot)->Draw("hist");
  TH1D* h_NC = sEnergyNC.ToTH1(pot, kBlue); //sEnergyNC.ToTH1(pot, kBlue)->Draw("hist same");
  TH1D* h_ccnue = sEnergy_nue.ToTH1(pot, kGreen); 
 lg1->AddEntry(h_ccnumu, "CC-#nu_{#mu}", "l" );
 lg1->AddEntry(h_ccnue , "CC-#nu_{e}", "l" );
 lg1->AddEntry(h_NC, "NC", "l" );
  h_ccnumu->Scale(1,"width");
  h_NC->Scale(1,"width");
  h_ccnue->Scale(1,"width");
  h_ccnumu->GetYaxis()->SetTitle("Events/.25 Gev"); 
  h_ccnumu->SetMaximum(h_ccnumu->GetMaximum()* 1.45); 
  h_ccnumu->GetYaxis()->SetNoExponent(false);
  h_ccnumu->GetYaxis()->SetMaxDigits(3);
  h_ccnumu->Draw("hist");
  h_NC->Draw("hist same");
  h_ccnue->Draw("hist same");
  lg1->Draw("same");
  c1 -> Print(pdf_title);
  /////////////////////////////////////
  TH1D* h_NDFHC_numu =sEnergy_NDFHC.ToTH1(pot);
  TH1D* h_NDRHC_numu =sEnergy_NDRHC.ToTH1(pot, kBlue);
  lg3->AddEntry(h_NDFHC_numu, "CC-#nu_{#mu} (ND)FHC", "l" );
  lg3->AddEntry(h_NDRHC_numu , "CC-#nu_{#mu}(ND) RHC", "l" );
  //h_NDFHC_numu->Scale(1,"width");
  //h_NDRHC_numu->Scale(1,"width");
  h_NDFHC_numu->Draw("hist");
  h_NDRHC_numu->Draw("hist same");
  lg3->Draw("same");
  c1 -> Print(pdf_title);
  
    /////////////////////////////////////
      
  TH1D* h_FDFHC_numu =sEnergy_FDNUMU.ToTH1(pot);
  TH1D* h_FD_nue =sEnergy_FDNUE.ToTH1(pot, kBlue);
  lg4->AddEntry(h_FDFHC_numu, "CC-#nu_{#mu} (FD)", "l" );
  lg4->AddEntry(h_FD_nue , "CC-#nu_{e}(FD) ", "l" );
  //h_FDFHC_numu->Scale(1,"width");
  //h_FD_nue->Scale(1,"width");
  h_FDFHC_numu->Draw("hist");
  h_FD_nue->Draw("hist same");
  lg4->Draw("same");
  c1 -> Print(pdf_title);
  //sCVN.ToTH1(pot)->Draw("hist");
  //sCVNNC.ToTH1(pot, kBlue)->Draw("hist same");

  int myGreen = TColor::GetColor("#578A4A");
  
  TH1D* h_sCVN = sCVN.ToTH1(pot, kBlue);
  //TH1D* h_sCVNnue = sCVNnue.ToTH1(pot,kGreen);
  TH1D* h_sCVNNC = sCVNNC.ToTH1(pot);
  TH1D* h_sCVN_tau = sCVNNC_tau.ToTH1(pot, myGreen);
  h_sCVN->SetMaximum(h_sCVN->GetMaximum()* 1.45); 
  h_sCVN->SetMinimum(3);
  h_sCVN->GetYaxis()->SetTitle("Events");
  h_sCVN->GetXaxis()->SetTitle("CVN  #nu_{#mu} Probability");
  h_sCVN->Draw("hist");
  //h_sCVNnue->Draw("hist same");
  h_sCVNNC->Draw("hist same");
  h_sCVN_tau->Draw("hist same");
  lg2->AddEntry(h_sCVN, "CC-#nu_{#mu} signal", "l" );
  //lg2->AddEntry(h_sCVNnue, "CVN_{e}", "l" );
  lg2->AddEntry(h_sCVN_tau, "#nu_{#tau} background", "l" );
  lg2->AddEntry(h_sCVNNC, "NC #nu background", "l" );
  lg2->Draw("same");
  c1 -> Print(pdf_title);
  ////////////////////////////////////////
  
  
  TH1D* hsCVN_nue_nue=sCVNnue_nue.ToTH1(pot, kRed);
  TH1D* hsCVN_nue_numuCC=sCVNnue_numuCC.ToTH1(pot, kBlue);
  TH1D* hsCVN_nue_NC=sCVNnue_NC.ToTH1(pot);
  
  lg5->AddEntry(hsCVN_nue_nue, "#nu_{e} signal", "l" );
  lg5->AddEntry(hsCVN_nue_numuCC, "#nu_{#mu} background", "l" );
  lg5->AddEntry(hsCVN_nue_NC, "NC  background", "l" );
  hsCVN_nue_nue->Draw("hist ");
  hsCVN_nue_numuCC->Draw("hist same");
  hsCVN_nue_NC->Draw("hist same");
    lg5->Draw("same");
  
  c1 -> Print(pdf_title);
  
 TH1D* h_LepE_FD =  sEnergy_FDNUMU_LepE_FD.ToTH1(pot);

 h_LepE_FD->Draw("hist");
 h_LepE_FD->Scale(1,"width");
 h_LepE_FD->GetYaxis()->SetTitle("Events/.25 Gev"); 
 c1 -> Print(pdf_title);
 
 TH1D* h_LepAngle_FD =  sEnergy_FDNUMU_LepAngle_FD.ToTH1(pot);
  h_LepAngle_FD->Draw("hist");
  h_LepE_FD->Scale(1,"width");
 c1 -> Print(pdf_title);
// 
 TH1D* h_LepEhad_FD =  sEnergy_FDNUMU_LepEhad_FD.ToTH1(pot);
  h_LepEhad_FD->Draw("hist");
 c1 -> Print(pdf_title);
  
  
  
 sprintf(pdf_title, "%s.pdf)", Pdf_name.c_str());
 c1 -> Print(pdf_title);
  
}
