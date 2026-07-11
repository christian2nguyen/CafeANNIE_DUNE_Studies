// Make a simple contour
// cafe demo3.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Analysis/Calcs.h"
//#include "CAFAna/Analysis/TDRLoaders.h"
#include "OscLib/func/OscCalculatorPMNSOpt.h"
#include "StandardRecord/StandardRecord.h"
#include "TCanvas.h"
#include "TH1.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Vars/FitVars.h"

// New includes
#include "CAFAna/Analysis/Surface.h"
#include "CAFAna/Experiment/MultiExperiment.h"


#include "TPad.h"
#include <string>
#include "TAxis.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"
#include <vector>


// std
#include <memory>
#include <vector>

// ROOT core
#include "TROOT.h"
#include "TObject.h"
#include "TString.h"
#include "TClass.h"

// Legend & containers
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TList.h"
#include "TObjArray.h"
#include "TLatex.h"
// Style traits
#include "TAttLine.h"
#include "TAttMarker.h"
#include "TAttFill.h"
#include "TLine.h"

// The actual drawables you touch
#include "TH2.h"       // fixes: incomplete type 'TH2'
#include "TH2D.h"
#include "TGraph.h"    // fixes: incomplete type 'TGraph'
#include "TPolyLine.h" // fixes: undeclared identifier 'TPolyLine'
#include "TBox.h"
#include "TMath.h"


using namespace ana;


// =========================================================
// AddAnyLegendEntryStyled: legend entry with explicit style
// - Works for TH2 (contours), graphs, histos, custom classes
// - Uses a clone/proxy so legend style won't change your plot
// =========================================================

struct LegendDetect {
  enum Kind { kPoint, kLine, kFill, kUnknown } kind = kUnknown;
  int lineColor = kBlack, lineStyle = 1, lineWidth = 1;
  int markerStyle = 0, markerColor = kBlack; double markerSize = 0;
  int fillColor = 0, fillStyle = 0;
};

template <class T> static void _copyLine(const TObject* o, LegendDetect& d) {
  if (auto t = dynamic_cast<const T*>(o)) {
    d.lineColor = t->GetLineColor();
    d.lineStyle = t->GetLineStyle();
    d.lineWidth = t->GetLineWidth();
  }
}
template <class T> static void _copyMarker(const TObject* o, LegendDetect& d) {
  if (auto t = dynamic_cast<const T*>(o)) {
    d.markerStyle = t->GetMarkerStyle();
    d.markerSize  = t->GetMarkerSize();
    d.markerColor = t->GetMarkerColor();
  }
}
template <class T> static void _copyFill(const TObject* o, LegendDetect& d) {
  if (auto t = dynamic_cast<const T*>(o)) {
    d.fillColor = t->GetFillColor();
    d.fillStyle = t->GetFillStyle();
  }
}

static LegendDetect DetectLegendStyle(const TObject* obj) {
  LegendDetect d;
  if (!obj) return d;
  _copyLine<TAttLine>(obj, d);
  _copyMarker<TAttMarker>(obj, d);
  _copyFill<TAttFill>(obj, d);

  const bool hasLine   = (d.lineWidth > 0);
  const bool hasMarker = (d.markerStyle != 0 && d.markerSize > 0);
  const bool hasFill   = (d.fillStyle != 0 && d.fillColor != 0);

  if (hasFill && !hasLine && !hasMarker) d.kind = LegendDetect::kFill;
  else if (hasMarker && !hasLine)        d.kind = LegendDetect::kPoint;
  else if (hasLine && !hasMarker)        d.kind = LegendDetect::kLine;
  else if (hasLine && hasMarker)         d.kind = LegendDetect::kLine;
  else                                   d.kind = LegendDetect::kUnknown;

  return d;
}

// --- find first contour segment from ROOT "contours" store (after CONT LIST) ---
static TGraph* _FirstContourSeg() {
  auto* conts = dynamic_cast<TObjArray*>(gROOT->GetListOfSpecials()->FindObject("contours"));
  if (!conts) return nullptr;
  for (int iLev = 0; iLev < conts->GetEntries(); ++iLev) {
    if (auto* level = dynamic_cast<TList*>(conts->At(iLev))) {
      if (!level->IsEmpty()) return dynamic_cast<TGraph*>(level->First());
    }
  }
  return nullptr;
}

// keep legend proxy objects alive
static std::vector<std::unique_ptr<TObject>>& _LegendKeep() {
  static std::vector<std::unique_ptr<TObject>> keep;
  return keep;
}

// build a small sample object (line / point / fill) and apply style
static TObject* _MakeSample(LegendDetect::Kind kind) {
  if (kind == LegendDetect::kFill) {
    auto b = std::make_unique<TBox>(0,0,1,1);
    _LegendKeep().push_back(std::move(b));
    return _LegendKeep().back().get();
  }
  if (kind == LegendDetect::kPoint) {
    auto g = std::make_unique<TGraph>(1);
    g->SetPoint(0, 0., 0.);
    _LegendKeep().push_back(std::move(g));
    return _LegendKeep().back().get();
  }
  auto pl = std::make_unique<TPolyLine>(2);
  pl->SetPoint(0,0.,0.); pl->SetPoint(1,1.,0.);
  _LegendKeep().push_back(std::move(pl));
  return _LegendKeep().back().get();
}

static void _ApplyStyle(TObject* o,
                        int lColor, int lStyle, int lWidth,
                        int mStyle, double mSize, int mColor,
                        int fStyle, int fColor) {
  if (auto L = dynamic_cast<TAttLine*>(o)) {
    if (lColor >= 0) L->SetLineColor(lColor);
    if (lStyle >= 0) L->SetLineStyle(lStyle);
    if (lWidth >= 0) L->SetLineWidth(lWidth);
  }
  if (auto M = dynamic_cast<TAttMarker*>(o)) {
    if (mStyle >= 0) M->SetMarkerStyle(mStyle);
    if (mSize  >= 0) M->SetMarkerSize(mSize);
    if (mColor >= 0) M->SetMarkerColor(mColor);
  }
  if (auto F = dynamic_cast<TAttFill*>(o)) {
    if (fStyle >= 0) F->SetFillStyle(fStyle);
    if (fColor >= 0) F->SetFillColor(fColor);
  }
}

/**
 * AddAnyLegendEntryStyled
 * Lets you *explicitly* choose line style/color/width etc for the legend sample
 * without changing your plotted object.
 *
 * kindHint: "contour" | "line" | "point" | "fill" | nullptr (auto)
 * optOverride: legend text option "l"|"p"|"f"|"lp" or nullptr (auto)
 *
 * For overrides: pass -1 to keep detected defaults. Any >=0 value is applied.
 */
static TLegendEntry* AddAnyLegendEntryStyled(
  TLegend* leg, TObject* obj, const char* label,
  const char* kindHint = nullptr,
  const char* optOverride = nullptr,
  int  legendLineStyle  = -1,
  int  legendLineWidth  = -1,
  int  legendLineColor  = -1,
  int  legendMarkerStyle = -1,
  double legendMarkerSize = -1,
  int  legendMarkerColor = -1,
  int  legendFillStyle  = -1,
  int  legendFillColor  = -1
) {
  if (!leg || !obj) return nullptr;

  // Decide desired glyph ("l"/"p"/"f")
  LegendDetect det = DetectLegendStyle(obj);
  TString hint = kindHint ? kindHint : "";
  hint.ToLower();
  LegendDetect::Kind k =
      hint == "point"   ? LegendDetect::kPoint :
      hint == "fill"    ? LegendDetect::kFill  :
      LegendDetect::kLine; // default/contour/line

  if (hint == "" ) k = det.kind == LegendDetect::kUnknown ? LegendDetect::kLine : det.kind;
  TString opt = optOverride && *optOverride ? optOverride :
                (k == LegendDetect::kPoint ? "p" : k == LegendDetect::kFill ? "f" : "l");

  TObject* sample = nullptr;

  // Special handling: TH2 contours -> clone first contour segment
  if (obj->InheritsFrom(TH2::Class()) && (hint == "contour" || k == LegendDetect::kLine)) {
    if (auto* seg = _FirstContourSeg()) {
      auto* clone = (TGraph*)seg->Clone(Form("legseg_%p", seg));
      _LegendKeep().emplace_back(clone); // own it
      sample = clone;
    }
  }

  // If no TH2 contour seg, try to use the object itself if it matches the glyph
  if (!sample) {
    bool matches =
      (opt == "l" && obj->InheritsFrom(TAttLine::Class())) ||
      (opt == "p" && obj->InheritsFrom(TAttMarker::Class())) ||
      (opt == "f" && obj->InheritsFrom(TAttFill::Class()));
    if (matches && legendLineStyle < 0 && legendLineWidth < 0 &&
        legendLineColor < 0 && legendMarkerStyle < 0 &&
        legendMarkerSize < 0 && legendMarkerColor < 0 &&
        legendFillStyle < 0 && legendFillColor < 0) {
      // no overrides requested → use object directly
      return leg->AddEntry(obj, label, opt);
    }
    // Otherwise build a small proxy sample
    sample = _MakeSample(k);
    // initialize proxy from detected style to look similar
    _ApplyStyle(sample,
                det.lineColor, det.lineStyle, det.lineWidth,
                det.markerStyle, det.markerSize, det.markerColor,
                det.fillStyle, det.fillColor);
  }

  // Apply user overrides to the sample (clone/proxy), not to the real object
  _ApplyStyle(sample,
              legendLineColor, legendLineStyle, legendLineWidth,
              legendMarkerStyle, legendMarkerSize, legendMarkerColor,
              legendFillStyle, legendFillColor);

  return leg->AddEntry(sample, label, opt);
}



struct OscParams {
  double dmsq21 = 7.53e-5;
  double dmsq32_nh = 2.453e-3;   // Normal ordering
  double dmsq32_ih = -2.536e-3;  // Inverted ordering  
  double th12 = 0.5843;
  double th13 = 0.150;
  double th23 = 0.739;           // maximal mixing
  double dcp = 0.0;
};



// Setup oscillation calculator with standard parameters
//osc::IOscCalcAdjustable* ConfigureCalc(const OscParams& p, bool normalOrder = true) {
//  osc::IOscCalcAdjustable* calc = DefaultOscCalc();
//  calc->SetL(kBaseline);
//  calc->SetRho(kDensity);
//  calc->SetDmsq21(p.dmsq21);
//  calc->SetDmsq32(normalOrder ? p.dmsq32_nh : p.dmsq32_ih);
//  calc->SetTh12(p.th12);
//  calc->SetTh13(p.th13);
//  calc->SetTh23(p.th23);
//  calc->SetdCP(p.dcp);
//  return calc;
//}

void checkDefaultOscCalc()
{
  osc::IOscCalculatorAdjustable* calc = DefaultOscCalc();
  
  std::cout << "=== Default Oscillation Calculator Settings ===" << std::endl;
  std::cout << "Baseline L: " << calc->GetL() << " km" << std::endl;
  std::cout << "Density rho: " << calc->GetRho() << " g/cm^3" << std::endl;
  std::cout << std::endl;
  
  std::cout << "=== Oscillation Parameters ===" << std::endl;
  std::cout << "theta12: " << calc->GetTh12() << " rad" << std::endl;
  std::cout << "theta13: " << calc->GetTh13() << " rad" << std::endl;
  std::cout << "theta23: " << calc->GetTh23() << " rad" << std::endl;
  std::cout << "sin^2(theta23): " << TMath::Sin(calc->GetTh23())*TMath::Sin(calc->GetTh23()) << std::endl;
  std::cout << std::endl;
  
  std::cout << "deltaCP: " << calc->GetdCP() << " rad" << std::endl;
  std::cout << "deltaCP/pi: " << calc->GetdCP()/TMath::Pi() << std::endl;
  std::cout << std::endl;
  
  std::cout << "Dm21^2: " << calc->GetDmsq21() << " eV^2" << std::endl;
  std::cout << "Dm32^2: " << calc->GetDmsq32() << " eV^2" << std::endl;
  std::cout << "Dm32^2 (scaled, x10^-3): " << calc->GetDmsq32()*1000 << std::endl;
  std::cout << std::endl;
}

void demo3_v1_1D_v2()
{
  
  TLegend* lg1 = new TLegend( 0.7, 0.7, 0.85, 0.88 );
    lg1->SetNColumns(1);
    lg1->SetBorderSize(0);
    lg1->SetTextSize(.03); //  

 char pdf_title[1024];
 std::string Pdf_name = "DUNEexample_delta_cp_v6";
  TCanvas *c1 = new TCanvas("c1");
   c1->SetMargin(0.12, 0.03, 0.12, 0.06); // (left, right, bottom, top)
 sprintf(pdf_title, "%s.pdf(", Pdf_name.c_str());
 c1 -> Print(pdf_title);
 sprintf(pdf_title, "%s.pdf", Pdf_name.c_str());


  // Repeat most of demo2.C
  //
  // Except we're introducing "Loaders", which knows the latest CAF locations
  // and packages the three types of loader together into one handy class.
  // Didn't have access to the loader included removed and using the Spectrum loader class 
  
  const std::string fnameNonSwap = "/pnfs/dune/persistent/users/cnguyen/cafana_example/FD_FHC_nonswap.root";
  const std::string fnameNueSwap = "/pnfs/dune/persistent/users/cnguyen/cafana_example/FD_FHC_nueswap.root";
  const std::string fnameTauSwap = "/pnfs/dune/persistent/users/cnguyen/cafana_example/FD_FHC_tauswap.root";
  SpectrumLoader loaderNonSwap(fnameNonSwap);
  SpectrumLoader loaderNueSwap(fnameNueSwap);
  SpectrumLoader loaderTauSwap(fnameTauSwap);
  
  
  
  const std::string fnameNonSwap_RHC = "/pnfs/dune/persistent/users/cnguyen/cafana_example/FD_RHC_nonswap.root";
  const std::string fnameNueSwap_RHC = "/pnfs/dune/persistent/users/cnguyen/cafana_example/FD_RHC_nueswap.root";
  const std::string fnameTauSwap_RHC = "/pnfs/dune/persistent/users/cnguyen/cafana_example/FD_RHC_tauswap.root";
  
  SpectrumLoader loaderNonSwap_RHC(fnameNonSwap_RHC);
  SpectrumLoader loaderNueSwap_RHC(fnameNueSwap_RHC);
  SpectrumLoader loaderTauSwap_RHC(fnameTauSwap_RHC);
    
   checkDefaultOscCalc();
    
    
  //TDRLoaders loaders(Loaders::kFHC);
  const Var kRecoEnergy = SIMPLEVAR(dune.Ev_reco_numu);
  const Binning binsEnergy = Binning::Simple(40, 0, 10);
  const HistAxis axEnergy("Reco energy (GeV)", binsEnergy, kRecoEnergy);
  const double pot = 3.5 * 1.47e21 * 40/1.13;
  const Cut kPassesCVN = SIMPLEVAR(dune.cvnnumu) > .5;
  //PredictionNoExtrap pred(loaders, axEnergy, kPassesCVN);
   // Create predictions
  PredictionNoExtrap pred(loaderNonSwap, loaderNueSwap, loaderTauSwap, axEnergy, kPassesCVN);
  loaderNonSwap.Go();
  loaderNueSwap.Go();
  loaderTauSwap.Go();
  
  PredictionNoExtrap predRHC(loaderNonSwap_RHC, loaderNueSwap_RHC, loaderTauSwap_RHC, axEnergy, kPassesCVN);
  loaderNonSwap_RHC.Go();
  loaderNueSwap_RHC.Go();
  loaderTauSwap_RHC.Go();

  // Set POT for staged running (3.5 years FHC + 3.5 years RHC)
  const double pot_staged = 3.5 * 1.47e21 * 40/1.13;
  
  // Setup for CP scan
  const int nPoints = 100; // -1 to 1 in units of pi
  std::vector<double> deltaCPValues;
  std::vector<double> sigmaValues;
  std::vector<double> chisqValues;
  
  // TRUE value of delta_CP for Asimov dataset (e.g., -pi/2)
  const double trueDeltaCP = -TMath::Pi()/2.0;
  
  //const double trueDeltaCP = 0.0;
  // Generate Asimov data at true value
  osc::IOscCalculatorAdjustable* calcTrue = DefaultOscCalc();
  

  
  //OscParams inputos_parameters; 
  
  //calcTrue->SetL(kBaseline);
  //calcTrue->SetRho(kDensity);
  //calcTrue->SetDmsq21(inputos_parameters.dmsq21);
  //calcTrue->SetDmsq32(inputos_parameters.dmsq32_nh);
  //calcTrue->SetTh12(inputos_parameters.th12);
  //calcTrue->SetTh13(inputos_parameters.th13);
  //calcTrue->SetTh23(inputos_parameters.th23);
  //calcTrue->SetdCP(inputos_parameters.dcp);
  
  //calcTrue->SetdCP(trueDeltaCP);
  
  const Spectrum dataFHC = pred.Predict(calcTrue).MockData(pot_staged);
  const Spectrum dataRHC = predRHC.Predict(calcTrue).MockData(pot_staged);
  
  SingleSampleExperiment exptFHC(&pred, dataFHC);
  SingleSampleExperiment exptRHC(&predRHC, dataRHC);
  MultiExperiment exptMulti({&exptFHC, &exptRHC});
  
  
  
    // First, find the best fit (minimum chi-square with all parameters free)
  osc::IOscCalculatorAdjustable* calcBestFit = DefaultOscCalc();
  
  
   calcBestFit->SetTh12(0.5903);    // theta12 in radians
  calcBestFit->SetTh13(0.150);      // theta13 in radians  
  calcBestFit->SetTh23(0.866); // maximal mixing
  
  //calc->SetdCP(-TMath::Pi()/2.0);  // delta_CP = -90 degrees
  
  calcBestFit->SetDmsq21(7.39e-5);   // Dm21^2 in eV^2
  calcBestFit->SetDmsq32(2.451e-3); // Dm32^2 in eV^2 (normal ordering)
  
  
  
  Fitter fitBestFit(&exptMulti, {&kFitSinSqTheta23, &kFitDmSq32Scaled, &kFitDeltaInPiUnits});
  double minChisq = fitBestFit.Fit(calcBestFit);
  
  std::cout << "Best fit delta_CP = " << calcBestFit->GetdCP() << std::endl;
  std::cout << "Minimum chi-square = " << minChisq << std::endl;
  
  
  // Scan over test values of delta_CP
  for(int i = 0; i < nPoints; i++) {
    double testDeltaCP = -TMath::Pi() + 2.0 * TMath::Pi() * i / (nPoints - 1);
    
    // Create test calculator
    osc::IOscCalculatorAdjustable* calcTest = DefaultOscCalc();
    calcTest->SetdCP(testDeltaCP);
    
    // For a proper sensitivity study, you should minimize over other parameters
    // This is a simplified version - ideally use a Fitter to profile over nuisance parameters
    Fitter fit(&exptMulti, {&kFitSinSqTheta23, &kFitDmSq32Scaled});

    
    double chisq = fit.Fit(calcTest);
    
    deltaCPValues.push_back(testDeltaCP / TMath::Pi()); // Convert to units of pi
    chisqValues.push_back(chisq);
    
    if(i % 10 == 0) {
    double test = testDeltaCP / TMath::Pi() ; 
        double sigma = (test > 0) ? sqrt(test) : 0;
      std::cout << "delta_CP/pi = " << test 
                << ", sigma = " << sigma << std::endl;
    }
    
    
  }
  
  /*
  for(int i = 0; i < nPoints; i++) {
    double testDeltaCP = -TMath::Pi() + 2.0 * TMath::Pi() * i / (nPoints - 1);
    
    osc::IOscCalculatorAdjustable* calcTest = DefaultOscCalc();
    calcTest->SetdCP(testDeltaCP);
    
    // Fit with delta_CP FIXED, but profile over theta23 and dm32
    Fitter fitProfile(&exptMulti, {&kFitSinSqTheta23, &kFitDmSq32Scaled});
    double chisq = fitProfile.Fit(calcTest);
    
    double deltaChisq = chisq - minChisq;
    double sigma = (deltaChisq > 0) ? sqrt(deltaChisq) : 0;
    
    deltaCPValues.push_back(testDeltaCP / TMath::Pi());
    sigmaValues.push_back(sigma);
    
    if(i % 10 == 0) {
      std::cout << "delta_CP/pi = " << testDeltaCP/TMath::Pi() 
                << ", sigma = " << sigma << std::endl;
    }
  }
  */
  
  
  // Find minimum chi-square
  double minChisq1 = *std::min_element(chisqValues.begin(), chisqValues.end());
  
  // Create graph of Delta chi-square
  TGraph* grSensitivity = new TGraph(nPoints);
  for(int i = 0; i < nPoints; i++) {
    grSensitivity->SetPoint(i, deltaCPValues[i], sqrt(chisqValues[i] - minChisq1));
  }
  
  // Styling
  grSensitivity->SetLineColor(kBlue);
  grSensitivity->SetLineWidth(3);
  grSensitivity->SetTitle("CP Violation Sensitivity");
  grSensitivity->GetXaxis()->SetTitle("#delta_{CP}/#pi");
  grSensitivity->GetYaxis()->SetTitle("#sigma = #sqrt{#Delta#chi^{2}}");
  grSensitivity->GetXaxis()->SetTitleSize(0.05);
  grSensitivity->GetYaxis()->SetTitleSize(0.05);
  grSensitivity->GetXaxis()->SetLabelSize(0.045);
  grSensitivity->GetYaxis()->SetLabelSize(0.045);
  grSensitivity->GetYaxis()->SetRangeUser(0, 10);
  grSensitivity->GetXaxis()->SetLimits(-1, 1);
  
  grSensitivity->Draw("AL");
  
  // Add significance lines
  TLine* line3sig = new TLine(-1, 3, 1, 3);
  line3sig->SetLineStyle(2);
  line3sig->SetLineWidth(2);
  line3sig->Draw("SAME");
  
  TLine* line5sig = new TLine(-1, 5, 1, 5);
  line5sig->SetLineStyle(4);
  line5sig->SetLineWidth(2);
  line5sig->Draw("SAME");
  
  // Add text labels
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.035);
  latex->SetTextAlign(12);
  latex->DrawLatexNDC(0.15, 0.85, "DUNE Sensitivity");
  //latex->DrawLatexNDC(0.15, 0.80, "All Systematics");
  latex->DrawLatexNDC(0.15, 0.75, "Normal Ordering");
  latex->DrawLatexNDC(0.15, 0.70, "sin^{2}2#theta_{13} = 0.088 #pm 0.003");
  latex->DrawLatexNDC(0.15, 0.65, "0.4 < sin^{2}#theta_{23} < 0.6");
  
  // Legend

  lg1->AddEntry(grSensitivity, "7 years (staged)", "l");
  lg1->AddEntry(line3sig, "3#sigma", "l");
  lg1->AddEntry(line5sig, "5#sigma", "l");
  
  lg1->Draw();
  
  c1->Print(pdf_title);
  

 
  sprintf(pdf_title, "%s.pdf)", Pdf_name.c_str());
  c1 -> Print(pdf_title);
}

