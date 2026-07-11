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



void demo3_v1_1D()
{
  
  TLegend* lg1 = new TLegend( 0.15, 0.7, 0.8, 0.88 );
    lg1->SetNColumns(2);
    lg1->SetBorderSize(0);
    lg1->SetTextSize(.035); //  

 char pdf_title[1024];
 std::string Pdf_name = "DUNEexample_delta_cp";
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
  
  //TDRLoaders loaders(Loaders::kFHC);
  const Var kRecoEnergy = SIMPLEVAR(dune.Ev_reco_numu);
  const Binning binsEnergy = Binning::Simple(40, 0, 10);
  const HistAxis axEnergy("Reco energy (GeV)", binsEnergy, kRecoEnergy);
  const double pot = 3.5 * 1.47e21 * 40/1.13;
  const Cut kPassesCVN = SIMPLEVAR(dune.cvnnumu) > .5;
  //PredictionNoExtrap pred(loaders, axEnergy, kPassesCVN);
  
  PredictionNoExtrap pred(loaderNonSwap, loaderNueSwap, loaderTauSwap, axEnergy, kPassesCVN);
  loaderNonSwap.Go();
  loaderNueSwap.Go();
  loaderTauSwap.Go();
  
  //loaders.Go();
  osc::IOscCalculatorAdjustable* calc = DefaultOscCalc();
  const Spectrum data = pred.Predict(calc).MockData(pot);
  SingleSampleExperiment expt(&pred, data);
  
  // A Surface evaluates the experiment's chisq across a grid
  //Surface surf(&expt, calc,
  //             &kFitSinSqTheta23, 100, 0.455, 0.575,
  //             &kFitDmSq32Scaled, 100, 2.4, 2.5);

 double pi = TMath::Pi();



//kFitTheta13
//kFitSinSq2Theta13
//kFitDeltaInPiUnits 
//kFitTheta23
//kFitSinSqTheta23
//kFitSinSqTheta23LowerOctant
//kFitSinSqTheta23UpperOctant
//kFitSinSqTheta23BelowSymmetry
//kFitSinSqTheta23AboveSymmetry
//kFitSinSq2Theta23
//kFitDmSq32
//kFitDmSq32Scaled
//kFitDmSq32NHScaled
//kFitDmSq32IHScaled
//kFitTanSqTheta12
//kFitSinSq2Theta12
//kFitDmSq21
//kFitDmSq21Scaled
//kFitRho

 const int nPoints = 100;
  TGraph* grFHC = new TGraph(nPoints);
  TGraph* grCombined = new TGraph(nPoints);
  
  // Scan over delta_CP values
  for(int i = 0; i < nPoints; i++) {
    double deltaCP = -1.0 + 2.0 * i / (nPoints - 1); // -1 to +1 (in units of pi)
    
    calc = DefaultOscCalc(); // Reset calculator
    calc->SetdCP(deltaCP * TMath::Pi()); // Set delta_CP
    
    // Calculate chi-square for FHC only
    double chisqFHC = expt.ChiSq(calc);
    grFHC->SetPoint(i, deltaCP, chisqFHC);
  }
  
  
  
    PredictionNoExtrap predRHC(loaderNonSwap_RHC, loaderNueSwap_RHC, loaderTauSwap_RHC, axEnergy, kPassesCVN);
  loaderNonSwap_RHC.Go();
  loaderNueSwap_RHC.Go();
  loaderTauSwap_RHC.Go();
  
  
  calc = DefaultOscCalc(); // Remember to reset, since fits modified it
  const Spectrum dataRHC = predRHC.Predict(calc).MockData(pot);
  SingleSampleExperiment exptRHC(&predRHC, dataRHC);
  MultiExperiment exptMulti({&expt, &exptRHC});
  
  
  
    for(int i = 0; i < nPoints; i++) {
    double deltaCP = -1.0 + 2.0 * i / (nPoints - 1);
    
    calc = DefaultOscCalc();
    calc->SetdCP(deltaCP * TMath::Pi());
    
    double chisqCombined = exptMulti.ChiSq(calc);
    grCombined->SetPoint(i, deltaCP, chisqCombined);
  }
  
  // Find minimum chi-square to calculate Delta chi-square
  double minChisqFHC = TMath::MinElement(grFHC->GetN(), grFHC->GetY());
  double minChisqCombined = TMath::MinElement(grCombined->GetN(), grCombined->GetY());
  
    // Scan for combined FHC+RHC
  for(int i = 0; i < nPoints; i++) {
    double deltaCP = -1.0 + 2.0 * i / (nPoints - 1);
    
    calc = DefaultOscCalc();
    calc->SetdCP(deltaCP * TMath::Pi());
    
    double chisqCombined = exptMulti.ChiSq(calc);
    grCombined->SetPoint(i, deltaCP, chisqCombined);
  }
  
  // Find minimum chi-square to calculate Delta chi-square

  
  // Convert to Delta chi-square
  for(int i = 0; i < nPoints; i++) {
    grFHC->GetY()[i] -= minChisqFHC;
    grCombined->GetY()[i] -= minChisqCombined;
  }
  
  // Plot
  grFHC->SetLineColor(kBlue);
  grFHC->SetLineWidth(2);
  grFHC->GetXaxis()->SetTitle("#delta_{CP} / #pi");
  grFHC->GetYaxis()->SetTitle("#Delta#chi^{2}");
  grFHC->GetYaxis()->SetRangeUser(0, 20);
  grFHC->Draw("AL");
  
  grCombined->SetLineColor(kRed);
  grCombined->SetLineWidth(2);
  grCombined->Draw("L SAME");
  
  // Add significance lines
  TLine* line1sig = new TLine(-1, 1, 1, 1);
  line1sig->SetLineStyle(2);
  line1sig->Draw("SAME");
  
  TLine* line3sig = new TLine(-1, 9, 1, 9);
  line3sig->SetLineStyle(2);
  line3sig->Draw("SAME");
  
  TLine* line5sig = new TLine(-1, 25, 1, 25);
  line5sig->SetLineStyle(2);
  line5sig->Draw("SAME");
  
  // Update legend

  lg1->AddEntry(grFHC, "FHC only", "l");
  lg1->AddEntry(grCombined, "FHC + RHC", "l");
  lg1->Draw();
  
  c1->Print(pdf_title);
  

  
  sprintf(pdf_title, "%s.pdf)", Pdf_name.c_str());
  c1 -> Print(pdf_title);
}

