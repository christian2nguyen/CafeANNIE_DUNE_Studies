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



void demo3_v1()
{
  
  TLegend* lg1 = new TLegend( 0.15, 0.7, 0.8, 0.88 );
    lg1->SetNColumns(2);
    lg1->SetBorderSize(0);
    lg1->SetTextSize(.035); //  

 char pdf_title[1024];
 std::string Pdf_name = "DUNEexample_demo3_v5_kFitDmSq32Scaled_kFitDeltaInPiUnits_v1";
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

  Surface surf(&expt, calc,
               &kFitDeltaInPiUnits, 40, -1, 1,
               &kFitDmSq32Scaled, 100, 2.4, 2.5 // I think scaled is in units of 10^-5
               );
//&kFitSinSq2Theta13, 80, 0.035,0.205
  //Surface surf(&expt, calc,
  //             &kFitSinSqTheta23, 40, 0.455, 0.575,
  //             &kFitDeltaInPiUnits, 40, -1, 1);

  //Surface surf(&expt, calc,
  //             &kFitDeltaInPiUnits, 100, -1, 1,
  //             &kFitSinSqTheta23, 100, 0.455, 0.6);


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




  //surf.Draw();
  surf.DrawBestFit(kBlue);

  // In a full Feldman-Cousins analysis you need to provide a critical value
  // surface to be able to draw a contour. But we provide these helper
  // functions to use the gaussian up-values.
  TH2* crit1sig = Gaussian68Percent2D(surf);
  TH2* crit2sig = Gaussian2Sigma2D(surf);

  surf.DrawContour(crit1sig, 7, kBlue);
  surf.DrawContour(crit2sig, kSolid, kBlue);
  lg1->Draw("same");
 c1 -> Print(pdf_title);

 TH2D* crit1sig_clone =  (TH2D*)crit1sig->Clone("crit1sig_clone");
 TH2D* crit2sig_clone =  (TH2D*)crit2sig->Clone("crit2sig_clone");

AddAnyLegendEntryStyled(lg1, crit1sig, "(FHC) Mock Data(Best fit)", "point", "p",
                        /*lineStyle*/-1, /*lineWidth*/-1, /*lineColor*/-1,
                        /*markerStyle*/20, /*markerSize*/1.2, /*markerColor*/kBlue);
                        
AddAnyLegendEntryStyled(lg1, crit2sig, "(RHC) Mock Data(Best fit)", "point", "p",
/*lineStyle*/-1, /*lineWidth*/-1, /*lineColor*/-1,
/*markerStyle*/20, /*markerSize*/1.2, /*markerColor*/kRed);
                        
AddAnyLegendEntryStyled(lg1, crit2sig_clone, "(FHC) 2#pm#sigma", "contour", "l", /*lineStyle*/1, /*lineWidth*/2, /*lineColor*/kBlue);
AddAnyLegendEntryStyled(lg1, crit1sig_clone, "(FHC) 68% C.L.", "contour", "l", /*lineStyle*/2, /*lineWidth*/2, /*lineColor*/kBlue);


  // Let's try to add in the effect of 3.5yrs of RHC data too
  //TDRLoaders loadersRHC(Loaders::kRHC);
  //PredictionNoExtrap predRHC(loadersRHC, axEnergy, kPassesCVN);
  //loadersRHC.Go();
  
    PredictionNoExtrap predRHC(loaderNonSwap_RHC, loaderNueSwap_RHC, loaderTauSwap_RHC, axEnergy, kPassesCVN);
  loaderNonSwap_RHC.Go();
  loaderNueSwap_RHC.Go();
  loaderTauSwap_RHC.Go();
  
  
  calc = DefaultOscCalc(); // Remember to reset, since fits modified it
  const Spectrum dataRHC = predRHC.Predict(calc).MockData(pot);
  SingleSampleExperiment exptRHC(&predRHC, dataRHC);

  // A MultiExperiment gets its chisq just by adding together its component
  // parts. Use to implement joint fits
  MultiExperiment exptMulti({&expt, &exptRHC});

  //Surface surfMulti(&exptMulti, calc,
  //                  &kFitSinSqTheta23, 100, 0.455, 0.575,
  //                  &kFitDmSq32Scaled, 100, 2.4, 2.5);
  
  
  Surface surfMulti(&exptMulti, calc,
                &kFitDeltaInPiUnits, 40, -1, 1,
              &kFitDmSq32Scaled, 100, 2.4, 2.5 // I think scaled is in units of 10^-5
               );
  
  

  surfMulti.DrawBestFit(kRed);
  surfMulti.DrawContour(crit1sig, 7, kRed);
  surfMulti.DrawContour(crit2sig, kSolid, kRed);
 
  
TH2D* crit1sig_cloneII =  (TH2D*)crit1sig->Clone("crit1sig_cloneIIII");
TH2D* crit2sig_cloneII =  (TH2D*)crit2sig->Clone("crit2sig_cloneIIII");

AddAnyLegendEntryStyled(lg1, crit2sig_cloneII, "(RHC) 2#pm#sigma", "contour", "l", /*lineStyle*/1, /*lineWidth*/2, /*lineColor*/kRed);
AddAnyLegendEntryStyled(lg1, crit1sig_cloneII, "(RHC) 68% C.L.", "contour", "l", /*lineStyle*/2, /*lineWidth*/2, /*lineColor*/kRed);
lg1->Draw("same");
  
   c1 -> Print(pdf_title);
  
  sprintf(pdf_title, "%s.pdf)", Pdf_name.c_str());
  c1 -> Print(pdf_title);
}

