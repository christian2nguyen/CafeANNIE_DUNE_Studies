#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLine.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
  std::string BaseName(const std::string& path)
  {
    const size_t pos = path.find_last_of("/\\");
    if(pos == std::string::npos) return path;
    return path.substr(pos + 1);
  }

  void PrintHelp()
  {
    std::cout << "Usage:\n"
              << "  root -l -q 'lblpwgtools/CAFAna/tute/compare_sens_cpv_nh.C(\"file1.root\",\"Legend 1\",\"file2.root\",\"Legend 2\")'\n"
              << "\n"
              << "Optional arguments:\n"
              << "  2nd argument: legend 1 name (default: basename of file1)\n"
              << "  4th argument: legend 2 name (default: basename of file2)\n"
              << "  5th argument: object name   (default: sens_cpv_nh)\n"
              << "  6th argument: output prefix (default: sens_cpv_nh_compare)\n"
              << "\n"
              << "Example:\n"
              << "  root -l -q 'lblpwgtools/CAFAna/tute/compare_sens_cpv_nh.C(\"first.root\",\"Nominal\",\"second.root\",\"Shifted\",\"sens_cpv_nh\",\"my_compare\")'\n"
              << "\n"
              << "What it does:\n"
              << "  - Opens two ROOT files\n"
              << "  - Reads the object named sens_cpv_nh (or your chosen object)\n"
              << "  - Compares the two objects point-by-point\n"
              << "  - Prints Chi2, Chi2/NDF, mean |difference|, RMS difference, and max |difference|\n"
              << "  - Saves overlay and difference plots as PDF and PNG\n"
              << "\n"
              << "Note:\n"
              << "  In this codebase, sens_cpv_nh is usually written as a TGraph, not a TH1.\n"
              << "  This macro supports both TGraph and TH1 inputs.\n";
  }

  struct PointData
  {
    double x = 0.0;
    double y = 0.0;
    double err = 1.0;
  };

  struct ComparisonSummary
  {
    double chi2 = 0.0;
    int ndf = 0;
    double meanAbsDiff = 0.0;
    double rmsDiff = 0.0;
    double maxAbsDiff = 0.0;
    bool usedUnitErrors = false;
  };

  TObject* GetRequiredObject(TFile& file, const char* objectName)
  {
    TObject* obj = file.Get(objectName);
    if(!obj){
      throw std::runtime_error(std::string("Could not find object '") + objectName +
                               "' in file " + file.GetName());
    }
    return obj;
  }

  double GetGraphYError(const TGraph* graph, int idx, bool& usedUnitErrors)
  {
    if(const auto* asymm = dynamic_cast<const TGraphAsymmErrors*>(graph)){
      const double errLo = asymm->GetErrorYlow(idx);
      const double errHi = asymm->GetErrorYhigh(idx);
      const double err = 0.5*(std::abs(errLo) + std::abs(errHi));
      if(err > 0.0) return err;
    }

    if(const auto* errs = dynamic_cast<const TGraphErrors*>(graph)){
      const double err = std::abs(errs->GetErrorY(idx));
      if(err > 0.0) return err;
    }

    usedUnitErrors = true;
    return 1.0;
  }

  std::vector<PointData> ExtractPoints(const TObject* obj, bool& usedUnitErrors)
  {
    std::vector<PointData> points;

    if(const auto* hist = dynamic_cast<const TH1*>(obj)){
      points.reserve(hist->GetNbinsX());
      for(int bin = 1; bin <= hist->GetNbinsX(); ++bin){
        PointData point;
        point.x = hist->GetXaxis()->GetBinCenter(bin);
        point.y = hist->GetBinContent(bin);
        point.err = hist->GetBinError(bin);
        if(point.err <= 0.0){
          point.err = 1.0;
          usedUnitErrors = true;
        }
        points.push_back(point);
      }
      return points;
    }

    if(const auto* graph = dynamic_cast<const TGraph*>(obj)){
      points.reserve(graph->GetN());
      for(int idx = 0; idx < graph->GetN(); ++idx){
        PointData point;
        graph->GetPoint(idx, point.x, point.y);
        point.err = GetGraphYError(graph, idx, usedUnitErrors);
        points.push_back(point);
      }
      return points;
    }

    throw std::runtime_error(std::string("Object '") + obj->GetName() +
                             "' is neither TH1 nor TGraph.");
  }

  void ValidateComparable(const std::vector<PointData>& a,
                          const std::vector<PointData>& b,
                          double xTolerance = 1e-9)
  {
    if(a.size() != b.size()){
      throw std::runtime_error("Objects do not have the same number of bins/points.");
    }

    for(size_t idx = 0; idx < a.size(); ++idx){
      if(std::abs(a[idx].x - b[idx].x) > xTolerance){
        throw std::runtime_error("Objects do not share the same x coordinates/bin centers.");
      }
    }
  }

  ComparisonSummary ComparePoints(const std::vector<PointData>& a,
                                  const std::vector<PointData>& b)
  {
    ComparisonSummary summary;
    summary.ndf = static_cast<int>(a.size());

    if(a.empty()) return summary;

    double sumAbs = 0.0;
    double sumSqDiff = 0.0;

    for(size_t idx = 0; idx < a.size(); ++idx){
      const double diff = a[idx].y - b[idx].y;
      const double sigma2 = a[idx].err*a[idx].err + b[idx].err*b[idx].err;
      const double safeSigma2 = sigma2 > 0.0 ? sigma2 : 1.0;

      summary.chi2 += diff*diff / safeSigma2;
      sumAbs += std::abs(diff);
      sumSqDiff += diff*diff;
      summary.maxAbsDiff = std::max(summary.maxAbsDiff, std::abs(diff));
    }

    summary.meanAbsDiff = sumAbs / a.size();
    summary.rmsDiff = std::sqrt(sumSqDiff / a.size());
    return summary;
  }

  TGraph* MakeGraphFromPoints(const std::vector<PointData>& points,
                              const char* name,
                              const char* title)
  {
    auto* graph = new TGraph(static_cast<int>(points.size()));
    graph->SetName(name);
    graph->SetTitle(title);
    for(size_t idx = 0; idx < points.size(); ++idx){
      graph->SetPoint(static_cast<int>(idx), points[idx].x, points[idx].y);
    }
    return graph;
  }

  TGraph* MakeDifferenceGraph(const std::vector<PointData>& a,
                              const std::vector<PointData>& b,
                              const char* name,
                              const char* title)
  {
    auto* graph = new TGraph(static_cast<int>(a.size()));
    graph->SetName(name);
    graph->SetTitle(title);
    for(size_t idx = 0; idx < a.size(); ++idx){
      graph->SetPoint(static_cast<int>(idx), a[idx].x, a[idx].y - b[idx].y);
    }
    return graph;
  }

  void StyleGraph(TGraph* graph, int color, int markerStyle)
  {
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(1.0);
    graph->SetLineWidth(2);
  }
}

void compare_sens_cpv_nh()
{
  PrintHelp();
}

void compare_sens_cpv_nh(const char* file1,
                         const char* legendLabel1 = "",
                         const char* file2 = "",
                         const char* legendLabel2 = "",
                         const char* objectName = "sens_cpv_nh",
                         const char* outputPrefix = "sens_cpv_nh_compare")
{
  if(!file1 || !file2 || std::string(file1).empty() || std::string(file2).empty()){
    PrintHelp();
    return;
  }

  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TFile in1(file1, "READ");
  TFile in2(file2, "READ");

  if(in1.IsZombie()){
    throw std::runtime_error(std::string("Could not open file: ") + file1);
  }
  if(in2.IsZombie()){
    throw std::runtime_error(std::string("Could not open file: ") + file2);
  }

  TObject* obj1 = GetRequiredObject(in1, objectName);
  TObject* obj2 = GetRequiredObject(in2, objectName);

  bool usedUnitErrors1 = false;
  bool usedUnitErrors2 = false;
  const std::vector<PointData> points1 = ExtractPoints(obj1, usedUnitErrors1);
  const std::vector<PointData> points2 = ExtractPoints(obj2, usedUnitErrors2);

  ValidateComparable(points1, points2);

  ComparisonSummary summary = ComparePoints(points1, points2);
  summary.usedUnitErrors = usedUnitErrors1 || usedUnitErrors2;

  const std::string label1 = (legendLabel1 && std::string(legendLabel1).size())
                           ? legendLabel1
                           : BaseName(file1);
  const std::string label2 = (legendLabel2 && std::string(legendLabel2).size())
                           ? legendLabel2
                           : BaseName(file2);

  std::unique_ptr<TGraph> graph1(MakeGraphFromPoints(points1, "graph1",
      ";#delta_{CP}/#pi;Sensitivity"));
  std::unique_ptr<TGraph> graph2(MakeGraphFromPoints(points2, "graph2",
      ";#delta_{CP}/#pi;Sensitivity"));
  std::unique_ptr<TGraph> diffGraph(MakeDifferenceGraph(points1, points2, "diffGraph",
      ";#delta_{CP}/#pi;Difference (file1 - file2)"));

  StyleGraph(graph1.get(), kBlue + 1, 20);
  StyleGraph(graph2.get(), kRed + 1, 21);
  StyleGraph(diffGraph.get(), kBlack, 20);

  TCanvas canvas("canvas", "sens_cpv_nh comparison", 900, 900);
  canvas.Divide(1, 2);

  canvas.cd(1);
  gPad->SetBottomMargin(0.12);
  gPad->SetTopMargin(0.08);
  gPad->SetRightMargin(0.06);
  graph1->GetYaxis()->SetRangeUser(0.0, 9.0);
  graph1->Draw("ALP");
  graph1->GetXaxis()->SetTitle("#delta_{CP}/#pi");
  graph1->GetYaxis()->SetTitle("Sensitivity");
  graph2->Draw("LP SAME");
  auto* legend = new TLegend(0.14, 0.76, 0.46, 0.90);
  legend->SetFillColor(0);
  legend->SetBorderSize(1);
  legend->SetTextFont(42);
  legend->SetTextSize(0.032);
  legend->AddEntry(graph1.get(), label1.c_str(), "lp");
  legend->AddEntry(graph2.get(), label2.c_str(), "lp");
  legend->Draw();

  auto* summaryBox = new TPaveText(0.58, 0.66, 0.92, 0.90, "NDC");
  summaryBox->SetFillColor(0);
  summaryBox->SetBorderSize(1);
  summaryBox->SetTextAlign(12);
  summaryBox->SetTextFont(42);
  summaryBox->SetTextSize(0.032);
  summaryBox->AddText(Form("Chi2 = %.6g", summary.chi2));
  if(summary.ndf > 0){
    summaryBox->AddText(Form("Chi2/NDF = %.6g", summary.chi2 / summary.ndf));
  }
  summaryBox->AddText(Form("Mean |diff| = %.6g", summary.meanAbsDiff));
  summaryBox->AddText(Form("RMS diff = %.6g", summary.rmsDiff));
  summaryBox->AddText(Form("Max |diff| = %.6g", summary.maxAbsDiff));
  summaryBox->Draw();

  canvas.cd(2);
  gPad->SetBottomMargin(0.12);
  diffGraph->Draw("ALP");
  diffGraph->GetXaxis()->SetTitle("#delta_{CP}/#pi");
  diffGraph->GetYaxis()->SetTitle("Sensitivity difference");
  TLine zeroLine(points1.front().x, 0.0, points1.back().x, 0.0);
  zeroLine.SetLineStyle(2);
  zeroLine.Draw();

  canvas.SaveAs((std::string(outputPrefix) + ".pdf").c_str());
  canvas.SaveAs((std::string(outputPrefix) + ".png").c_str());

  std::cout << "Compared object: " << objectName << "\n";
  std::cout << "File 1: " << file1 << "\n";
  std::cout << "File 2: " << file2 << "\n";
  std::cout << "Points compared: " << summary.ndf << "\n";
  std::cout << "Chi2 = " << summary.chi2 << "\n";
  if(summary.ndf > 0){
    std::cout << "Chi2/NDF = " << summary.chi2 / summary.ndf << "\n";
  }
  std::cout << "Mean |difference| = " << summary.meanAbsDiff << "\n";
  std::cout << "RMS difference = " << summary.rmsDiff << "\n";
  std::cout << "Max |difference| = " << summary.maxAbsDiff << "\n";
  if(summary.usedUnitErrors){
    std::cout << "Note: one or both inputs had no y-errors, so unit uncertainties were assumed where needed.\n";
  }
  std::cout << "Saved plots to " << outputPrefix << ".pdf and "
            << outputPrefix << ".png\n";
}
