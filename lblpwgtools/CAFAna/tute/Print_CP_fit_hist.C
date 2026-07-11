#include "TCanvas.h"
#include "TCollection.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  void SaveCanvas(TCanvas& can,
                  const std::string& outBase,
                  const std::string& combinedPdf,
                  bool& combinedPdfOpen,
                  bool saveExtraFiles);

  void RenderHist(TH1* hist,
                  const std::string& outBase,
                  const std::string& relPath,
                  const std::string& contextLabel,
                  const std::string& combinedPdf,
                  bool& combinedPdfOpen,
                  bool saveExtraFiles);

  void RenderGraph(TGraph* graph,
                   const std::string& outBase,
                   const std::string& relPath,
                   const std::string& combinedPdf,
                   bool& combinedPdfOpen,
                   bool saveExtraFiles);

  std::string Sanitize(const std::string& name)
  {
    const std::regex badChars("[^A-Za-z0-9._-]+");
    std::string cleaned = std::regex_replace(name, badChars, "_");
    while(!cleaned.empty() && cleaned.front() == '_') cleaned.erase(cleaned.begin());
    while(!cleaned.empty() && cleaned.back() == '_') cleaned.pop_back();
    return cleaned.empty() ? "root" : cleaned;
  }

  std::string JoinPath(const std::string& left, const std::string& right)
  {
    if(left.empty()) return right;
    if(right.empty()) return left;
    if(left.back() == '/') return left + right;
    return left + "/" + right;
  }

  void EnsureDir(const std::string& path)
  {
    if(path.empty()) return;
    gSystem->mkdir(path.c_str(), true);
  }

  bool HasBinLabels(TAxis* axis)
  {
    if(!axis) return false;
    for(int idx = 1; idx <= axis->GetNbins(); ++idx){
      const char* label = axis->GetBinLabel(idx);
      if(label && label[0] != '\0') return true;
    }
    return false;
  }

  bool AxisHasTitle(TH1* hist, char axisName)
  {
    if(!hist) return false;
    TAxis* axis = nullptr;
    if(axisName == 'X') axis = hist->GetXaxis();
    if(axisName == 'Y') axis = hist->GetYaxis();
    if(axisName == 'Z') axis = hist->GetZaxis();
    return axis && axis->GetTitle() && axis->GetTitle()[0] != '\0';
  }

  std::string FormatFloat(double value, double scale = 1.0)
  {
    if(!std::isfinite(value)) return "n/a";
    char buffer[128];
    std::snprintf(buffer, sizeof(buffer), "%.6g", value/scale);
    return buffer;
  }

  std::string Trim(const std::string& text)
  {
    const size_t first = text.find_first_not_of(" \t");
    if(first == std::string::npos) return "";
    const size_t last = text.find_last_not_of(" \t");
    return text.substr(first, last-first+1);
  }

  std::set<int> ParseIndexSpec(const std::string& indexSpec)
  {
    std::set<int> indices;
    if(Trim(indexSpec).empty()) return indices;

    std::stringstream ss(indexSpec);
    std::string token;
    while(std::getline(ss, token, ',')){
      token = Trim(token);
      if(token.empty()) continue;

      const size_t dash = token.find('-');
      if(dash == std::string::npos){
        indices.insert(std::stoi(token));
        continue;
      }

      const int start = std::stoi(Trim(token.substr(0, dash)));
      const int stop = std::stoi(Trim(token.substr(dash+1)));
      const int lo = std::min(start, stop);
      const int hi = std::max(start, stop);
      for(int idx = lo; idx <= hi; ++idx) indices.insert(idx);
    }

    return indices;
  }

  bool ExtractIndexFromPath(const std::string& relPath, int& index)
  {
    static const std::regex indexRegex("(^|[^A-Za-z0-9])index_([0-9]+)([^A-Za-z0-9]|$)");
    std::smatch match;
    if(!std::regex_search(relPath, match, indexRegex)) return false;
    index = std::stoi(match[2].str());
    return true;
  }

  bool ShouldRenderPath(const std::string& relPath, const std::set<int>& allowedIndices)
  {
    if(allowedIndices.empty()) return true;

    int index = -1;
    if(ExtractIndexFromPath(relPath, index)){
      return allowedIndices.count(index);
    }

    // Keep top-level summary objects such as the CPV graph when filtering.
    return relPath.find('/') == std::string::npos;
  }

  bool IsPriorityCovarianceName(const std::string& name)
  {
    static const std::set<std::string> priorityNames = {
      "corr",
      "covar",
      "covar_xsec",
      "covar_xsec_corr",
      "covar_other",
      "covar_other_corr",
      "xsec_corr_strength_matrix_corr"
    };
    return priorityNames.count(name);
  }

  bool StartsWith(const std::string& text, const std::string& prefix)
  {
    return text.rfind(prefix, 0) == 0;
  }

  bool IsEventDistributionPath(const std::string& relPath, const TH1* hist)
  {
    if(!hist) return false;
    if(hist->InheritsFrom(TH2::Class())) return false;

    const std::string name = hist->GetName();
    return StartsWith(name, "data_")
        || StartsWith(name, "prefit_")
        || StartsWith(name, "postfit_");
  }

  bool IsPullTermHist(const TH1* hist)
  {
    if(!hist) return false;
    const std::string name = hist->GetName();
    return name.find("_pull_from_central_prefit") != std::string::npos
        || name.find("_residual_to_fake_prefit") != std::string::npos;
  }

  std::string GetDiagnosticYAxisTitle(const std::string& histName)
  {
    if(histName.find("_pull_from_central_prefit") != std::string::npos){
      return "(postfit - central) / #sigma_{prefit}";
    }
    if(histName.find("_residual_to_fake_prefit") != std::string::npos){
      return "(postfit - fake) / #sigma_{prefit}";
    }
    if(histName.find("_postfit_minus_central") != std::string::npos){
      return "postfit - central";
    }
    if(histName.find("_postfit_minus_fake") != std::string::npos){
      return "postfit - fake";
    }
    if(histName.find("_postfit_error") != std::string::npos){
      return "#sigma_{postfit}";
    }
    if(histName.find("_prefit_error") != std::string::npos){
      return "#sigma_{prefit}";
    }
    if(histName.find("_postfit_value") != std::string::npos){
      return "postfit value";
    }
    if(histName.find("_central_value") != std::string::npos){
      return "central value";
    }
    if(histName.find("_fake_value") != std::string::npos){
      return "fake value";
    }
    return "Value";
  }

  bool ReadFitContext(TDirectory* dir, double& trueDcp)
  {
    if(!dir) return false;
    TTree* fitContext = static_cast<TTree*>(dir->Get("fit_context"));
    if(!fitContext || fitContext->GetEntries() < 1) return false;
    fitContext->SetBranchAddress("true_dcp", &trueDcp);
    fitContext->GetEntry(0);
    return true;
  }

  std::string BuildContextLabel(TDirectory* dir, const std::string& relPath)
  {
    int index = -1;
    const bool hasIndex = ExtractIndexFromPath(relPath, index);

    double trueDcp = 0.0;
    const bool hasDcp = ReadFitContext(dir, trueDcp);

    if(!hasIndex && !hasDcp) return "";

    std::string label;
    if(hasIndex) label += "scan index = " + std::to_string(index);
    if(hasDcp){
      if(!label.empty()) label += ", ";
      label += "true #delta_{CP}/#pi = " + FormatFloat(trueDcp, TMath::Pi());
    }
    return label;
  }

  struct ScanPullPoint
  {
    int index = -1;
    double trueDcp = 0.0;
    std::vector<double> pulls;
    std::vector<std::string> labels;
  };

  bool LoadScanPullPoint(TDirectory* dir, ScanPullPoint& point)
  {
    if(!dir) return false;

    TTree* fitContext = static_cast<TTree*>(dir->Get("fit_context"));
    TH1* pullHist = static_cast<TH1*>(dir->Get("xsec_pull_from_central_prefit"));
    if(!fitContext || !pullHist) return false;

    double true_dcp = 0.0;
    fitContext->SetBranchAddress("true_dcp", &true_dcp);
    if(fitContext->GetEntries() < 1) return false;
    fitContext->GetEntry(0);

    point.trueDcp = true_dcp;
    point.pulls.clear();
    point.labels.clear();
    for(int ibin = 1; ibin <= pullHist->GetNbinsX(); ++ibin){
      point.pulls.push_back(pullHist->GetBinContent(ibin));
      point.labels.push_back(pullHist->GetXaxis()->GetBinLabel(ibin));
    }
    return true;
  }

  std::vector<ScanPullPoint> CollectScanPullPoints(TDirectory* topDir,
                                                   const std::set<int>& allowedIndices)
  {
    std::vector<ScanPullPoint> points;
    if(!topDir) return points;

    static const std::regex scanDirRegex("^scan_true_dcp_index_([0-9]+)_bestfit$");

    TIter next(topDir->GetListOfKeys());
    while(TKey* key = static_cast<TKey*>(next())){
      if(std::string(key->GetClassName()) != "TDirectoryFile") continue;

      const std::string dirName = key->GetName();
      std::smatch match;
      if(!std::regex_match(dirName, match, scanDirRegex)) continue;

      const int index = std::stoi(match[1].str());
      if(!allowedIndices.empty() && !allowedIndices.count(index)) continue;

      std::unique_ptr<TObject> obj(key->ReadObj());
      TDirectory* dir = dynamic_cast<TDirectory*>(obj.get());
      if(!dir) continue;

      ScanPullPoint point;
      point.index = index;
      if(LoadScanPullPoint(dir, point)) points.push_back(point);
    }

    std::sort(points.begin(), points.end(),
              [](const ScanPullPoint& a, const ScanPullPoint& b){ return a.index < b.index; });
    return points;
  }

  void RenderReadableScanHeatmap(TH2D* heatmap,
                                 const std::string& outBase,
                                 const std::string& relPath,
                                 const std::string& contextLabel,
                                 const std::string& combinedPdf,
                                 bool& combinedPdfOpen,
                                 bool saveExtraFiles)
  {
    if(!heatmap) return;

    const int nY = heatmap->GetNbinsY();
    const int width = 3200;
    const int height = std::max(1500, std::min(3000, 950 + 20*nY));
    TCanvas can(("can_" + Sanitize(relPath)).c_str(), "", width, height);
    can.SetLeftMargin(0.36);
    can.SetRightMargin(0.10);
    can.SetBottomMargin(0.12);
    can.SetTopMargin(0.06);

    heatmap->SetContour(100);
    heatmap->GetXaxis()->SetLabelSize(0.020);
    heatmap->GetXaxis()->SetTitleSize(0.028);
    heatmap->GetXaxis()->SetTitleOffset(1.4);
    heatmap->GetYaxis()->SetLabelSize(std::max(0.013, std::min(0.020, 0.85/std::max(1, nY))));
    heatmap->GetYaxis()->SetTitleSize(0.028);
    heatmap->GetYaxis()->SetTitleOffset(2.8);
    heatmap->GetZaxis()->SetLabelSize(0.020);
    heatmap->GetZaxis()->SetTitleSize(0.028);
    heatmap->GetZaxis()->SetTitleOffset(1.1);
    heatmap->GetZaxis()->SetTitle("pull");
    heatmap->SetTitle("Xsec pull vs true #delta_{CP} index (rows sorted by max |pull|)");
    heatmap->Draw("COLZ");

    SaveCanvas(can, outBase, combinedPdf, combinedPdfOpen, saveExtraFiles);
  }

  void RenderXSecPullScanPlots(TDirectory* topDir,
                               const std::string& outputDir,
                               const std::string& combinedPdf,
                               bool& combinedPdfOpen,
                               bool saveExtraFiles,
                               std::vector<std::string>& saved)
  {
    const std::set<int> noIndexFilter;
    const std::vector<ScanPullPoint> points = CollectScanPullPoints(topDir, noIndexFilter);
    if(points.empty()) return;

    const int nX = points.size();
    const int nY = points.front().pulls.size();
    if(nY == 0) return;

    std::vector<double> rowMaxAbs(nY, 0.0);
    for(int iy = 0; iy < nY; ++iy){
      for(int ix = 0; ix < nX; ++ix){
        rowMaxAbs[iy] = std::max(rowMaxAbs[iy], std::abs(points[ix].pulls[iy]));
      }
    }

    std::vector<int> rowOrder(nY);
    for(int iy = 0; iy < nY; ++iy) rowOrder[iy] = iy;
    std::sort(rowOrder.begin(), rowOrder.end(),
              [&](int a, int b){
                if(rowMaxAbs[a] != rowMaxAbs[b]) return rowMaxAbs[a] > rowMaxAbs[b];
                return points.front().labels[a] < points.front().labels[b];
              });

    TH2D heatmap("xsec_pull_scan_heatmap",
                 "Xsec pull vs true #delta_{CP} index;scan index;Xsec parameter",
                 nX, 0, nX, nY, 0, nY);

    TGraph maxAbsGraph;
    maxAbsGraph.SetNameTitle("xsec_pull_maxabs_vs_dcp",
                             "Max |xsec pull| vs true #delta_{CP};true #delta_{CP}/#pi;max |pull|");

    TGraph rmsGraph;
    rmsGraph.SetNameTitle("xsec_pull_rms_vs_dcp",
                          "RMS xsec pull vs true #delta_{CP};true #delta_{CP}/#pi;RMS pull");

    for(int rank = 0; rank < nY; ++rank){
      const int src = rowOrder[rank];
      const int ybin = nY - rank;
      heatmap.GetYaxis()->SetBinLabel(ybin, points.front().labels[src].c_str());
    }

    double globalMaxAbs = 0.0;

    for(int ix = 0; ix < nX; ++ix){
      const ScanPullPoint& point = points[ix];
      heatmap.GetXaxis()->SetBinLabel(ix+1, std::to_string(point.index).c_str());

      double maxAbs = 0.0;
      double sumSq = 0.0;
      for(int rank = 0; rank < nY; ++rank){
        const int src = rowOrder[rank];
        const int ybin = nY - rank;
        const double pull = point.pulls[src];
        heatmap.SetBinContent(ix+1, ybin, pull);
        maxAbs = std::max(maxAbs, std::abs(pull));
        sumSq += pull*pull;
      }
      globalMaxAbs = std::max(globalMaxAbs, maxAbs);

      const double dcpOverPi = point.trueDcp / TMath::Pi();
      const double rms = std::sqrt(sumSq / std::max(1, nY));
      maxAbsGraph.SetPoint(maxAbsGraph.GetN(), dcpOverPi, maxAbs);
      rmsGraph.SetPoint(rmsGraph.GetN(), dcpOverPi, rms);
    }

    if(globalMaxAbs > 0.0){
      heatmap.SetMinimum(-globalMaxAbs);
      heatmap.SetMaximum(globalMaxAbs);
    }

    RenderReadableScanHeatmap(&heatmap,
                              JoinPath(outputDir, "xsec_pull_scan_heatmap"),
                              "summary/xsec_pull_scan_heatmap",
                              "",
                              combinedPdf,
                              combinedPdfOpen,
                              saveExtraFiles);
    saved.push_back("summary/xsec_pull_scan_heatmap");

    RenderGraph(&maxAbsGraph,
                JoinPath(outputDir, "xsec_pull_maxabs_vs_dcp"),
                "summary/xsec_pull_maxabs_vs_dcp",
                combinedPdf,
                combinedPdfOpen,
                saveExtraFiles);
    saved.push_back("summary/xsec_pull_maxabs_vs_dcp");

    RenderGraph(&rmsGraph,
                JoinPath(outputDir, "xsec_pull_rms_vs_dcp"),
                "summary/xsec_pull_rms_vs_dcp",
                combinedPdf,
                combinedPdfOpen,
                saveExtraFiles);
    saved.push_back("summary/xsec_pull_rms_vs_dcp");
  }

  void SaveCanvas(TCanvas& can,
                  const std::string& outBase,
                  const std::string& combinedPdf,
                  bool& combinedPdfOpen,
                  bool saveExtraFiles)
  {
    can.Modified();
    can.Update();
    if(saveExtraFiles){
      can.SaveAs((outBase + ".png").c_str());
      can.SaveAs((outBase + ".pdf").c_str());
    }
    if(!combinedPdf.empty()){
      if(!combinedPdfOpen){
        can.Print((combinedPdf + "(").c_str(), "pdf");
        combinedPdfOpen = true;
      }
      else{
        can.Print(combinedPdf.c_str(), "pdf");
      }
    }
  }

  void RenderHist(TH1* hist,
                  const std::string& outBase,
                  const std::string& relPath,
                  const std::string& contextLabel,
                  const std::string& combinedPdf,
                  bool& combinedPdfOpen,
                  bool saveExtraFiles)
  {
    if(!hist) return;

    if(IsEventDistributionPath(relPath, hist)){
      hist->Scale(1.0, "width");
    }

    const bool labeledX = HasBinLabels(hist->GetXaxis());
    const int width = (labeledX || hist->GetNbinsX() > 20) ? 1800 : 1100;
    const int height = labeledX ? 900 : 800;
    TCanvas can(("can_" + Sanitize(relPath)).c_str(), "", width, height);

    if(hist->InheritsFrom(TH2::Class())){
      can.SetRightMargin(0.16);
      can.SetLeftMargin(0.13);
      can.SetBottomMargin(0.12);
      can.SetTopMargin(0.08);
      hist->SetContour(80);
      hist->Draw("COLZ");
    }
    else{
      can.SetLeftMargin(0.12);
      can.SetRightMargin(0.04);
      can.SetTopMargin(0.08);
      can.SetBottomMargin(labeledX ? 0.28 : 0.12);
      hist->SetLineWidth(2);
      hist->SetMarkerStyle(20);
      hist->SetMarkerSize(0.9);
      hist->Draw("E1");
      if(labeledX) hist->LabelsOption("v", "X");
    }

    if(std::string(hist->GetTitle()).empty()) hist->SetTitle(relPath.c_str());
    if(!contextLabel.empty() && IsPullTermHist(hist)){
      std::string title = hist->GetTitle();
      if(!title.empty()) title += " | ";
      title += contextLabel;
      hist->SetTitle(title.c_str());
    }
    if(IsEventDistributionPath(relPath, hist)){
      hist->GetXaxis()->SetTitle("E_{#nu} [GeV]");
      hist->GetYaxis()->SetTitle("Events / GeV");
    }
    else{
      if(!AxisHasTitle(hist, 'X')) hist->GetXaxis()->SetTitle("Bin");
      hist->GetYaxis()->SetTitle(GetDiagnosticYAxisTitle(hist->GetName()).c_str());
    }

    SaveCanvas(can, outBase, combinedPdf, combinedPdfOpen, saveExtraFiles);
  }

  void RenderGraph(TGraph* graph,
                   const std::string& outBase,
                   const std::string& relPath,
                   const std::string& combinedPdf,
                   bool& combinedPdfOpen,
                   bool saveExtraFiles)
  {
    if(!graph) return;

    TCanvas can(("can_" + Sanitize(relPath)).c_str(), "", 1100, 800);
    can.SetLeftMargin(0.12);
    can.SetRightMargin(0.04);
    can.SetTopMargin(0.08);
    can.SetBottomMargin(0.12);

    graph->SetLineWidth(2);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.0);

    std::string drawOpt = "ALP";
    if(graph->InheritsFrom(TGraphErrors::Class()) ||
       graph->InheritsFrom(TGraphAsymmErrors::Class())){
      drawOpt = "AP";
    }

    graph->Draw(drawOpt.c_str());

    TH1* frame = graph->GetHistogram();
    if(frame){
      if(std::string(graph->GetTitle()).empty()){
        graph->SetTitle(relPath.c_str());
        frame->SetTitle(relPath.c_str());
      }
      if(!AxisHasTitle(frame, 'X')) frame->GetXaxis()->SetTitle("x");
      if(!AxisHasTitle(frame, 'Y')) frame->GetYaxis()->SetTitle("y");
    }

    SaveCanvas(can, outBase, combinedPdf, combinedPdfOpen, saveExtraFiles);
  }

  void RenderTextCard(const std::vector<std::string>& lines,
                      const std::string& outBase,
                      const std::string& title,
                      const std::string& combinedPdf,
                      bool& combinedPdfOpen,
                      bool saveExtraFiles)
  {
    TCanvas can(("can_" + Sanitize(title + "_" + outBase)).c_str(), "", 1000, 700);
    can.SetLeftMargin(0.04);
    can.SetRightMargin(0.04);
    can.SetTopMargin(0.04);
    can.SetBottomMargin(0.04);

    TPaveText pave(0.03, 0.03, 0.97, 0.97, "NDC");
    pave.SetFillColor(0);
    pave.SetBorderSize(1);
    pave.SetTextAlign(12);
    pave.SetTextFont(42);
    pave.SetTextSize(0.032);
    pave.AddText(title.c_str());
    pave.AddText("");
    for(const std::string& line: lines) pave.AddText(line.c_str());
    pave.Draw();

    SaveCanvas(can, outBase, combinedPdf, combinedPdfOpen, saveExtraFiles);
  }

  bool RenderFitContextTree(TTree* tree,
                            const std::string& outBase,
                            const std::string& relPath,
                            const std::string& combinedPdf,
                            bool& combinedPdfOpen,
                            bool saveExtraFiles)
  {
    int true_hie = 0;
    int test_hie = 0;
    double true_dcp = 0;
    double test_dcp = 0;
    double chisq = 0;

    tree->SetBranchAddress("true_hie", &true_hie);
    tree->SetBranchAddress("test_hie", &test_hie);
    tree->SetBranchAddress("true_dcp", &true_dcp);
    tree->SetBranchAddress("test_dcp", &test_dcp);
    tree->SetBranchAddress("chisq", &chisq);
    if(tree->GetEntries() > 0) tree->GetEntry(0);

    RenderTextCard({
        "path: " + relPath,
        "true_hie: " + std::to_string(true_hie),
        "test_hie: " + std::to_string(test_hie),
        "true_dcp/pi: " + FormatFloat(true_dcp, TMath::Pi()),
        "test_dcp/pi: " + FormatFloat(test_dcp, TMath::Pi()),
        "chisq: " + FormatFloat(chisq)
      }, outBase, "fit_context", combinedPdf, combinedPdfOpen, saveExtraFiles);
    return true;
  }

  bool RenderFitInfoTree(TTree* tree,
                         const std::string& outBase,
                         const std::string& relPath,
                         const std::string& combinedPdf,
                         bool& combinedPdfOpen,
                         bool saveExtraFiles)
  {
    double chisq = 0;
    double resMem = 0;
    double virtMem = 0;
    double nfcn = 0;
    double edm = 0;
    unsigned nSeconds = 0;
    unsigned nOscSeeds = 0;
    bool isValid = false;

    tree->SetBranchAddress("chisq", &chisq);
    tree->SetBranchAddress("ResMemUsage", &resMem);
    tree->SetBranchAddress("VirtMemUsage", &virtMem);
    tree->SetBranchAddress("NFCN", &nfcn);
    tree->SetBranchAddress("EDM", &edm);
    tree->SetBranchAddress("NSeconds", &nSeconds);
    tree->SetBranchAddress("NOscSeeds", &nOscSeeds);
    tree->SetBranchAddress("IsValid", &isValid);
    if(tree->GetEntries() > 0) tree->GetEntry(0);

    RenderTextCard({
        "path: " + relPath,
        "chisq: " + FormatFloat(chisq),
        "NFCN: " + FormatFloat(nfcn),
        "EDM: " + FormatFloat(edm),
        std::string("IsValid: ") + (isValid ? "1" : "0"),
        "NSeconds: " + std::to_string(nSeconds),
        "NOscSeeds: " + std::to_string(nOscSeeds),
        "ResMemUsage (kB): " + FormatFloat(resMem),
        "VirtMemUsage (kB): " + FormatFloat(virtMem)
      }, outBase, "fit_info", combinedPdf, combinedPdfOpen, saveExtraFiles);
    return true;
  }

  bool RenderMetaTree(TTree* tree,
                      const std::string& outBase,
                      const std::string& relPath,
                      const std::string& combinedPdf,
                      bool& combinedPdfOpen,
                      bool saveExtraFiles)
  {
    std::vector<std::string>* paramNames = nullptr;
    std::vector<std::string>* envNames = nullptr;
    tree->SetBranchAddress("fParamNames", &paramNames);
    tree->SetBranchAddress("fEnvVarNames", &envNames);
    if(tree->GetEntries() > 0) tree->GetEntry(0);

    std::vector<std::string> lines;
    lines.push_back("path: " + relPath);
    lines.push_back("n_parameters: " + std::to_string(paramNames ? paramNames->size() : 0));
    lines.push_back("n_env_vars: " + std::to_string(envNames ? envNames->size() : 0));

    if(paramNames && !paramNames->empty()){
      std::string preview = "param preview: ";
      const size_t nPreview = std::min<size_t>(10, paramNames->size());
      for(size_t i = 0; i < nPreview; ++i){
        if(i) preview += ", ";
        preview += paramNames->at(i);
      }
      if(paramNames->size() > nPreview) preview += ", ...";
      lines.push_back(preview);
    }

    RenderTextCard(lines, outBase, "meta_tree", combinedPdf, combinedPdfOpen, saveExtraFiles);
    return true;
  }

  bool RenderFitParameterSummaryTree(TTree* tree,
                                     const std::string& outBase,
                                     const std::string& relPath,
                                     const std::string& combinedPdf,
                                     bool& combinedPdfOpen,
                                     bool saveExtraFiles)
  {
    int is_xsec = 0;
    int is_osc = 0;
    int is_other_nuisance = 0;
    tree->SetBranchAddress("is_xsec", &is_xsec);
    tree->SetBranchAddress("is_osc", &is_osc);
    tree->SetBranchAddress("is_other_nuisance", &is_other_nuisance);

    int nXSec = 0;
    int nOsc = 0;
    int nOther = 0;
    for(Long64_t i = 0; i < tree->GetEntries(); ++i){
      tree->GetEntry(i);
      if(is_xsec) ++nXSec;
      else if(is_osc) ++nOsc;
      else if(is_other_nuisance) ++nOther;
      else ++nOther;
    }

    RenderTextCard({
        "path: " + relPath,
        "entries: " + std::to_string(tree->GetEntries()),
        "xsec parameters: " + std::to_string(nXSec),
        "osc parameters: " + std::to_string(nOsc),
        "other nuisance parameters: " + std::to_string(nOther)
      }, outBase, "fit_parameter_summary", combinedPdf, combinedPdfOpen, saveExtraFiles);
    return true;
  }

  bool RenderTree(TTree* tree,
                  const std::string& outBase,
                  const std::string& relPath,
                  const std::string& combinedPdf,
                  bool& combinedPdfOpen,
                  bool saveExtraFiles)
  {
    if(!tree) return false;
    const std::string name = tree->GetName();
    if(name == "fit_context") return RenderFitContextTree(tree, outBase, relPath, combinedPdf, combinedPdfOpen, saveExtraFiles);
    if(name == "fit_info") return RenderFitInfoTree(tree, outBase, relPath, combinedPdf, combinedPdfOpen, saveExtraFiles);
    if(name == "meta_tree") return RenderMetaTree(tree, outBase, relPath, combinedPdf, combinedPdfOpen, saveExtraFiles);
    if(name == "fit_parameter_summary") return RenderFitParameterSummaryTree(tree, outBase, relPath, combinedPdf, combinedPdfOpen, saveExtraFiles);
    return false;
  }

  void WalkDirectory(TDirectory* dir,
                     const std::string& outDir,
                     const std::string& relDir,
                     const std::string& combinedPdf,
                     bool& combinedPdfOpen,
                     const std::set<int>& allowedIndices,
                     bool saveExtraFiles,
                     std::vector<std::string>& saved,
                     std::vector<std::string>& skipped)
  {
    if(!dir) return;
    EnsureDir(outDir);

    std::vector<TKey*> priorityKeys;
    std::vector<TKey*> otherKeys;
    TIter next(dir->GetListOfKeys());
    while(TKey* key = static_cast<TKey*>(next())){
      if(IsPriorityCovarianceName(key->GetName())) priorityKeys.push_back(key);
      else otherKeys.push_back(key);
    }

    std::vector<TKey*> orderedKeys;
    orderedKeys.insert(orderedKeys.end(), priorityKeys.begin(), priorityKeys.end());
    orderedKeys.insert(orderedKeys.end(), otherKeys.begin(), otherKeys.end());

    for(TKey* key: orderedKeys){
      TObject* obj = key->ReadObj();
      if(!obj) continue;

      const std::string name = key->GetName();
      const std::string relPath = relDir.empty() ? name : relDir + "/" + name;
      const std::string outBase = JoinPath(outDir, Sanitize(name));
      const std::string contextLabel = BuildContextLabel(dir, relPath);

      if(obj->InheritsFrom(TDirectory::Class())){
        WalkDirectory(static_cast<TDirectory*>(obj),
                      JoinPath(outDir, Sanitize(name)),
                      relPath, combinedPdf, combinedPdfOpen, allowedIndices, saveExtraFiles, saved, skipped);
        delete obj;
        continue;
      }

      if(!ShouldRenderPath(relPath, allowedIndices)){
        delete obj;
        continue;
      }

      if(obj->InheritsFrom(TH1::Class())){
        RenderHist(static_cast<TH1*>(obj), outBase, relPath, contextLabel, combinedPdf, combinedPdfOpen, saveExtraFiles);
        saved.push_back(relPath);
        delete obj;
        continue;
      }

      if(obj->InheritsFrom(TGraph::Class())){
        RenderGraph(static_cast<TGraph*>(obj), outBase, relPath, combinedPdf, combinedPdfOpen, saveExtraFiles);
        saved.push_back(relPath);
        delete obj;
        continue;
      }

      if(obj->InheritsFrom(TTree::Class())){
        if(RenderTree(static_cast<TTree*>(obj), outBase, relPath, combinedPdf, combinedPdfOpen, saveExtraFiles)) saved.push_back(relPath);
        else skipped.push_back(relPath + " [TTree]");
        delete obj;
        continue;
      }

      skipped.push_back(relPath + " [" + std::string(obj->ClassName()) + "]");
      delete obj;
    }
  }
}

void Print_CP_fit_hist(std::string inputRoot = "",
                       std::string outputDir = "",
                       std::string indexSpec = "",
                       std::string outputPdf = "",
                       bool saveExtraFiles = false)
{
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.045, "XY");
  gStyle->SetLabelSize(0.035, "XY");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  if(inputRoot.empty()){
    std::cout << "Usage:" << std::endl;
    std::cout << "root -l -b -q 'lblpwgtools/CAFAna/tute/Print_CP_fit_hist.C(\"/path/to/output.root\", \"\", \"0,1,4-7\", \"my_plots.pdf\", true)'" << std::endl;
    return;
  }

  if(gSystem->AccessPathName(inputRoot.c_str())){
    std::cerr << "Input file not found: " << inputRoot << std::endl;
    return;
  }

  if(outputDir.empty()){
    std::string base = gSystem->BaseName(inputRoot.c_str());
    const size_t dot = base.rfind('.');
    if(dot != std::string::npos) base = base.substr(0, dot);
    outputDir = JoinPath(gSystem->DirName(inputRoot.c_str()), base + "_plots");
  }
  EnsureDir(outputDir);
  std::string combinedPdf = outputPdf;
  if(combinedPdf.empty()){
    combinedPdf = JoinPath(gSystem->WorkingDirectory(), "all_plots.pdf");
  }
  else if(combinedPdf.find('/') == std::string::npos){
    combinedPdf = JoinPath(gSystem->WorkingDirectory(), combinedPdf);
  }
  const std::set<int> allowedIndices = ParseIndexSpec(indexSpec);

  TFile* fin = TFile::Open(inputRoot.c_str(), "READ");
  if(!fin || fin->IsZombie()){
    std::cerr << "Could not open ROOT file: " << inputRoot << std::endl;
    if(fin) fin->Close();
    delete fin;
    return;
  }

  std::vector<std::string> saved;
  std::vector<std::string> skipped;
  bool combinedPdfOpen = false;
  WalkDirectory(fin, outputDir, "", combinedPdf, combinedPdfOpen, allowedIndices, saveExtraFiles, saved, skipped);
  RenderXSecPullScanPlots(fin, outputDir, combinedPdf, combinedPdfOpen, saveExtraFiles, saved);
  fin->Close();
  delete fin;

  if(combinedPdfOpen){
    TCanvas closeCan("closeCan", "", 10, 10);
    closeCan.Print((combinedPdf + ")").c_str(), "pdf");
  }

  const std::string manifest = JoinPath(outputDir, "plot_manifest.txt");
  std::ofstream out(manifest.c_str());
  out << "Input ROOT file: " << inputRoot << "\n";
  out << "Index filter: " << (indexSpec.empty() ? "all" : indexSpec) << "\n";
  out << "Save extra files: " << (saveExtraFiles ? "true" : "false") << "\n";
  out << "Plots written to: " << outputDir << "\n\n";
  out << "Combined PDF: " << combinedPdf << "\n\n";
  out << "Saved objects:\n";
  for(const std::string& item: saved) out << "  " << item << "\n";
  out << "\nSkipped objects:\n";
  for(const std::string& item: skipped) out << "  " << item << "\n";
  out.close();

  std::cout << "Saved " << saved.size() << " plotted objects to " << outputDir << std::endl;
  std::cout << "Combined PDF written to " << combinedPdf << std::endl;
  std::cout << "Manifest written to " << manifest << std::endl;
  if(!skipped.empty()) std::cout << "Skipped " << skipped.size() << " unsupported objects" << std::endl;
}
