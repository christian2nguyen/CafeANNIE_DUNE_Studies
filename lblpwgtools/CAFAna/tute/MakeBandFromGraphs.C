#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>
#include <cmath>

#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include <TGraphErrors.h>
#include <stdexcept>
#include <cmath>
#include "TKey.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TH1F.h"

/**
 * Struct to hold the band and median graph together
 */
struct BandResult {
  std::unique_ptr<TGraphAsymmErrors> band;
  std::unique_ptr<TGraph> median;
  int nThrows;
  
  bool IsValid() const { return band && median; }
};

/**
 * Configuration for quantile calculation
 */
struct QuantileConfig {
  double lower = 0.16;   // 1-sigma lower (~16%)
  double median = 0.50;  // median
  double upper = 0.84;   // 1-sigma upper (~84%)
  
  // For 2-sigma bands
  static QuantileConfig TwoSigma() {
    QuantileConfig config;
    config.lower = 0.025;  // 2.5%
    config.upper = 0.975;  // 97.5%
    return config;
  }
};

/**
 * Helper function to get all ROOT files from a directory
 */
std::vector<std::string> GetRootFiles(const char* dirpath) {
  std::vector<std::string> rootfiles;
  
  TSystemDirectory dir("dir", dirpath);
  TList *files = dir.GetListOfFiles();
  
  if (!files) {
    std::cerr << "Error: Directory not found: " << dirpath << "\n";
    return rootfiles;
  }

  TIter next(files);
  while (TSystemFile *f = (TSystemFile*)next()) {
    if (f->IsDirectory()) continue;
    
    TString fname = f->GetName();
    if (fname.EndsWith(".root")) {
      rootfiles.emplace_back(std::string(dirpath) + "/" + fname.Data());
    }
  }

  return rootfiles;
}

/**
 * Helper to read a TGraph from a file (assumes first key is the graph)
 */
TGraph* ReadFirstGraph(const char* filename) {
  TFile f(filename, "READ");
  if (f.IsZombie()) {
    std::cerr << "Error: Cannot open " << filename << "\n";
    return nullptr;
  }

  TKey *key = (TKey*)f.GetListOfKeys()->At(0);
  if (!key) {
    std::cerr << "Error: No keys in " << filename << "\n";
    return nullptr;
  }

  TGraph *g = dynamic_cast<TGraph*>(key->ReadObj());
  if (!g) {
    std::cerr << "Error: First object in " << filename << " is not a TGraph\n";
    return nullptr;
  }

  // Detach from file so it survives f.Close()
//  g->SetDirectory(nullptr);
  return g;
}

/**
 * Validate that x-grids match (with tolerance)
 */
bool ValidateXGrid(const std::vector<double>& x_ref, TGraph* g, 
                   double tolerance = 1e-6) {
  if (g->GetN() != (int)x_ref.size()) {
    return false;
  }
  
  for (int i = 0; i < g->GetN(); ++i) {
    double xi, yi;
    g->GetPoint(i, xi, yi);
    if (std::abs(xi - x_ref[i]) > tolerance) {
      std::cerr << "Warning: X-grid mismatch at point " << i 
                << " (expected " << x_ref[i] << ", got " << xi << ")\n";
      return false;
    }
  }
  return true;
}

/**
 * Calculate quantile from sorted vector
 */
double GetQuantile(const std::vector<double>& sorted_data, double quantile) {
  if (sorted_data.empty()) return 0.0;
  
  double index = quantile * (sorted_data.size() - 1);
  int lower_idx = (int)std::floor(index);
  int upper_idx = (int)std::ceil(index);
  
  // Linear interpolation between indices
  if (lower_idx == upper_idx) {
    return sorted_data[lower_idx];
  }
  
  double weight = index - lower_idx;
  return sorted_data[lower_idx] * (1.0 - weight) + 
         sorted_data[upper_idx] * weight;
}

/**
 * Main function: Creates a band from multiple TGraphs
 * 
 * dirpath : directory containing .root files with TGraphs
 * config  : quantile configuration (default 1-sigma band)
 * color   : base color for the band/line
 * label   : legend label for the median line
 *
 * Returns: BandResult containing band and median graph
 */
 /*
BandResult MakeBandFromGraphs(const char* dirpath, 
                               const QuantileConfig& config = QuantileConfig(),
                               Color_t color = kAzure+1,
                               const char* label = "Median of throws")
{
  BandResult result;

  // 1. Get all ROOT files
  auto rootfiles = GetRootFiles(dirpath);
  if (rootfiles.empty()) {
    std::cerr << "Error: No .root files found in " << dirpath << "\n";
    return result;
  }
  
  std::cout << "Found " << rootfiles.size() << " root files in " << dirpath << "\n";

  // 2. Read first graph to establish x-grid
  TGraph *g0 = ReadFirstGraph(rootfiles[0].c_str());
  if (!g0) return result;

  int nPoints = g0->GetN();
  std::vector<double> x(nPoints);
  for (int i = 0; i < nPoints; ++i) {
    double y_dummy;
    g0->GetPoint(i, x[i], y_dummy);
  }
  delete g0;

  // 3. Collect y-values from all files
  std::vector<std::vector<double>> yvals(nPoints);
  int nValidFiles = 0;

  for (const auto &fname : rootfiles) {
    TGraph *g = ReadFirstGraph(fname.c_str());
    if (!g) continue;
    
    if (!ValidateXGrid(x, g)) {
      std::cerr << "Warning: Skipping " << fname 
                << " due to x-grid mismatch\n";
      delete g;
      continue;
    }

    for (int i = 0; i < nPoints; ++i) {
      double xi, yi;
      g->GetPoint(i, xi, yi);
      yvals[i].push_back(yi);
    }
    
    delete g;
    ++nValidFiles;
  }

  if (nValidFiles == 0) {
    std::cerr << "Error: No valid graphs were read!\n";
    return result;
  }
  
  std::cout << "Successfully read " << nValidFiles << " graphs\n";

  // 4. Calculate quantiles for each point
  std::vector<double> y_med(nPoints), y_low(nPoints), y_high(nPoints);
  
  for (int i = 0; i < nPoints; ++i) {
    auto &vec = yvals[i];
    std::sort(vec.begin(), vec.end());

    y_med[i]  = GetQuantile(vec, config.median);
    y_low[i]  = GetQuantile(vec, config.lower);
    y_high[i] = GetQuantile(vec, config.upper);
  }

  // 5. Create output graphs
  result.median = std::make_unique<TGraph>(nPoints, x.data(), y_med.data());
  result.median->SetLineColor(color + 2);
  result.median->SetLineWidth(3);
  result.median->SetTitle(label);

  result.band = std::make_unique<TGraphAsymmErrors>(nPoints);
  for (int i = 0; i < nPoints; ++i) {
    double yerr_low  = y_med[i] - y_low[i];
    double yerr_high = y_high[i] - y_med[i];
    result.band->SetPoint(i, x[i], y_med[i]);
    result.band->SetPointError(i, 0.0, 0.0, yerr_low, yerr_high);
   
  }

  result.band->SetFillColorAlpha(color, 0.4);
  result.band->SetLineColor(color);
  result.band->SetLineWidth(1);
  result.nThrows = nValidFiles;

  return result;
}
*/

BandResult MakeBandFromGraphs(const char* dirpath, 
                               const QuantileConfig& config = QuantileConfig(),
                               Color_t color = kAzure+1,
                               const char* label = "Mean of throws")
{
  BandResult result;

  // 1. Get all ROOT files
  auto rootfiles = GetRootFiles(dirpath);
  if (rootfiles.empty()) {
    std::cerr << "Error: No .root files found in " << dirpath << "\n";
    return result;
  }
  
  std::cout << "Found " << rootfiles.size() << " root files in " << dirpath << "\n";

  // 2. Read first graph to establish x-grid
  TGraph *g0 = ReadFirstGraph(rootfiles[0].c_str());
  if (!g0) return result;

  int nPoints = g0->GetN();
  std::vector<double> x(nPoints);
  for (int i = 0; i < nPoints; ++i) {
    double y_dummy;
    g0->GetPoint(i, x[i], y_dummy);
  }
  delete g0;

  // 3. Collect y-values from all files
  std::vector<std::vector<double>> yvals(nPoints);
  int nValidFiles = 0;

  for (const auto &fname : rootfiles) {
    TGraph *g = ReadFirstGraph(fname.c_str());
    if (!g) continue;
    
    if (!ValidateXGrid(x, g)) {
      std::cerr << "Warning: Skipping " << fname 
                << " due to x-grid mismatch\n";
      delete g;
      continue;
    }

    for (int i = 0; i < nPoints; ++i) {
      double xi, yi;
      g->GetPoint(i, xi, yi);
      yvals[i].push_back(yi);
    }
    
    delete g;
    ++nValidFiles;
  }

  if (nValidFiles == 0) {
    std::cerr << "Error: No valid graphs were read!\n";
    return result;
  }
  
  std::cout << "Successfully read " << nValidFiles << " graphs\n";

  // 4. Calculate mean and quantiles for each point
  std::vector<double> y_mean(nPoints), y_low(nPoints), y_high(nPoints);
  
  for (int i = 0; i < nPoints; ++i) {
    auto &vec = yvals[i];
    
    // Calculate mean (average)
    double sum = 0.0;
    for (double val : vec) {
      sum += val;
    }
    y_mean[i] = sum / vec.size();
    
    // Sort for quantiles
    std::sort(vec.begin(), vec.end());
    y_low[i]  = GetQuantile(vec, config.lower);
    y_high[i] = GetQuantile(vec, config.upper);
  }

  // 5. Create output graphs
  result.median = std::make_unique<TGraph>(nPoints, x.data(), y_mean.data());
  result.median->SetLineColor(color + 2);
  result.median->SetLineWidth(1);
  result.median->SetTitle(label);

  result.band = std::make_unique<TGraphAsymmErrors>(nPoints);
  for (int i = 0; i < nPoints; ++i) {
    double yerr_low  = y_mean[i] - y_low[i];
    double yerr_high = y_high[i] - y_mean[i];
    result.band->SetPoint(i, x[i], y_mean[i]);
    result.band->SetPointError(i, 0.0, 0.0, yerr_low, yerr_high);
  }

  result.band->SetFillColorAlpha(color, 0.4);
  result.band->SetLineColor(color);
  result.band->SetLineWidth(1);
  result.nThrows = nValidFiles;

  return result;
}
/**
 * Convenience overload without QuantileConfig
 */
BandResult MakeBandFromGraphs(const char* dirpath, Color_t color = kAzure+1,
                               const char* label = "Median of throws")
{
  return MakeBandFromGraphs(dirpath, QuantileConfig(), color, label);
}

BandResult MakeBandFromGraphs_new(const char* dirpath, 
                               const char* median_file = nullptr,  // NEW: path to specific median file
                               Color_t color = kAzure+1,
                               const char* label = "Median of throws",const QuantileConfig& config = QuantileConfig())
{
  BandResult result;

  // 1. Get all ROOT files
  auto rootfiles = GetRootFiles(dirpath);
  if (rootfiles.empty()) {
    std::cerr << "Error: No .root files found in " << dirpath << "\n";
    return result;
  }
  
  std::cout << "Found " << rootfiles.size() << " root files in " << dirpath << "\n";

  // 2. Read first graph to establish x-grid
  TGraph *g0 = ReadFirstGraph(rootfiles[0].c_str());
  if (!g0) return result;

  int nPoints = g0->GetN();
  std::vector<double> x(nPoints);
  for (int i = 0; i < nPoints; ++i) {
    double y_dummy;
    g0->GetPoint(i, x[i], y_dummy);
  }
  delete g0;

  // 3. Collect y-values from all files
  std::vector<std::vector<double>> yvals(nPoints);
  int nValidFiles = 0;

  for (const auto &fname : rootfiles) {
    TGraph *g = ReadFirstGraph(fname.c_str());
    if (!g) continue;
    
    if (!ValidateXGrid(x, g)) {
      std::cerr << "Warning: Skipping " << fname 
                << " due to x-grid mismatch\n";
      delete g;
      continue;
    }

    for (int i = 0; i < nPoints; ++i) {
      double xi, yi;
      g->GetPoint(i, xi, yi);
      yvals[i].push_back(yi);
    }
    
    delete g;
    ++nValidFiles;
  }

  if (nValidFiles == 0) {
    std::cerr << "Error: No valid graphs were read!\n";
    return result;
  }
  
  std::cout << "Successfully read " << nValidFiles << " graphs\n";

  // 4. Calculate quantiles for bands OR read median from file
  std::vector<double> y_med(nPoints), y_low(nPoints), y_high(nPoints);
  
  if (median_file != nullptr) {
    // Use specified file for median
    std::cout << "Using specified median file: " << median_file << "\n";
    TGraph *g_median = ReadFirstGraph(median_file);
    if (!g_median) {
      std::cerr << "Error: Failed to read median file " << median_file << "\n";
      return result;
    }
    
    if (!ValidateXGrid(x, g_median)) {
      std::cerr << "Error: Median file x-grid doesn't match!\n";
      delete g_median;
      return result;
    }
    
    // Extract median y-values from the file
    for (int i = 0; i < nPoints; ++i) {
      double xi, yi;
      g_median->GetPoint(i, xi, yi);
      y_med[i] = yi;
    }
    delete g_median;
    
    // Still calculate quantiles for the band (using throws)
    for (int i = 0; i < nPoints; ++i) {
      auto &vec = yvals[i];
      std::sort(vec.begin(), vec.end());
      y_low[i]  = GetQuantile(vec, config.lower);
      y_high[i] = GetQuantile(vec, config.upper);
    }
  } else {
    // Original behavior: use median from throws
    for (int i = 0; i < nPoints; ++i) {
      auto &vec = yvals[i];
      std::sort(vec.begin(), vec.end());
      y_med[i]  = GetQuantile(vec, config.median);
      y_low[i]  = GetQuantile(vec, config.lower);
      y_high[i] = GetQuantile(vec, config.upper);
    }
  }

  // 5. Create output graphs
  result.median = std::make_unique<TGraph>(nPoints, x.data(), y_med.data());
  result.median->SetLineColor(color + 2);
  result.median->SetLineWidth(3);
  result.median->SetTitle(label);

  result.band = std::make_unique<TGraphAsymmErrors>(nPoints);
  for (int i = 0; i < nPoints; ++i) {
    double yerr_low  = y_med[i] - y_low[i];
    double yerr_high = y_high[i] - y_med[i];
    result.band->SetPoint(i, x[i], y_med[i]);
    result.band->SetPointError(i, 0.0, 0.0, yerr_low, yerr_high);
  }

  result.band->SetFillColorAlpha(color, 0.4);
  result.band->SetLineColor(color);
  result.band->SetLineWidth(1);
  result.nThrows = nValidFiles;

  return result;
}



/**
 * Draw all individual throws from a directory
 */
void DrawAllThrows(const char* dirpath, Color_t color, 
                   const char* title, TCanvas* c, bool draw_axes = true) {
  auto rootfiles = GetRootFiles(dirpath);
  if (rootfiles.empty()) {
    std::cerr << "No files in " << dirpath << "\n";
    return;
  }

  std::cout << "Drawing " << rootfiles.size() << " throws from " << dirpath << "\n";

  bool first = true;
  for (const auto &fname : rootfiles) {
    TGraph *g = ReadFirstGraph(fname.c_str());
    if (!g) continue;

    g->SetLineColor(color);
    g->SetLineWidth(1);
    g->SetLineStyle(1);
    g->SetMarkerSize(0);
    
    if (first && draw_axes) {
      g->SetTitle(title);
      g->GetXaxis()->SetTitle("#delta_{CP}/#pi");
      g->GetYaxis()->SetTitle("#sqrt{#Delta#chi^{2}}");
      g->GetXaxis()->SetLimits(-1.0, 1.0);
      g->SetMinimum(0.0);
      g->SetMaximum(8.5);
      g->Draw("AL");
      first = false;
    } else {
      g->Draw("L same");
    }
  }
}

 std::vector<TGraph*> Get_figures(const char* dirpath){

 std::vector<TGraph*> output_tgraphs;

  auto rootfiles = GetRootFiles(dirpath);
  if (rootfiles.empty()) {
    std::cerr << "No files in " << dirpath << "\n";
    return output_tgraphs;
  }
  
  for (const auto &fname : rootfiles) {
    TGraph *g = ReadFirstGraph(fname.c_str());
    output_tgraphs.push_back(g);
  }


return output_tgraphs;

}

#include <string>
#include <algorithm>

// Helper function to extract the key name from filename
std::string ExtractKeyName(const std::string& filename) {
  // Expected format: StateFeb4_OUTPUT_<KEY_NAME>_624ktmwyr_...
  
  // Find "OUTPUT_"
  size_t start = filename.find("OUTPUT_");
  if (start == std::string::npos) {
    return filename; // Return full name if pattern not found
  }
  start += 7; // Move past "OUTPUT_"
  
  // Find the next underscore after OUTPUT_
  size_t end = filename.find("_", start);
  if (end == std::string::npos) {
    return filename.substr(start); // Return rest of string if no underscore found
  }
  
  return filename.substr(start, end - start);
}

void DrawAllThrows(const char* dirpath, Color_t base_color, 
                   const char* title, TCanvas* c, TLegend* legend,
                   bool draw_axes = true) {
  auto rootfiles = GetRootFiles(dirpath);
  if (rootfiles.empty()) {
    std::cerr << "No files in " << dirpath << "\n";
    return;
  }

  std::cout << "Drawing " << rootfiles.size() << " throws from " << dirpath << "\n";

  // Define color palette and line styles
  std::vector<Color_t> colors = {kRed, kBlue, kGreen+2, kMagenta, kCyan+1, 
                                  kOrange+7, kViolet, kTeal-5, kSpring+3, kAzure+7,
                                  kPink-3, kYellow+2, kRed+2, kBlue+2, kGreen-2};
  
  std::vector<Int_t> line_styles = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};  // Solid, dashed, dotted, etc.

  bool first = true;
  int file_index = 0;
  
  for (const auto &fname : rootfiles) {
    TGraph *g = ReadFirstGraph(fname.c_str());
    if (!g) continue;

    // Cycle through colors and line styles
    Color_t color = colors[file_index % colors.size()];
    Int_t line_style = line_styles[file_index % line_styles.size()];
    
    g->SetLineColor(color);
    g->SetLineWidth(1);  // Make lines a bit thicker for visibility
    g->SetLineStyle(line_style);
    g->SetMarkerSize(0);
    
    if (first && draw_axes) {
      g->SetTitle(title);
      g->GetXaxis()->SetTitle("#delta_{CP}/#pi");
      g->GetYaxis()->SetTitle("#sqrt{#Delta#chi^{2}}");
      g->GetXaxis()->SetLimits(-1.0, 1.0);
      g->SetMinimum(0.0);
      g->SetMaximum(10.5);
      g->Draw("AL");
      first = false;
    } else {
      g->Draw("L same");
    }
    
    // Add entry to legend for each graph
    if (legend) {
      // Extract key name from filename
      std::string basename = fname.substr(fname.find_last_of("/\\") + 1);
      std::string keyname = ExtractKeyName(basename);
      
      // Create a unique label with file index
      char label[256];
      sprintf(label, "%s (#%d)", keyname.c_str(), file_index + 1);
      legend->AddEntry(g, label, "l");
    }
    
    file_index++;
  }
}

std::vector<TGraph*> Get_figures(const char* dirpath, TLegend* legend) {

  std::vector<TGraph*> output_tgraphs;

  auto rootfiles = GetRootFiles(dirpath);
  if (rootfiles.empty()) {
    std::cerr << "No files in " << dirpath << "\n";
    return output_tgraphs;
  }
  
  bool first = true;
  for (const auto &fname : rootfiles) {
    TGraph *g = ReadFirstGraph(fname.c_str());
    if (g) {
      output_tgraphs.push_back(g);
      
      // Add to legend if provided (only from first file)
      if (legend && first) {
        // Extract key name from filename
        std::string basename = fname.substr(fname.find_last_of("/\\") + 1);
        std::string keyname = ExtractKeyName(basename);
        legend->AddEntry(g, keyname.c_str(), "l");
        first = false;
      }
    }
  }

  return output_tgraphs;
}

TGraph* AverageGraphs(const std::vector<TGraph*>& graphs) {
    if (graphs.empty()) {
        throw std::runtime_error("AverageGraphs: no graphs provided.");
    }

    // All graphs must have the same number of points
    int nPoints = graphs[0]->GetN();
    for (size_t i = 1; i < graphs.size(); ++i) {
        if (graphs[i]->GetN() != nPoints) {
            throw std::runtime_error("AverageGraphs: graphs have different numbers of points.");
        }
    }

    // Prepare arrays for averages
    std::vector<double> xAvg(nPoints, 0.0);
    std::vector<double> yAvg(nPoints, 0.0);

    // Loop over graphs and accumulate
    for (TGraph* g : graphs) {
        const double* x = g->GetX();
        const double* y = g->GetY();

        for (int i = 0; i < nPoints; ++i) {
            xAvg[i] += x[i];
            yAvg[i] += y[i];
        }
    }

    // Divide by number of graphs
    int N = graphs.size();
    for (int i = 0; i < nPoints; ++i) {
        xAvg[i] /= N;
        yAvg[i] /= N;
    }

    // Create a new averaged TGraph
    TGraph* gAvg = new TGraph(nPoints, xAvg.data(), yAvg.data());
    gAvg->SetName("gAverage");
    gAvg->SetTitle("Average of graphs");

    return gAvg;
}



TGraphErrors* MakeAverageGraphWithUncertainty(const std::vector<TGraph*>& graphs)
{
    if (graphs.empty()) {
        throw std::runtime_error("MakeAverageGraphWithUncertainty: no graphs provided.");
    }

    // All graphs must have the same number of points
    int nPoints = graphs[0]->GetN();
    for (size_t ig = 1; ig < graphs.size(); ++ig) {
        if (!graphs[ig]) {
            throw std::runtime_error("MakeAverageGraphWithUncertainty: null TGraph* in input.");
        }
        if (graphs[ig]->GetN() != nPoints) {
            throw std::runtime_error("MakeAverageGraphWithUncertainty: graphs have different numbers of points.");
        }
    }

    const int nGraphs = graphs.size();

    std::vector<double> x(nPoints, 0.0);
    std::vector<double> yMean(nPoints, 0.0);
    std::vector<double> ySigma(nPoints, 0.0); // uncertainty (std dev)

    // Use X from the first graph
    const double* x0 = graphs[0]->GetX();
    for (int i = 0; i < nPoints; ++i) {
        x[i] = x0[i];
    }

    // --- First pass: compute mean Y at each point ---
    for (TGraph* g : graphs) {
        const double* y = g->GetY();
        for (int i = 0; i < nPoints; ++i) {
            yMean[i] += y[i];
        }
    }
    for (int i = 0; i < nPoints; ++i) {
        yMean[i] /= static_cast<double>(nGraphs);
    }

    // --- Second pass: compute standard deviation at each point ---
    for (TGraph* g : graphs) {
        const double* y = g->GetY();
        for (int i = 0; i < nPoints; ++i) {
            double diff = y[i] - yMean[i];
            ySigma[i] += diff * diff;
        }
    }

    // Use sample standard deviation: sqrt( sum (y - mean)^2 / (N - 1) )
    // If you prefer population std dev, divide by N instead of (N-1).
    for (int i = 0; i < nPoints; ++i) {
        if (nGraphs > 1) {
            ySigma[i] = std::sqrt(ySigma[i] / static_cast<double>(nGraphs - 1));
        } else {
            ySigma[i] = 0.0; // only one graph -> no spread
        }

        // If you wanted the error on the mean instead of the spread:
        // ySigma[i] /= std::sqrt(static_cast<double>(nGraphs));
    }

    // X-errors are 0 here
    std::vector<double> ex(nPoints, 0.0);

    // Build TGraphErrors
    TGraphErrors* gAvgErr = new TGraphErrors(
        nPoints,
        x.data(),
        yMean.data(),
        ex.data(),
        ySigma.data()
    );

    gAvgErr->SetName("gAverageWithUnc");
    gAvgErr->SetTitle("Average graph with spread as uncertainty");

    return gAvgErr;
}




/**
 * Draw CP violation sensitivity bands with customization
 * Also creates plots showing all individual throws
 */
void DrawCPBandExample(const char* dir7, const char* dir10,
                       const char* output_pdf = "cp_sensitivity.pdf")
{
  auto band7  = MakeBandFromGraphs(dir7,  kAzure+1, "7 years (median)");
  auto band10 = MakeBandFromGraphs(dir10, kOrange+7, "10 years (median)");

  if (!band7.IsValid() || !band10.IsValid()) {
    std::cerr << "Error: Failed to create bands\n";
    return;
  }

  std::cout << "7-year band uses " << band7.nThrows << " throws\n";
  std::cout << "10-year band uses " << band10.nThrows << " throws\n";

  // ===== PLOT 1: Band summary =====
  TCanvas *c1 = new TCanvas("c_cp_bands", "CP sensitivity bands", 800, 600);
  c1->SetMargin(0.12, 0.04, 0.12, 0.06);

  // Set up axes
  band10.median->SetTitle(";#delta_{CP}/#pi;#sqrt{#Delta#chi^{2}}");

 
  band10.median->GetXaxis()->SetTitleSize(0.045);
  band10.median->GetYaxis()->SetTitleSize(0.045);
  
  band7.median->SetLineColor(kAzure+2);
  band7.median->SetLineWidth(3);
  //band7.median->Draw("L same");
  
  band10.median->SetLineColor(kOrange+7);
  band10.median->SetLineWidth(3);
  //band10.median->Draw("L same");
  // Draw bands and medians
  
  band10.band->SetFillColor(kOrange+7);
  band10.band->SetFillStyle(3005);
  band10.median->GetXaxis()->SetLimits(-1.0, 1.0);
  band10.median->GetYaxis()->SetLimits(0.0, 10.5);
  band10.median->SetMinimum(0.0);
  band10.median->SetMaximum(9.5);

  band10.median->Draw("AL");
  band10.median->GetXaxis()->SetLimits(-1.0, 1.0);
  band10.median->GetYaxis()->SetLimits(0.0, 10.5);
     c1->Update();
  band10.band->Draw("a4 same");
  band10.band->GetXaxis()->SetLimits(-1.0, 1.0);
  band10.band->GetYaxis()->SetLimits(0.0, 10.5);
    
  band10.band->SetMinimum(0.0);
  band10.band->SetMaximum(10.5);

  band10.median->Draw("L same");
  
  band7.band->SetFillColor(kAzure+2);
  band7.band->SetFillStyle(3005);
  band7.band->Draw("a4 same");
  

  
  band7.median->Draw("L same");

  
  // Reference lines for 3σ and 5σ
  TLine *l3 = new TLine(-1.0, 3.0, 1.0, 3.0);
  l3->SetLineStyle(2);
  l3->SetLineColor(kGray+2);
  l3->Draw("same");
  
  TLine *l5 = new TLine(-1.0, 5.0, 1.0, 5.0);
  l5->SetLineStyle(2);
  l5->SetLineColor(kGray+2);
  l5->Draw("same");

  // Legend
  auto leg1 = new TLegend(0.55, 0.65, 0.88, 0.88);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(band7.median.get(),  "7 years (staged)",  "l");
  leg1->AddEntry(band7.band.get(),    "7 yr 1#sigma band", "f");
  leg1->AddEntry(band10.median.get(), "10 years (staged)", "l");
  leg1->AddEntry(band10.band.get(),   "10 yr 1#sigma band","f");
  leg1->AddEntry(l3, "3#sigma", "l");
  leg1->AddEntry(l5, "5#sigma", "l");
  leg1->Draw("same");

 
  c1->SaveAs(output_pdf);
  std::cout << "Band plot saved to " << output_pdf << "\n";

  // ===== PLOT 2: All 7-year throws =====
  TCanvas *c2 = new TCanvas("c_cp_7yr_all", "CP sensitivity - 7 year throws", 800, 600);
  c2->SetMargin(0.12, 0.04, 0.12, 0.06);
  
  DrawAllThrows(dir7, kAzure-9, ";#delta_{CP}/#pi;#sqrt{#Delta#chi^{2}}", c2, true);
  
  // Overlay median on top
  band7.median->SetLineColor(kBlack);
  band7.median->SetLineWidth(2);
  band7.median->Draw("L same");
  
  auto leg2 = new TLegend(0.65, 0.75, 0.88, 0.88);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->AddEntry((TObject*)0, "7 year throws", "");
  leg2->AddEntry(band7.median.get(), "Median", "l");
  leg2->Draw("SAME");
  
  c2->Update();
  TString pdf7 = output_pdf;
  pdf7.ReplaceAll(".pdf", "_7yr_all.pdf");
  c2->SaveAs(pdf7);
  std::cout << "7-year throws plot saved to " << pdf7 << "\n";

  // ===== PLOT 3: All 10-year throws =====
  TCanvas *c3 = new TCanvas("c_cp_10yr_all", "CP sensitivity - 10 year throws", 800, 600);
  c3->SetMargin(0.12, 0.04, 0.12, 0.06);
  
  DrawAllThrows(dir10, kOrange-9, ";#delta_{CP}/#pi;#sqrt{#Delta#chi^{2}}", c3, true);
  
  // Overlay median on top
  band10.median->SetLineColor(kBlack);
  band10.median->SetLineWidth(3);
  band10.median->Draw("L same");
  
  auto leg3 = new TLegend(0.65, 0.75, 0.88, 0.88);
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->AddEntry((TObject*)0, "10 year throws", "");
  leg3->AddEntry(band10.median.get(), "Median", "l");
  leg3->Draw("SAME");
  
  c3->Update();
  TString pdf10 = output_pdf;
  pdf10.ReplaceAll(".pdf", "_10yr_all.pdf");
  c3->SaveAs(pdf10);
  std::cout << "10-year throws plot saved to " << pdf10 << "\n";

  // ===== PLOT 4: Both overlaid =====
  TCanvas *c4 = new TCanvas("c_cp_both_all", "CP sensitivity - All throws", 800, 600);
  c4->SetMargin(0.12, 0.04, 0.12, 0.06);
  
  DrawAllThrows(dir7, kAzure-9, ";#delta_{CP}/#pi;#sqrt{#Delta#chi^{2}}", c4, true);
  DrawAllThrows(dir10, kOrange-9, "", c4, false);
  
  // Draw medians on top
  band7.median->SetLineColor(kAzure+2);
  band7.median->SetLineWidth(2);
  band7.median->Draw("L same");
  
  band10.median->SetLineColor(kOrange+7);
  band10.median->SetLineWidth(2);
  band10.median->Draw("L same");
  
  auto leg4 = new TLegend(0.60, 0.70, 0.88, 0.88);
  leg4->SetBorderSize(0);
  leg4->SetFillStyle(0);
  leg4->AddEntry((TObject*)0, "Individual throws:", "");
  //leg4->AddEntry((TGraph*)0, "  7 years", "l")->SetLineColor(kAzure-9);
  //leg4->AddEntry((TGraph*)0, "  10 years", "l")->SetLineColor(kOrange-9);
  leg4->AddEntry(band7.median.get(), "7 yr median", "l");
  leg4->AddEntry(band10.median.get(), "10 yr median", "l");
  leg4->Draw("SAME");
  
  c4->Update();
  TString pdf_both = output_pdf;
  pdf_both.ReplaceAll(".pdf", "_both_all.pdf");
  c4->SaveAs(pdf_both);
  std::cout << "Combined throws plot saved to " << pdf_both << "\n";
}


void DrawCPBandExampleII(const char* dir7, const char* dir10,
                       const char* output_pdf = "cp_sensitivity_II.pdf")
{
  //auto band7  = MakeBandFromGraphs(dir7,  kAzure+1, "7 years (median)");
  //auto band10 = MakeBandFromGraphs(dir10, kOrange+7, "10 years (median)");


   std::vector<TGraph*> dir10_figures = Get_figures(dir10);
   TGraph* dir10_figures_aver = AverageGraphs(dir10_figures);
   
   
   TGraphErrors* dir10_figures_uncer =  MakeAverageGraphWithUncertainty(dir10_figures);
  //std::cout << "7-year band uses " << band7.nThrows << " throws\n";
  //std::cout << "10-year band uses " << band10.nThrows << " throws\n";

  // ===== PLOT 1: Band summary =====
  TCanvas *c1 = new TCanvas("c_cp_bands", "CP sensitivity bands", 800, 600);
  c1->SetMargin(0.12, 0.04, 0.12, 0.06);
  bool first = true; 
  int color = 2; 
  int next=0; 
  for(TGraph* g : dir10_figures){
  g->SetLineWidth(3);
  std::cout<<"index:"<<next<< std::endl;
  g->SetLineColor(color);
  if (first ) {
     //g->SetTitle(title);
      g->GetXaxis()->SetTitle("#delta_{CP}/#pi");
      g->GetYaxis()->SetTitle("#sqrt{#Delta#chi^{2}}");
      g->GetXaxis()->SetLimits(-1.0, 1.0);
      g->SetMinimum(0.0);
      g->SetMaximum(10.5);
      g->Draw("AL");
      first = false;
    } else {
      g->Draw("L same");
    }
   color ++; 
   next ++; 
  }
  dir10_figures_uncer->SetFillColor(kAzure+2);
  dir10_figures_uncer->SetFillStyle(3010);

  dir10_figures_uncer->Draw("a4 same");

  dir10_figures_aver->Draw("L same");



 /*
  // Set up axes
  band10.median->SetTitle(";#delta_{CP}/#pi;#sqrt{#Delta#chi^{2}}");

 
  band10.median->GetXaxis()->SetTitleSize(0.045);
  band10.median->GetYaxis()->SetTitleSize(0.045);
  
  band7.median->SetLineColor(kAzure+2);
  band7.median->SetLineWidth(3);
  //band7.median->Draw("L same");
  
  band10.median->SetLineColor(kOrange+7);
  band10.median->SetLineWidth(3);
  //band10.median->Draw("L same");
  // Draw bands and medians
  
  band10.band->SetFillColor(kOrange+7);
  band10.band->SetFillStyle(3005);
  band10.median->GetXaxis()->SetLimits(-1.0, 1.0);
  band10.median->GetYaxis()->SetLimits(0.0, 10.5);
  band10.median->SetMinimum(0.0);
  band10.median->SetMaximum(10.5);

  band10.median->Draw("AL");
  band10.median->GetXaxis()->SetLimits(-1.0, 1.0);
  band10.median->GetYaxis()->SetLimits(0.0, 10.5);
     c1->Update();
  band10.band->Draw("a4 same");
  band10.band->GetXaxis()->SetLimits(-1.0, 1.0);
  band10.band->GetYaxis()->SetLimits(0.0, 10.5);
    
  band10.band->SetMinimum(0.0);
  band10.band->SetMaximum(10.5);

  band10.median->Draw("L same");
  
  band7.band->SetFillColor(kAzure+2);
  band7.band->SetFillStyle(3005);
  band7.band->Draw("a4 same");
  

  
  band7.median->Draw("L same");

  */
  // Reference lines for 3σ and 5σ
  TLine *l3 = new TLine(-1.0, 3.0, 1.0, 3.0);
  l3->SetLineStyle(2);
  l3->SetLineColor(kGray+2);
  l3->Draw("same");
  
  TLine *l5 = new TLine(-1.0, 5.0, 1.0, 5.0);
  l5->SetLineStyle(2);
  l5->SetLineColor(kGray+2);
  l5->Draw("same");

  // Legend
  auto leg1 = new TLegend(0.55, 0.65, 0.88, 0.88);
  //leg1->SetBorderSize(0);
  //leg1->SetFillStyle(0);
  //leg1->AddEntry(band7.median.get(),  "7 years (staged)",  "l");
  //leg1->AddEntry(band7.band.get(),    "7 yr 1#sigma band", "f");
  //leg1->AddEntry(band10.median.get(), "10 years (staged)", "l");
  //leg1->AddEntry(band10.band.get(),   "10 yr 1#sigma band","f");
  //leg1->AddEntry(l3, "3#sigma", "l");
  //leg1->AddEntry(l5, "5#sigma", "l");
  //leg1->Draw("same");

 
  c1->SaveAs(output_pdf);
  std::cout << "Band plot saved to " << output_pdf << "\n";

  
  //std::cout << "Combined throws plot saved to " << pdf_both << "\n";
}

void DrawCPBandExampleIII(const char* dir7, const char* dir10,const char* output_pdf = "cp_sensitivity" )
{
  auto band7  = MakeBandFromGraphs(dir7,  kAzure+1, "7 years (median)");
  auto band10 = MakeBandFromGraphs(dir10, kOrange+7, "10 years (median)");
   char output_pdf_name[1024];
  if (!band7.IsValid() || !band10.IsValid()) {
    std::cerr << "Error: Failed to create bands\n";
    return;
  }

  std::cout << "7-year band uses " << band7.nThrows << " throws\n";
  std::cout << "10-year band uses " << band10.nThrows << " throws\n";

  // ===== PLOT 1: Band summary =====
  TCanvas *c1 = new TCanvas("c_cp_bands", "CP sensitivity bands", 800, 600);
  c1->SetMargin(0.11, 0.04, 0.12, 0.06);

  // Create a dummy histogram for axis setup
  TH1F *hframe = c1->DrawFrame(-1.0, 0.0, 1.0, 9.);
  hframe->SetTitle(";#delta_{CP}/#pi;#sqrt{#Delta#chi^{2}}");
  hframe->GetXaxis()->SetTitleSize(0.045);
  hframe->GetYaxis()->SetTitleSize(0.045);
  hframe->GetXaxis()->SetTitleOffset(1.1);
  hframe->GetYaxis()->SetTitleOffset(1.2);

  // Draw 10-year band with smooth shading
  band10.band->SetFillColorAlpha(kOrange+1, 0.35);  // Semi-transparent
  band10.band->SetFillStyle(1001);  // Solid fill
  band10.band->SetLineWidth(0);     // No border on band
  band10.band->SetMaximum(12.0);
  band10.band->Draw("E3");      // Draw filled area only
  
  band10.median->SetLineColor(kOrange+7);
  band10.median->SetLineWidth(2);
  band10.median->SetLineStyle(1);
  band10.median->SetMaximum(11.0);
  band10.median->Draw("L same");
  
  // Draw 7-year band with smooth shading
  band7.band->SetFillColorAlpha(kAzure+1, 0.35);   // Semi-transparent
  band7.band->SetFillStyle(1001);   // Solid fill
  band7.band->SetLineWidth(0);      // No border on band
  band7.band->Draw("E3 same");       // Draw filled area only
  
  band7.median->SetLineColor(kAzure+2);
  band7.median->SetLineWidth(2);
  band7.median->SetLineStyle(1);
  band7.median->Draw("L same");

  // Reference lines for 3σ and 5σ
  TLine *l3 = new TLine(-1.0, 3.0, 1.0, 3.0);
  l3->SetLineStyle(3);
  l3->SetLineWidth(2);
  l3->SetLineColor(kGray+2);
  l3->Draw("same");
  
  TLine *l5 = new TLine(-1.0, 5.0, 1.0, 5.0);
  l5->SetLineStyle(4);
  l5->SetLineWidth(2);
  l5->SetLineColor(kGray+2);
  l5->Draw("same");

  // Legend with improved styling
  auto leg1 = new TLegend(0.12, 0.73, 0.92, 0.93);
  leg1->SetBorderSize(1);
  leg1->SetNColumns(1);
  leg1->SetTextSize(.039);
  leg1->SetBorderSize(0);
  //leg1->SetFillColorAlpha(kWhite, 0.9);
  TH1D* clone1 =(TH1D*)band7.band.get()->Clone("1");
  TH1D* clone2 =(TH1D*)band10.band.get()->Clone("2");
  clone2->SetLineColor(kOrange+7);
  clone2->SetLineWidth(2);
  clone2->SetLineStyle(1);
  clone2->SetFillColorAlpha(kOrange+1, 0.35);  // Semi-transparent
  clone2->SetFillStyle(1001);
  
  clone1->SetLineColor(kAzure+2);
  clone1->SetLineWidth(2);
  clone1->SetLineStyle(1);
  clone1->SetFillColorAlpha(kAzure+1, 0.35);  // Semi-transparent
  clone1->SetFillStyle(1001);
  
  leg1->AddEntry(clone1,  "7 years[syst:10Flux+xsec+det]",  "lf");
  //leg1->AddEntry(band7.band.get(),    "7 yr band", "f");
  leg1->AddEntry(clone2, "10 years[syst:10Flux+xsec+det]", "lf");
  //leg1->AddEntry(band10.band.get(),   "10 yr band","f");
  leg1->AddEntry(l3, "3#sigma", "l");
  leg1->AddEntry(l5, "5#sigma", "l");
  leg1->Draw("same");

  c1->Update();
  sprintf(output_pdf_name,"%s.pdf(",output_pdf);
  c1->Print(output_pdf_name);
  
  std::cout << "Band plot saved to " << output_pdf << "\n";

  // ===== PLOT 2: All 7-year throws =====
  TCanvas *c2 = new TCanvas("c_cp_7yr_all", "CP sensitivity - 7 year throws", 800, 600);
  c2->SetMargin(0.12, 0.04, 0.12, 0.06);
auto leg2 = new TLegend(0.35, 0.4, 0.88, 0.88);
    leg2->SetNColumns(2);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->AddEntry((TObject*)0, "7 year throws", "");
  leg2->AddEntry(band7.median.get(), "Mean", "l");
  
  DrawAllThrows(dir7, kAzure-9, ";#delta_{CP}/#pi;#sqrt{#Delta#chi^{2}}", c2,leg2, true);
  
  // Overlay median on top
  band7.median->SetLineColor(kBlack);
  band7.median->SetLineWidth(2);
  band7.median->Draw("L same");
  

  leg2->Draw("SAME");
  
  c2->Update();
  sprintf(output_pdf_name,"%s.pdf",output_pdf);
  c2->Print(output_pdf_name);
  


  // ===== PLOT 3: All 10-year throws =====
  TCanvas *c3 = new TCanvas("c_cp_10yr_all", "CP sensitivity - 10 year throws", 800, 600);
  c3->SetMargin(0.12, 0.04, 0.12, 0.06);
    auto leg3 = new TLegend(0.35, 0.3, 0.88, 0.88);
      leg3->SetNColumns(2);
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->AddEntry((TObject*)0, "10 year throws", "");
  leg3->AddEntry(band10.median.get(), "Median", "l");
  
  DrawAllThrows(dir10, kOrange-9, ";#delta_{CP}/#pi;#sqrt{#Delta#chi^{2}}", c3,leg3, true);
  
  // Overlay median on top
  band10.median->SetLineColor(kBlack);
  band10.median->SetLineWidth(3);
  band10.median->Draw("L same");
  

  leg3->Draw("SAME");
  
  c3->Update();
  sprintf(output_pdf_name,"%s.pdf",output_pdf);
  c3->Print(output_pdf_name);
  


  // ===== PLOT 4: Both overlaid =====
  TCanvas *c4 = new TCanvas("c_cp_both_all", "CP sensitivity - All throws", 800, 600);
  c4->SetMargin(0.12, 0.04, 0.12, 0.06);
   auto leg4 = new TLegend(0.2, 0.45, 0.88, 0.88);
     leg4->SetBorderSize(0);
  leg4->SetFillStyle(0);
  leg4->SetNColumns(2);
  leg4->AddEntry((TObject*)0, "Individual throws:", "");
  leg4->AddEntry(band7.median.get(), "7 yr mean", "l");
  leg4->AddEntry(band10.median.get(), "10 yr mean", "l");
   
   
  DrawAllThrows(dir7, kAzure-9, ";#delta_{CP}/#pi;#sqrt{#Delta#chi^{2}}", c4,leg4, true);
  DrawAllThrows(dir10, kOrange-9, "", c4,leg4, false);
  
  // Draw medians on top
  band7.median->SetLineColor(kAzure+2);
  band7.median->SetLineWidth(2);
  band7.median->Draw("L same");
  
  band10.median->SetLineColor(kOrange+7);
  band10.median->SetLineWidth(2);
  band10.median->Draw("L same");
  
 

  leg4->Draw("SAME");
  
  c4->Update();
  sprintf(output_pdf_name,"%s.pdf)",output_pdf);
  c4->Print(output_pdf_name);

}