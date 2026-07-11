//check
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TKey.h"
#include "TObject.h"
#include "TMatrixT.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

std::string CleanName(const std::string& name)
{
  std::string out = name;
  std::replace(out.begin(), out.end(), ' ', '_');
  std::replace(out.begin(), out.end(), '/', '_');
  std::replace(out.begin(), out.end(), '#', '_');
  std::replace(out.begin(), out.end(), ':', '_');
  std::replace(out.begin(), out.end(), '<', '_');
  std::replace(out.begin(), out.end(), '>', '_');
  return out;
}

TH2D* MatrixToHist(const TMatrixT<double>& mat,
                   const char* hname,
                   const char* htitle)
{
  const int nRows = mat.GetNrows();
  const int nCols = mat.GetNcols();

  TH2D* h = new TH2D(hname, htitle, nCols, 0, nCols, nRows, 0, nRows);

  for (int i = 0; i < nRows; ++i) {
    for (int j = 0; j < nCols; ++j) {
      h->SetBinContent(j + 1, i + 1, mat(i, j));
    }
  }

  return h;
}

TMatrixT<double> CovToCorr(const TMatrixT<double>& cov)
{
  const int nRows = cov.GetNrows();
  const int nCols = cov.GetNcols();

  TMatrixT<double> corr(nRows, nCols);
  corr.Zero();

  if (nRows != nCols) {
    std::cerr << "Warning: matrix is not square, cannot build correlation matrix."
              << std::endl;
    return corr;
  }

  for (int i = 0; i < nRows; ++i) {
    for (int j = 0; j < nCols; ++j) {
      const double di = cov(i, i);
      const double dj = cov(j, j);

      if (di > 0.0 && dj > 0.0) {
        corr(i, j) = cov(i, j) / std::sqrt(di * dj);
      } else {
        corr(i, j) = 0.0;
      }
    }
  }

  return corr;
}

void StyleHist(TH2D* h, bool isCorr)
{
  h->SetStats(0);
  h->GetXaxis()->SetTitle("Bin j");
  h->GetYaxis()->SetTitle("Bin i");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetZaxis()->CenterTitle();

  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetZaxis()->SetTitleSize(0.045);

  h->GetXaxis()->SetLabelSize(0.035);
  h->GetYaxis()->SetLabelSize(0.035);
  h->GetZaxis()->SetLabelSize(0.035);

  if (isCorr) {
    h->SetMinimum(-1.0);
    h->SetMaximum(1.0);
  }
}

void DrawCovAndCorrToImages(const char* filename = "input.root",
                            const char* outdir   = "matrix_plots",
                            const char* imgext   = "png")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".2f");

  gSystem->mkdir(outdir, kTRUE);

  TFile* file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: could not open file " << filename << std::endl;
    return;
  }

  TIter next(file->GetListOfKeys());
  TKey* key = nullptr;

  int matrixCount = 0;

  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();
    if (!obj) continue;

    if (!obj->InheritsFrom("TMatrixT<double>")) {
      delete obj;
      continue;
    }

    TMatrixT<double>* cov = (TMatrixT<double>*)obj;

    // Use the ROOT key name, not obj->GetName()
    std::string rootName = key->GetName();
    std::string cleanName = CleanName(rootName);

    std::cout << "Processing matrix: " << rootName
              << " (" << cov->GetNrows() << " x " << cov->GetNcols() << ")"
              << std::endl;

    TH2D* hCov = MatrixToHist(
      *cov,
      Form("hCov_%s", cleanName.c_str()),
      Form("%s;Bin j;Bin i;Covariance", rootName.c_str())
    );
    StyleHist(hCov, false);

    TMatrixT<double> corr = CovToCorr(*cov);
    TH2D* hCorr = MatrixToHist(
      corr,
      Form("hCorr_%s", cleanName.c_str()),
      Form("%s;Bin j;Bin i;Correlation", rootName.c_str())
    );
    StyleHist(hCorr, true);

    TCanvas* cCov = new TCanvas(Form("cCov_%s", cleanName.c_str()), "", 900, 700);
    cCov->SetRightMargin(0.15);
    cCov->SetLeftMargin(0.12);
    cCov->SetBottomMargin(0.12);
    hCov->Draw("COLZ");
    cCov->SaveAs(Form("%s/%s_cov.%s", outdir, cleanName.c_str(), imgext));

    TCanvas* cCorr = new TCanvas(Form("cCorr_%s", cleanName.c_str()), "", 900, 700);
    cCorr->SetRightMargin(0.15);
    cCorr->SetLeftMargin(0.12);
    cCorr->SetBottomMargin(0.12);
    hCorr->Draw("COLZ");
    cCorr->SaveAs(Form("%s/%s_corr.%s", outdir, cleanName.c_str(), imgext));

    delete cCov;
    delete cCorr;
    delete hCov;
    delete hCorr;
    delete obj;

    ++matrixCount;
  }

  file->Close();
  delete file;

  std::cout << "Saved " << matrixCount << " covariance/correlation plot pairs in: "
            << outdir << std::endl;
}