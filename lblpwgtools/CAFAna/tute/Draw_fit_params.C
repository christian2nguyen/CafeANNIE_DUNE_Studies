//Checking
#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TClass.h"
#include "TROOT.h"
#include "TStyle.h"
#include <iostream>
#include <vector>
#include <map>
#include <string>

void DrawHistogramsFromDirectory(TFile* file,
                                 const char* dirName,
                                 const char* pdfOutput )
{
    if (!file || file->IsZombie()) {
        std::cerr << "Error: invalid or zombie TFile." << std::endl;
        return;
    }

    TDirectory* dir = file->GetDirectory(dirName);
    if (!dir) {
        std::cerr << "Error: directory '" << dirName << "' not found." << std::endl;
        return;
    }

    std::cout << "Reading directory: " << dirName << std::endl;

    // -----------------------------------------------------------------------
    // Collect histograms
    // -----------------------------------------------------------------------
    std::map<std::string, std::vector<TH1D*>> th1dMap;
    std::vector<std::pair<std::string, TH2D*>> th2dList;

    TIter  next(dir->GetListOfKeys());
    TKey*  key = nullptr;

    while ((key = (TKey*)next())) {
        TClass* cl = gROOT->GetClass(key->GetClassName());
        if (!cl) continue;

        std::string objName = key->GetName();
        int         cycle   = key->GetCycle();

        if (cl->InheritsFrom("TH2D")) {
            TH2D* h = dynamic_cast<TH2D*>(key->ReadObj());
            if (h) {
                h->SetDirectory(nullptr);
                th2dList.push_back({objName + ";" + std::to_string(cycle), h});
                std::cout << "  TH2D: " << objName << ";" << cycle << std::endl;
            }
        }
        else if (cl->InheritsFrom("TH1D")) {
            TH1D* h = dynamic_cast<TH1D*>(key->ReadObj());
            if (h) {
                h->SetDirectory(nullptr);
                th1dMap[objName].push_back(h);
                std::cout << "  TH1D: " << objName << ";" << cycle << std::endl;
            }
        }
    }

    // -----------------------------------------------------------------------
    // Style
    // -----------------------------------------------------------------------
    gStyle->SetOptStat(0);   // turn off the stat box — it would overlap when overlaying
    gStyle->SetOptTitle(1);

    const int kColors[4] = { kBlue+1, kRed+1, kGreen+2, kOrange+7 };
    const int kFills [4] = { 3004,    3005,    3006,     3007      };

    // -----------------------------------------------------------------------
    // Open multi-page PDF
    // -----------------------------------------------------------------------
    TCanvas* dummy = new TCanvas("_pdfDriver", "", 1, 1);
    std::string openCmd  = std::string(pdfOutput) + "(";
    std::string closeCmd = std::string(pdfOutput) + ")";
    dummy->Print(openCmd.c_str());

    // -----------------------------------------------------------------------
    // TH1D: one page per base name, all cycles overlaid on the same pad
    // -----------------------------------------------------------------------
  for (auto& kv1 : th1dMap) {
    const std::string&        baseName = kv1.first;
    std::vector<TH1D*>&       histVec  = kv1.second;

        int nHists = (int)histVec.size();

        TCanvas* c = new TCanvas(("c1D_" + baseName).c_str(),
                                 ("TH1D: " + baseName).c_str(),
                                 800, 600);
        c->SetLeftMargin(0.12);
        c->SetBottomMargin(0.12);

        // Find the global Y maximum across all cycles so the first Draw
        // sets an axis range that accommodates every histogram
        double yMax = 0.0;
        for (TH1D* h : histVec){
            h->Scale(1,"width");
            yMax = std::max(yMax, h->GetMaximum());
		}
        // Legend — placed top-right, sized to the number of cycles
        TLegend* leg = new TLegend(0.35, 0.88 - nHists * 0.07, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(1001);
        leg->SetTextSize(0.035);

        for (int i = 0; i < nHists; ++i) {
            TH1D* h = histVec[i];
            h->SetLineColor(kColors[i % 4]);
           // h->SetFillColor(kColors[i % 4]);
            h->SetFillStyle(kFills [i % 4]);
            h->SetLineWidth(2);

            // Let the first histogram define the axes; the rest use "SAME"
            if (i == 0) {
                h->SetTitle(baseName.c_str());
                 std::string xTitle = (h->GetXaxis()->GetXmax() <= 10.0) ? "E_{#nu}[GeV]" : "Bin Number";
                h->GetXaxis()->SetTitle(xTitle.c_str());
                h->GetYaxis()->SetTitle("Events / Bin Width");
                h->GetYaxis()->SetRangeUser(0, yMax * 1.7);  // 15% headroom
                h->Draw("HIST");
            } else {
                h->Draw("HIST SAME");
            }
            
            std::string legLabel = std::string(h->GetTitle()) + "  cycle: " + std::to_string(i + 1);
             leg->AddEntry(h, legLabel.c_str(), "lf");
        }

        leg->Draw();
        c->Update();
        c->Print(pdfOutput);
        delete c;
    }

    // -----------------------------------------------------------------------
    // TH2D: one page per histogram (unchanged)
    // -----------------------------------------------------------------------
    for (auto& kv2 : th2dList) {
    const std::string& label = kv2.first;
    TH2D*              h     = kv2.second;

        std::string cName = "c2D_" + label;
        for (char& ch : cName) if (ch == ';') ch = '_';

        TCanvas* c = new TCanvas(cName.c_str(),
                                 ("TH2D: " + label).c_str(),
                                 800, 650);
        c->SetLeftMargin(0.12);
        c->SetRightMargin(0.16);
        c->SetBottomMargin(0.12);
        
         h->SetTitle(("TH2D: " + label).c_str());
         //std::cout<<"h->GetXaxis()->GetXmax() = "<< h->GetXaxis()->GetXmax() << std::endl;
        std::string xTitle = (h->GetXaxis()->GetXmax() == 10.0) ? "E_{#nu} [GeV]" : "Bin Number";
        h->GetXaxis()->SetTitle(xTitle.c_str());
        
        std::string yTitle = (h->GetYaxis()->GetXmax() <= 1.0) ? "RECO Bjorken y (inelasticity)" : "Bin Number";
        h->GetYaxis()->SetTitle(yTitle.c_str());

       
        h->Draw("COLZ");

        TLegend* leg = new TLegend(0.13, 0.88, 0.55, 0.96);
        leg->SetBorderSize(1);
        leg->SetFillStyle(1001);
        leg->SetTextSize(0.032);
        leg->AddEntry(h, label.c_str(), "f");
        //leg->Draw();

        c->Update();
        c->Print(pdfOutput);
        delete c;
    }

    // -----------------------------------------------------------------------
    // Close PDF
    // -----------------------------------------------------------------------
    dummy->Print(closeCmd.c_str());
    delete dummy;

    std::cout << "\nSaved: " << pdfOutput << std::endl;
    std::cout << "  TH1D groups : " << th1dMap.size()  << " page(s)" << std::endl;
    std::cout << "  TH2D hists  : " << th2dList.size() << " page(s)" << std::endl;
}

void Draw_fit_params(){


std::cout<<"drawing fit params "<< std::endl;


std::string MC_file3 =  " /pnfs/dune/persistent/users/cnguyen/cafana_example/StateFeb10_OUTPUT_kNuFitTh23LoNH_624ktmwyr_ndfd:10year_det:xsec:nflux=30_nopen_asimov_set_1.root";
  char inputName[1024];
  
  
  sprintf(inputName, "%s", MC_file3.c_str());
  std::cout<<"InputFileName MC = "<< inputName << std::endl;
  TFile *TFile_MC_fit = new TFile(inputName);

   DrawHistogramsFromDirectory(TFile_MC_fit,
                                 "fit_0.000000",
                                 "histograms_fit_0.000000.pdf");

}