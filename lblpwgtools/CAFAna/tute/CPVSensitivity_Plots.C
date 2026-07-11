#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>
#include <map>
#include <string>
#include <vector>



#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"

#include "CAFAna/Cuts/TruthCuts.h"

#include "StandardRecord/StandardRecord.h"

#include "TH1.h"
#include "TPad.h"
#include <string>
#include "TAxis.h"
#include "THStack.h"
#include "TLegend.h"


void CPVSensitivity_Plots() {
    // Set ROOT style
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    double YMax = 10;
    
    double XtitleSize = 0.05;
    double YtitleSize = 0.05;
    char pdf_title[1024];
  std::string Pdf_name = "DUNE_DeltaCP";
  //TCanvas *c1 = new TCanvas("c1");
  TCanvas *c1 = new TCanvas("c1", "Original Cases", 1000, 600);
   c1->SetMargin(0.12, 0.03, 0.12, 0.06); // (left, right, bottom, top)
   sprintf(pdf_title, "%s.pdf(", Pdf_name.c_str());
   c1 -> Print(pdf_title);
   sprintf(pdf_title, "%s.pdf", Pdf_name.c_str());

    
    
    
    // Original file set
    std::map<std::string, std::string> original_filenames = {
        {"All Syst", "State_allsyst_624ktmwyr_2nd.root"},
        {"No Syst", "State_nosyst_624ktmwyr_2nd.root"},
        {"Flux only", "State_noxsec:nodet_624ktmwyr_2nd.root"},
        {"Flux + Det", "State_noxsec_624ktmwyr_2nd.root"}
    };
    
    // My own file set
    std::map<std::string, std::string> myown_filenames = {
        {"All Syst", "myown_all_allsyst_624ktmwyr_2nd.root"},
        {"No Syst", "myown_all_nosyst_624ktmwyr_2nd.root"},
        {"Flux only", "myown_all_noxsec:nodet_624ktmwyr_2nd.root"},
        {"Flux + Det", "myown_all_noxsec_624ktmwyr_2nd.root"}
    };
    
    std::string Path_files = "/pnfs/dune/persistent/users/cnguyen/cafana_example/";
    // Case colors
    std::vector<int> case_colors = {kGreen, kRed, kCyan, kMagenta};
    std::vector<std::string> labels = {"All Syst", "No Syst", "Flux only", "Flux + Det"};
    
    // ========== Plot all original cases together ==========
    //TCanvas *c1 = new TCanvas("c1", "Original Cases", 1000, 600);
    TLegend *leg1 = new TLegend(0.15, 0.7, 0.7, 0.88);
    leg1->SetTextSize(0.035);
    leg1->SetNColumns(2);
    bool first = true;
    for (size_t i = 0; i < labels.size(); i++) {
    std::string fullpath = Path_files + original_filenames[labels[i]]; 
        TFile *f = TFile::Open(fullpath.c_str());
        if (!f || f->IsZombie()) continue;
        
        TGraph *g = (TGraph*)f->Get("sens_cpv_nh");
        if (!g) {
            f->Close();
            continue;
        }
        
        g->SetLineColor(case_colors[i]);
        g->SetLineWidth(2);
        g->SetTitle("CP Violation Sensitivity: All Original Cases");
        g->GetXaxis()->SetTitle("#delta_{CP}/#pi");
        g->GetYaxis()->SetTitle("Sensitivity");
        g->GetXaxis()->SetTitleSize(XtitleSize);
        
        g->GetYaxis()->SetTitleSize(YtitleSize);
        g->SetMaximum(YMax);
        if (first) {
            g->Draw("AL");
            first = false;
        } else {
            g->Draw("L SAME");
        }
        
        leg1->AddEntry(g, ("Original " + labels[i]).c_str(), "l");
        f->Close();
    }
    
    leg1->Draw();
    c1->SetGrid();
    c1 -> Print(pdf_title);
    //c1->SaveAs("cpv_sensitivity_all_original_cases.png");
    
    // ========== Plot all myown cases together ==========
    //TCanvas *c2 = new TCanvas("c2", "My Own Cases", 1000, 600);
    TLegend *leg2 = new TLegend(0.15, 0.7, 0.7, 0.88);
    leg2->SetTextSize(0.035);
    leg2->SetNColumns(2);
    first = true;
    for (size_t i = 0; i < labels.size(); i++) {
     std::string fullpath2 = Path_files + myown_filenames[labels[i]]; 
        TFile *f = TFile::Open(fullpath2.c_str());
        if (!f || f->IsZombie()) continue;
        
        TGraph *g = (TGraph*)f->Get("sens_cpv_nh");
        if (!g) {
            f->Close();
            continue;
        }
        
        g->SetLineColor(case_colors[i]);
        g->SetLineWidth(2);
        g->SetTitle("CP Violation Sensitivity: All My Own Cases");
        g->GetXaxis()->SetTitle("#delta_{CP}/#pi");
        g->GetYaxis()->SetTitle("Sensitivity");
        g->GetXaxis()->SetTitleSize(XtitleSize);
        g->GetYaxis()->SetTitleSize(YtitleSize);
        g->SetMaximum(YMax);
        if (first) {
            g->Draw("AL");
            first = false;
        } else {
            g->Draw("L SAME");
        }
        
        leg2->AddEntry(g, ("My Own " + labels[i]).c_str(), "l");
        f->Close();
    }
    
    leg2->Draw();
    c1->SetGrid();
    c1 -> Print(pdf_title);
    //c2->SaveAs("cpv_sensitivity_all_myown_cases.png");
    
    // ========== Individual comparison plots ==========
    for (size_t i = 0; i < labels.size(); i++) {
        std::string label = labels[i];
        std::string canvas_name = "c_" + std::to_string(i+3);
        //TCanvas *c = new TCanvas(canvas_name.c_str(), label.c_str(), 1000, 600);
        TLegend *leg = new TLegend(0.15, 0.75, 0.45, 0.88);
        leg->SetTextSize(0.035);
        
        // Load original file
        std::string fullpath3 = Path_files + original_filenames[labels[i]]; 
        TFile *f_orig = TFile::Open(fullpath3.c_str());
        if (!f_orig || f_orig->IsZombie()) continue;
        
        TGraph *g_orig = (TGraph*)f_orig->Get("sens_cpv_nh");
        if (!g_orig) {
            f_orig->Close();
            continue;
        }
        
        // Load myown file
        std::string fullpath4 = Path_files + myown_filenames[labels[i]]; 
        TFile *f_myown = TFile::Open(fullpath4.c_str());
        if (!f_myown || f_myown->IsZombie()) {
            f_orig->Close();
            continue;
        }
        
        TGraph *g_myown = (TGraph*)f_myown->Get("sens_cpv_nh");
        if (!g_myown) {
            f_orig->Close();
            f_myown->Close();
            continue;
        }
        
        // Configure original graph
        g_orig->SetLineColor(kRed);
        g_orig->SetLineWidth(2);
        g_orig->SetTitle(("CP Violation Sensitivity: Original vs My Own - " + label).c_str());
        g_orig->GetXaxis()->SetTitle("#delta_{CP}/#pi");
        g_orig->GetYaxis()->SetTitle("Sensitivity");
        g_orig->GetXaxis()->SetTitleSize(XtitleSize);
        g_orig->GetYaxis()->SetTitleSize(YtitleSize);
        g_orig->SetMaximum(YMax);
        // Configure myown graph
        g_myown->SetLineColor(kBlue);
        g_myown->SetLineWidth(2);
        g_myown->SetLineStyle(2); // dashed
        
        // Draw graphs
        g_orig->Draw("AL");
        g_myown->Draw("L SAME");
        
        // Add legend entries
        leg->AddEntry(g_orig, ("Original " + label).c_str(), "l");
        leg->AddEntry(g_myown, ("My Own " + label).c_str(), "l");
        leg->Draw();
        c1 -> Print(pdf_title);
        
        // Save with sanitized filename
        std::string filename = "cpv_sensitivity_comparison_" + label;
        std::replace(filename.begin(), filename.end(), ' ', '_');
        //filename += ".png";
        //c->SaveAs(filename.c_str());
        c1 -> Print(pdf_title);
        f_orig->Close();
        f_myown->Close();
    }
    
    
   sprintf(pdf_title, "%s.pdf)", Pdf_name.c_str());
   c1 -> Print(pdf_title);
    
    
}// ENd of Function 

/*
int main() {
    CPVSensitivity_Plots();
    return 0;
}
*/