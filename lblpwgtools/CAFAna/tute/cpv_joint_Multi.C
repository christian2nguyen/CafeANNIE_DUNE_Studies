#include "CAFAna/Analysis/common_fit_definitions.h"


#include "TPad.h"
#include <string>
#include "TAxis.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"
#include <vector>

#include "TCanvas.h"
#include "TH1.h"

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


#include <algorithm>

using namespace ana;

std::string GetAsimovName(int asimov_set, int hie)
{
    if (asimov_set == 1 && hie == +1) return "kNuFitTh23LoNH";
    if (asimov_set == 2 && hie == +1) return "kNuFitTh23HiNH";
    if (asimov_set == 1 && hie == -1) return "kNuFitTh23LoIH";
    if (asimov_set == 2 && hie == -1) return "kNuFitTh23HiIH";

    if (asimov_set == 3) return "kNuFitTh23MM";

    if (asimov_set == 4 && hie == +1) return "kNuFitTh23MinNH";
    if (asimov_set == 4 && hie == -1) return "kNuFitTh23MinIH";

    if (asimov_set == 5 && hie == +1) return "kNuFitTh23MaxNH";
    if (asimov_set == 5 && hie == -1) return "kNuFitTh23MaxIH";

    if (asimov_set == 6 && hie == +1) return "kNuFitTh13MinNH";
    if (asimov_set == 6 && hie == -1) return "kNuFitTh13MinIH";

    if (asimov_set == 7 && hie == +1) return "kNuFitTh13MaxNH";
    if (asimov_set == 7 && hie == -1) return "kNuFitTh13MaxIH";

    if (asimov_set == 8 && hie == +1) return "kNuFitDmsq32MinNH";
    if (asimov_set == 8 && hie == -1) return "kNuFitDmsq32MinIH";

    if (asimov_set == 9 && hie == +1) return "kNuFitDmsq32MaxNH";
    if (asimov_set == 9 && hie == -1) return "kNuFitDmsq32MaxIH";

    return "UNKNOWN";
}


struct VectorCompareResult {
    bool same;
    std::vector<std::string> onlyInA;
    std::vector<std::string> onlyInB;
};

VectorCompareResult compareVectorsAsSets(const std::vector<std::string>& a,
                                         const std::vector<std::string>& b)
{
    VectorCompareResult result;
    
    for (const auto& s : a)
        if (std::find(b.begin(), b.end(), s) == b.end())
            result.onlyInA.push_back(s);

    for (const auto& s : b)
        if (std::find(a.begin(), a.end(), s) == a.end())
            result.onlyInB.push_back(s);

    result.same = result.onlyInA.empty() && result.onlyInB.empty();
    return result;
}





void cpv_joint(int hie=9, std::string sampleString="ndfd:10year", 
std::string outputFname_input="State_",
std::string stateFname  = "/pnfs/dune/persistent/users/cnguyen/cafana_example/State"/*="/pnfs/dune/persistent/users/cnguyen/cafana_example/toy_2perEnergy_test_FD_FHC_v1.root"*//*"common_state_mcc11v3.root"*/,
	       std::string systSet = "allsyst"/*allsyst,nosyst*/, 
	       std::string penaltyString="nopen",  std::string asimov_set="0",
	       int isetCP = 0,
	       std::string fakeDataShift = "", int fitBias = 0){

std::string paraset_name = GetAsimovName(hie, 1);
std::string outputFname =  outputFname_input + "State_"+ paraset_name + "_624ktmwyr_" + sampleString + "_" + systSet  + "_" + penaltyString + ".root";
 std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;
std::cout<<"make File: "<< outputFname<< std::endl;
 std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;
 char pdf_title[1024];
 std::string Pdf_name = "DUNEexample_cpv_joint_test_FHC";
 // TCanvas *c1 = new TCanvas("c1");
   //c1->SetMargin(0.12, 0.03, 0.12, 0.06); // (left, right, bottom, top)
 sprintf(pdf_title, "%s.pdf(", Pdf_name.c_str());
 //c1 -> Print(pdf_title);
 sprintf(pdf_title, "%s.pdf", Pdf_name.c_str());



  gROOT->SetBatch(1);
  gRandom->SetSeed(0);

  // Allow a fake data bias
  SystShifts trueSyst = GetFakeDataSystShift(fakeDataShift);

  // Get the systematics to use
  std::vector<const ISyst*> systlist = GetListOfSysts(systSet);
  std::vector<const ISyst*> systlist_flux = GetListOfSysts(systSet);
  
  std::vector<const ISyst *> systlist_det = GetListOfSysts(false, false,
                                          true, true, true,
                                          false, false,
                                          false, false,
                                          false);
  
    std::vector<const ISyst *> systlist_xsec = GetListOfSysts(false, true,
                                          false, false, false,
                                          false, false,
                                          false, false,
                                          false);
                                          
  std::vector<std::string> syslist_vector{"flux","det","xsec"};
                                          
  std::map<std::string,std::vector<const ISyst *>> ISYSYT_Map;
  ISYSYT_Map[syslist_vector[0]]=systlist_flux;
  ISYSYT_Map[syslist_vector[1]]=systlist_det;
  ISYSYT_Map[syslist_vector[2]]=systlist_xsec;
                                          
  //std::vector<const ISyst *> GetListOfSysts(bool fluxsyst_Nov17, bool xsecsyst,
  //                                        bool detsyst, bool useND, bool useFD,
  //                                        bool useNueOnE, bool useFakeDataDials,
  //                                        bool fluxsyst_CDR, int NFluxSysts,
  //                                        bool removeFDNonFitDials)
  //
  
  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;
  std::cout<<"~~~~~~~~~~~~~~~~~Starting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;
  std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;
  std::cout<<"~~~~~~~systlist.size() = "<<systlist.size()<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;
   std::vector<std::string> non_flux_systs;
   std::vector<std::string> systs_check;
   std::vector<std::string> systs_complete_check;
 
  std::cout<<"~~~~~~~~~~~~~~~~~ALL  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;
  for(auto syst:systlist){
   std::cout<<"syst name : "<<syst->ShortName() << std::endl;
      systs_complete_check.push_back(syst->ShortName());
   if (syst->ShortName().find("flux_Nov17_") == 0){}
   else{non_flux_systs.push_back(syst->ShortName()); }
  }
 
  RemoveSysts(systlist_flux, non_flux_systs);
  
    std::cout<<"~~~~~~~~~~~~~~~~~ flux  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;
    for(auto syst:systlist_flux){
   std::cout<<"syst name : "<<syst->ShortName() << std::endl;
   systs_check.push_back(syst->ShortName());
    }
     std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;
    
     std::cout<<"~~~~~~~~~~~~~~~~~systlist_det~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;
       for(auto syst:systlist_det){
   std::cout<<"syst name : "<<syst->ShortName() << std::endl;
   systs_check.push_back(syst->ShortName());
    }
       std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;
              std::cout<<"~~~~~~~~~~xsec~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;
      for(auto syst:systlist_xsec){
   std::cout<<"syst name : "<<syst->ShortName() << std::endl;
   systs_check.push_back(syst->ShortName());
    }
       std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;

 auto cmp = compareVectorsAsSets(systs_check, systs_complete_check);

 if (cmp.same) {
    std::cout << "Same!\n";
 } else {
    std::cout << "Only in systs_check:\n";
    for (auto& s : cmp.onlyInA) std::cout << "  " << s << "\n";

    std::cout << "Only in systs_complete_check:\n";
    for (auto& s : cmp.onlyInB) std::cout << "  " << s << "\n";
 }


  // Should the fake data dials be removed from the systlist?
  if (!fitBias){

    std::vector<std::string> bias_syst_names;
    // Loop over all systs set to a non-nominal value and remove
    for (auto syst : trueSyst.ActiveSysts()){
      std::cout << "Removing " << syst->ShortName() <<std::endl;
      bias_syst_names.push_back(syst->ShortName());
    }
    
    RemoveSysts(systlist, bias_syst_names);
  }



  TFile* fout = new TFile(outputFname.c_str(), "RECREATE");
  fout->cd();

  // First find the minimum for dcp = 0
  std::map<const IFitVar*, std::vector<double>> oscSeeds = {};
  oscSeeds[&kFitSinSqTheta23] = {0.4, 0.6};

  osc::IOscCalculatorAdjustable* trueOscGlob = NuFitOscCalc(hie, 1, asimov_set);
  trueOscGlob->SetdCP(0.0);
  std::map<std::string,double>glob_chisqmin_sys_map;
  
  double glob_chisqmin = 99999;
  double thischisq,thischisq_sep;

  // Now loop over all true values
  int nsteps = 36;
  double dcpstep = 2*TMath::Pi()/nsteps;
  TGraph* gCPV = new TGraph();
  // After the main loop, before fout->Close()
  std::map<std::string, TGraph*> gCPV_per_syst; // want to save indivdiual 
  std::cout<<"~~~~~~~~~~~~~Trying to loop throught univeres ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;
  
  for(auto syst_name: syslist_vector){
   gCPV_per_syst[syst_name] = new TGraph();
   glob_chisqmin_sys_map[syst_name]= 99999;
  }
  
//////////////////////////////////////////////////////////
// Have to loop because I need to change the testOsc dcp value
  for (int idcp = 0; idcp < 2; ++idcp) {
  //for (int idcp = 0; idcp < 1; ++idcp) {
    for(int ihie = -1; ihie <= +1; ihie += 2) {
      double dcptest = idcp*TMath::Pi();
      //double dcptest = -TMath::Pi() + isetCP*dcpstep;

      std::vector<const IFitVar*> oscVars = GetOscVars("th23:th13:dmsq32", ihie);      
      osc::IOscCalculatorAdjustable* testOsc = NuFitOscCalc(ihie, 1, asimov_set);
      testOsc->SetdCP(dcptest);
      
      IExperiment *penalty = GetPenalty(ihie, 1, penaltyString, asimov_set);
      SystShifts testSyst = kNoShift;
      
      thischisq = RunFitPoint(stateFname, sampleString,
			      trueOscGlob, trueSyst, false,
			      oscVars, systlist,
			      testOsc, testSyst,
			      oscSeeds, penalty, 
			      Fitter::kNormal, nullptr);
      
      glob_chisqmin = TMath::Min(thischisq,glob_chisqmin);
      
      for(auto syslist_name : syslist_vector){
      
            thischisq_sep = RunFitPoint(stateFname, sampleString,
			      trueOscGlob, trueSyst, false,
			      oscVars, ISYSYT_Map[syslist_name],
			      testOsc, testSyst,
			      oscSeeds, penalty, 
			      Fitter::kNormal, nullptr);
      
           glob_chisqmin_sys_map[syslist_name] = TMath::Min(thischisq_sep,glob_chisqmin_sys_map[syslist_name]);
      
      }
      
      delete penalty;
      delete testOsc;
    }
  }
  
  delete trueOscGlob;
  std::cout << "Global chi2 at dCP=0,pi: " << glob_chisqmin << std::endl;


  for(double idcp = 0; idcp < nsteps+1; ++idcp) {
    std::cout << "Trying idcp = " << idcp << std::endl;
    
    double thisdcp = -TMath::Pi() + idcp*dcpstep;
    osc::IOscCalculatorAdjustable* trueOsc = NuFitOscCalc(hie, 1, asimov_set);
    trueOsc->SetdCP(thisdcp);

    double chisqmin = 99999.0;
    std::map<std::string,double > chisqmin_map;
    chisqmin_map[syslist_vector[0]]= 99999.0;
    chisqmin_map[syslist_vector[1]]= 99999.0;
    chisqmin_map[syslist_vector[2]]= 99999.0;
    
    // Still need to loop over dcp choices
    for (int idcp = 0; idcp < 2; ++idcp) {
    //for (int idcp = 0; idcp < 1; ++idcp) {
      double dcptest = idcp*TMath::Pi();
      //double dcptest = -TMath::Pi() + isetCP*dcpstep;
      std::cout<< "idcp (index) = "<<idcp << " ---- the true and test dcp values: "<<thisdcp<<" "<<dcptest<<std::endl;

      for(int ihie = -1; ihie <= +1; ihie += 2) {

	std::vector<const IFitVar*> oscVars = GetOscVars("th23:th13:dmsq32", ihie);	
	osc::IOscCalculatorAdjustable* testOsc = NuFitOscCalc(ihie, 1, asimov_set);
	testOsc->SetdCP(dcptest);
	
	IExperiment *penalty = GetPenalty(ihie, 1, penaltyString, asimov_set);
	SystShifts testSyst = kNoShift;
	
	thischisq = RunFitPoint(stateFname, sampleString,
				trueOsc, trueSyst, false,
				oscVars, systlist,
				testOsc, testSyst,
				oscSeeds, penalty, Fitter::kNormal, nullptr);
	
	chisqmin = TMath::Min(thischisq,chisqmin);
	
	
	 for(auto syslist_name :syslist_vector){
	 
	 
	 	thischisq_sep = RunFitPoint(stateFname, sampleString,
				trueOsc, trueSyst, false,
				oscVars, ISYSYT_Map[syslist_name],
				testOsc, testSyst,
				oscSeeds, penalty, Fitter::kNormal, nullptr);
	
	    chisqmin_map[syslist_name] = TMath::Min(thischisq_sep,chisqmin_map[syslist_name]);
	 
	 
	 
	 }
	
	delete penalty;
	delete testOsc;
	
	

	
      }
    }// end of idcp inner 
    
    
    
    
    delete trueOsc;
    std::cout<<" ---- current chi2 "<<chisqmin<<std::endl;   
 
    chisqmin = TMath::Max(chisqmin,1e-6);
    double diff = chisqmin-glob_chisqmin;
    diff = TMath::Max(diff, 0.);
    gCPV->SetPoint(gCPV->GetN(),thisdcp/TMath::Pi(),sqrt(diff));
    
    for(auto syslist_name :syslist_vector){
        chisqmin_map[syslist_name] = TMath::Max(chisqmin_map[syslist_name],1e-6);
       double diff = chisqmin_map[syslist_name]-glob_chisqmin_sys_map[syslist_name];
        diff = TMath::Max(diff, 0.);
        gCPV_per_syst[syslist_name]->SetPoint(gCPV_per_syst[syslist_name]->GetN(),thisdcp/TMath::Pi(),sqrt(diff));
    }
    
    
    
  }
  /////////////// end of idcp ouuter 

  fout->cd();
  gCPV->Draw("ALP");
  //c1 -> Print(pdf_title);
  gCPV->Write(hie > 0 ? "sens_cpv_nh" : "sens_cpv_ih");
  
  for(auto syslist_name :syslist_vector){
  std::string name = hie > 0 ? "sens_cpv_nh" : "sens_cpv_ih";
  std::string final_name = name  + "_" + syslist_name;
  gCPV_per_syst[syslist_name]->Draw("ALP");
  //c1 -> Print(pdf_title);
  gCPV_per_syst[syslist_name]->Write(final_name.c_str());
  
  }
  
  
  fout->Close();



 sprintf(pdf_title, "%s.pdf)", Pdf_name.c_str());
 //c1 -> Print(pdf_title);
 //c1->Close(); 

}

/*int main(int argc, char* argv[]){
  //const char* stateFname = argv[1];
  //const char* outputFname = argv[2];
  //const char* systSet = argv[3];
  //int icp = atoi(argv[4]);
  
  
  std::string input_stateFname = "/pnfs/dune/persistent/users/cnguyen/cafana_example/State_";
  std::string input_outputFname = "State_kNuFitDmsq32MaxNH_nosyst_624ktmwyr_test_v2.root";
  
  int icp_input = 9 ; 
   std::string systSet = "nosyst";
  // if this value is (pos) neg means (normal)inverted mass ordering 
  cpv_joint(input_stateFname, input_outputFname,
            systSet, "ndfd:10year",
            "nopen", 1, "0", icp_input);
    // icp_input - hits this if statments 
    //if (asimov_set == 1 && hie == +1) ret->SetTh23(kNuFitTh23LoNH);
    //if (asimov_set == 2 && hie == +1) ret->SetTh23(kNuFitTh23HiNH);
    //if (asimov_set == 1 && hie == -1) ret->SetTh23(kNuFitTh23LoIH);
    //if (asimov_set == 2 && hie == -1) ret->SetTh23(kNuFitTh23HiIH);
    //if (asimov_set == 3) ret->SetTh23(kNuFitTh23MM);
    //if (asimov_set == 4 && hie == +1) ret->SetTh23(kNuFitTh23MinNH);
    //if (asimov_set == 4 && hie == -1) ret->SetTh23(kNuFitTh23MinIH);
    //if (asimov_set == 5 && hie == +1) ret->SetTh23(kNuFitTh23MaxNH);
    //if (asimov_set == 5 && hie == -1) ret->SetTh23(kNuFitTh23MaxIH);
    //if (asimov_set == 6 && hie == +1) ret->SetTh13(kNuFitTh13MinNH);
    //if (asimov_set == 6 && hie == -1) ret->SetTh13(kNuFitTh13MinIH);
    //if (asimov_set == 7 && hie == +1) ret->SetTh13(kNuFitTh13MaxNH);
    //if (asimov_set == 7 && hie == -1) ret->SetTh13(kNuFitTh13MaxIH);
    //if (asimov_set == 8 && hie == +1) ret->SetDmsq32(kNuFitDmsq32MinNH);
    //if (asimov_set == 8 && hie == -1) ret->SetDmsq32(kNuFitDmsq32MinIH);
    //if (asimov_set == 9 && hie == +1) ret->SetDmsq32(kNuFitDmsq32MaxNH);
    //if (asimov_set == 9 && hie == -1) ret->SetDmsq32(kNuFitDmsq32MaxIH);
    //if (asimov_set == 10) ret->SetdCP(0);
    //if (asimov_set == 11) ret->SetdCP(-TMath::Pi()/2); 
            
            
}

*/

void cpv_joint_Multi(){
std::string outputname = "/pnfs/dune/persistent/users/cnguyen/cafana_example/State_";
std::string sampleString_10="ndfd:10year";
std::string sampleString_7="ndfd:7year";
std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;


std::cout<<"Starting Loop "<< std::endl;
for(int i = 1; i < 10; i++){
std::cout<<"Index:"<<i<<std::endl;
cpv_joint(i, sampleString_10, outputname);
}

for(int i = 1; i < 10; i++){

cpv_joint(i, sampleString_7, outputname);
}


}
