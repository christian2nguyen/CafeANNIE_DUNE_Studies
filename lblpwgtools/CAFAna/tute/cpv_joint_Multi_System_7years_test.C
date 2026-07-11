#include "CAFAna/Analysis/common_fit_definitions.h"
#include "CAFAna/Analysis/XSecSystList.h"


#include "TPad.h"
#include <string>
#include "TAxis.h"
#include "TDirectory.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTree.h"
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
#include <cmath>
#include <set>

using namespace ana;

std::string GetAsimovName(int asimov_set, int hie)
{   if (asimov_set == 0 && hie == +1) return "kNuFitCV";  
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

namespace
{
  std::set<std::string> GetOscParamNameSet(const std::vector<const IFitVar*>& oscVars)
  {
    std::set<std::string> names;
    for(const auto* var: oscVars){
      if(var) names.insert(var->ShortName());
    }
    return names;
  }

  bool IsXSecParamName(const std::string& name)
  {
    static const std::set<std::string> xsecNames = [](){
      std::set<std::string> names;
      for(const auto& systName: ana::GetAllXSecSystNames()){
        names.insert(systName);
      }
      return names;
    }();

    return xsecNames.count(name);
  }

  double GetVecVal(const std::vector<double>* values, size_t idx)
  {
    if(!values) return 0.0;
    if(idx >= values->size()) return 0.0;
    return values->at(idx);
  }

  double SafeDivide(double numerator, double denominator)
  {
    if(!std::isfinite(numerator) || !std::isfinite(denominator)) return 0.0;
    if(std::abs(denominator) < 1e-12) return 0.0;
    return numerator/denominator;
  }

  void WriteFitContext(TDirectory* outDir,
                       int true_hie,
                       int test_hie,
                       double true_dcp,
                       double test_dcp,
                       double chisq)
  {
    if(!outDir) return;

    outDir->cd();

    TTree fitContext("fit_context", "fit_context");
    int trueHierarchy = true_hie;
    int testHierarchy = test_hie;
    double trueDcp = true_dcp;
    double testDcp = test_dcp;
    double chiSq = chisq;

    fitContext.Branch("true_hie", &trueHierarchy);
    fitContext.Branch("test_hie", &testHierarchy);
    fitContext.Branch("true_dcp", &trueDcp);
    fitContext.Branch("test_dcp", &testDcp);
    fitContext.Branch("chisq", &chiSq);
    fitContext.Fill();
    fitContext.Write();
  }

  void WriteParameterSummaryTree(TDirectory* outDir,
                                 const FitTreeBlob& fitBlob,
                                 const std::vector<const IFitVar*>& oscVars)
  {
    if(!outDir || !fitBlob.fParamNames) return;

    const auto oscNames = GetOscParamNameSet(oscVars);
    const auto& names = *fitBlob.fParamNames;

    outDir->cd();

    TTree summary("fit_parameter_summary", "fit_parameter_summary");
    int index = -1;
    int isXSec = 0;
    int isOsc = 0;
    int isOtherNuisance = 0;
    std::string name;
    double prefitValue = 0.0;
    double prefitError = 0.0;
    double postfitValue = 0.0;
    double postfitError = 0.0;
    double centralValue = 0.0;
    double fakeValue = 0.0;
    double postfitMinusFake = 0.0;
    double postfitMinusCentral = 0.0;
    double pullFromCentralPrefit = 0.0;
    double residualToFakePrefit = 0.0;

    summary.Branch("index", &index);
    summary.Branch("name", &name);
    summary.Branch("is_xsec", &isXSec);
    summary.Branch("is_osc", &isOsc);
    summary.Branch("is_other_nuisance", &isOtherNuisance);
    summary.Branch("prefit_value", &prefitValue);
    summary.Branch("prefit_error", &prefitError);
    summary.Branch("postfit_value", &postfitValue);
    summary.Branch("postfit_error", &postfitError);
    summary.Branch("central_value", &centralValue);
    summary.Branch("fake_value", &fakeValue);
    summary.Branch("postfit_minus_fake", &postfitMinusFake);
    summary.Branch("postfit_minus_central", &postfitMinusCentral);
    summary.Branch("pull_from_central_prefit", &pullFromCentralPrefit);
    summary.Branch("residual_to_fake_prefit", &residualToFakePrefit);

    for(size_t idx = 0; idx < names.size(); ++idx){
      name = names[idx];
      index = idx;
      isXSec = IsXSecParamName(name) ? 1 : 0;
      isOsc = oscNames.count(name) ? 1 : 0;
      isOtherNuisance = (!isXSec && !isOsc) ? 1 : 0;
      prefitValue = GetVecVal(fitBlob.fPreFitValues, idx);
      prefitError = GetVecVal(fitBlob.fPreFitErrors, idx);
      postfitValue = GetVecVal(fitBlob.fPostFitValues, idx);
      postfitError = GetVecVal(fitBlob.fPostFitErrors, idx);
      centralValue = GetVecVal(fitBlob.fCentralValues, idx);
      fakeValue = GetVecVal(fitBlob.fFakeDataVals, idx);
      postfitMinusFake = postfitValue - fakeValue;
      postfitMinusCentral = postfitValue - centralValue;
      pullFromCentralPrefit = SafeDivide(postfitMinusCentral, prefitError);
      residualToFakePrefit = SafeDivide(postfitMinusFake, prefitError);
      summary.Fill();
    }

    summary.Write();
  }

  void WriteSubsetSummary(TDirectory* outDir,
                          const FitTreeBlob& fitBlob,
                          const std::vector<int>& indices,
                          const std::string& prefix)
  {
    if(!outDir || !fitBlob.fParamNames || indices.empty()) return;

    const auto& names = *fitBlob.fParamNames;
    const int nbins = indices.size();

    outDir->cd();

    TH1D postfitValue((prefix + "_postfit_value").c_str(),
                      (prefix + " postfit value").c_str(), nbins, 0, nbins);
    TH1D postfitError((prefix + "_postfit_error").c_str(),
                      (prefix + " postfit error").c_str(), nbins, 0, nbins);
    TH1D prefitError((prefix + "_prefit_error").c_str(),
                     (prefix + " prefit error").c_str(), nbins, 0, nbins);
    TH1D centralValue((prefix + "_central_value").c_str(),
                      (prefix + " central value").c_str(), nbins, 0, nbins);
    TH1D fakeValue((prefix + "_fake_value").c_str(),
                   (prefix + " fake value").c_str(), nbins, 0, nbins);
    TH1D postfitMinusFake((prefix + "_postfit_minus_fake").c_str(),
                          (prefix + " postfit minus fake").c_str(), nbins, 0, nbins);
    TH1D postfitMinusCentral((prefix + "_postfit_minus_central").c_str(),
                             (prefix + " postfit minus central").c_str(), nbins, 0, nbins);
    TH1D pullFromCentralPrefit((prefix + "_pull_from_central_prefit").c_str(),
                               (prefix + " pull from central / prefit error").c_str(),
                               nbins, 0, nbins);
    TH1D residualToFakePrefit((prefix + "_residual_to_fake_prefit").c_str(),
                              (prefix + " residual to fake / prefit error").c_str(),
                              nbins, 0, nbins);

    for(int ibin = 0; ibin < nbins; ++ibin){
      const int idx = indices[ibin];
      const char* label = names[idx].c_str();
      const double postfit = GetVecVal(fitBlob.fPostFitValues, idx);
      const double prefitErr = GetVecVal(fitBlob.fPreFitErrors, idx);
      const double postfitErr = GetVecVal(fitBlob.fPostFitErrors, idx);
      const double central = GetVecVal(fitBlob.fCentralValues, idx);
      const double fake = GetVecVal(fitBlob.fFakeDataVals, idx);

      postfitValue.GetXaxis()->SetBinLabel(ibin + 1, label);
      postfitError.GetXaxis()->SetBinLabel(ibin + 1, label);
      prefitError.GetXaxis()->SetBinLabel(ibin + 1, label);
      centralValue.GetXaxis()->SetBinLabel(ibin + 1, label);
      fakeValue.GetXaxis()->SetBinLabel(ibin + 1, label);
      postfitMinusFake.GetXaxis()->SetBinLabel(ibin + 1, label);
      postfitMinusCentral.GetXaxis()->SetBinLabel(ibin + 1, label);
      pullFromCentralPrefit.GetXaxis()->SetBinLabel(ibin + 1, label);
      residualToFakePrefit.GetXaxis()->SetBinLabel(ibin + 1, label);

      postfitValue.SetBinContent(ibin + 1, postfit);
      postfitError.SetBinContent(ibin + 1, postfitErr);
      prefitError.SetBinContent(ibin + 1, prefitErr);
      centralValue.SetBinContent(ibin + 1, central);
      fakeValue.SetBinContent(ibin + 1, fake);
      postfitMinusFake.SetBinContent(ibin + 1, postfit - fake);
      postfitMinusCentral.SetBinContent(ibin + 1, postfit - central);
      pullFromCentralPrefit.SetBinContent(ibin + 1, SafeDivide(postfit - central, prefitErr));
      residualToFakePrefit.SetBinContent(ibin + 1, SafeDivide(postfit - fake, prefitErr));
    }

    postfitValue.Write();
    postfitError.Write();
    prefitError.Write();
    centralValue.Write();
    fakeValue.Write();
    postfitMinusFake.Write();
    postfitMinusCentral.Write();
    pullFromCentralPrefit.Write();
    residualToFakePrefit.Write();
  }

  void WriteFitSubsetSummaries(TDirectory* outDir,
                               const FitTreeBlob& fitBlob,
                               const std::vector<const IFitVar*>& oscVars)
  {
    if(!outDir || !fitBlob.fParamNames) return;

    const auto oscNames = GetOscParamNameSet(oscVars);
    const auto& names = *fitBlob.fParamNames;

    std::vector<int> xsecIndices;
    std::vector<int> otherNuisanceIndices;
    std::vector<int> oscIndices;

    for(size_t idx = 0; idx < names.size(); ++idx){
      if(IsXSecParamName(names[idx])){
        xsecIndices.push_back(idx);
      }
      else if(oscNames.count(names[idx])){
        oscIndices.push_back(idx);
      }
      else{
        otherNuisanceIndices.push_back(idx);
      }
    }

    WriteSubsetSummary(outDir, fitBlob, xsecIndices, "xsec");
    WriteSubsetSummary(outDir, fitBlob, otherNuisanceIndices, "other_nuisance");
    WriteSubsetSummary(outDir, fitBlob, oscIndices, "osc");
  }

  void FinalizeFitDiagnostics(TDirectory* outDir,
                              FitTreeBlob& fitBlob,
                              const std::vector<const IFitVar*>& oscVars,
                              int true_hie,
                              int test_hie,
                              double true_dcp,
                              double test_dcp,
                              double chisq)
  {
    if(!outDir) return;

    fitBlob.SetDirectory(outDir);
    fitBlob.Write();
    WriteFitContext(outDir, true_hie, test_hie, true_dcp, test_dcp, chisq);
    WriteParameterSummaryTree(outDir, fitBlob, oscVars);
    WriteFitSubsetSummaries(outDir, fitBlob, oscVars);
  }
}




void cpv_joint(int hie_input=9, std::string sampleString="ndfd:10year", 
std::string outputFname_input="State_Feb11_",
std::string stateFname  = "/pnfs/dune/persistent/users/cnguyen/cafana_example/State"/*="/pnfs/dune/persistent/users/cnguyen/cafana_example/toy_2perEnergy_test_FD_FHC_v1.root"*//*"common_state_mcc11v3.root"*/,
	       std::string systSet = "det:xsec:nflux=30"/*allsyst,nosyst, det:xsec:nflux=20*/, 
	       std::string penaltyString="nopen",  std::string asimov_set="0",
	       int isetCP = 0,
	       std::string fakeDataShift = "", int fitBias = 0){
	       
	int NormMass_hie = 1; 
	int invertedNormMass_hie = -1; 
 asimov_set = std::to_string(hie_input);
 int asimov_int = hie_input;
std::string paraset_name = GetAsimovName(hie_input, 1);
std::string outputFname =  stateFname + "Feb10_OUTPUT_"+ paraset_name + "_624ktmwyr_" + sampleString + "_" + systSet  + "_" + penaltyString +"_asimov_set_" + asimov_set +  ".root";
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

std::vector<std::string> flux_universe_remove{
"flux_Nov17_11",
"flux_Nov17_12",
"flux_Nov17_13",
"flux_Nov17_14",
"flux_Nov17_15",
"flux_Nov17_16",
"flux_Nov17_17",
"flux_Nov17_18",
"flux_Nov17_19",
"flux_Nov17_20",
"flux_Nov17_21",
"flux_Nov17_22",
"flux_Nov17_23",
"flux_Nov17_24",
"flux_Nov17_25",
"flux_Nov17_26",
"flux_Nov17_27",
"flux_Nov17_28",
"flux_Nov17_29"

};

  gROOT->SetBatch(1);
  gRandom->SetSeed(0);

  // Allow a fake data bias
  SystShifts trueSyst = GetFakeDataSystShift(fakeDataShift);

  // Get the systematics to use
  std::vector<const ISyst*> systlist = GetListOfSysts(systSet);
                                          
                                          
  std::map<std::string,std::vector<const ISyst *>> ISYSYT_Map;

                                          
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

  // removing all flux universes only keeping 10 
  //RemoveSysts(systlist, flux_universe_remove);




  TFile* fout = new TFile(outputFname.c_str(), "RECREATE");
  fout->cd();

  // First find the minimum for dcp = 0
  std::map<const IFitVar*, std::vector<double>> oscSeeds = {};
  oscSeeds[&kFitSinSqTheta23] = {0.4, 0.6};

  osc::IOscCalculatorAdjustable* trueOscGlob = NuFitOscCalc(NormMass_hie, 1, asimov_set);
  trueOscGlob->SetdCP(0.0);
  //std::map<std::string,double>glob_chisqmin_sys_map;
  
  double glob_chisqmin = 99999;
  double thischisq;

  // Now loop over all true values
  //
 // 
  //
  // After the main loop, before fout->Close()
  //std::map<std::string, TGraph*> gCPV_per_syst; // want to save indivdiual 
  std::cout<<"~~~~~~~~~~~~~Trying to loop throught univeres ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<< std::endl;
  
//////////////////////////////////////////////////////////
// Have to loop because I need to change the testOsc dcp value
  for (int idcp = 0; idcp < 2; ++idcp) {  
    for(int ihie = -1; ihie <= +1; ihie += 2) {
      double dcptest = idcp*TMath::Pi();
      //double dcptest = -TMath::Pi() + isetCP*dcpstep;

      std::vector<const IFitVar*> oscVars = GetOscVars("th23:th13:dmsq32", ihie);      
      osc::IOscCalculatorAdjustable* testOsc = NuFitOscCalc(ihie, 1, asimov_set);
      testOsc->SetdCP(dcptest);
      
      IExperiment *penalty = GetPenalty(ihie, 1, penaltyString, asimov_set);
      SystShifts testSyst = kNoShift;
      std::string fitDirName = "global_prefit_seed_testhie_" + std::to_string(ihie)
                             + "_testdcp_index_" + std::to_string(idcp);
      TDirectory* fitDir = fout->mkdir(fitDirName.c_str());
      FitTreeBlob fitBlob("fit_info", "meta_tree");
      
      thischisq = RunFitPoint(stateFname, sampleString,
			      trueOscGlob, trueSyst, false,
			      oscVars, systlist,
			      testOsc, testSyst,
			      oscSeeds, penalty, 
			      Fitter::kNormal, fitDir, &fitBlob, nullptr, ana::junkShifts);

      FinalizeFitDiagnostics(fitDir, fitBlob, oscVars,
                             NormMass_hie, ihie, 0.0, dcptest, thischisq);
      
      glob_chisqmin = TMath::Min(thischisq,glob_chisqmin);
      
      delete penalty;
      delete testOsc;
    }
  }
  
  delete trueOscGlob;
  std::cout << "Global chi2 at dCP=0,pi: " << glob_chisqmin << std::endl;
  
  int nsteps = 36;
  double dcpstep = 2*TMath::Pi()/nsteps;
  TGraph* gCPV = new TGraph();

////////////////////////////////////////////////////////////

  for(int idcp = 0; idcp < nsteps+1; ++idcp) {
    std::cout << "Trying idcp = " << idcp << std::endl;
    
    double thisdcp = -TMath::Pi() + idcp*dcpstep;
    osc::IOscCalculatorAdjustable* trueOsc = NuFitOscCalc(NormMass_hie, 1, asimov_set);
    trueOsc->SetdCP(thisdcp);

    double chisqmin = 99999.0;
    double bestTestDcp = 0.0;
    double bestChiSq = 99999.0;
    int bestTestHie = NormMass_hie;
    //std::map<std::string,double > chisqmin_map;

    // Still need to loop over dcp choices
    for (int itestdcp = 0; itestdcp < 2; ++itestdcp) {

      double dcptest = itestdcp*TMath::Pi();
      //double dcptest = -TMath::Pi() + isetCP*dcpstep;
      //std::cout<< "idcp (index) = "<<idcp << " ---- the true and test dcp values: "<<thisdcp<<" "<<dcptest<<std::endl;

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
        if(thischisq < bestChiSq){
          bestChiSq = thischisq;
          bestTestDcp = dcptest;
          bestTestHie = ihie;
        }
	
	delete penalty;
	delete testOsc;	
      }
    }// end of idcp inner 

    {
      std::vector<const IFitVar*> bestOscVars = GetOscVars("th23:th13:dmsq32", bestTestHie);
      osc::IOscCalculatorAdjustable* bestTestOsc = NuFitOscCalc(bestTestHie, 1, asimov_set);
      bestTestOsc->SetdCP(bestTestDcp);
      IExperiment *bestPenalty = GetPenalty(bestTestHie, 1, penaltyString, asimov_set);
      SystShifts bestTestSyst = kNoShift;
      std::string fitDirName = "scan_true_dcp_index_" + std::to_string(idcp)
                             + "_bestfit";
      TDirectory* fitDir = fout->mkdir(fitDirName.c_str());
      FitTreeBlob fitBlob("fit_info", "meta_tree");

      bestChiSq = RunFitPoint(stateFname, sampleString,
                              trueOsc, trueSyst, false,
                              bestOscVars, systlist,
                              bestTestOsc, bestTestSyst,
                              oscSeeds, bestPenalty,
                              Fitter::kNormal, fitDir, &fitBlob, nullptr, ana::junkShifts);

      FinalizeFitDiagnostics(fitDir, fitBlob, bestOscVars,
                             NormMass_hie, bestTestHie, thisdcp, bestTestDcp, bestChiSq);

      chisqmin = TMath::Min(chisqmin, bestChiSq);

      delete bestPenalty;
      delete bestTestOsc;
    }
    
    
    delete trueOsc;
    std::cout<<" ---- current chi2 "<<chisqmin<<std::endl;   
 
    chisqmin = TMath::Max(chisqmin,1e-6);
    double diff = chisqmin-glob_chisqmin;
    diff = TMath::Max(diff, 0.);
    gCPV->SetPoint(gCPV->GetN(),thisdcp/TMath::Pi(),sqrt(diff));    
    
  }
  /////////////// end of idcp ouuter 

  fout->cd();
  gCPV->Draw("ALP");
  //c1 -> Print(pdf_title);
  gCPV->Write(hie_input > 0 ? "sens_cpv_nh" : "sens_cpv_ih");
    
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

void cpv_joint_Multi_System_7years_test(){
std::string outputname = "/pnfs/dune/persistent/users/cnguyen/cafana_example/State_withsystem_only10flux_";
std::string sampleString_10="ndfd:10year";
std::string sampleString_7="ndfd:7year";
std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;


std::cout<<"Starting Loop 10 YEAR  "<< std::endl;
//for(int i = 0; i < 10; i++){
//std::cout<<"Index:"<<i<<std::endl;
//cpv_joint(i, sampleString_10, outputname);
//}

std::vector<int> Makeuplist{9};

std::cout<<"Starting Loop 7 YEAR  "<< std::endl;
//for(int i = 0; i < 10; i++){
//cpv_joint(i, sampleString_7, outputname);
//}

for(int i = 0; i < Makeuplist.size(); i++){
cpv_joint(Makeuplist.at(i), sampleString_7, outputname);
}


//cpv_joint(7, sampleString_7, outputname);

}
