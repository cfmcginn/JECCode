#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TDatime.h"
#include "TLorentzVector.h"

#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <fstream>

#include "include/getLogBins.h"
#include "include/getLinBins.h"
#include "include/etaPhiFunc.h"
#include "include/getBkgEstimate.h"

#include "include/returnFileList.h"
#include "include/returnRootFileContentsList.h"
#include "include/getResidualJetCorr.h"

const std::string storeStr = "/store";
const std::string xrootdStr = "root://xrootd.unl.edu/";


const Bool_t debugMode = false;

const Bool_t doGetBkg = false;

const Bool_t deriveResidualMode = false;
const Bool_t applyResidualMode = false;
const Bool_t applyResidualModeFull = false;

const Int_t nJtCat = 3;
const std::string jtCat[nJtCat] = {"Eta2", "Eta1", "Dijet"};

const Int_t nMeanFit = 2;
const std::string meanFit[nMeanFit] = {"", "Fit"};

const Int_t nQG = 3;
const std::string qg[nQG] = {"Inc", "Q", "G"};

const Int_t nMaxJets = 500;

const Int_t nCentBins = 2;
const Int_t centBins[nCentBins+1] = {200, 60, 0};
const Float_t centBins2[nCentBins+1] = {0.001, 30, 99.999};
//const Int_t nCentBins = 4;
//const Int_t centBins[nCentBins+1] = {200, 100, 60, 20, 0};
//const Float_t centBins2[nCentBins+1] = {0.001, 10, 30, 50, 99.999};
//const Int_t nCentBins = 8;
//const Int_t centBins[nCentBins+1] = {200, 140, 100, 80, 60, 40, 20, 10, 0};

const Int_t nMaxGen = 100000;

const Int_t nJetRecoCuts = 6;
const Int_t jetRecoCuts[nJetRecoCuts] = {5, 10, 15, 20, 25, 30};

const Float_t eleMass = 0.000510998;
const Float_t muMass = .105658371;


float findNcoll(int hiBin) {
  const int nbins = 200;
  const float Ncoll[nbins] = {1976.95, 1944.02, 1927.29, 1891.9, 1845.3, 1807.2, 1760.45, 1729.18, 1674.8, 1630.3, 1590.52, 1561.72, 1516.1, 1486.5, 1444.68, 1410.88, 1376.4, 1347.32, 1309.71, 1279.98, 1255.31, 1219.89, 1195.13, 1165.96, 1138.92, 1113.37, 1082.26, 1062.42, 1030.6, 1009.96, 980.229, 955.443, 936.501, 915.97, 892.063, 871.289, 847.364, 825.127, 806.584, 789.163, 765.42, 751.187, 733.001, 708.31, 690.972, 677.711, 660.682, 640.431, 623.839, 607.456, 593.307, 576.364, 560.967, 548.909, 530.475, 519.575, 505.105, 490.027, 478.133, 462.372, 451.115, 442.642, 425.76, 416.364, 405.154, 392.688, 380.565, 371.167, 360.28, 348.239, 340.587, 328.746, 320.268, 311.752, 300.742, 292.172, 281.361, 274.249, 267.025, 258.625, 249.931, 240.497, 235.423, 228.63, 219.854, 214.004, 205.425, 199.114, 193.618, 185.644, 180.923, 174.289, 169.641, 161.016, 157.398, 152.151, 147.425, 140.933, 135.924, 132.365, 127.017, 122.127, 117.817, 113.076, 109.055, 105.16, 101.323, 98.098, 95.0548, 90.729, 87.6495, 84.0899, 80.2237, 77.2201, 74.8848, 71.3554, 68.7745, 65.9911, 63.4136, 61.3859, 58.1903, 56.4155, 53.8486, 52.0196, 49.2921, 47.0735, 45.4345, 43.8434, 41.7181, 39.8988, 38.2262, 36.4435, 34.8984, 33.4664, 31.8056, 30.351, 29.2074, 27.6924, 26.7754, 25.4965, 24.2802, 22.9651, 22.0059, 21.0915, 19.9129, 19.1041, 18.1487, 17.3218, 16.5957, 15.5323, 14.8035, 14.2514, 13.3782, 12.8667, 12.2891, 11.61, 11.0026, 10.3747, 9.90294, 9.42648, 8.85324, 8.50121, 7.89834, 7.65197, 7.22768, 6.7755, 6.34855, 5.98336, 5.76555, 5.38056, 5.11024, 4.7748, 4.59117, 4.23247, 4.00814, 3.79607, 3.68702, 3.3767, 3.16309, 2.98282, 2.8095, 2.65875, 2.50561, 2.32516, 2.16357, 2.03235, 1.84061, 1.72628, 1.62305, 1.48916, 1.38784, 1.28366, 1.24693, 1.18552, 1.16085, 1.12596, 1.09298, 1.07402, 1.06105, 1.02954};
  return Ncoll[hiBin];
}


Bool_t isGoodCaloJet(TString algo, Float_t hcalSum, Float_t ecalSum)
{
  if(algo.Index("Calo") >= 0){
    Float_t totSum = hcalSum+ecalSum;
    if(hcalSum/totSum < .1) return false;
    else if(ecalSum/totSum < .05) return false;
    else return true;
  }
  else return false;
}

Bool_t isGoodPFJet(TString algo, Int_t chargedN, Float_t chargedSum, Float_t neutralSum, Float_t rawPt)
{
  if(algo.Index("PF") >= 0){
    if(chargedN <= 0) return false;
    //    else if(chargedSum/rawPt >= .98) return false;
    //    else if(neutralSum/rawPt >= .98) return false;
    else return true;

    //following lines just to get around unused variable temporrarily
    if(chargedSum/rawPt >= .98) return false;
    if(neutralSum/rawPt >= .98) return false;
  }
  else return false;
}


void genSort(Int_t nGenJt, Float_t genJtPt[], Float_t genJtPhi[], Float_t genJtEta[])
{
  for(Int_t iter = 0; iter < nGenJt-1; iter++){
    for(Int_t iter2 = iter+1; iter2 < nGenJt; iter2++){
      if(genJtPt[iter] < genJtPt[iter2]){
	Float_t tempPt = genJtPt[iter];
	Float_t tempPhi = genJtPhi[iter];
	Float_t tempEta = genJtEta[iter];

	genJtPt[iter] = genJtPt[iter2];
	genJtPhi[iter] = genJtPhi[iter2];
	genJtEta[iter] = genJtEta[iter2];

	genJtPt[iter2] = tempPt;
	genJtPhi[iter2] = tempPhi;
	genJtEta[iter2] = tempEta;
      }
      else continue;
    }
  }

  return;
}


void FitGauss(TH1F* hist_p, Bool_t isPbPb, Float_t& mean, Float_t& meanErr, Float_t& res, Float_t& resErr)
{
  if(hist_p->Integral() == 0) return;
  if(hist_p->GetEntries() == 0) return;

  TF1* f1_p = new TF1("f1_p", "gaus", hist_p->GetXaxis()->GetXmin(), hist_p->GetXaxis()->GetXmax());

  hist_p->Fit("f1_p", "Q M");

  mean = f1_p->GetParameter(1);
  res = f1_p->GetParameter(2);

  meanErr = f1_p->GetParError(1);
  resErr = f1_p->GetParError(2);

  //  if(TMath::Abs(mean - 1.0) < 0.01) return;
  if(f1_p->GetProb() > .01) return;
  if(!isPbPb) return;

  for(Int_t fitIter = 0; fitIter < 1; fitIter++){
    Float_t max = -1;
    Int_t maxBin = -1;

    for(Int_t iter = 0; iter < hist_p->GetNbinsX(); iter++){
      if(hist_p->GetBinContent(iter+1) > max){
	max = hist_p->GetBinContent(iter+1);
	maxBin = iter+1;
      }
    }
    
    Int_t nBins = 1;

    Float_t fitLow = -1;
    Float_t fitHi = -1;
    
    for(Int_t iter = 1; iter < hist_p->GetNbinsX(); iter++){
      Int_t tempBin = maxBin - iter;
      if(tempBin < 1) tempBin = 1;
      
      nBins += 2;
      
      if(hist_p->Integral(tempBin, maxBin+iter) > (.95-.05*fitIter)*hist_p->Integral() || tempBin == 1){
	fitLow = hist_p->GetBinCenter(tempBin);
	fitHi = hist_p->GetBinCenter(maxBin+iter);
	break;
      }
    }
    
    if(nBins >= 7){
      hist_p->Fit("f1_p", "Q M", "", fitLow, fitHi);
      
      mean = f1_p->GetParameter(1);
      res = f1_p->GetParameter(2);
      
      meanErr = f1_p->GetParError(1);
      resErr = f1_p->GetParError(2);
      
      //  if(TMath::Abs(mean - 1.0) < 0.01) return;
      if(f1_p->GetProb() > .01) return;
      return;
    }
    nBins = 1;
    
    Int_t meanBin = hist_p->FindBin(hist_p->GetMean());
    fitLow = -1;
    fitHi = -1;
    
    for(Int_t iter = 1; iter < hist_p->GetNbinsX(); iter++){
      Int_t tempBin = meanBin - iter;
      if(tempBin < 1) tempBin = 1;
      
      nBins += 2;
      
      if(hist_p->Integral(tempBin, meanBin+iter) > (.95-.05*fitIter)*hist_p->Integral() || tempBin == 1){
	fitLow = hist_p->GetBinCenter(tempBin);
	fitHi = hist_p->GetBinCenter(meanBin+iter);
	break;
      }
    }
    
    
    if(nBins >= 7){
      hist_p->Fit("f1_p", "Q M", "", fitLow, fitHi);
      
      mean = f1_p->GetParameter(1);
      res = f1_p->GetParameter(2);
      
      meanErr = f1_p->GetParError(1);
      resErr = f1_p->GetParError(2);
    }
  }

  delete f1_p;

  return;
}


void FitCSN(TH1F* hist_p, Bool_t isPbPb, const Int_t jtAlgoR = -1, const std::string inPPFileName = "")
{
  if(hist_p->GetEntries() == 0) return;
  if(hist_p->Integral() == 0) return;

  TF1* f1_p;

  if(isPbPb){
    if(jtAlgoR == -1 || !strcmp(inPPFileName.c_str(), "")) return;

    TFile* ppFile_p = new TFile(inPPFileName.c_str(), "READ");
    TH1F* ppHist_p = (TH1F*)ppFile_p->Get(Form("ak%dPF/jtRecoOverGenVPt_Inc_FitRes_ak%dPF_PP_h", jtAlgoR, jtAlgoR));
    TF1* f1PP_p = (TF1*)ppHist_p->GetFunction("f1PP_p");

    Float_t par0 = f1PP_p->GetParameter(0);
    Float_t par1 = f1PP_p->GetParameter(1);

    //    par0 = .0246;
    //    par1 = 1.213;

    Float_t parErr0 = f1PP_p->GetParError(0);
    Float_t parErr1 = f1PP_p->GetParError(1);

    ppFile_p->Close();
    delete ppFile_p;

    f1_p = new TF1("f1_p", Form("TMath::Sqrt([0]*[0] + [1]*[1]/(x) + [2]*[2]/(x*x))"), 30/*hist_p->GetXaxis()->GetXmin()*/, hist_p->GetXaxis()->GetXmax());
    f1_p->SetParameter(0, par0); 
    std::cout << "par0, par1: " << par0 << ", " << par1 << std::endl;
    f1_p->SetParLimits(0, par0-parErr0, par0+parErr0); 
    f1_p->SetParameter(1, par1); 
    f1_p->SetParLimits(1, par1-parErr1, par1+parErr1); 
    f1_p->SetParameter(2, .1);
    //    f1_p->SetParLimits(0, 0, 1000000000);
  }
  else{
    f1_p = new TF1("f1PP_p", "TMath::Sqrt([0]*[0] + [1]*[1]/(x) + [2]*[2]/(x*x))", 30/*hist_p->GetXaxis()->GetXmin()*/, hist_p->GetXaxis()->GetXmax());
    f1_p->SetParameter(0, .03);
    f1_p->SetParameter(1, 1.2);
    f1_p->SetParameter(2, .001);
    //    f1_p->SetParLimits(0, 0, 1000000000);
    //  f1_p->SetParLimits(1, 0, 1000000000);
    f1_p->SetParLimits(2, -.1, .1);
  }

  if(isPbPb) hist_p->Fit("f1_p", "Q M", "", 30, hist_p->GetXaxis()->GetXmax());
  else hist_p->Fit("f1PP_p", "Q M", "", 30, hist_p->GetXaxis()->GetXmax());

  delete f1_p;

  return;
}


void FitCSN2(TH1F* hist_p, Bool_t isPerph, TH1F* perphHist_p)
{
  if(hist_p->GetEntries() == 0) return;
  if(hist_p->Integral() == 0) return;

  TF1* f1_p;

  if(!isPerph){
    //    perphHist_p->Print("ALL");
    TF1* f1Perph_p = (TF1*)perphHist_p->GetFunction("f1Perph_p");

    Float_t par0 = f1Perph_p->GetParameter(0);
    Float_t par1 = f1Perph_p->GetParameter(1);

    //    par0 = .0246;
    //    par1 = 1.213;

    Float_t parErr0 = f1Perph_p->GetParError(0);
    Float_t parErr1 = f1Perph_p->GetParError(1);

    f1_p = new TF1("f1_p", Form("TMath::Sqrt([0]*[0] + [1]*[1]/(x) + [2]*[2]/(x*x))"), 30/*hist_p->GetXaxis()->GetXmin()*/, 200/*hist_p->GetXaxis()->GetXmax()*/);
    f1_p->SetParameter(0, par0); 
    std::cout << "par0, par1: " << par0 << ", " << par1 << std::endl;
    f1_p->SetParLimits(0, par0-parErr0, par0+parErr0); 
    f1_p->SetParameter(1, par1); 
    f1_p->SetParLimits(1, par1-parErr1, par1+parErr1); 
    f1_p->SetParameter(2, .1);
    //    f1_p->SetParLimits(0, 0, 1000000000);
  }
  else{
    f1_p = new TF1("f1Perph_p", "TMath::Sqrt([0]*[0] + [1]*[1]/(x) + [2]*[2]/(x*x))", 30/*hist_p->GetXaxis()->GetXmin()*/, 300/*hist_p->GetXaxis()->GetXmax()*/);
    f1_p->SetParameter(0, .03);
    f1_p->SetParameter(1, 1.2);
    f1_p->SetParameter(2, .001);
    //    f1_p->SetParLimits(0, 0, 1000000000);
    //  f1_p->SetParLimits(1, 0, 1000000000);
    f1_p->SetParLimits(2, -.1, .1);
  }

  if(!isPerph) hist_p->Fit("f1_p", "Q M", "", 30, 200/*hist_p->GetXaxis()->GetXmax()*/);
  else hist_p->Fit("f1Perph_p", "Q M", "", 30, 300/*hist_p->GetXaxis()->GetXmax()*/);

  delete f1_p;

  return;
}


void FitCSNSimple(TH1F* hist_p)
{
  if(hist_p->GetEntries() == 0) return;
  if(hist_p->Integral() == 0) return;

  TF1* f1_p = new TF1("f1_p", "TMath::Sqrt([0]*[0] + [1]*[1]/(x) + [2]*[2]/(x*x))", 30/*hist_p->GetXaxis()->GetXmin()*/, hist_p->GetXaxis()->GetXmax());
  f1_p->SetParameter(0, .05);
  f1_p->SetParameter(1, 1);
  f1_p->SetParameter(2, .05);

  hist_p->Fit("f1_p", "Q M", "", 30, hist_p->GetXaxis()->GetXmax());

  delete f1_p;

  return;
}


void FitErf(TH1F* hist_p)
{
  if(hist_p->GetEntries() == 0) return;
  if(hist_p->Integral() == 0) return;

  TF1* f1_p = new TF1("f1_p", ".5*[0]*(1+TMath::Erf([1]+[2]*x))", hist_p->GetXaxis()->GetXmin(), hist_p->GetXaxis()->GetXmax());
  f1_p->SetParameter(0, .99);
  f1_p->SetParameter(1, -.17);
  f1_p->SetParameter(2, .06);

  hist_p->Fit("f1_p", "Q M", "", hist_p->GetXaxis()->GetXmin(), hist_p->GetXaxis()->GetXmax());

  delete f1_p;

  return;
}


void makeJECHist(const std::string inFileName15, const std::string inFileName30, const std::string inFileName50, const std::string inFileName80, const Bool_t isPbPb, const std::string inPPFileName = "", const std::string hibinWeightFileName = "")
{
  /*  
  if(deriveResidualMode && applyResidualMode){
    std::cout << "Both derive and applyResidual mode active. Please choose one. Return." << std::endl;
    return;
  }
  */

  if(debugMode) std::cout << __LINE__ << std::endl;
  if(applyResidualMode){
    if(debugMode) std::cout << __LINE__ << std::endl;

    initGetResidualJetCorr("outputDir/merged_dgulhan-Pythia8_Dijet30_pp_TuneCUETP8M1_Hydjet_MinBias_5020GeV_RECODEBUG_758_PrivMC_forest_v28_0_20160512_pthat_30_RESIDUALCORR.root");

    //    initGetResidualJetCorr("outputDir/merged_dgulhan-Pythia8_Dijet50_pp_TuneCUETP8M1_Hydjet_MinBias_5020GeV_RECODEBUG_758_PrivMC_forest_v28_0_20160506_pthat_50_RESIDUALCORR_TEMPPREV.root");

    //initGetResidualJetCorr("merged_dgulhan-Pythia8_Dijet30_pp_TuneCUETP8M1_Hydjet_MinBias_5020GeV_RECODEBUG_758_PrivMC_forest_v28_0_20160507_pthat_30_RESIDUALCORR.root");

    if(debugMode) std::cout << __LINE__ << std::endl;
    if(debugMode) std::cout << __LINE__ << std::endl;
  }

  const Int_t nFiles = 4;
  std::string inFileNames[nFiles] = {inFileName15, inFileName30, inFileName50, inFileName80};
  const std::string inFilePtHats[nFiles] = {"pthat15", "pthat30", "pthat50", "pthat80"};
  const std::string inFilePtHats2[nFiles] = {"15", "30", "50", "80"};
  const Int_t pthatVals[nFiles+1] = {15, 30, 50, 80, 9999999};
  Float_t pthatVals2[nFiles+1] = {15, 30, 50, 80, 9999999};
  const Float_t pthatWeights[nFiles+1] = {.5335, .03378, .003778, .0004412, .00000001};
  Bool_t isFile[nFiles] = {false, false, false, false};

  for(Int_t iter = 0; iter < nFiles; iter++){
    if(strcmp(inFileNames[iter].c_str(), "") != 0) isFile[iter] = true;
  }


  if(debugMode) std::cout << __LINE__ << std::endl;

  std::vector<std::string> listOfFiles[nFiles];

  if(debugMode) std::cout << __LINE__ << std::endl;

  for(Int_t iter = 0; iter < nFiles; iter++){
    if(isFile[iter]){

      if(inFileNames[iter].substr(0, storeStr.size()).find(storeStr) != std::string::npos) inFileNames[iter] = xrootdStr + inFileNames[iter];

      if(inFileNames[iter].substr(inFileNames[iter].size()-5, 5).find(".root") != std::string::npos) listOfFiles[iter].push_back(inFileNames[iter]);
      else listOfFiles[iter] = returnFileList(inFileNames[iter], "HiForest", listOfFiles[iter].size());
    }
  }



  if(debugMode) std::cout << __LINE__ << std::endl;

  Bool_t isThereSomeFiles = false;

  for(Int_t iter = 0; iter < nFiles; iter++){
    if(listOfFiles[iter].size() != 0){
      isThereSomeFiles = true;
      break;
    }
  }

  if(debugMode) std::cout << __LINE__ << std::endl;


  if(!isThereSomeFiles){
    std::cout << "No input files given. Return." << std::endl;
    return;
  }

  if(debugMode) std::cout << __LINE__ << std::endl;


  for(Int_t iter = 0; iter < nFiles; iter++){
    int fileFailIter = 0;
    while(fileFailIter < (int)listOfFiles[iter].size()){
      if(listOfFiles[iter].at(fileFailIter).find("/failed/") != std::string::npos) listOfFiles[iter].erase(listOfFiles[iter].begin() + fileFailIter);
      else fileFailIter++;
    }
  }

  if(debugMode) std::cout << __LINE__ << std::endl;

  std::string outName;
  for(Int_t iter = 0; iter < nFiles; iter++){
    if(isFile[iter]){
      outName = listOfFiles[iter].at(0).c_str();
      break;
    }
  }

  std::string tempJetFileName = outName;

  if(debugMode) std::cout << __LINE__ << std::endl;
  const std::string inString = ".root";
  const std::string inString2 = "*";
  TDatime* date = new TDatime();
  std::string inHatString = "pthat";
  if(debugMode) std::cout << __LINE__ << std::endl;
  for(Int_t iter = 0; iter < nFiles; iter++){
    if(isFile[iter]) inHatString = inHatString + "_" + inFilePtHats2[iter];
  }

  std::string residualInsertStr = "";
  if(deriveResidualMode && applyResidualMode) residualInsertStr = "RESIDUAL2";
  else if(deriveResidualMode) residualInsertStr = "RESIDUAL";
  else if(applyResidualMode) residualInsertStr = "APPLYRES";

  const std::string outString = Form("_%d_%s_%sHIST.root", date->GetDate(), inHatString.c_str(), residualInsertStr.c_str());
  std::size_t strIndex = 0;
  if(debugMode) std::cout << __LINE__ << std::endl;
  std::string tempOutName = outName;
  Int_t outNum = 0;
  Int_t outIter = 0;
  Int_t totOutNum = 0;

  if(debugMode) std::cout << __LINE__ << std::endl;

  while(tempOutName.find("/") != std::string::npos){
    tempOutName.replace(0, tempOutName.find("/")+1, "");
    if(tempOutName.size() > 40) outNum++;
    totOutNum++;
  }

  while(totOutNum - outNum < 2) outNum--;

  if(debugMode) std::cout << __LINE__ << std::endl;

  while(outName.find("/") != std::string::npos){
    if(outIter < outNum){
      outIter++;
      outName.replace(0, outName.find("/")+1, "");
    }
    else outName.replace(outName.find("/"), 1, "_");
  }


  if(debugMode) std::cout << __LINE__ << std::endl;

  strIndex = outName.find(inString);
  if(!(strIndex == std::string::npos)){
    outName.replace(strIndex, inString.length(), outString); 
  }

  strIndex = outName.find(inString2);
  if(!(strIndex == std::string::npos)){
    outName.replace(strIndex, inString2.length(), outString);
  }

  TFile* outFile_p = new TFile(outName.c_str(), "UPDATE");
  
  TFile* tempJetInFile_p = TFile::Open(tempJetFileName.c_str(), "READ");
  const std::vector<std::string> jetAlgoInFile = returnRootFileContentsList(tempJetInFile_p, "TTree", "3PFJetAnalyzer");

  const Int_t nJetAlgo = (Int_t)jetAlgoInFile.size();
  std::cout << nJetAlgo << std::endl;
  std::vector<std::string> jetAlgo;
  Int_t jetAlgoR[nJetAlgo];

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    jetAlgo.push_back(jetAlgoInFile.at(iter).substr(0, jetAlgoInFile.at(iter).find("JetAnalyzer")));

    for(Int_t rIter = 0; rIter < 10; rIter++){
      if(jetAlgo.at(iter).find(std::to_string(rIter)) != std::string::npos){
	jetAlgoR[iter] = rIter;
	break;
      }
    }
  }
  

  if(debugMode) std::cout << __LINE__ << std::endl;


  if(debugMode) std::cout << __LINE__ << std::endl;

  std::ofstream myfile;
  myfile.open ("example.txt");

  
  if(debugMode) std::cout << __LINE__ << std::endl;
  
  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    std::cout << "Algos: " << jetAlgo.at(iter) << std::endl;
  }

  Int_t hiBin_ = -1;
  Float_t vz_;

  UInt_t run_;
  ULong64_t evt_;

  std::vector<float>* genPt_p = 0;
  std::vector<float>* genEta_p = 0;
  std::vector<float>* genPhi_p = 0;
  std::vector<int>* genPDG_p = 0;

  std::vector<double>* rho_p = 0;

  Int_t nJt_[nJetAlgo]; 
  Float_t jtPt_[nJetAlgo][nMaxJets];
  Float_t jtRawPt_[nJetAlgo][nMaxJets];
  Float_t jtEta_[nJetAlgo][nMaxJets];
  Float_t jtPhi_[nJetAlgo][nMaxJets];
  Float_t jtHCalSum_[nJetAlgo][nMaxJets];
  Float_t jtECalSum_[nJetAlgo][nMaxJets];
  Float_t jtNeutralSum_[nJetAlgo][nMaxJets];
  Float_t jtChargedSum_[nJetAlgo][nMaxJets];
  Float_t jtPhotonSum_[nJetAlgo][nMaxJets];
  Float_t jtESum_[nJetAlgo][nMaxJets];
  Int_t jtNeutralN_[nJetAlgo][nMaxJets];
  Int_t jtEN_[nJetAlgo][nMaxJets];
  Int_t jtMuN_[nJetAlgo][nMaxJets];
  Int_t jtPhotonN_[nJetAlgo][nMaxJets];
  Int_t jtPhotonHardN_[nJetAlgo][nMaxJets];
  Int_t jtTrackN_[nJetAlgo][nMaxJets];
  Float_t jtPU_[nJetAlgo][nMaxJets];
  Int_t jtChargedN_[nJetAlgo][nMaxJets];
  Float_t refPt_[nJetAlgo][nMaxJets];
  Float_t refEta_[nJetAlgo][nMaxJets];
  Float_t refPhi_[nJetAlgo][nMaxJets];
  Int_t refSubID_[nJetAlgo][nMaxJets];
  Int_t refPartFlav_[nJetAlgo][nMaxJets];
  Float_t ptHat_[nJetAlgo];

  Float_t jtTrackMax_[nJetAlgo][nMaxJets];
  Float_t jtNeutralMax_[nJetAlgo][nMaxJets];
  Float_t jtPhotonMax_[nJetAlgo][nMaxJets];
  Float_t jtChargedMax_[nJetAlgo][nMaxJets];

  Int_t nGenJt_[nJetAlgo]; 
  Float_t genJtPt_[nJetAlgo][nMaxJets];
  Float_t genJtEta_[nJetAlgo][nMaxJets];
  Float_t genJtPhi_[nJetAlgo][nMaxJets];
  Int_t genJtMatchIndex_[nJetAlgo][nMaxJets];
  Int_t genSubId_[nJetAlgo][nMaxJets];

  
  //  const Int_t nJtPtBins = 30;
  const Int_t nJtPtBins = 20;
  const Float_t jtPtLow = 28;
  const Float_t jtPtHi = 125;

  const Int_t nJtPtBinsPERP = 25;
  const Float_t jtPtLowPERP = 50;
  const Float_t jtPtHiPERP = 200;

  //const Float_t jtPtLow = 60;
  //const Float_t jtPtHi = 120;
  //const Float_t jtPtLow = 100;
  //const Float_t jtPtHi = 200;
  Double_t jtPtBins[nJtPtBins+1];
  getLogBins(jtPtLow, jtPtHi, nJtPtBins, jtPtBins);

  Double_t jtPtBinsPERP[nJtPtBinsPERP+1];
  getLogBins(jtPtLowPERP, jtPtHiPERP, nJtPtBinsPERP, jtPtBinsPERP);
  
  //DERIVE
  /*
  const Int_t nJtPtBins = 30;
  //const Int_t nJtPtBins = 15;
  const Float_t jtPtLow = 25;
  const Float_t jtPtHi = 300;
  Double_t jtPtBins[nJtPtBins+1];
  getLogBins(jtPtLow, jtPtHi, nJtPtBins, jtPtBins);
  */

  for(Int_t pthatIter = 1; pthatIter < nFiles; pthatIter++){
    for(Int_t iter = 0; iter < nJtPtBins+1; iter++){
      if(jtPtBins[iter] > pthatVals2[pthatIter]){
	pthatVals2[pthatIter] = jtPtBins[iter-1];
	break;
      }
    }
  }

  for(Int_t pthatIter = 0; pthatIter < nFiles+1; pthatIter++){
    std::cout << "Pthat remix: " << pthatVals2[pthatIter] << std::endl;
  }

  const Int_t nPtHatBins = 40;
  const Float_t ptHatLow = 7.999;
  const Float_t ptHatHi = 1000.001;
  Double_t ptHatBins[nPtHatBins+1];
  getLogBins(ptHatLow, ptHatHi, nPtHatBins, ptHatBins);

  const Int_t nJtEtaBins = 8;
  const Float_t jtEtaLow = -1.6;
  const Float_t jtEtaHi = 1.6;
  Double_t jtEtaBins[nJtEtaBins+1];
  getLinBins(jtEtaLow, jtEtaHi, nJtEtaBins, jtEtaBins);

  const Int_t nJtPtBins2 = 5;
  const Float_t jtPtLow2 = 30;
  const Float_t jtPtHi2 = 300;
  Double_t jtPtBins2[nJtPtBins2+1];
  getLogBins(jtPtLow2, jtPtHi2, nJtPtBins2, jtPtBins2);

  const Int_t nJtEtaBins2 = 3;
  const Float_t jtEtaLow2 = -1.6;
  const Float_t jtEtaHi2 = 1.6;
  Double_t jtEtaBins2[nJtEtaBins2+1] = {-2.1, -1.3, 1.3, 2.1};
  //  getLinBins(jtEtaLow2, jtEtaHi2, nJtEtaBins2, jtEtaBins2);


  const Int_t nRhoBins = 10;
  const Float_t rhoLow = 1;
  const Float_t rhoHi = 300;
  Double_t rhoBins[nRhoBins+1];
  getLogBins(rhoLow, rhoHi, nRhoBins, rhoBins);

  const Int_t nRhoBins2 = 10;
  Double_t rhoBins2[nRhoBins2+1];
  getLogBins(rhoLow, rhoHi, nRhoBins2, rhoBins2);

  const Int_t nRhoBins3 = 100;
  Double_t rhoBins3[nRhoBins3+1];
  getLogBins(rhoLow, rhoHi, nRhoBins3, rhoBins3);

  for(Int_t iter = 0; iter < nJtPtBins+1; iter++){
    std::cout << "Pt Bin, val: " << iter << ", " << jtPtBins[iter] << std::endl;
  }

  for(Int_t iter = 0; iter < nJtEtaBins+1; iter++){
    std::cout << "Eta BFin, val: " << iter << ", " << jtEtaBins[iter] << std::endl;
  }

  Int_t nCentBinsTemp = nCentBins;
  if(!isPbPb) nCentBinsTemp = 1;

  const Int_t nCentBins2 = nCentBinsTemp;


  if(debugMode) std::cout << __LINE__ << std::endl;

  const Int_t nRhoCentBins = 25;
  const Float_t rhoCentLow = 0;
  const Float_t rhoCentHi = 100;
  Double_t rhoCentBins[nRhoCentBins+1];
  getLinBins(rhoCentLow, rhoCentHi, nRhoCentBins, rhoCentBins);

  TH2F* centVRho_p = new TH2F("centVRho_h", ";Centrality;#rho", nRhoCentBins, rhoCentBins, nRhoBins2, rhoBins2);
  TH1F* centVRho_Proj_p[nRhoCentBins];
  for(Int_t rhoCentIter = 0; rhoCentIter < nRhoCentBins; rhoCentIter++){
    centVRho_Proj_p[rhoCentIter] = new TH1F(Form("centVRho_Proj_Cent%dto%d_h", (Int_t)rhoCentBins[rhoCentIter], (Int_t)rhoCentBins[rhoCentIter+1]), ";#rho;Events", nRhoBins3, rhoBins3);
  }

  TH1F* pthat_Unweighted_p[nJetAlgo];
  TH1F* pthat_Weighted_p[nJetAlgo];
  TH1F* pthat_Unweighted_Sample_p[nJetAlgo][nFiles];
  TH1F* jtPerPthat_p[nJetAlgo][nFiles];


  TH1F* bkgEstimate_p[nJetAlgo][nCentBins2];
  TH1F* areaEstimate_p[nJetAlgo][nCentBins2];
  TH2F* bkgAreaEstimate_p[nJetAlgo][nCentBins2];

  TH1F* jtPt_p[nJetAlgo][nCentBins2][nQG];
  TH1F* genJtPt_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtEta_p[nJetAlgo][nCentBins2][nQG];

  TH2F* jtRecoGenDRVPt_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoGenDRVPt_Mean_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoGenDRVPt_Res_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoGenDRVPt_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];
  TH1F* jtRecoGenDPhiVPt_Mean_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoGenDPhiVPt_Res_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoGenDPhiVPt_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];
  TH1F* jtRecoGenDEtaVPt_Mean_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoGenDEtaVPt_Res_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoGenDEtaVPt_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];

  TH1F* jtRecoGenDRVEta_Mean_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoGenDRVEta_Res_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoGenDRVEta_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtEtaBins];
  TH1F* jtRecoGenDPhiVEta_Mean_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoGenDPhiVEta_Res_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoGenDPhiVEta_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtEtaBins];
  TH1F* jtRecoGenDEtaVEta_Mean_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoGenDEtaVEta_Res_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoGenDEtaVEta_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtEtaBins];

  TH2F* jtRecoVGen_p[nJetAlgo][nCentBins2][nQG];
  TH2F* jtRecoOverGenVPt_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoOverGenVPt_Mean_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRecoOverGenVPt_Res_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRecoOverGenVPt_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];

  TH1F* jtRecoOverGenVPt_PERP_Res_p[nJetAlgo][nQG][nMeanFit];
  TH1F* jtRecoOverGenVPt_PERP_MeanResPts_p[nJetAlgo][nQG][nJtPtBinsPERP];

  TH1F* jtRecoOverGenVPt_Rho_Res_p[nJetAlgo][nRhoBins][nQG][nMeanFit];
  TH1F* jtRecoOverGenVPt_Rho_MeanResPts_p[nJetAlgo][nRhoBins][nQG][nJtPtBins];

  TH1F* jtRecoOverGenVRecoPt_Mean_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRecoOverGenVRecoPt_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];

  TH1F* jtRecoOverGenVCent_Mean_p[nJetAlgo][nJtPtBins2][nJtEtaBins2][nQG][nMeanFit];
  TH1F* jtRecoOverGenVCent_Res_p[nJetAlgo][nJtPtBins2][nJtEtaBins2][nQG][nMeanFit];
  TH1F* jtRecoOverGenVCent_MeanResPts_p[nJetAlgo][nJtPtBins2][nJtEtaBins2][nQG][nCentBins2];

  TH2F* jtRecoOverRawVPt_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoOverRawVPt_Mean_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRecoOverRawVPt_Res_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRecoOverRawVPt_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];

  TH2F* jtRawOverGenVPt_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRawOverGenVPt_Mean_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRawOverGenVPt_Res_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRawOverGenVPt_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];

  TH2F* jtEffVPt_p[nJetAlgo][nCentBins2][nQG][nJetRecoCuts];
  TH1F* jtEffVPt_Mean_p[nJetAlgo][nCentBins2][nQG][nJetRecoCuts];
  TH1F* jtEffVPt_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJetRecoCuts][nJtPtBins];

  TH2F* jtFakeVPt_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtFakeVPt_Mean_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtFakeVPt_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];

  TH2F* jtEffVPtEta_p[nJetAlgo][nCentBins2][nQG][nJetRecoCuts];
  TH2F* jtEffVPtEta_Denom_p[nJetAlgo][nCentBins2][nQG][nJetRecoCuts];

  TH2F* jtRecoOverGenVEta_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoOverGenVEta_Mean_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRecoOverGenVEta_Res_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRecoOverGenVEta_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtEtaBins];

  TH2F* jtRecoOverRawVEta_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoOverRawVEta_Mean_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRecoOverRawVEta_Res_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRecoOverRawVEta_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtEtaBins];

  TH2F* jtRawOverGenVEta_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRawOverGenVEta_Mean_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRawOverGenVEta_Res_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRawOverGenVEta_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtEtaBins];

  TH2F* jtEffVEta_p[nJetAlgo][nCentBins2][nQG][nJetRecoCuts];
  TH1F* jtEffVEta_Mean_p[nJetAlgo][nCentBins2][nQG][nJetRecoCuts];
  TH1F* jtEffVEta_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJetRecoCuts][nJtEtaBins];

  TH2F* jtFakeVEta_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtFakeVEta_Mean_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtFakeVEta_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtEtaBins];

  Int_t nNoLeadEle = 0;
  Int_t nNoSubleadEle = 0;

  Int_t nNoLeadMu = 0;
  Int_t nNoSubleadMu = 0;

  Int_t nNoLeadPho = 0;

  Int_t etaFills[nJetAlgo][nCentBins2];
  Int_t testFill1[nJetAlgo][nCentBins2];
  Int_t testFill2[nJetAlgo][nCentBins2];
  Int_t testFill3[nJetAlgo][nCentBins2];


  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    for(Int_t iter2 = 0; iter2 < nCentBins2; iter2++){
      etaFills[iter][iter2] = 0;

      testFill1[iter][iter2] = 0;
      testFill2[iter][iter2] = 0;
      testFill3[iter][iter2] = 0;
    }
  }

  for(Int_t iter = 0; iter < nJetAlgo; iter++){

    pthat_Unweighted_p[iter] = new TH1F(Form("pthat_Unweighted_%s_h", jetAlgo.at(iter).c_str()), ";pthat;Events", nPtHatBins, ptHatBins);
    pthat_Weighted_p[iter] = new TH1F(Form("pthat_Weighted_%s_h", jetAlgo.at(iter).c_str()), ";pthat;Events", nPtHatBins, ptHatBins);
    pthat_Weighted_p[iter]->Sumw2();

    for(Int_t pthatIter = 0; pthatIter < nFiles; pthatIter++){
      pthat_Unweighted_Sample_p[iter][pthatIter] = new TH1F(Form("pthat_Unweighted_%s_%d_h", jetAlgo.at(iter).c_str(), pthatVals[pthatIter]), Form(";pthat;Events (min pthat = %d", pthatVals[pthatIter]), nPtHatBins, ptHatBins);

      jtPerPthat_p[iter][pthatIter] = new TH1F(Form("jtPerPthat_%s_%d_h", jetAlgo.at(iter).c_str(), pthatVals[pthatIter]), Form(";pthat;Events (min pthat = %d)", pthatVals[pthatIter]), nJtPtBins, jtPtBins);
    }

    for(Int_t ptIter = 0; ptIter < nJtPtBins2; ptIter++){
      Int_t ptLowInt = std::trunc(jtPtBins2[ptIter]);
      Int_t ptHiInt = std::trunc(jtPtBins2[ptIter+1]);
      
      Int_t ptLowDec = std::trunc(jtPtBins2[ptIter]*10 - ptLowInt*10);
      Int_t ptHiDec = std::trunc(jtPtBins2[ptIter+1]*10 - ptHiInt*10);

      std::string ptStr1 = Form("Pt%dp%dTo%dp%d", ptLowInt, ptLowDec, ptHiInt, ptHiDec);
      std::string ptStr2 = Form("(%d.%d<p_{T,Gen.}<%d.%d)", ptLowInt, ptLowDec, ptHiInt, ptHiDec);

      for(Int_t etaIter = 0; etaIter < nJtEtaBins2; etaIter++){
	Int_t etaLowInt = std::trunc(jtEtaBins2[etaIter]);
	Int_t etaHiInt = std::trunc(jtEtaBins2[etaIter+1]);

	Int_t etaLowDec = std::trunc(jtEtaBins2[etaIter]*10 - etaLowInt*10);
	Int_t etaHiDec = std::trunc(jtEtaBins2[etaIter+1]*10 - etaHiInt*10);

	std::string etaStr1 = Form("Eta%dp%dTo%dp%d", etaLowInt, etaLowDec, etaHiInt, etaHiDec);
	std::string etaStr2 = Form("(%d.%d<#eta_{Gen.}<%d.%d)", etaLowInt, etaLowDec, etaHiInt, etaHiDec);

	for(Int_t qgIter = 0; qgIter < nQG; qgIter++){

	  for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	    jtRecoOverGenVCent_Mean_p[iter][ptIter][etaIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVCent_%s_%sMean_%s_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), ptStr1.c_str(), etaStr1.c_str()), Form(";Centrality;#mu_{Reco./Gen.} (%s)", jetAlgo.at(iter).c_str()), nCentBins2, centBins2);
	    jtRecoOverGenVCent_Res_p[iter][ptIter][etaIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVCent_%s_%sRes_%s_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), ptStr1.c_str(), etaStr1.c_str()), Form(";Centrality;#sigma_{Reco./Gen.} (%s)", jetAlgo.at(iter).c_str()), nCentBins2, centBins2);
	  }

          for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
	    std::string centStr = "PP";
	    if(isPbPb) centStr = Form("cent%dto%d", centBins[centIter+1]/2, centBins[centIter]/2);

	    jtRecoOverGenVCent_MeanResPts_p[iter][ptIter][etaIter][qgIter][centIter] = new TH1F(Form("jtRecoOverGenVCent_%s_MeanResPts_%s_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), ptStr1.c_str(), etaStr1.c_str(), centStr.c_str()), Form(";(Reco. %s Jet p_{T})/(Gen. Jet p_{T});Events (%s)", jetAlgo.at(iter).c_str(), centStr.c_str()), 30, 0, 3);
	  }
	}
      }
    }

    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      std::string centStr = "PP";
      if(isPbPb) centStr = Form("cent%dto%d", centBins[centIter+1]/2, centBins[centIter]/2);

      
      bkgEstimate_p[iter][centIter] = new TH1F(Form("bkgEstimate_%s_%s_h", jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Bkg. Estimate;Events"), 100, 0, 50);

      areaEstimate_p[iter][centIter] = new TH1F(Form("areaEstimate_%s_%s_h", jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Area;Events"), 100, 0, 4*2*TMath::Pi());

      bkgAreaEstimate_p[iter][centIter] = new TH2F(Form("bkgAreaEstimate_%s_%s_h", jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Area;Bkg."), 20, 0, 4*2*TMath::Pi(), 20, 0, 50);

      for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	jtPt_p[iter][centIter][qgIter] = new TH1F(Form("jtPt_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";p_{T};Events (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);
	genJtPt_p[iter][centIter][qgIter] = new TH1F(Form("genJtPt_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. p_{T};Events (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);

	jtEta_p[iter][centIter][qgIter] = new TH1F(Form("jtEta_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";#eta;Events (%s)", jetAlgo.at(iter).c_str()), 100, -3., 3.);

	jtRecoGenDRVPt_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoGenDRVPt_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; #DeltaR_{Reco. %s, Gen.}", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins, 20, 0, 0.4);

	jtRecoGenDRVPt_Mean_p[iter][centIter][qgIter] = new TH1F(Form("jtRecoGenDRVPt_Mean_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; <#DeltaR_{Reco. %s, Gen.}>", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);

	jtRecoGenDRVPt_Res_p[iter][centIter][qgIter] = new TH1F(Form("jtRecoGenDRVPt_Res_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; #sigma(#DeltaR_{Reco. %s, Gen.})", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);

	jtRecoGenDPhiVPt_Mean_p[iter][centIter][qgIter] = new TH1F(Form("jtRecoGenDPhiVPt_Mean_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; <#Delta#phi_{Reco. %s, Gen.}>", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);

	jtRecoGenDPhiVPt_Res_p[iter][centIter][qgIter] = new TH1F(Form("jtRecoGenDPhiVPt_Res_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; #sigma(#Delta#phi_{Reco. %s, Gen.})", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);

	jtRecoGenDEtaVPt_Mean_p[iter][centIter][qgIter] = new TH1F(Form("jtRecoGenDEtaVPt_Mean_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; <#Delta#eta_{Reco. %s, Gen.}>", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);

	jtRecoGenDEtaVPt_Res_p[iter][centIter][qgIter] = new TH1F(Form("jtRecoGenDEtaVPt_Res_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; #sigma(#Delta#eta_{Reco. %s, Gen.})", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);

	for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	  Int_t ptLowInt = std::trunc(jtPtBins[jtIter]);
	  Int_t ptHiInt = std::trunc(jtPtBins[jtIter+1]);
	  
	  Int_t ptLowDec = std::trunc(jtPtBins[jtIter]*10 - ptLowInt*10);
	  Int_t ptHiDec = std::trunc(jtPtBins[jtIter+1]*10 - ptHiInt*10);
	  
	  jtRecoGenDRVPt_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoGenDRVPt_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";#DeltaR;Events (%d.%d<p_{T,Gen.}<%d.%d)", ptLowInt, ptLowDec, ptHiInt, ptHiDec), 30, 0, 0.4);

	  jtRecoGenDPhiVPt_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoGenDPhiVPt_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";#Delta#phi;Events (%d.%d<p_{T,Gen.}<%d.%d)", ptLowInt, ptLowDec, ptHiInt, ptHiDec), 30, -0.4, 0.4);

	  jtRecoGenDEtaVPt_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoGenDEtaVPt_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";#Delta#eta;Events (%d.%d<p_{T,Gen.}<%d.%d)", ptLowInt, ptLowDec, ptHiInt, ptHiDec), 30, -0.4, 0.4);
	}

	jtRecoGenDRVEta_Mean_p[iter][centIter][qgIter] = new TH1F(Form("jtRecoGenDRVEta_Mean_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta; <#DeltaR_{Reco. %s, Gen.}>", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);

	jtRecoGenDRVEta_Res_p[iter][centIter][qgIter] = new TH1F(Form("jtRecoGenDRVEta_Res_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta; #sigma(#DeltaR_{Reco. %s, Gen.})", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);

	jtRecoGenDPhiVEta_Mean_p[iter][centIter][qgIter] = new TH1F(Form("jtRecoGenDPhiVEta_Mean_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta; <#Delta#phi_{Reco. %s, Gen.}>", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);

	jtRecoGenDPhiVEta_Res_p[iter][centIter][qgIter] = new TH1F(Form("jtRecoGenDPhiVEta_Res_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta; #sigma(#Delta#phi_{Reco. %s, Gen.})", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);

	jtRecoGenDEtaVEta_Mean_p[iter][centIter][qgIter] = new TH1F(Form("jtRecoGenDEtaVEta_Mean_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta; <#Delta#eta_{Reco. %s, Gen.}>", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);

	jtRecoGenDEtaVEta_Res_p[iter][centIter][qgIter] = new TH1F(Form("jtRecoGenDEtaVEta_Res_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta; #sigma(#Delta#eta_{Reco. %s, Gen.})", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);

	for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	  Int_t etaLowInt = std::trunc(jtEtaBins[jtIter]);
	  Int_t etaHiInt = std::trunc(jtEtaBins[jtIter+1]);
	  
	  Int_t etaLowDec = std::trunc(jtEtaBins[jtIter]*10 - etaLowInt*10);
	  Int_t etaHiDec = std::trunc(jtEtaBins[jtIter+1]*10 - etaHiInt*10);


	  jtRecoGenDRVEta_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoGenDRVEta_%s_MeanResPts_Eta%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";#DeltaR;Events (%d.%d<#eta_{Gen.}<%d.%d)", etaLowInt, etaLowDec, etaHiInt, etaHiDec), 30, 0, 0.4);

	  jtRecoGenDPhiVEta_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoGenDPhiVEta_%s_MeanResPts_Eta%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";#Delta#phi;Events (%d.%d<#eta_{Gen.}<%d.%d)", etaLowInt, etaLowDec, etaHiInt, etaHiDec), 30, -0.4, 0.4);

	  jtRecoGenDEtaVEta_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoGenDEtaVEta_%s_MeanResPts_Eta%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";#Delta#eta;Events (%d.%d<#eta_{Gen.}<%d.%d)", etaLowInt, etaLowDec, etaHiInt, etaHiDec), 30, -0.4, 0.4);
	}



	jtRecoVGen_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoVGen_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; Reco. %s Jet p_{T}", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins, nJtPtBins, jtPtBins);

	jtRecoOverGenVPt_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoOverGenVPt_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; (Reco. %s Jet p_{T})/(Gen. Jet p_{T})", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins, 30, 0, 3);

	jtRecoOverRawVPt_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoOverRawVPt_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Reco. Jet p_{T}; (Reco %s Jet p_{T})/(Raw Jet p_{T})", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins, 30, 0, 3);

	jtRawOverGenVPt_p[iter][centIter][qgIter] = new TH2F(Form("jtRawOverGenVPt_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; (Raw %s Jet p_{T})/(Gen Jet p_{T})", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins, 30, 0, 3);

	for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	  jtEffVPt_p[iter][centIter][qgIter][recoCutIter] = new TH2F(Form("jtEffVPt_Reco%d_%s_%s_%s_h", jetRecoCuts[recoCutIter], qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};Eff. (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins, 2, -0.5, 1.5);
	}

	jtFakeVPt_p[iter][centIter][qgIter] = new TH2F(Form("jtFakeVPt_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Reco. Jet p_{T};Fake (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins, 2, -0.5, 1.5);

	for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	  jtEffVPtEta_p[iter][centIter][qgIter][recoCutIter] = new TH2F(Form("jtEffVPtEta_Reco%d_%s_%s_%s_h", jetRecoCuts[recoCutIter], qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta;Gen. Jet p_{T} (%s)", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins, nJtPtBins, jtPtBins);
	  jtEffVPtEta_Denom_p[iter][centIter][qgIter][recoCutIter] = new TH2F(Form("jtEffVPtEta_Reco%d_%s_Denom_%s_%s_h", jetRecoCuts[recoCutIter], qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta;Gen. Jet p_{T} (%s)", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins, nJtPtBins, jtPtBins);
	}

	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  jtRecoOverGenVPt_Mean_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVPt_%s_%sMean_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#mu_{Reco./Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);

	  jtRecoOverGenVRecoPt_Mean_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVRecoPt_%s_%sMean_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Reco. Jet p_{T};#mu_{Reco./Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);

	  jtRecoOverRawVPt_Mean_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverRawVPt_%s_%sMean_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Reco. Jet p_{T};#mu_{Reco./Raw} (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);

	  jtRawOverGenVPt_Mean_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRawOverGenVPt_%s_%sMean_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#mu_{Raw/Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);
	}

        for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	  jtEffVPt_Mean_p[iter][centIter][qgIter][recoCutIter] = new TH1F(Form("jtEffVPt_Reco%d_%s_Mean_%s_%s_h", jetRecoCuts[recoCutIter], qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};Eff. (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);
	}

	jtFakeVPt_Mean_p[iter][centIter][qgIter] = new TH1F(Form("jtFakeVPt_%s_Mean_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Reco. Jet p_{T};Fake (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);
	
	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVPt_%s_%sRes_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#sigma_{Reco./Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);
	  
	  jtRecoOverRawVPt_Res_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverRawVPt_%s_%sRes_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Reco. Jet p_{T};#sigma_{Reco./Raw} (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);
	  
	  jtRawOverGenVPt_Res_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRawOverGenVPt_%s_%sRes_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#sigma_{Raw/Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);


	  if(centIter == 0){
	    jtRecoOverGenVPt_PERP_Res_p[iter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVPt_PERP_%s_%sRes_%s__h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str()), Form(";Gen. Jet p_{T};#sigma_{Reco./Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtPtBinsPERP, jtPtBinsPERP);

	     for(Int_t rhoIter = 0; rhoIter < nRhoBins; rhoIter++){
	       Int_t rhoLowInt = std::trunc(rhoBins[rhoIter]);
	       Int_t rhoHiInt = std::trunc(rhoBins[rhoIter+1]);
	       
	       Int_t rhoLowDec = std::trunc(rhoBins[rhoIter]*10 - rhoLowInt*10);
	       Int_t rhoHiDec = std::trunc(rhoBins[rhoIter+1]*10 - rhoHiInt*10);
	       
	       jtRecoOverGenVPt_Rho_Res_p[iter][rhoIter][qgIter][mIter] = new TH1F(Form("jtRecoRho%dOverGenVPt_Rho%dp%dto%dp%d_%s_%sRes_%s_h", rhoIter, rhoLowInt, rhoLowDec, rhoHiInt, rhoHiDec, qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str()), Form(";Gen. Jet p_{T};#sigma_{Reco./Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);
	     }
	  }

	}
	
	jtRecoOverGenVEta_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoOverGenVEta_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta; (Reco. %s Jet p_{T})/(Gen. Jet p_{T})", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins, 30, 0, 3);

	jtRecoOverRawVEta_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoOverRawVEta_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Reco. Jet #eta; (Reco %s Jet p_{T})/(Raw Jet p_{T})", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins, 30, 0, 3);

	jtRawOverGenVEta_p[iter][centIter][qgIter] = new TH2F(Form("jtRawOverGenVEta_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta; (Raw %s Jet p_{T})/(Gen Jet p_{T})", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins, 30, 0, 3);

        for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	  jtEffVEta_p[iter][centIter][qgIter][recoCutIter] = new TH2F(Form("jtEffVEta_Reco%d_%s_%s_%s_h", jetRecoCuts[recoCutIter], qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta;Eff. (%s)", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins, 2, -0.5, 1.5);
	}

	jtFakeVEta_p[iter][centIter][qgIter] = new TH2F(Form("jtFakeVEta_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Reco. Jet #eta;Fake (%s)", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins, 2, -0.5, 1.5);

	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  jtRecoOverGenVEta_Mean_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVEta_%s_%sMean_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta;#mu_{Reco./Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);

	  jtRecoOverRawVEta_Mean_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverRawVEta_%s_%sMean_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Reco. Jet #eta;#mu_{Reco./Raw} (%s)", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);
	  
	  jtRawOverGenVEta_Mean_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRawOverGenVEta_%s_%sMean_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta;#mu_{Raw/Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);
	}
	
        for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	  jtEffVEta_Mean_p[iter][centIter][qgIter][recoCutIter] = new TH1F(Form("jtEffVEta_Reco%d_%s_Mean_%s_%s_h", jetRecoCuts[recoCutIter], qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta;Eff. (%s)", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);
	}

	jtFakeVEta_Mean_p[iter][centIter][qgIter] = new TH1F(Form("jtFakeVEta_%s_Mean_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Reco. Jet #eta;Fake (%s)", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);
	
	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  jtRecoOverGenVEta_Res_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVEta_%s_%sRes_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta;#sigma_{Reco./Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);
	  
	  jtRecoOverRawVEta_Res_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverRawVEta_%s_%sRes_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Reco. Jet #eta;#sigma_{Reco./Raw} (%s)", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);
	  
	  jtRawOverGenVEta_Res_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRawOverGenVEta_%s_%sRes_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet #eta;#sigma_{Raw/Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);
	}
	

	if(centIter == 0){
	  for(Int_t jtIter = 0; jtIter < nJtPtBinsPERP; jtIter++){
	    Int_t ptLowInt = std::trunc(jtPtBinsPERP[jtIter]);
	    Int_t ptHiInt = std::trunc(jtPtBinsPERP[jtIter+1]);
	    
	    Int_t ptLowDec = std::trunc(jtPtBinsPERP[jtIter]*10 - ptLowInt*10);
	    Int_t ptHiDec = std::trunc(jtPtBinsPERP[jtIter+1]*10 - ptHiInt*10);

	    jtRecoOverGenVPt_PERP_MeanResPts_p[iter][qgIter][jtIter] = new TH1F(Form("jtRecoOverGenVPt_PERP_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo.at(iter).c_str()), Form(";(Reco. %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo.at(iter).c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 45, 0, 3);
	  }
	}


	for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	  Int_t ptLowInt = std::trunc(jtPtBins[jtIter]);
	  Int_t ptHiInt = std::trunc(jtPtBins[jtIter+1]);
	  
	  Int_t ptLowDec = std::trunc(jtPtBins[jtIter]*10 - ptLowInt*10);
	  Int_t ptHiDec = std::trunc(jtPtBins[jtIter+1]*10 - ptHiInt*10);
	  
	  jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoOverGenVPt_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";(Reco. %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo.at(iter).c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 45, 0, 3);
	  jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Sumw2();

	  jtRecoOverGenVRecoPt_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoOverGenVRecoPt_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";(Reco. %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<p_{T,Reco.}<%d.%d)", jetAlgo.at(iter).c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 45, 0, 3);
	  jtRecoOverGenVRecoPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Sumw2();
	  
	  if(centIter == 0){
	    //	    jtRecoOverGenVPt_PERP_MeanResPts_p[iter][qgIter][jtIter] = new TH1F(Form("jtRecoOverGenVPt_PERP_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo.at(iter).c_str()), Form(";(Reco. %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo.at(iter).c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 45, 0, 3);

	    for(Int_t rhoIter = 0; rhoIter < nRhoBins; rhoIter++){
	      Int_t rhoLowInt = std::trunc(rhoBins[rhoIter]);
	      Int_t rhoHiInt = std::trunc(rhoBins[rhoIter+1]);
	      
	      Int_t rhoLowDec = std::trunc(rhoBins[rhoIter]*10 - rhoLowInt*10);
	      Int_t rhoHiDec = std::trunc(rhoBins[rhoIter+1]*10 - rhoHiInt*10);
	      
	      jtRecoOverGenVPt_Rho_MeanResPts_p[iter][rhoIter][qgIter][jtIter] = new TH1F(Form("jtRecoRho%dOverGenVPt_Rho%dp%dto%dp%d_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_h", rhoIter, rhoLowInt, rhoLowDec, rhoHiInt, rhoHiDec, qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo.at(iter).c_str()), Form(";(Reco. %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo.at(iter).c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 30, 0, 3);
	    }
	  }
	  
	  jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoOverRawVPt_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";(Reco. %s Jet p_{T})/(Raw. Jet p_{T});Events (%d.%d<p_{T,Reco.}<%d.%d)", jetAlgo.at(iter).c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 30, 0, 3);
	  
	  jtRawOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRawOverGenVPt_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";(Raw %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo.at(iter).c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 30, 0, 3);
	  
	  for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	    jtEffVPt_MeanResPts_p[iter][centIter][qgIter][recoCutIter][jtIter] = new TH1F(Form("jtEffVPt_Reco%d_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", jetRecoCuts[recoCutIter], qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";(Eff. (%s));Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo.at(iter).c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 2, -0.5, 1.5);
	  }

	  jtFakeVPt_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtFakeVPt_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";(Fake (%s));Events (%d.%d<p_{T,Reco.}<%d.%d)", jetAlgo.at(iter).c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 2, -0.5, 1.5);
	}
	
	
	for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	  Int_t etaLowInt = std::trunc(jtEtaBins[jtIter]);
	  Int_t etaHiInt = std::trunc(jtEtaBins[jtIter+1]);
	  
	  Int_t etaLowDec = std::trunc(jtEtaBins[jtIter]*10 - etaLowInt*10);
	  Int_t etaHiDec = std::trunc(jtEtaBins[jtIter+1]*10 - etaHiInt*10);
	  
	  jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoOverGenVEta_%s_MeanResPts_Eta%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";(Reco. %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<#eta_{Gen.}<%d.%d)", jetAlgo.at(iter).c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec), 30, 0, 3);
	  
	  jtRecoOverRawVEta_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoOverRawVEta_%s_MeanResPts_Eta%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";(Reco. %s Jet p_{T})/(Raw. Jet p_{T});Events (%d.%d<#eta_{Reco.}<%d.%d)", jetAlgo.at(iter).c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec), 30, 0, 3);
	  
	  jtRawOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRawOverGenVEta_%s_MeanResPts_Eta%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";(Raw %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<#eta_{Gen.}<%d.%d)", jetAlgo.at(iter).c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec), 30, 0, 3);
	  
	  for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	    jtEffVEta_MeanResPts_p[iter][centIter][qgIter][recoCutIter][jtIter] = new TH1F(Form("jtEffVEta_Reco%d_%s_MeanResPts_Eta%dp%dTo%dp%d_%s_%s_h", jetRecoCuts[recoCutIter], qg[qgIter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";(Eff. (%s));Events (%d.%d<#eta_{Gen.}<%d.%d)", jetAlgo.at(iter).c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec), 2, -0.5, 1.5);
	  }

	  jtFakeVEta_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtFakeVEta_%s_MeanResPts_Eta%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";(Fake (%s));Events (%d.%d<#eta_{Reco.}<%d.%d)", jetAlgo.at(iter).c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec), 2, -0.5, 1.5);
	}
      }
    }
  }


  if(debugMode) std::cout << __LINE__ << std::endl;


  Int_t nPerPtHat[nFiles];
  for(Int_t iter = 0; iter < nFiles; iter++){
    nPerPtHat[iter] = 0;
  }

  TH1F* localHibin_p = new TH1F("localHibin_h", ";hiBin;Events", 200, -0.5, 199.5);
  TH1F* hibinWeighted_p = new TH1F("hibinWeighted_h", ";hiBin;Events", 200, -0.5, 199.5);

  for(Int_t pthatIter = 0; pthatIter < nFiles; pthatIter++){
    if(!isFile[pthatIter]) continue;

    const Int_t numberOfFiles = (Int_t)listOfFiles[pthatIter].size();
    Int_t fileDiv = ((Int_t)(numberOfFiles/10));
    if(fileDiv < 1) fileDiv = 1;
    
    for(Int_t fileIter = 0; fileIter < numberOfFiles; fileIter++){
      TFile* inFile_p = TFile::Open(listOfFiles[pthatIter].at(fileIter).c_str(), "READ");
      if(inFile_p == NULL) continue;
      if(inFile_p->GetSize() < 1000){
	std::cout << "File " << listOfFiles[pthatIter].at(fileIter) << " less than 1 kb. Continue" << std::endl;
        continue;
      }
      else if(inFile_p->IsZombie()){
	std::cout << "File " << listOfFiles[pthatIter].at(fileIter) << " is zombie. Continue" << std::endl;
        continue;
      }

      if(debugMode) std::cout << __LINE__ << std::endl;

      TTree* jetTree_p = (TTree*)inFile_p->Get(jetAlgoInFile[0].c_str());
      TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
      if(debugMode) std::cout << __LINE__ << std::endl;

      jetTree_p->SetBranchStatus("*", 0);
      jetTree_p->SetBranchStatus("pthat", 1);
      jetTree_p->SetBranchAddress("pthat", &ptHat_[0]);

      if(debugMode) std::cout << __LINE__ << std::endl;

      hiTree_p->SetBranchStatus("*", 0);
      hiTree_p->SetBranchStatus("hiBin", 1);

      hiTree_p->SetBranchAddress("hiBin", &hiBin_);

      if(debugMode) std::cout << __LINE__ << std::endl;

      Int_t tempStartPos = 0;
      Int_t tempNEntries = jetTree_p->GetEntries();
      if(deriveResidualMode) tempNEntries /= 2;
      else if(applyResidualMode) tempStartPos = tempNEntries/2;

      if(deriveResidualMode && applyResidualMode) tempStartPos = 0;

      if(applyResidualModeFull) tempStartPos = 0;

      const Int_t startPos = tempStartPos;
      const Int_t nEntries = tempNEntries;
      Int_t entryDiv = ((Int_t)((nEntries-startPos)/10));

      if(debugMode) std::cout << __LINE__ << std::endl;

      for(Int_t entry = startPos; entry < nEntries; entry++){
	jetTree_p->GetEntry(entry);
	hiTree_p->GetEntry(entry);

        Int_t pthatLowerPos = -1;
        for(Int_t pthatIter2 = nFiles-1; pthatIter2 >= 0; pthatIter2--){
          if(ptHat_[0] > pthatVals[pthatIter2] && isFile[pthatIter2]){
            pthatLowerPos = pthatIter2;
            break;
          }
        }
        if(pthatLowerPos == -1){
	  //std::cout << "Error: pthat == " << ptHat_[0] << " unassigned. Continue." << std::endl;
	  //	  continue;
	  pthatLowerPos = 0;
        }

	localHibin_p->Fill(hiBin_);

	nPerPtHat[pthatLowerPos]++;
      }

      inFile_p->Close();
    }
  }

  Bool_t doHibinWeight = false;
  if(strcmp(hibinWeightFileName.c_str(), "") != 0){
    doHibinWeight = true;

    localHibin_p->Sumw2();

    std::cout << "DOES REWEIGHT" << std::endl;

    TFile* hibinWeightFile_p = new TFile(hibinWeightFileName.c_str(), "READ");
    TH1F* hibinWeight_p = (TH1F*)hibinWeightFile_p->Get("zJetHiBin_h");
    hibinWeight_p->Sumw2();
    hibinWeight_p->Scale(1./hibinWeight_p->Integral());

    localHibin_p->Scale(1./localHibin_p->Integral());
    localHibin_p->Multiply(hibinWeight_p);

    hibinWeightFile_p->Close();
    delete hibinWeightFile_p;
  }

  if(debugMode) std::cout << __LINE__ << std::endl;

  for(Int_t iter = 0; iter < nFiles; iter++){
    std::cout << "pthat bin " << iter << ": " << nPerPtHat[iter] << std::endl;
  }

  //  TFile* ratioFile_p = new TFile("merged_dgulhan-Pythia8_Dijet30_pp_TuneCUETP8M1_Hydjet_MinBias_5020GeV_RECODEBUG_758_PrivMC_forest_v28_0_20160506_QGFRACTIONHIST.root", "READ");
  //  TFile* ratioFile_p = new TFile("outputDir/merged_dgulhan-Pythia8_Dijet30_pp_TuneCUETP8M1_Hydjet_MinBias_5020GeV_RECODEBUG_758_PrivMC_forest_v28_0_20160512_QGFRACTIONHIST.root", "READ");
  //  TH1F* qRatio_p = (TH1F*)ratioFile_p->Get("qRatioZOverDijet_h");
  //  TH1F* gRatio_p = (TH1F*)ratioFile_p->Get("gRatioZOverDijet_h");


    //  fileDiv = 1;
    
    //  std::cout << "Number of files: " << numberOfFiles << std::endl;
   

  for(Int_t pthatIter = 0; pthatIter < nFiles; pthatIter++){
    if(!isFile[pthatIter]) continue;
    const Int_t numberOfFiles = (Int_t)listOfFiles[pthatIter].size();
    Int_t fileDiv = ((Int_t)(numberOfFiles/10));
    if(fileDiv < 1) fileDiv = 1;
 
    for(Int_t fileIter = 0; fileIter < numberOfFiles; fileIter++){
      if(fileIter%fileDiv == 0) std::cout << "File # " << fileIter << "/" << numberOfFiles << std::endl;
      
      TFile* inFile_p = TFile::Open(listOfFiles[pthatIter].at(fileIter).c_str(), "READ");
      
      
      if(inFile_p == NULL) continue;
      if(inFile_p->GetSize() < 1000){
	std::cout << "File " << listOfFiles[pthatIter].at(fileIter) << " less than 1 kb. Continue" << std::endl;
	continue;
      }
      else if(inFile_p->IsZombie()){
	std::cout << "File " << listOfFiles[pthatIter].at(fileIter) << " is zombie. Continue" << std::endl;
	continue;
      }
      
      //    std::cout << listOfFiles.at(fileIter).c_str() << std::endl;

      TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
      TTree* genTree_p = (TTree*)inFile_p->Get("HiGenParticleAna/hi");
      TTree* rhoTree_p = (TTree*)inFile_p->Get("hiFJRhoAnalyzer/t");
      if(debugMode) std::cout << __LINE__ << std::endl;
      
      Bool_t isRho = true;
      if(rhoTree_p == NULL) isRho = false;
      
      TTree* jetTree_p[nJetAlgo];
      
      for(Int_t iter = 0; iter < nJetAlgo; iter++){
	jetTree_p[iter] = (TTree*)inFile_p->Get(Form("%sJetAnalyzer/t", jetAlgo.at(iter).c_str()));
      }

      //      std::cout << "pthatPos, UpperPos, weight hi, events: " << pthatIter << ", " << pthatUpperPos << ", " << pthatWeights[pthatUpperPos] << ", " << jetTree_p[0]->GetEntries() << std::endl;
      
      hiTree_p->SetBranchStatus("*", 0);
      if(isPbPb) hiTree_p->SetBranchStatus("hiBin", 1);
      hiTree_p->SetBranchStatus("vz", 1);
      hiTree_p->SetBranchStatus("run", 1);
      hiTree_p->SetBranchStatus("evt", 1);
    
      if(isPbPb) hiTree_p->SetBranchAddress("hiBin", &hiBin_);
      hiTree_p->SetBranchAddress("vz", &vz_);
      hiTree_p->SetBranchAddress("run", &run_);
      hiTree_p->SetBranchAddress("evt", &evt_);
      
      genTree_p->SetBranchStatus("*", 0);
      genTree_p->SetBranchStatus("pt", 1);
      genTree_p->SetBranchStatus("phi", 1);
      genTree_p->SetBranchStatus("eta", 1);
      genTree_p->SetBranchStatus("pdg", 1);
      
      genTree_p->SetBranchAddress("pt", &genPt_p);
      genTree_p->SetBranchAddress("phi", &genPhi_p);
      genTree_p->SetBranchAddress("eta", &genEta_p);
      genTree_p->SetBranchAddress("pdg", &genPDG_p);
      
      
      if(debugMode) std::cout << __LINE__ << std::endl;
      
      if(isPbPb && isRho){
	rhoTree_p->SetBranchStatus("*", 0);
	rhoTree_p->SetBranchStatus("rho", 1);
	
	rhoTree_p->SetBranchAddress("rho", &rho_p);
      }
      
      if(debugMode) std::cout << __LINE__ << std::endl;
      
      
      //    std::cout << "Gets Here A" << std::endl;
      
      for(Int_t iter = 0; iter < nJetAlgo; iter++){
	jetTree_p[iter]->SetBranchStatus("*", 0);
	jetTree_p[iter]->SetBranchStatus("nref", 1);
	jetTree_p[iter]->SetBranchStatus("jtpt", 1);
	jetTree_p[iter]->SetBranchStatus("rawpt", 1);
	jetTree_p[iter]->SetBranchStatus("jteta", 1);
	jetTree_p[iter]->SetBranchStatus("jtphi", 1);
	jetTree_p[iter]->SetBranchStatus("hcalSum", 1);
	jetTree_p[iter]->SetBranchStatus("ecalSum", 1);
	jetTree_p[iter]->SetBranchStatus("chargedSum", 1);
	jetTree_p[iter]->SetBranchStatus("neutralSum", 1);
	jetTree_p[iter]->SetBranchStatus("photonSum", 1);
	jetTree_p[iter]->SetBranchStatus("eSum", 1);
	jetTree_p[iter]->SetBranchStatus("chargedN", 1);
	jetTree_p[iter]->SetBranchStatus("neutralN", 1);
	jetTree_p[iter]->SetBranchStatus("photonHardN", 1);
	jetTree_p[iter]->SetBranchStatus("photonN", 1);
	jetTree_p[iter]->SetBranchStatus("eN", 1);
	jetTree_p[iter]->SetBranchStatus("muN", 1);
	jetTree_p[iter]->SetBranchStatus("trackN", 1);
	jetTree_p[iter]->SetBranchStatus("refpt", 1);
	jetTree_p[iter]->SetBranchStatus("refeta", 1);
	jetTree_p[iter]->SetBranchStatus("refphi", 1);
	jetTree_p[iter]->SetBranchStatus("subid", 1);
	jetTree_p[iter]->SetBranchStatus("refparton_flavor", 1);
	jetTree_p[iter]->SetBranchStatus("pthat", 1);
	
	jetTree_p[iter]->SetBranchStatus("ngen", 1);
	jetTree_p[iter]->SetBranchStatus("genpt", 1);
	jetTree_p[iter]->SetBranchStatus("geneta", 1);
	jetTree_p[iter]->SetBranchStatus("genphi", 1);
	jetTree_p[iter]->SetBranchStatus("genmatchindex", 1);
	jetTree_p[iter]->SetBranchStatus("gensubid", 1);

	jetTree_p[iter]->SetBranchStatus("jtpu", 1);
	jetTree_p[iter]->SetBranchStatus("trackMax", 1);
	jetTree_p[iter]->SetBranchStatus("neutralMax", 1);
	jetTree_p[iter]->SetBranchStatus("photonMax", 1);
	jetTree_p[iter]->SetBranchStatus("chargedMax", 1);
	
	jetTree_p[iter]->SetBranchAddress("nref", &nJt_[iter]);
	jetTree_p[iter]->SetBranchAddress("jtpt", jtPt_[iter]);
	jetTree_p[iter]->SetBranchAddress("rawpt", jtRawPt_[iter]);
	jetTree_p[iter]->SetBranchAddress("jteta", jtEta_[iter]);
	jetTree_p[iter]->SetBranchAddress("jtphi", jtPhi_[iter]);
	jetTree_p[iter]->SetBranchAddress("hcalSum", jtHCalSum_[iter]);
	jetTree_p[iter]->SetBranchAddress("ecalSum", jtECalSum_[iter]);
	jetTree_p[iter]->SetBranchAddress("neutralSum", jtNeutralSum_[iter]);
	jetTree_p[iter]->SetBranchAddress("photonSum", jtPhotonSum_[iter]);
	jetTree_p[iter]->SetBranchAddress("chargedSum", jtChargedSum_[iter]);
	jetTree_p[iter]->SetBranchAddress("eSum", jtESum_[iter]);
	jetTree_p[iter]->SetBranchAddress("chargedN", jtChargedN_[iter]);
	jetTree_p[iter]->SetBranchAddress("neutralN", &jtNeutralN_[iter]);
	jetTree_p[iter]->SetBranchAddress("photonN", &jtPhotonN_[iter]);
	jetTree_p[iter]->SetBranchAddress("eN", &jtEN_[iter]);
	jetTree_p[iter]->SetBranchAddress("muN", &jtMuN_[iter]);
	jetTree_p[iter]->SetBranchAddress("photonHardN", &jtPhotonHardN_[iter]);
	jetTree_p[iter]->SetBranchAddress("trackN", &jtTrackN_[iter]);
	jetTree_p[iter]->SetBranchAddress("jtpu", jtPU_[iter]);
	jetTree_p[iter]->SetBranchAddress("refpt", refPt_[iter]);
	jetTree_p[iter]->SetBranchAddress("refeta", refEta_[iter]);
	jetTree_p[iter]->SetBranchAddress("refphi", refPhi_[iter]);
	jetTree_p[iter]->SetBranchAddress("subid", refSubID_[iter]);
	jetTree_p[iter]->SetBranchAddress("refparton_flavor", refPartFlav_[iter]);
	jetTree_p[iter]->SetBranchAddress("pthat", &ptHat_[iter]);
	jetTree_p[iter]->SetBranchAddress("ngen", &nGenJt_[iter]);
	jetTree_p[iter]->SetBranchAddress("genpt", genJtPt_[iter]);
	jetTree_p[iter]->SetBranchAddress("geneta", genJtEta_[iter]);
	jetTree_p[iter]->SetBranchAddress("genphi", genJtPhi_[iter]);
	jetTree_p[iter]->SetBranchAddress("genmatchindex", genJtMatchIndex_[iter]);
	jetTree_p[iter]->SetBranchAddress("gensubid", genSubId_[iter]);

	jetTree_p[iter]->SetBranchAddress("trackMax", jtTrackMax_[iter]);
	jetTree_p[iter]->SetBranchAddress("neutralMax", jtNeutralMax_[iter]);
	jetTree_p[iter]->SetBranchAddress("photonMax", jtPhotonMax_[iter]);
	jetTree_p[iter]->SetBranchAddress("chargedMax", jtChargedMax_[iter]);

      }
      
      if(debugMode) std::cout << __LINE__ << std::endl;


      //std::cout << "Gets Here B" << std::endl;
      Int_t tempStartPos = 0;
      Int_t tempNEntries = jetTree_p[0]->GetEntries();
      if(deriveResidualMode) tempNEntries /= 2;
      else if(applyResidualMode) tempStartPos = tempNEntries/2;

      if(deriveResidualMode && applyResidualMode) tempStartPos= 0;

      if(applyResidualModeFull) tempStartPos = 0;

      const Int_t startPos = tempStartPos;
      const Int_t nEntries = tempNEntries;
      Int_t entryDiv = TMath::Max(1, ((Int_t)((nEntries-startPos)/200)));
      
      //    std::cout << "Gets Here C" << std::endl;
      
      for(Int_t entry = startPos; entry < nEntries; entry++){
	//      std::cout << "Gets here d" << std::endl;
	
	if(entry%entryDiv == 0 && nEntries >= 10000) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
	//      std::cout << "Gets here e" << std::endl;
	
	hiTree_p->GetEntry(entry);
	genTree_p->GetEntry(entry);
	if(debugMode) std::cout << __LINE__ << std::endl;
	if(isPbPb && isRho) rhoTree_p->GetEntry(entry);
	
	if(debugMode) std::cout << __LINE__ << std::endl;
	
	Int_t centPos = -1;

	if(TMath::Abs(vz_) > 15) continue;
	
	Float_t rhoVal = -1;
	if(isPbPb && isRho) rhoVal = rho_p->at(3);
	
	if(isPbPb && isRho) centVRho_p->Fill(hiBin_/2, rhoVal);
	
	for(Int_t rhoIter = 0; rhoIter < nRhoCentBins; rhoIter++){
	  if(hiBin_/2 > rhoCentBins[rhoIter] && hiBin_/2 < rhoCentBins[rhoIter+1]){
	    centVRho_Proj_p[rhoIter]->Fill(rhoVal);
	    break;
	  }
	}
	
	Float_t hibinWeight = 1;
	if(doHibinWeight){
	  //  hibinWeight *= localHibin_p->GetBinContent(localHibin_p->FindBin(hiBin_));

	  hibinWeight = findNcoll(hiBin_);
	}

	hibinWeighted_p->Fill(hiBin_, hibinWeight);

	if(isPbPb)
	  for(Int_t iter = 0; iter < nCentBins2; iter++){
	    if(hiBin_ < centBins[iter] && hiBin_ >= centBins[iter+1]){
	      centPos = iter;
	      break;
	    }
	  }
	else centPos = 0;
	
	Int_t rhoPos = -1;
	if(isPbPb){
	  for(Int_t rhoIter = 0; rhoIter < nRhoBins; rhoIter++){
	    if(rhoVal > rhoBins[rhoIter] && rhoVal < rhoBins[rhoIter+1]){
	      rhoPos = rhoIter;
	      break;
	    }
	  }
	}
	else rhoPos = 0;
	
	if(debugMode) std::cout << __LINE__ << std::endl;
	
	Float_t maxElePt = -1;
	Float_t maxEleEta = -100;
	Float_t maxElePhi = -100;
	
	Float_t twoElePt = -1;
	Float_t twoEleEta = -100;
	Float_t twoElePhi = -100;
	
	Float_t maxMuPt = -1;
	Float_t maxMuEta = -100;
	Float_t maxMuPhi = -100;
	
	Float_t twoMuPt = -1;
	Float_t twoMuEta = -100;
	Float_t twoMuPhi = -100;
	
	Float_t maxPhoPt = -1;
	Float_t maxPhoEta = -100;
	Float_t maxPhoPhi = -100;
	
	const Int_t nMult_ = genPt_p->size();
	
	for(Int_t iter = 0; iter < nMult_; iter++){
	  if(TMath::Abs(genPDG_p->at(iter)) != 11 && TMath::Abs(genPDG_p->at(iter)) != 13 && TMath::Abs(genPDG_p->at(iter)) != 22) continue;
	  
	  if(TMath::Abs(genPDG_p->at(iter)) == 11){
	    if(genPt_p->at(iter) > maxElePt){
	      twoElePt = maxElePt;
	      twoEleEta = maxEleEta;
	      twoElePhi = maxElePhi;
	      
	      maxElePt = genPt_p->at(iter);
	      maxEleEta = genEta_p->at(iter);
	      maxElePhi = genPhi_p->at(iter);
	    }
	    else if(genPt_p->at(iter) > twoElePt){
	      twoElePt = genPt_p->at(iter);
	      twoEleEta = genEta_p->at(iter);
	      twoElePhi = genPhi_p->at(iter);
	    }
	  }
	  
	  if(TMath::Abs(genPDG_p->at(iter)) == 13){
	    if(genPt_p->at(iter) > maxMuPt){
	      twoMuPt = maxMuPt;
	      twoMuEta = maxMuEta;
	      twoMuPhi = maxMuPhi;
	      
	      maxMuPt = genPt_p->at(iter);
	      maxMuEta = genEta_p->at(iter);
	      maxMuPhi = genPhi_p->at(iter);
	    }
	    else if(genPt_p->at(iter) > twoMuPt){
	      twoMuPt = genPt_p->at(iter);
	      twoMuEta = genEta_p->at(iter);
	      twoMuPhi = genPhi_p->at(iter);
	    }
	  }

	  if(TMath::Abs(genPDG_p->at(iter)) == 22){
	    if(genPt_p->at(iter) > maxPhoPt){
	      maxPhoPt = genPt_p->at(iter);
	      maxPhoEta = genEta_p->at(iter);
	      maxPhoPhi = genPhi_p->at(iter);
	    }
	  }
	}
	
	if(debugMode) std::cout << __LINE__ << std::endl;
	
	
	if(maxElePt < 10) nNoLeadEle++;
	if(twoElePt < 5) nNoSubleadEle++;
	
	if(maxMuPt < 10) nNoLeadMu++;
	if(twoMuPt < 5) nNoSubleadMu++;
	
	if(maxPhoPt < 10) nNoLeadPho++;

	TLorentzVector zEE, zMuMu;
	if(maxElePt > 10 && twoElePt > 5){
	  TLorentzVector tempLeadEle, tempSubLeadEle;
	  tempLeadEle.SetPtEtaPhiM(maxElePt, maxEleEta, maxElePhi, eleMass);
	  tempSubLeadEle.SetPtEtaPhiM(twoElePt, twoEleEta, twoElePhi, eleMass);
	  zEE = tempLeadEle+tempSubLeadEle;
	}
	if(maxMuPt > 10 && twoMuPt > 5){
	  TLorentzVector tempLeadMu, tempSubLeadMu;
	  tempLeadMu.SetPtEtaPhiM(maxMuPt, maxMuEta, maxMuPhi, muMass);
	  tempSubLeadMu.SetPtEtaPhiM(twoMuPt, twoMuEta, twoMuPhi, muMass);
	  zMuMu = tempLeadMu+tempSubLeadMu;
	}
	
	
	TLorentzVector zFinal;
	if(zEE.Pt() > zMuMu.Pt()) zFinal = zEE;
	else zFinal = zMuMu;
	
	for(Int_t iter = 0; iter < nJetAlgo; iter++){
	  jetTree_p[iter]->GetEntry(entry);
	}

	//	if(ptHat_[0] < 50) continue;

	Int_t pthatLowerPos = -1;
	for(Int_t pthatIter2 = nFiles-1; pthatIter2 >= 0; pthatIter2--){
	  if(ptHat_[0] > pthatVals[pthatIter2] && isFile[pthatIter2]){
	    pthatLowerPos = pthatIter2;
	    break;
	  }
	}
	if(pthatLowerPos == -1){
	  //std::cout << "Error: pthat == " << ptHat_[0] << " unassigned. Continue." << std::endl;
	  //continue;

	  pthatLowerPos = 0;
	}
        
	Int_t pthatUpperPos = -1;
	for(Int_t pthatIter2 = pthatLowerPos+1; pthatIter2 < nFiles; pthatIter2++){
	  if(isFile[pthatIter2]){
	    pthatUpperPos = pthatIter2;
	    break;
	  }
	}
	if(pthatUpperPos == -1) pthatUpperPos = nFiles;
	
	
	if(debugMode) std::cout << __LINE__ << std::endl;
	
	Float_t pthatLowerCut = pthatVals[pthatLowerPos];
	Float_t pthatUpperCut = pthatVals[pthatUpperPos];
	
	Float_t pthatWeight = (pthatWeights[pthatLowerPos] - pthatWeights[pthatUpperPos])/(nPerPtHat[pthatLowerPos]);

	
	//      std::cout << "ALPHA" << std::endl;
	
	Bool_t skipEvent = false;
	
	for(Int_t iter = 0; iter < nGenJt_[0]; iter++){
	  if(genJtPt_[0][iter] < 15) continue;
	  /*
	    if(TMath::Abs(genJtEta_[0][iter]) > 3.0){
	    skipEvent = true;
	    break;
	    }
	  */
	}
	
	if(skipEvent) continue;
	
	if(debugMode) std::cout << __LINE__ << std::endl;	  
      
	if(applyResidualMode){
	  for(Int_t algoIter = 0; algoIter < nJetAlgo; algoIter++){
	    for(Int_t jtIter = 0; jtIter < nJt_[algoIter]; jtIter++){
	      jtPt_[algoIter][jtIter] = jtPt_[algoIter][jtIter]*getResCorrJetPt(jtPt_[algoIter][jtIter], hiBin_);
	      //	      if(refPt_[algoIter][jtIter] > 0) jtPt_[algoIter][jtIter] = jtPt_[algoIter][jtIter]*getResCorrJetPt(refPt_[algoIter][jtIter], hiBin_);
	    }
	  }
	}
	/*
	for(Int_t algoIter = 0; algoIter < nJetAlgo; algoIter++){
	  for(Int_t jtIter = 0; jtIter < nJt_[algoIter]; jtIter++){
	    jtPt_[algoIter][jtIter] = .98*jtPt_[algoIter][jtIter];
	  }
	  }*/

	if(debugMode) std::cout << __LINE__ << std::endl;	  

	for(Int_t algoIter = 0; algoIter < nJetAlgo; algoIter++){
	  //	  if(ptHat_[algoIter] > pthatUpperCut) continue;

	  for(Int_t jtIter = 0; jtIter < nJt_[algoIter]; jtIter++){
	    if(TMath::Abs(jtEta_[algoIter][jtIter]) > 1.6) continue;
	    if(jtPt_[algoIter][jtIter] < 5.0) continue;
	    

	    if(isPbPb && !isGoodCaloJet(jetAlgo.at(algoIter).c_str(), jtHCalSum_[algoIter][jtIter], jtECalSum_[algoIter][jtIter]) && !isGoodPFJet(jetAlgo.at(algoIter).c_str(), jtChargedN_[algoIter][jtIter], jtChargedSum_[algoIter][jtIter], jtNeutralSum_[algoIter][jtIter], jtRawPt_[algoIter][jtIter])) continue;
	    
	    
	    
	    //	if(refSubID_[algoIter][jtIter] != 0 && refSubID_[algoIter][jtIter] != -1) continue;
	    
	    if(maxElePt > 10)
	      if(getDR(jtEta_[algoIter][jtIter], jtPhi_[algoIter][jtIter], maxEleEta, maxElePhi) < 0.4) continue;
	    
	    if(twoElePt > 5)
	      if(getDR(jtEta_[algoIter][jtIter], jtPhi_[algoIter][jtIter], twoEleEta, twoElePhi) < 0.4) continue;
	    
	    if(maxMuPt > 10)
	      if(getDR(jtEta_[algoIter][jtIter], jtPhi_[algoIter][jtIter], maxMuEta, maxMuPhi) < 0.4) continue;
	    
	    if(twoMuPt > 5)
	      if(getDR(jtEta_[algoIter][jtIter], jtPhi_[algoIter][jtIter], twoMuEta, twoMuPhi) < 0.4) continue;
	    
	    if(maxPhoPt > 10)
	      if(getDR(jtEta_[algoIter][jtIter], jtPhi_[algoIter][jtIter], maxPhoEta, maxPhoPhi) < 0.4) continue;

	    Int_t qgPos[2] = {0, -1};
	    if(TMath::Abs(refPartFlav_[algoIter][jtIter]) < 9) qgPos[1] = 1;
	    else if(TMath::Abs(refPartFlav_[algoIter][jtIter]) == 21) qgPos[1] = 2;
	    
	    for(Int_t qgIter = 0; qgIter < 2; qgIter++){
	      if(qgPos[qgIter] == -1) continue;
	      
	      jtRecoOverRawVPt_p[algoIter][centPos][qgPos[qgIter]]->Fill(jtPt_[algoIter][jtIter], jtPt_[algoIter][jtIter]/jtRawPt_[algoIter][jtIter]);
	      
	      for(Int_t jtIter2 = 0; jtIter2 < nJtPtBins; jtIter2++){
		if(jtPt_[algoIter][jtIter] > jtPtBins[jtIter2] && jtPt_[algoIter][jtIter] < jtPtBins[jtIter2+1]){
		  jtRecoOverRawVPt_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(jtPt_[algoIter][jtIter]/jtRawPt_[algoIter][jtIter]);
		  break;
		}
	      }
	      
	      if(jtPt_[algoIter][jtIter] > 30){
		jtRecoOverRawVEta_p[algoIter][centPos][qgPos[qgIter]]->Fill(jtEta_[algoIter][jtIter], jtPt_[algoIter][jtIter]/jtRawPt_[algoIter][jtIter]);
		
		for(Int_t jtIter2 = 0; jtIter2 < nJtEtaBins; jtIter2++){
		  if(jtEta_[algoIter][jtIter] > jtEtaBins[jtIter2] && jtEta_[algoIter][jtIter] < jtEtaBins[jtIter2+1]){
		    jtRecoOverRawVEta_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(jtPt_[algoIter][jtIter]/jtRawPt_[algoIter][jtIter]);
		    break;
		  }
		}
	      }
	    }
	  }
	  
	  //	std::cout << "ALPHA2" << std::endl;
	  
	  genSort(nGenJt_[algoIter], genJtPt_[algoIter], genJtPhi_[algoIter], genJtEta_[algoIter]);
	  
	  Float_t evtBkg = 0;
	  Float_t evtArea = 0;
	  
	  std::vector<Bool_t>* isBkg_p = new std::vector<Bool_t>;
	  	  
	  if(evtBkg != 0){
	    bkgEstimate_p[algoIter][centPos]->Fill(evtBkg);
	    areaEstimate_p[algoIter][centPos]->Fill(evtArea);
	    bkgAreaEstimate_p[algoIter][centPos]->Fill(evtArea, evtBkg);
	  }
	  
	  if(evtBkg > 5){
	    if(evtBkg > 20){
	      std::cout << "evtBkg, run, evt: " << evtBkg << ", " << run_ << ", " << evt_ << std::endl;
	      const Int_t bkgNum = ((Int_t)(isBkg_p->size()));
	      
	      for(Int_t bkgIter = 0; bkgIter < bkgNum; bkgIter++){
		//	      if(isBkg_p->at(bkgIter)) std::cout << "  particle pt, phi, eta: " << genPt_p->at(bkgIter) << ", " << genPhi_p->at(bkgIter) << ", " << genEta_p->at(bkgIter) << std::endl;
		//	    else std::cout << "skipping: " <<  genPt_p->at(bkgIter) << ", " << genPhi_p->at(bkgIter) << ", " << genEta_p->at(bkgIter) << std::endl;
	      }
	    }
	    evtBkg = 5;
	  }
	  
	  isBkg_p->clear();
	  delete isBkg_p;
	  
	  for(Int_t genIter = 0; genIter < nGenJt_[algoIter]; genIter++){
	    genJtPt_[algoIter][genIter] -= evtBkg;
	  }
	  
	  const Int_t boolSize = nJt_[algoIter];
	  Bool_t isUsedRecoJet[boolSize];
	  Int_t genPosJet[boolSize];
	  const Int_t genJtSize = nGenJt_[algoIter];
	  Int_t recoPosJet[genJtSize];
	  for(Int_t boolIter = 0; boolIter < boolSize; boolIter++){
	    isUsedRecoJet[boolIter] = false;
	  }
	  for(Int_t gIter = 0; gIter < genJtSize; gIter++){
	    recoPosJet[gIter] = -1;
	  }
	  
	  //      if(nGen_[algoIter] != 0) std::cout << "NGENJT: " << nGen_[algoIter] << std::endl;
	  
	  
	  Bool_t isLeadDone = false;
	  Bool_t isSubleadDone = false;	  
	  for(Int_t jtIter = 0; jtIter < nGenJt_[algoIter]; jtIter++){
	    if(TMath::Abs(genJtEta_[algoIter][jtIter]) > 1.6) continue;
	    if(genJtPt_[algoIter][jtIter] < 5.0) break;
	    //	  if(genSubId_[algoIter][jtIter] != 0) continue;
	    
	    if(maxElePt > 10)
	      if(getDR(genJtEta_[algoIter][jtIter], genJtPhi_[algoIter][jtIter], maxEleEta, maxElePhi) < 0.4) continue;
	    
	    if(twoElePt > 5)
	      if(getDR(genJtEta_[algoIter][jtIter], genJtPhi_[algoIter][jtIter], twoEleEta, twoElePhi) < 0.4) continue;


	    Int_t qgPos[2] = {0, -1};
            if(TMath::Abs(refPartFlav_[algoIter][jtIter/*tempPos*/]) < 9) qgPos[1] = 1;
            else if(TMath::Abs(refPartFlav_[algoIter][jtIter/*tempPos*/]) == 21) qgPos[1] = 2;

	    //	    if(qgPos[1] != 1 && qgPos[1] != 2) continue;

	    Float_t qgWeight = 1;
	    
	    /*
	    if((applyResidualMode || deriveResidualMode)){
	      if(qgPos[1] == 1) qgWeight *= qRatio_p->GetBinContent(qRatio_p->FindBin(refPt_[algoIter][jtIter]));
	      else if(qgPos[1] == 2) qgWeight *= gRatio_p->GetBinContent(gRatio_p->FindBin(refPt_[algoIter][jtIter]));
	    }
	    */
	    
	    
	    /*
	    if(isLeadDone){
	      isLeadDone = false;
	      break;
	    }
	    else isLeadDone = true;
	    */
	    /*
	    if(isLeadDone){
	      if(isSubleadDone){
		isLeadDone = false;
		isSubleadDone = false;
		break;
	      }
	      else isSubleadDone = true;
	    }
	    else isLeadDone = true;
	    */

	    /*	    
	    if(pthatLowerPos + 1 < nFiles){
	      if(genJtPt_[algoIter][jtIter] > pthatVals2[pthatLowerPos+1]) continue;
	    }
	    */

	    //	    if(genJtPt_[algoIter][jtIter] < pthatLowerCut*1.2 || genJtPt_[algoIter][jtIter] > pthatUpperCut*1.2) continue;
	    
	    
	    Int_t fillVal = 0;
	    Float_t tempDR = 999;
	    Int_t tempPos = -1;
	    
	    for(Int_t jtIter2 = 0; jtIter2 < nJt_[algoIter]; jtIter2++){
	      if(isUsedRecoJet[jtIter2]) continue;
	      
	      //	    if(isPbPb && !isGoodCaloJet(jetAlgoTemp[algoIter].c_str(), jtHCalSum_[algoIter][jtIter2], jtECalSum_[algoIter][jtIter2]) && !isGoodPFJet(jetAlgoTemp[algoIter].c_str(), jtChargedN_[algoIter][jtIter2], jtChargedSum_[algoIter][jtIter2], jtNeutralSum_[algoIter][jtIter2], jtRawPt_[algoIter][jtIter2])) continue;
	    
	      if(getDR(genJtEta_[algoIter][jtIter], genJtPhi_[algoIter][jtIter], jtEta_[algoIter][jtIter2], jtPhi_[algoIter][jtIter2]) < tempDR){
		tempPos = jtIter2;
		tempDR = getDR(genJtEta_[algoIter][jtIter], genJtPhi_[algoIter][jtIter], jtEta_[algoIter][jtIter2], jtPhi_[algoIter][jtIter2]);
	      }
	    }
	    
	    if(tempDR < 0.4){
	      fillVal = 1;
	      isUsedRecoJet[tempPos] = true;
	      genPosJet[tempPos] = jtIter;
	      recoPosJet[jtIter] = tempPos;
	    }
	    

	    for(Int_t jtIter2 = 0; jtIter2 < nJtPtBins; jtIter2++){
	      if(genJtPt_[algoIter][jtIter] > jtPtBins[jtIter2] && genJtPt_[algoIter][jtIter] < jtPtBins[jtIter2+1]){
		for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
		  Int_t tempFillVal = 1;
		  
		  if(recoPosJet[jtIter] < 0){
		    tempFillVal = 0;
		    if(genJtPt_[algoIter][jtIter] > 150 && recoCutIter == 0) std::cout << "missed: " << entry << ", " << genJtPt_[algoIter][jtIter] << ", " << jetAlgo.at(algoIter) << std::endl;
		  }
		  else if(jtPt_[algoIter][recoPosJet[jtIter]] < jetRecoCuts[recoCutIter]) tempFillVal = 0;
		  
		  jtEffVPt_MeanResPts_p[algoIter][centPos][0][recoCutIter][jtIter2]->Fill(tempFillVal, qgWeight*hibinWeight);
		}
		
		break;
	      }
	    }
	  }
	  
	  isSubleadDone = false;
	  isLeadDone = false;

	  pthat_Unweighted_p[algoIter]->Fill(ptHat_[algoIter]);
	  pthat_Weighted_p[algoIter]->Fill(ptHat_[algoIter], pthatWeight);
	  
	  pthat_Unweighted_Sample_p[algoIter][pthatIter]->Fill(ptHat_[algoIter]);

	  //MODDING FOR REFPT USAGE
	  for(Int_t jtIter = 0; jtIter < nJt_[algoIter]/*nGenJt_[algoIter]*/; jtIter++){
	    if(TMath::Abs(refEta_[algoIter][jtIter]) > 1.6) continue;
	    if(refPt_[algoIter][jtIter] < 5.0) continue;
	    if(refSubID_[algoIter][jtIter] != 0) continue;
	    /*
	    if((jtNeutralMax_[algoIter][jtIter]/jtRawPt_[algoIter][jtIter]*(0.085)+
	       jtPhotonMax_[algoIter][jtIter]/jtRawPt_[algoIter][jtIter]*(-0.337)+
	       jtChargedMax_[algoIter][jtIter]/jtRawPt_[algoIter][jtIter]*0.584+
	       jtNeutralSum_[algoIter][jtIter]/jtRawPt_[algoIter][jtIter]*(-0.454)+
	       jtPhotonSum_[algoIter][jtIter]/jtRawPt_[algoIter][jtIter]*(-0.127)+
	       jtChargedSum_[algoIter][jtIter]/jtRawPt_[algoIter][jtIter]*(-0.239)+
	       jtPU_[algoIter][jtIter]/jtRawPt_[algoIter][jtIter]*(-0.184)+0.173 )< -0.45) continue;
	    */

	    //	    if(jtNeutralN_[algoIter][jtIter] == 0 && jtPhotonHardN_[algoIter][jtIter] == 0 && jtTrackN_[algoIter][jtIter] == 1) continue;
	    
	    //	    if(TMath::Abs(jtTrackMax_[algoIter][jtIter]/jtPt_[algoIter][jtIter])>1.2) continue;
	    /*	    
	    Bool_t keepBool = false;
	    if(jtNeutralSum/jtRawPt < 0.9
               &&  jtChargedSum/jtRawPt > 0.01
               && (jtChargedN + jtPhotonN + jtNeutralN + eN + muN) > 0
               &&  jtChargedMax/jtRawPt < 0.99
               &&  jtPhotonSum/jtRawPt  < 0.99
               &&  eSum/jtRawPt       < 0.99) keepBool = true;

	    if(!keepBool) continue;
	    */
	    testFill1[algoIter][centPos]++;
	    
	    if(maxElePt > 10)
	      if(getDR(refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter], maxEleEta, maxElePhi) < 0.4) continue;
	    
	    if(twoElePt > 5)
	      if(getDR(refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter], twoEleEta, twoElePhi) < 0.4) continue;
	    
	    if(maxMuPt > 10)
	      if(getDR(refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter], maxMuEta, maxMuPhi) < 0.4) continue;
	    
	    if(twoMuPt > 5)
	      if(getDR(refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter], twoMuEta, twoMuPhi) < 0.4) continue;

	    if(maxPhoPt > 10)
	      if(getDR(refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter], maxPhoEta, maxPhoPhi) < 0.4) continue;
	    //	    testFill2[algoIter][centPos]++;



	    //	    if(TMath::Abs(getDPHI(zFinal.Phi(), refPhi_[algoIter][jtIter])) < 7*TMath::Pi()/8) continue;
	    
	    Int_t fillVal = 0;
	    Float_t tempDR = 999;
	    Int_t tempPos = -1;

	    Int_t qgPos[2] = {0, -1};
	    if(TMath::Abs(refPartFlav_[algoIter][jtIter/*tempPos*/]) < 9) qgPos[1] = 1;
	    else if(TMath::Abs(refPartFlav_[algoIter][jtIter/*tempPos*/]) == 21) qgPos[1] = 2;
	    

	    Float_t qgWeight = 1;
	    
	    /*
	    if((applyResidualMode || deriveResidualMode)){
	      if(qgPos[1] == 1) qgWeight *= qRatio_p->GetBinContent(qRatio_p->FindBin(refPt_[algoIter][jtIter]));
	      else if(qgPos[1] == 2) qgWeight *= gRatio_p->GetBinContent(gRatio_p->FindBin(refPt_[algoIter][jtIter]));
	    }
	    */
	    
	    if(refPt_[algoIter][jtIter] < 25 || refPt_[algoIter][jtIter]  > 100) testFill3[algoIter][centPos]++;
	    

	    //	    if(qgPos[1] == 1 && (applyResidualMode || deriveResidualMode)) qgWeight *= qgRatio_p->GetBinContent(qgRatio_p->FindBin(refPt_[algoIter][jtIter]));
	    //	    if(qgPos[1] == 1 && (applyResidualMode/* || deriveResidualMode*/)) qgWeight /= qgRatio_p->GetBinContent(qgRatio_p->FindBin(refPt_[algoIter][jtIter]));

	    //	    if(qgPos[1] != 1 && qgPos[1] != 2) continue;

	    for(Int_t qgIter = 0; qgIter < 2; qgIter++){
	      if(qgPos[qgIter] == -1) continue;
	      
	      if(true/*recoPosJet[jtIter] != -1*/){
		jtEta_p[algoIter][centPos][qgPos[qgIter]]->Fill(jtEta_[algoIter][jtIter]);
		
		jtRecoVGen_p[algoIter][centPos][qgPos[qgIter]]->Fill(/*genJtPt_*/refPt_[algoIter][jtIter], jtPt_[algoIter][jtIter/*recoPosJet[jtIter]*/]);
		
		jtRecoGenDRVPt_p[algoIter][centPos][qgPos[qgIter]]->Fill(/*genJtPt_*/refPt_[algoIter][jtIter], getDR(jtEta_[algoIter][jtIter/*recoPosJet[jtIter]*/], jtPhi_[algoIter][jtIter/*recoPosJet[jtIter]*/], /*genJtEta_*/refEta_[algoIter][jtIter], /*genJtPhi_*/refPhi_[algoIter][jtIter]));
		
		jtRecoOverGenVPt_p[algoIter][centPos][qgPos[qgIter]]->Fill(/*genJtPt_*/refPt_[algoIter][jtIter], jtPt_[algoIter][jtIter/*recoPosJet[jtIter]*/]/*genJtPt_*//refPt_[algoIter][jtIter]);
		jtRawOverGenVPt_p[algoIter][centPos][qgPos[qgIter]]->Fill(/*genJtPt_*/refPt_[algoIter][jtIter], jtRawPt_[algoIter][jtIter/*recoPosJet[jtIter]*/]/*genJtPt_*//refPt_[algoIter][jtIter]);
	      }
	      
	      for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
		Int_t tempFillVal = /*fillVal*/1;
		if(jtPt_[algoIter][jtIter/*recoPosJet[jtIter]*/] < jetRecoCuts[recoCutIter]) tempFillVal = 0;
		
		jtEffVPt_p[algoIter][centPos][qgPos[qgIter]][recoCutIter]->Fill(/*genJtPt_*/refPt_[algoIter][jtIter], tempFillVal);
		
		jtEffVPtEta_Denom_p[algoIter][centPos][qgPos[qgIter]][recoCutIter]->Fill(/*genJtEta_*/refEta_[algoIter][jtIter], /*genJtPt_*/refPt_[algoIter][jtIter]);
		if(tempFillVal) jtEffVPtEta_p[algoIter][centPos][qgPos[qgIter]][recoCutIter]->Fill(/*genJtEta_*/refEta_[algoIter][jtIter], /*genJtPt_*/refPt_[algoIter][jtIter]);
	      }	    
	      
	      for(Int_t jtIter2 = 0; jtIter2 < nJtPtBins2; jtIter2++){
		if(/*genJtPt_*/refPt_[algoIter][jtIter] > jtPtBins2[jtIter2] && /*genJtPt_*/refPt_[algoIter][jtIter] < jtPtBins2[jtIter2+1]){
		  for(Int_t etaIter2 = 0; etaIter2 < nJtEtaBins2; etaIter2++){
		    if(/*genJtEta_*/refEta_[algoIter][jtIter] > jtEtaBins2[etaIter2] && /*genJtEta_*/refEta_[algoIter][jtIter] < jtEtaBins2[etaIter2+1]){
		      jtRecoOverGenVCent_MeanResPts_p[algoIter][jtIter2][etaIter2][qgPos[qgIter]][centPos]->Fill(jtPt_[algoIter][jtIter/*recoPosJet[jtIter]*/]/*genJtPt_*//refPt_[algoIter][jtIter]);
		      
		      break;
		    }
		  }
		  break;
		}
	      }
	      
	      for(Int_t jtIter2 = 0; jtIter2 < nJtPtBins; jtIter2++){
		if(/*genJtPt_*/refPt_[algoIter][jtIter] > jtPtBins[jtIter2] && /*genJtPt_*/refPt_[algoIter][jtIter] < jtPtBins[jtIter2+1]){
		  
		  if(true/*recoPosJet[jtIter] != -1*/){
		    //		  if((entry == 0 || entry == 4 || entry == 5 || entry == 10 || entry == 15 || entry == 26 || entry == 31 || entry == 32 || entry == 37 || entry == 49 || entry == 54 || entry == 55) && algoIter == 0) std::cout << "entry, pt, refpt, phi, eta: " << entry << ", " << jtPt_[algoIter][jtIter] << ", " << refPt_[algoIter][jtIter] << ", " << jtPhi_[algoIter][jtIter] << ", " << jtEta_[algoIter][jtIter] << std::endl;
		    

		    Float_t tempDRFill = getDR(refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter], jtEta_[algoIter][jtIter], jtPhi_[algoIter][jtIter]);
		    Float_t tempDPhi = getDPHI(refPhi_[algoIter][jtIter], jtPhi_[algoIter][jtIter]);
		    Float_t tempDEta = refEta_[algoIter][jtIter] - jtEta_[algoIter][jtIter];
		    
		    jtPt_p[algoIter][centPos][qgPos[qgIter]]->Fill(jtPt_[algoIter][jtIter]);
		    genJtPt_p[algoIter][centPos][qgPos[qgIter]]->Fill(refPt_[algoIter][jtIter]);
		    

		    if(qgIter == 0) testFill2[algoIter][centPos]++;
		    hibinWeight = 1;
		    qgWeight = 1;
		    jtRecoOverGenVPt_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(jtPt_[algoIter][jtIter/*recoPosJet[jtIter]*/]/*genJtPt_*//refPt_[algoIter][jtIter], /*pthatWeight*/hibinWeight*qgWeight);


		    
		    //		    if(hiBin_ > 140) jtRecoOverGenVPt_PERP_MeanResPts_p[algoIter][qgPos[qgIter]][jtIter2]->Fill(jtPt_[algoIter][jtIter]/refPt_[algoIter][jtIter], hibinWeight*qgWeight);
		    
		    if(rhoPos != -1){
		      jtRecoOverGenVPt_Rho_MeanResPts_p[algoIter][rhoPos][qgPos[qgIter]][jtIter2]->Fill(jtPt_[algoIter][jtIter/*recoPosJet[jtIter]*/]/*genJtPt_*//refPt_[algoIter][jtIter]);
		    }
		    //		  else std::cout << "Out of range rho: " << rho_p->at(3) << std::endl;
		    
		    jtRawOverGenVPt_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(jtRawPt_[algoIter][jtIter/*recoPosJet[jtIter]*/]/*genJtPt_*//refPt_[algoIter][jtIter]);	
		    
		    jtRecoGenDRVPt_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(tempDRFill);
		    jtRecoGenDPhiVPt_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(tempDPhi);
		    jtRecoGenDEtaVPt_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(tempDEta);
		  }
		  /*
		    for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
		    Int_t tempFillVal = 1fillVal;
		    if(jtPt_[algoIter][jtIterrecoPosJet[jtIter]] < jetRecoCuts[recoCutIter]) tempFillVal = 0;
		    
		    jtEffVPt_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][recoCutIter][jtIter2]->Fill(tempFillVal);
		    }
		  */
		  break;
		}
	      }

	      if(hiBin_ > 140){
		for(Int_t jtIter2 = 0; jtIter2 < nJtPtBinsPERP; jtIter2++){
		  if(/*genJtPt_*/refPt_[algoIter][jtIter] > jtPtBinsPERP[jtIter2] && /*genJtPt_*/refPt_[algoIter][jtIter] < jtPtBinsPERP[jtIter2+1]){
		    jtRecoOverGenVPt_PERP_MeanResPts_p[algoIter][qgPos[qgIter]][jtIter2]->Fill(jtPt_[algoIter][jtIter]/refPt_[algoIter][jtIter], hibinWeight*qgWeight);
		    break;
		  }
		}
	      }

	      for(Int_t jtIter2 = 0; jtIter2 < nJtPtBins; jtIter2++){
		if(/*genJtPt_*/jtPt_[algoIter][jtIter] > jtPtBins[jtIter2] && /*genJtPt_*/jtPt_[algoIter][jtIter] < jtPtBins[jtIter2+1] && refPt_[algoIter][jtIter] > 0){
		  
		  if(true/*recoPosJet[jtIter] != -1*/){
		    jtRecoOverGenVRecoPt_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(jtPt_[algoIter][jtIter/*recoPosJet[jtIter]*/]/*genJtPt_*//refPt_[algoIter][jtIter], /*pthatWeight*/hibinWeight);
		  }
		  break;
		}
	      }
	    
	      //EDITING HERE
	      if(/*genJtPt_*/refPt_[algoIter][jtIter] > 30){
		if(true/*recoPosJet[jtIter] != -1*/){
		  
		  if(qgIter == 0) etaFills[algoIter][centPos]++;
		  
		  jtRecoOverGenVEta_p[algoIter][centPos][qgPos[qgIter]]->Fill(/*genJtEta_*/refEta_[algoIter][jtIter], jtPt_[algoIter][jtIter/*recoPosJet[jtIter]*/]/*genJtPt_*//refPt_[algoIter][jtIter]);
		  jtRawOverGenVEta_p[algoIter][centPos][qgPos[qgIter]]->Fill(/*genJtEta_*/refEta_[algoIter][jtIter], jtRawPt_[algoIter][jtIter/*recoPosJet[jtIter]*/]/*genJtPt_*//refPt_[algoIter][jtIter]);
		}	    
		
		for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
		  Int_t tempFillVal = 1/*fillVal*/;
		  if(jtPt_[algoIter][jtIter/*recoPosJet[jtIter]*/] < jetRecoCuts[recoCutIter]) tempFillVal = 0;
		  
		  jtEffVEta_p[algoIter][centPos][qgPos[qgIter]][recoCutIter]->Fill(/*genJtEta_*/refEta_[algoIter][jtIter], tempFillVal);
		}	      
		
		for(Int_t jtIter2 = 0; jtIter2 < nJtEtaBins; jtIter2++){
		  if(/*genJtEta_*/refEta_[algoIter][jtIter] > jtEtaBins[jtIter2] && /*genJtEta_*/refEta_[algoIter][jtIter] < jtEtaBins[jtIter2+1]){
		    if(true/*recoPosJet[jtIter] != -1*/){
		      
		      Float_t tempDRFill = getDR(refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter], jtEta_[algoIter][jtIter], jtPhi_[algoIter][jtIter]);
		      Float_t tempDPhi = getDPHI(refPhi_[algoIter][jtIter], jtPhi_[algoIter][jtIter]);
		      Float_t tempDEta = refEta_[algoIter][jtIter] - jtEta_[algoIter][jtIter];
		      
		      jtRecoOverGenVEta_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(jtPt_[algoIter][jtIter/*recoPosJet[jtIter]*/]/*genJtPt_*//refPt_[algoIter][jtIter]);
		      jtRawOverGenVEta_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(jtRawPt_[algoIter][jtIter/*recoPosJet[jtIter]*/]/*genJtPt_*//refPt_[algoIter][jtIter]);	
		      
		      jtRecoGenDRVEta_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(tempDRFill);
		      jtRecoGenDPhiVEta_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(tempDPhi);
		      jtRecoGenDEtaVEta_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(tempDEta);
		    }
		    
		    for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
		      Int_t tempFillVal = 1/*fillVal*/;
		      if(jtPt_[algoIter][jtIter/*recoPosJet[jtIter]*/] < jetRecoCuts[recoCutIter]) tempFillVal = 0;
		      
		      jtEffVEta_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][recoCutIter][jtIter2]->Fill(tempFillVal);
		    }
		    
		    break;
		  }
		}
	      }
	    }
	  }
	  
	  for(Int_t boolIter = 0; boolIter < boolSize; boolIter++){
	    if(TMath::Abs(jtEta_[algoIter][boolIter]) > 1.6) continue;
	    
	    if(jtPt_[algoIter][boolIter] < 5.0) continue;
	    
	    if(isPbPb && !isGoodCaloJet(jetAlgo.at(algoIter).c_str(), jtHCalSum_[algoIter][boolIter], jtECalSum_[algoIter][boolIter]) && !isGoodPFJet(jetAlgo.at(algoIter).c_str(), jtChargedN_[algoIter][boolIter], jtChargedSum_[algoIter][boolIter], jtNeutralSum_[algoIter][boolIter], jtRawPt_[algoIter][boolIter])) continue;
	    
	    Int_t fillVal = 1;
	    //	  if(isUsedRecoJet[boolIter]) fillVal = 0;
	    
	    Int_t qgPos[2] = {0, -1};
	    if(TMath::Abs(refPartFlav_[algoIter][boolIter]) < 9) qgPos[1] = 1;
	    else if(TMath::Abs(refPartFlav_[algoIter][boolIter]) == 21) qgPos[1] = 2;
	    
	    for(Int_t qgIter = 0; qgIter < 2; qgIter++){
	      if(qgPos[qgIter] == -1) continue;
	      
	      jtFakeVPt_p[algoIter][centPos][qgPos[qgIter]]->Fill(jtPt_[algoIter][boolIter], fillVal);
	      
	      for(Int_t jtIter2 = 0; jtIter2 < nJtPtBins; jtIter2++){
		if(jtPt_[algoIter][boolIter] > jtPtBins[jtIter2] && jtPt_[algoIter][boolIter] < jtPtBins[jtIter2+1]){
		  jtFakeVPt_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(fillVal);
		  break;
		}
	      }
	      
	      if(jtPt_[algoIter][boolIter] > 30){
		jtFakeVEta_p[algoIter][centPos][qgPos[qgIter]]->Fill(jtEta_[algoIter][boolIter], fillVal);
		
		for(Int_t jtIter2 = 0; jtIter2 < nJtEtaBins; jtIter2++){
		  if(jtEta_[algoIter][boolIter] > jtEtaBins[jtIter2] && jtEta_[algoIter][boolIter] < jtEtaBins[jtIter2+1]){
		    jtFakeVEta_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(fillVal);
		    break;
		  }
		}
	      }
	    }	   
	  }
      	}
      }
      inFile_p->Close();
    }
  }    

  myfile.close();

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    for(Int_t iter2 = 0; iter2 < nCentBins2; iter2++){
      std::cout << "Fills for " << jetAlgo.at(iter) << ", centrality " << centBins[iter2] << "-" << centBins[iter2+1]  << "%: "<< etaFills[iter][iter2] << std::endl;

      std::cout << testFill1[iter][iter2] << ", " << testFill2[iter][iter2] << ", " << testFill3[iter][iter2] << std::endl;

    }
  }

  std::cout << "#Events with no lead electron: " << nNoLeadEle << std::endl;
  std::cout << "#Events with no sublead electron: " << nNoSubleadEle << std::endl;

  std::cout << "#Events with no lead muon: " << nNoLeadMu << std::endl;
  std::cout << "#Events with no sublead muon: " << nNoSubleadMu << std::endl;

  std::cout << "#Events with no lead photon: " << nNoLeadPho << std::endl;

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    for(Int_t ptIter = 0; ptIter < nJtPtBins2; ptIter++){
      for(Int_t etaIter = 0; etaIter < nJtEtaBins2; etaIter++){
	for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	  for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
	    
	    Float_t tempMean[nMeanFit];
	    Float_t tempMeanErr[nMeanFit];
	    Float_t tempRes[nMeanFit];
	    Float_t tempResErr[nMeanFit];

	    tempMean[0] = jtRecoOverGenVCent_MeanResPts_p[iter][ptIter][etaIter][qgIter][centIter]->GetMean();
	    tempMeanErr[0] = jtRecoOverGenVCent_MeanResPts_p[iter][ptIter][etaIter][qgIter][centIter]->GetMeanError();
	    tempRes[0] = jtRecoOverGenVCent_MeanResPts_p[iter][ptIter][etaIter][qgIter][centIter]->GetStdDev();
	    tempResErr[0] = jtRecoOverGenVCent_MeanResPts_p[iter][ptIter][etaIter][qgIter][centIter]->GetStdDevError();

	    FitGauss(jtRecoOverGenVCent_MeanResPts_p[iter][ptIter][etaIter][qgIter][centIter], isPbPb, tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);

	    for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	      jtRecoOverGenVCent_Mean_p[iter][ptIter][etaIter][qgIter][mIter]->SetBinContent(nCentBins2 - centIter, tempMean[mIter]);
	      jtRecoOverGenVCent_Mean_p[iter][ptIter][etaIter][qgIter][mIter]->SetBinError(nCentBins2 - centIter, tempMeanErr[mIter]);
	      jtRecoOverGenVCent_Res_p[iter][ptIter][etaIter][qgIter][mIter]->SetBinContent(nCentBins2 - centIter, tempRes[mIter]);
	      jtRecoOverGenVCent_Res_p[iter][ptIter][etaIter][qgIter][mIter]->SetBinError(nCentBins2 - centIter, tempResErr[mIter]);
	    }
	    
	  }
	}
      }
    }

    if(debugMode) std::cout << __LINE__ << std::endl;	  

    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      for(Int_t qgIter = 0; qgIter < nQG; qgIter++){

	if(centIter == 0){
	  for(Int_t jtIter = 0; jtIter < nJtPtBinsPERP; jtIter++){

	    Float_t tempMean[nMeanFit];
	    Float_t tempMeanErr[nMeanFit];
	    Float_t tempRes[nMeanFit];
	    Float_t tempResErr[nMeanFit];

	    tempRes[0] = jtRecoOverGenVPt_PERP_MeanResPts_p[iter][qgIter][jtIter]->GetStdDev();
	    tempResErr[0] = jtRecoOverGenVPt_PERP_MeanResPts_p[iter][qgIter][jtIter]->GetStdDevError();
	  
	    FitGauss(jtRecoOverGenVPt_PERP_MeanResPts_p[iter][qgIter][jtIter], isPbPb, tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);
	
	    for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	      jtRecoOverGenVPt_PERP_Res_p[iter][qgIter][mIter]->SetBinContent(jtIter+1, tempRes[mIter]);
	      jtRecoOverGenVPt_PERP_Res_p[iter][qgIter][mIter]->SetBinError(jtIter+1, tempResErr[mIter]);
	    }
	  }
	}

	for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	  Float_t tempMean[nMeanFit];
	  Float_t tempMeanErr[nMeanFit];
	  Float_t tempRes[nMeanFit];
	  Float_t tempResErr[nMeanFit];

	  tempMean[0] = jtRecoGenDRVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRecoGenDRVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();	
	  tempRes[0] = jtRecoGenDRVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRecoGenDRVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();	

	  jtRecoGenDRVPt_Mean_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, tempMean[0]);
	  jtRecoGenDRVPt_Mean_p[iter][centIter][qgIter]->SetBinError(jtIter+1, tempMeanErr[0]);
	  jtRecoGenDRVPt_Res_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, tempRes[0]);
	  jtRecoGenDRVPt_Res_p[iter][centIter][qgIter]->SetBinError(jtIter+1, tempResErr[0]);

	  tempMean[0] = jtRecoGenDPhiVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRecoGenDPhiVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();
	  tempRes[0] = jtRecoGenDPhiVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRecoGenDPhiVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();

	  jtRecoGenDPhiVPt_Mean_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, tempMean[0]);
	  jtRecoGenDPhiVPt_Mean_p[iter][centIter][qgIter]->SetBinError(jtIter+1, tempMeanErr[0]);
	  jtRecoGenDPhiVPt_Res_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, tempRes[0]);
	  jtRecoGenDPhiVPt_Res_p[iter][centIter][qgIter]->SetBinError(jtIter+1, tempResErr[0]);

	  if(debugMode) std::cout << __LINE__ << std::endl;	  

	  tempMean[0] = jtRecoGenDEtaVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRecoGenDEtaVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();
	  tempRes[0] = jtRecoGenDEtaVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRecoGenDEtaVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();

	  jtRecoGenDEtaVPt_Mean_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, tempMean[0]);
	  jtRecoGenDEtaVPt_Mean_p[iter][centIter][qgIter]->SetBinError(jtIter+1, tempMeanErr[0]);
	  jtRecoGenDEtaVPt_Res_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, tempRes[0]);
	  jtRecoGenDEtaVPt_Res_p[iter][centIter][qgIter]->SetBinError(jtIter+1, tempResErr[0]);


	  tempMean[0] = jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();	
	  tempRes[0] = jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();
	  
	  FitGauss(jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter], isPbPb, tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);
	
	  for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	    jtRecoOverGenVPt_Mean_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempMean[mIter]);
	    jtRecoOverGenVPt_Mean_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempMeanErr[mIter]);
	    jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempRes[mIter]);
	    jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempResErr[mIter]);
	  }

	  tempMean[0] = jtRecoOverGenVRecoPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRecoOverGenVRecoPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();	
	  
	  FitGauss(jtRecoOverGenVRecoPt_MeanResPts_p[iter][centIter][qgIter][jtIter], isPbPb, tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);
	
	  for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	    jtRecoOverGenVRecoPt_Mean_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempMean[mIter]);
	    jtRecoOverGenVRecoPt_Mean_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempMeanErr[mIter]);
	  }


	  if(centIter == 0){
	    if(debugMode) std::cout << __LINE__ << std::endl;	  

	    for(Int_t rhoIter = 0; rhoIter < nRhoBins; rhoIter++){
	      tempMean[0] = jtRecoOverGenVPt_Rho_MeanResPts_p[iter][rhoIter][qgIter][jtIter]->GetMean();
	      tempMeanErr[0] = jtRecoOverGenVPt_Rho_MeanResPts_p[iter][rhoIter][qgIter][jtIter]->GetMeanError();	
	      tempRes[0] = jtRecoOverGenVPt_Rho_MeanResPts_p[iter][rhoIter][qgIter][jtIter]->GetStdDev();
	      tempResErr[0] = jtRecoOverGenVPt_Rho_MeanResPts_p[iter][rhoIter][qgIter][jtIter]->GetStdDevError();
	      
	      FitGauss(jtRecoOverGenVPt_Rho_MeanResPts_p[iter][rhoIter][qgIter][jtIter], isPbPb, tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);
	      
	      for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
		jtRecoOverGenVPt_Rho_Res_p[iter][rhoIter][qgIter][mIter]->SetBinContent(jtIter+1, tempRes[mIter]);
		jtRecoOverGenVPt_Rho_Res_p[iter][rhoIter][qgIter][mIter]->SetBinError(jtIter+1, tempResErr[mIter]);
	      }
	    }
	  }
	  
	  tempMean[0] = jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();
	  tempRes[0] = jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();

	  if(debugMode) std::cout << __LINE__ << std::endl;	  

	  FitGauss(jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter], isPbPb, tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);
	  
	  for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	    jtRecoOverRawVPt_Mean_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempMean[mIter]);
	    jtRecoOverRawVPt_Mean_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempMeanErr[mIter]);
	    jtRecoOverRawVPt_Res_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempRes[mIter]);
	    jtRecoOverRawVPt_Res_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempResErr[mIter]);
	  }
	  
	  tempMean[0] = jtRawOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRawOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();
	  tempRes[0] = jtRawOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRawOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();
	  
	  FitGauss(jtRawOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter], isPbPb, tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);

	  if(debugMode) std::cout << __LINE__ << std::endl;	  

	  for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	    jtRawOverGenVPt_Mean_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempMean[mIter]);
	    jtRawOverGenVPt_Mean_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempMeanErr[mIter]);
	    jtRawOverGenVPt_Res_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempRes[mIter]);
	    jtRawOverGenVPt_Res_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempResErr[mIter]);
	  }
	  
	  for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	    jtEffVPt_Mean_p[iter][centIter][qgIter][recoCutIter]->SetBinContent(jtIter+1, jtEffVPt_MeanResPts_p[iter][centIter][qgIter][recoCutIter][jtIter]->GetMean());
	    jtEffVPt_Mean_p[iter][centIter][qgIter][recoCutIter]->SetBinError(jtIter+1, jtEffVPt_MeanResPts_p[iter][centIter][qgIter][recoCutIter][jtIter]->GetMeanError());
	  }

	  jtFakeVPt_Mean_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, jtFakeVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean());
	  jtFakeVPt_Mean_p[iter][centIter][qgIter]->SetBinError(jtIter+1, jtFakeVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError());
	}

	for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	  Float_t tempMean[nMeanFit];
	  Float_t tempMeanErr[nMeanFit];
	  Float_t tempRes[nMeanFit];
	  Float_t tempResErr[nMeanFit];

	  if(debugMode) std::cout << __LINE__ << std::endl;	  


	  tempMean[0] = jtRecoGenDRVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRecoGenDRVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();
	  tempRes[0] = jtRecoGenDRVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRecoGenDRVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();

	  jtRecoGenDRVEta_Mean_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, tempMean[0]);
	  jtRecoGenDRVEta_Mean_p[iter][centIter][qgIter]->SetBinError(jtIter+1, tempMeanErr[0]);
	  jtRecoGenDRVEta_Res_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, tempRes[0]);
	  jtRecoGenDRVEta_Res_p[iter][centIter][qgIter]->SetBinError(jtIter+1, tempResErr[0]);

	  tempMean[0] = jtRecoGenDPhiVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRecoGenDPhiVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();	
	  tempRes[0] = jtRecoGenDPhiVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRecoGenDPhiVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();

	  jtRecoGenDPhiVEta_Mean_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, tempMean[0]);
	  jtRecoGenDPhiVEta_Mean_p[iter][centIter][qgIter]->SetBinError(jtIter+1, tempMeanErr[0]);
	  jtRecoGenDPhiVEta_Res_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, tempRes[0]);
	  jtRecoGenDPhiVEta_Res_p[iter][centIter][qgIter]->SetBinError(jtIter+1, tempResErr[0]);

	  tempMean[0] = jtRecoGenDEtaVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRecoGenDEtaVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();	

	  if(debugMode) std::cout << __LINE__ << std::endl;	  

	  tempRes[0] = jtRecoGenDEtaVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRecoGenDEtaVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();

	  jtRecoGenDEtaVEta_Mean_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, tempMean[0]);
	  jtRecoGenDEtaVEta_Mean_p[iter][centIter][qgIter]->SetBinError(jtIter+1, tempMeanErr[0]);
	  jtRecoGenDEtaVEta_Res_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, tempRes[0]);
	  jtRecoGenDEtaVEta_Res_p[iter][centIter][qgIter]->SetBinError(jtIter+1, tempResErr[0]);


	  tempMean[0] = jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();
	  tempRes[0] = jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();
	  
	  FitGauss(jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter], isPbPb, tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);

	  if(debugMode) std::cout << __LINE__ << std::endl;	  
	  
	  for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	    jtRecoOverGenVEta_Mean_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempMean[mIter]);
	    jtRecoOverGenVEta_Mean_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempMeanErr[mIter]);
	    jtRecoOverGenVEta_Res_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempRes[mIter]);
	    jtRecoOverGenVEta_Res_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempResErr[mIter]);
	  }
	  
	  tempMean[0] = jtRecoOverRawVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRecoOverRawVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();
	  tempRes[0] = jtRecoOverRawVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRecoOverRawVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();
	  
	  FitGauss(jtRecoOverRawVEta_MeanResPts_p[iter][centIter][qgIter][jtIter], isPbPb, tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);
	  
	  for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	    jtRecoOverRawVEta_Mean_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempMean[mIter]);
	    jtRecoOverRawVEta_Mean_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempMeanErr[mIter]);
	    jtRecoOverRawVEta_Res_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempRes[mIter]);
	    jtRecoOverRawVEta_Res_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempResErr[mIter]);
	  }
	  
	  tempMean[0] = jtRawOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRawOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();
	  tempRes[0] = jtRawOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRawOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();
	  
	  FitGauss(jtRawOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter], isPbPb, tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);

	  for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	    jtRawOverGenVEta_Mean_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempMean[mIter]);
	    jtRawOverGenVEta_Mean_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempMeanErr[mIter]);
	    jtRawOverGenVEta_Res_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempRes[mIter]);
	    jtRawOverGenVEta_Res_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempResErr[mIter]);
	  }
	  
	  for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	    jtEffVEta_Mean_p[iter][centIter][qgIter][recoCutIter]->SetBinContent(jtIter+1, jtEffVEta_MeanResPts_p[iter][centIter][qgIter][recoCutIter][jtIter]->GetMean());
	    jtEffVEta_Mean_p[iter][centIter][qgIter][recoCutIter]->SetBinError(jtIter+1, jtEffVEta_MeanResPts_p[iter][centIter][qgIter][recoCutIter][jtIter]->GetMeanError());
	  }

	  jtFakeVEta_Mean_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, jtFakeVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean());
	  jtFakeVEta_Mean_p[iter][centIter][qgIter]->SetBinError(jtIter+1, jtFakeVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError());
	}
      }
    }
  }

  outFile_p->cd();

  for(Int_t iter = 0; iter < nFiles; iter++){
    if(isFile[iter]){
      TNamed pathStr(Form("pathStr_%s", inFilePtHats[iter].c_str()), inFileNames[iter].c_str());
      pathStr.Write("", TObject::kOverwrite);
    }
  }

  localHibin_p->Write("", TObject::kOverwrite);
  hibinWeighted_p->Write("", TObject::kOverwrite);
  centVRho_p->Write("", TObject::kOverwrite);


  if(debugMode) std::cout << __LINE__ << std::endl;

  for(Int_t rhoIter = 0; rhoIter < nRhoCentBins; rhoIter++){
    centVRho_Proj_p[rhoIter]->Write("", TObject::kOverwrite);
  }

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    TDirectory* dir_p = outFile_p->GetDirectory(Form("%s", jetAlgo.at(iter).c_str()));
    if(dir_p){
      dir_p->cd();
    }
    else{
      dir_p = outFile_p->mkdir(Form("%s", jetAlgo.at(iter).c_str()));
      dir_p->cd();
    }

    pthat_Unweighted_p[iter]->Write("", TObject::kOverwrite);
    pthat_Weighted_p[iter]->Write("", TObject::kOverwrite);

    for(Int_t pthatIter = 0; pthatIter < nFiles; pthatIter++){
      pthat_Unweighted_Sample_p[iter][pthatIter]->Write("", TObject::kOverwrite);
      jtPerPthat_p[iter][pthatIter]->Write("", TObject::kOverwrite);
    }

    for(Int_t ptIter = 0; ptIter < nJtPtBins2; ptIter++){
      for(Int_t etaIter = 0; etaIter < nJtEtaBins2; etaIter++){
        for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
          for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
            jtRecoOverGenVCent_Mean_p[iter][ptIter][etaIter][qgIter][mIter]->Write("", TObject::kOverwrite);
            jtRecoOverGenVCent_Res_p[iter][ptIter][etaIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  }

	  for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
	    jtRecoOverGenVCent_MeanResPts_p[iter][ptIter][etaIter][qgIter][centIter]->Write("", TObject::kOverwrite);
	  }
	}
      }
    }


    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      bkgEstimate_p[iter][centIter]->Write("", TObject::kOverwrite);
      areaEstimate_p[iter][centIter]->Write("", TObject::kOverwrite);
      bkgAreaEstimate_p[iter][centIter]->Write("", TObject::kOverwrite);

      for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	jtPt_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	genJtPt_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtEta_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);

	jtRecoVGen_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoGenDRVPt_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoOverGenVPt_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoOverRawVPt_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRawOverGenVPt_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	
	jtRecoGenDRVPt_Mean_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoGenDEtaVPt_Mean_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoGenDPhiVPt_Mean_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);

	jtRecoGenDRVPt_Res_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoGenDEtaVPt_Res_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	
	if(debugMode) std::cout << __LINE__ << std::endl;

	FitCSNSimple(jtRecoGenDPhiVPt_Res_p[iter][centIter][qgIter]);
	dir_p->cd();
	jtRecoGenDPhiVPt_Res_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);

	jtRecoGenDRVEta_Mean_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoGenDEtaVEta_Mean_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoGenDPhiVEta_Mean_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);

	jtRecoGenDRVEta_Res_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoGenDEtaVEta_Res_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoGenDPhiVEta_Res_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);

	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  jtRecoOverGenVPt_Mean_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRecoOverGenVRecoPt_Mean_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);

	  if(centIter == 0){

	    if(debugMode) std::cout << __LINE__ << std::endl;

	    FitCSN2(jtRecoOverGenVPt_PERP_Res_p[iter][qgIter][mIter], 1, jtRecoOverGenVPt_PERP_Res_p[iter][qgIter][mIter]);
	    if(debugMode) std::cout << __LINE__ << std::endl;

	    dir_p->cd();
	    if(debugMode) std::cout << __LINE__ << std::endl;

	    jtRecoOverGenVPt_PERP_Res_p[iter][qgIter][mIter]->Write("", TObject::kOverwrite);
	
	    for(Int_t rhoIter = 0; rhoIter < nRhoBins; rhoIter++){
	      jtRecoOverGenVPt_Rho_Res_p[iter][rhoIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	    }
	  }

	  if(isPbPb) FitCSN2(jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter], 0, jtRecoOverGenVPt_PERP_Res_p[iter][qgIter][mIter]);
	  else FitCSN2(jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter], 1, jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter]);
	    
	  dir_p->cd();
	  jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);


	  jtRecoOverRawVPt_Mean_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRecoOverRawVPt_Res_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRawOverGenVPt_Mean_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRawOverGenVPt_Res_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	}
	
	for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	  jtEffVPt_p[iter][centIter][qgIter][recoCutIter]->Write("", TObject::kOverwrite);
	  FitErf(jtEffVPt_Mean_p[iter][centIter][qgIter][recoCutIter]);
	  dir_p->cd();
	  jtEffVPt_Mean_p[iter][centIter][qgIter][recoCutIter]->Write("", TObject::kOverwrite);
	}
	jtFakeVPt_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtFakeVPt_Mean_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	
      //      jtEffVPtEta_p[iter][centIter][qgIter]->Write(Form("%s_NUM",jtEffVPtEta_p[iter][centIter][qgIter]->GetName()), TObject::kOverwrite);
      //      jtEffVPtEta_Denom_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	  jtEffVPtEta_p[iter][centIter][qgIter][recoCutIter]->Divide(jtEffVPtEta_Denom_p[iter][centIter][qgIter][recoCutIter]);
	  jtEffVPtEta_p[iter][centIter][qgIter][recoCutIter]->Write("", TObject::kOverwrite);
	}

	if(centIter == 0){
	  for(Int_t jtIter = 0; jtIter < nJtPtBinsPERP; jtIter++){
	    jtRecoOverGenVPt_PERP_MeanResPts_p[iter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  }
	}

	for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	  jtRecoGenDRVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  jtRecoGenDPhiVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  jtRecoGenDEtaVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);

	  jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  jtRecoOverGenVRecoPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  if(centIter == 0){
	    //	    jtRecoOverGenVPt_PERP_MeanResPts_p[iter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	    for(Int_t rhoIter = 0; rhoIter < nRhoBins; rhoIter++){
	      jtRecoOverGenVPt_Rho_MeanResPts_p[iter][rhoIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	    }
	  }
	  jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  jtRawOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);

          for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	    jtEffVPt_MeanResPts_p[iter][centIter][qgIter][recoCutIter][jtIter]->Write("", TObject::kOverwrite);
	  }
	  jtFakeVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	}

	jtRecoOverGenVEta_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoOverRawVEta_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRawOverGenVEta_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	
	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  jtRecoOverGenVEta_Mean_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRecoOverGenVEta_Res_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRecoOverRawVEta_Mean_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRecoOverRawVEta_Res_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRawOverGenVEta_Mean_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRawOverGenVEta_Res_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	}
	
	for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	  jtEffVEta_p[iter][centIter][qgIter][recoCutIter]->Write("", TObject::kOverwrite);
	  jtEffVEta_Mean_p[iter][centIter][qgIter][recoCutIter]->Write("", TObject::kOverwrite);
	}
	jtFakeVEta_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtFakeVEta_Mean_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	
	for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	  jtRecoGenDRVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  jtRecoGenDPhiVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  jtRecoGenDEtaVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);

	  jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  jtRecoOverRawVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  jtRawOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);

          for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	    jtEffVEta_MeanResPts_p[iter][centIter][qgIter][recoCutIter][jtIter]->Write("", TObject::kOverwrite);
	  }
	  
	  jtFakeVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	}
      }
    }
  }

  delete centVRho_p;

  for(Int_t rhoIter = 0; rhoIter < nRhoCentBins; rhoIter++){
    delete centVRho_Proj_p[rhoIter];
  }


  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    delete pthat_Unweighted_p[iter];
    delete pthat_Weighted_p[iter];

    for(Int_t pthatIter = 0; pthatIter < nFiles; pthatIter++){
      delete pthat_Unweighted_Sample_p[iter][pthatIter];

      delete jtPerPthat_p[iter][pthatIter];
    }

    for(Int_t ptIter = 0; ptIter < nJtPtBins2; ptIter++){
      for(Int_t etaIter = 0; etaIter < nJtEtaBins2; etaIter++){
        for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
          for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
            delete jtRecoOverGenVCent_Mean_p[iter][ptIter][etaIter][qgIter][mIter];
            delete jtRecoOverGenVCent_Res_p[iter][ptIter][etaIter][qgIter][mIter];
          }

          for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
            delete jtRecoOverGenVCent_MeanResPts_p[iter][ptIter][etaIter][qgIter][centIter];
          }
        }
      }
    }

    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      delete bkgEstimate_p[iter][centIter];
      delete areaEstimate_p[iter][centIter];
      delete bkgAreaEstimate_p[iter][centIter];

      for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	delete jtEta_p[iter][centIter][qgIter];

	delete jtRecoVGen_p[iter][centIter][qgIter];
	delete jtRecoGenDRVPt_p[iter][centIter][qgIter];
	delete jtRecoOverGenVPt_p[iter][centIter][qgIter];
	delete jtRecoOverRawVPt_p[iter][centIter][qgIter];
	delete jtRawOverGenVPt_p[iter][centIter][qgIter];
	
	delete jtRecoGenDRVPt_Mean_p[iter][centIter][qgIter];
	delete jtRecoGenDPhiVPt_Mean_p[iter][centIter][qgIter];
	delete jtRecoGenDEtaVPt_Mean_p[iter][centIter][qgIter];

	delete jtRecoGenDRVPt_Res_p[iter][centIter][qgIter];
	delete jtRecoGenDPhiVPt_Res_p[iter][centIter][qgIter];
	delete jtRecoGenDEtaVPt_Res_p[iter][centIter][qgIter];

	delete jtRecoGenDRVEta_Mean_p[iter][centIter][qgIter];
	delete jtRecoGenDPhiVEta_Mean_p[iter][centIter][qgIter];
	delete jtRecoGenDEtaVEta_Mean_p[iter][centIter][qgIter];

	delete jtRecoGenDRVEta_Res_p[iter][centIter][qgIter];
	delete jtRecoGenDPhiVEta_Res_p[iter][centIter][qgIter];
	delete jtRecoGenDEtaVEta_Res_p[iter][centIter][qgIter];

	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  delete jtRecoOverGenVPt_Mean_p[iter][centIter][qgIter][mIter];
	  delete jtRecoOverGenVRecoPt_Mean_p[iter][centIter][qgIter][mIter];
	  delete jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter];

	  if(centIter == 0){
	    delete jtRecoOverGenVPt_PERP_Res_p[iter][qgIter][mIter];

	    for(Int_t rhoIter = 0; rhoIter < nRhoBins; rhoIter++){
	      delete jtRecoOverGenVPt_Rho_Res_p[iter][rhoIter][qgIter][mIter];
	    }
	  }

	  delete jtRecoOverRawVPt_Mean_p[iter][centIter][qgIter][mIter];
	  delete jtRecoOverRawVPt_Res_p[iter][centIter][qgIter][mIter];
	  delete jtRawOverGenVPt_Mean_p[iter][centIter][qgIter][mIter];
	  delete jtRawOverGenVPt_Res_p[iter][centIter][qgIter][mIter];
	}
	
	for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	  delete jtEffVPt_p[iter][centIter][qgIter][recoCutIter];
	  delete jtEffVPt_Mean_p[iter][centIter][qgIter][recoCutIter];
	}

	delete jtFakeVPt_p[iter][centIter][qgIter];
	delete jtFakeVPt_Mean_p[iter][centIter][qgIter];

	for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){	
	  delete jtEffVPtEta_p[iter][centIter][qgIter][recoCutIter];
	  delete jtEffVPtEta_Denom_p[iter][centIter][qgIter][recoCutIter];
	}


	if(centIter== 0){
	  for(Int_t jtIter = 0; jtIter < nJtPtBinsPERP; jtIter++){
	    delete jtRecoOverGenVPt_PERP_MeanResPts_p[iter][qgIter][jtIter];
	  }
	}

	for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	  delete jtRecoGenDRVPt_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtRecoGenDPhiVPt_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtRecoGenDEtaVPt_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtRecoOverGenVRecoPt_MeanResPts_p[iter][centIter][qgIter][jtIter];

	  if(centIter == 0){
	    for(Int_t rhoIter = 0; rhoIter < nRhoBins; rhoIter++){
	      delete jtRecoOverGenVPt_Rho_MeanResPts_p[iter][rhoIter][qgIter][jtIter];
	    } 
	  }

	  delete jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtRawOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter];

          for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	    delete jtEffVPt_MeanResPts_p[iter][centIter][qgIter][recoCutIter][jtIter];
	  }

	  delete jtFakeVPt_MeanResPts_p[iter][centIter][qgIter][jtIter];
	}
	
	delete jtRecoOverGenVEta_p[iter][centIter][qgIter];
	delete jtRecoOverRawVEta_p[iter][centIter][qgIter];
	delete jtRawOverGenVEta_p[iter][centIter][qgIter];
	
	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  delete jtRecoOverGenVEta_Mean_p[iter][centIter][qgIter][mIter];
	  delete jtRecoOverGenVEta_Res_p[iter][centIter][qgIter][mIter];
	  delete jtRecoOverRawVEta_Mean_p[iter][centIter][qgIter][mIter];
	  delete jtRecoOverRawVEta_Res_p[iter][centIter][qgIter][mIter];
	  delete jtRawOverGenVEta_Mean_p[iter][centIter][qgIter][mIter];
	  delete jtRawOverGenVEta_Res_p[iter][centIter][qgIter][mIter];
	}
	
	for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	  delete jtEffVEta_p[iter][centIter][qgIter][recoCutIter];
	  delete jtEffVEta_Mean_p[iter][centIter][qgIter][recoCutIter];
	}
	delete jtFakeVEta_p[iter][centIter][qgIter];
	delete jtFakeVEta_Mean_p[iter][centIter][qgIter];
	
	for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	  delete jtRecoGenDRVEta_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtRecoGenDPhiVEta_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtRecoGenDEtaVEta_MeanResPts_p[iter][centIter][qgIter][jtIter];


	  delete jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtRecoOverRawVEta_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtRawOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter];

          for(Int_t recoCutIter = 0; recoCutIter < nJetRecoCuts; recoCutIter++){
	    delete jtEffVEta_MeanResPts_p[iter][centIter][qgIter][recoCutIter][jtIter];
	  }

	  delete jtFakeVEta_MeanResPts_p[iter][centIter][qgIter][jtIter];
	}
      }
    }
  }
    
  outFile_p->Close();
  delete outFile_p;

  delete date;

  return;
}
