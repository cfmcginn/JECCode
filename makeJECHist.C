#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TDatime.h"

#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include "getLogBins.h"
#include "getLinBins.h"
#include "etaPhiFunc.h"
#include "getBkgEstimate.h"

#include "dirent.h"

const Bool_t doGetBkg = false;

//const Int_t nJetAlgo = 6;
//const std::string jetAlgoTemp[nJetAlgo] = {"akVs4Calo", "akPu4Calo", "akVs4PF", "akPu4PF", "akVs3PF", "akPu3PF"};
const Int_t nJetAlgo = 1;
const std::string jetAlgoTemp[nJetAlgo] = {"akCs4PF"};
//const std::string jetAlgoTemp[nJetAlgo] = {"ak4Calo", "ak3Calo", "ak4PF", "ak3PF"};
//const Int_t jetAlgoR[nJetAlgo] = {4, 3, 4, 3};
//const Int_t jetAlgoR[nJetAlgo] = {4, 4, 4, 4, 3, 3};
const Int_t jetAlgoR[nJetAlgo] = {4};

const Int_t nJtCat = 3;
const std::string jtCat[nJtCat] = {"Eta2", "Eta1", "Dijet"};

const Int_t nMeanFit = 2;
const std::string meanFit[nMeanFit] = {"", "Fit"};

const Int_t nQG = 3;
const std::string qg[nQG] = {"Inc", "Q", "G"};

const Int_t nMaxJets = 500;

const Int_t nCentBins = 4;
const Int_t centBins[nCentBins+1] = {200, 100, 60, 20, 0};
const Float_t centBins2[nCentBins+1] = {0.001, 10, 30, 50, 99.999};
//const Int_t nCentBins = 8;
//const Int_t centBins[nCentBins+1] = {200, 140, 100, 80, 60, 40, 20, 10, 0};

const Int_t nMaxGen = 100000;

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


void FitGauss(TH1F* hist_p, Float_t& mean, Float_t& meanErr, Float_t& res, Float_t& resErr)
{
  TF1* f1_p = new TF1("f1_p", "gaus", hist_p->GetXaxis()->GetXmin(), hist_p->GetXaxis()->GetXmax());

  hist_p->Fit("f1_p", "LL Q M");

  mean = f1_p->GetParameter(1);
  res = f1_p->GetParameter(2);

  meanErr = f1_p->GetParError(1);
  resErr = f1_p->GetParError(2);

  //  if(TMath::Abs(mean - 1.0) < 0.01) return;
  if(f1_p->GetProb() > .01) return;

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

    if(hist_p->Integral(tempBin, maxBin+iter) > .90*hist_p->Integral() || tempBin == 1){
      fitLow = hist_p->GetBinCenter(tempBin);
      fitHi = hist_p->GetBinCenter(maxBin+iter);
      break;
    }
  }

  if(nBins >= 7){
    hist_p->Fit("f1_p", "LL Q M", "", fitLow, fitHi);
    
    mean = f1_p->GetParameter(1);
    res = f1_p->GetParameter(2);
    
    meanErr = f1_p->GetParError(1);
    resErr = f1_p->GetParError(2);
    
    //  if(TMath::Abs(mean - 1.0) < 0.01) return;
    if(f1_p->GetProb() > .01) return;
  }
  nBins = 1;

  Int_t meanBin = hist_p->FindBin(hist_p->GetMean());
  fitLow = -1;
  fitHi = -1;

  for(Int_t iter = 1; iter < hist_p->GetNbinsX(); iter++){
    Int_t tempBin = meanBin - iter;
    if(tempBin < 1) tempBin = 1;
    
    nBins += 2;

    if(hist_p->Integral(tempBin, meanBin+iter) > .90*hist_p->Integral() || tempBin == 1){
      fitLow = hist_p->GetBinCenter(tempBin);
      fitHi = hist_p->GetBinCenter(meanBin+iter);
      break;
    }
  }


  if(nBins >= 7){
    hist_p->Fit("f1_p", "LL Q M", "", fitLow, fitHi);
    
    mean = f1_p->GetParameter(1);
    res = f1_p->GetParameter(2);
    
    meanErr = f1_p->GetParError(1);
    resErr = f1_p->GetParError(2);
  }

  delete f1_p;

  return;
}


void makeJECHist(const std::string inFileName, const Bool_t isPbPb)
{
  std::vector<std::string>* listOfFiles_p = new std::vector<std::string>;

  DIR *dpdf;
  struct dirent *epdf;

  if(!strcmp(&(inFileName.back()), "/")){
    dpdf = opendir(inFileName.c_str());
    
    if(dpdf != NULL){
      while(epdf = readdir(dpdf)){
	TString temp = epdf->d_name;
	
	if(temp.Index("HiForest") >= 0) listOfFiles_p->push_back(Form("%s%s", inFileName.c_str(), temp.Data()));
      }
    }
    else{
      std::cout << "NULL PATH" << std::endl;
      return;
    }
  }
  else listOfFiles_p->push_back(inFileName);

  std::string outName = listOfFiles_p->at(0).c_str();
  const std::string inString = ".root";
  const std::string inString2 = "*";
  TDatime* date = new TDatime();
  const std::string outString = Form("_%d_HIST.root", date->GetDate());
  std::size_t strIndex = 0;

  std::string tempOutName = outName;
  Int_t outNum = 0;
  Int_t outIter = 0;

  while(tempOutName.find("/") != std::string::npos){
    outNum++;
    tempOutName.replace(0, tempOutName.find("/")+1, "");
  }

  while(outName.find("/") != std::string::npos){
    if(outIter < outNum-2){
      outIter++;
      outName.replace(0, outName.find("/")+1, "");
    }
    else outName.replace(outName.find("/"), 1, "_");
  }

  strIndex = outName.find(inString);
  if(!(strIndex == std::string::npos)){
    outName.replace(strIndex, inString.length(), outString); 
  }

  strIndex = outName.find(inString2);
  if(!(strIndex == std::string::npos)){
    outName.replace(strIndex, inString2.length(), outString);
  }

  TFile* outFile_p = new TFile(outName.c_str(), "UPDATE");
  std::string jetAlgo[nJetAlgo];
  
  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    jetAlgo[iter] = jetAlgoTemp[iter];
    /*
      TString tempjt = jetAlgoTemp[iter].c_str();
      Int_t pos = tempjt.Index("Vs");
      
      if(!isPbPb && pos >= 0){
      jetAlgo[iter] = Form("%s%s", jetAlgoTemp[iter].substr(0, pos).c_str(), jetAlgoTemp[iter].substr(pos+2, jetAlgoTemp[iter].length() - pos - 2).c_str());
      }
      else jetAlgo[iter] = jetAlgoTemp[iter];
    */
  }
  
  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    std::cout << "Algos: " << jetAlgo[iter] << std::endl;
  }

  Int_t hiBin_ = -1;
  Float_t vz_;

  UInt_t run_;
  ULong64_t evt_;

  std::vector<float>* genPt_p = 0;
  std::vector<float>* genEta_p = 0;
  std::vector<float>* genPhi_p = 0;
  std::vector<int>* genPDG_p = 0;

  Int_t nJt_[nJetAlgo]; 
  Float_t jtPt_[nJetAlgo][nMaxJets];
  Float_t jtRawPt_[nJetAlgo][nMaxJets];
  Float_t jtEta_[nJetAlgo][nMaxJets];
  Float_t jtPhi_[nJetAlgo][nMaxJets];
  Float_t jtHCalSum_[nJetAlgo][nMaxJets];
  Float_t jtECalSum_[nJetAlgo][nMaxJets];
  Float_t jtNeutralSum_[nJetAlgo][nMaxJets];
  Float_t jtChargedSum_[nJetAlgo][nMaxJets];
  Int_t jtChargedN_[nJetAlgo][nMaxJets];
  Float_t refPt_[nJetAlgo][nMaxJets];
  Float_t refEta_[nJetAlgo][nMaxJets];
  Int_t refSubID_[nJetAlgo][nMaxJets];
  Int_t refPartFlav_[nJetAlgo][nMaxJets];

  Int_t nGenJt_[nJetAlgo]; 
  Float_t genJtPt_[nJetAlgo][nMaxJets];
  Float_t genJtEta_[nJetAlgo][nMaxJets];
  Float_t genJtPhi_[nJetAlgo][nMaxJets];
  Int_t genJtMatchIndex_[nJetAlgo][nMaxJets];
  Int_t genSubId_[nJetAlgo][nMaxJets];

  const Int_t nJtPtBins = 29;
  const Float_t jtPtLow = 7.999;
  const Float_t jtPtHi = 300.001;
  Double_t jtPtBins[nJtPtBins+1];
  getLogBins(jtPtLow, jtPtHi, nJtPtBins, jtPtBins);

  const Int_t nJtEtaBins = 20;
  const Float_t jtEtaLow = -2.0;
  const Float_t jtEtaHi = 2.0;
  Double_t jtEtaBins[nJtEtaBins+1];
  getLinBins(jtEtaLow, jtEtaHi, nJtEtaBins, jtEtaBins);

  const Int_t nJtPtBins2 = 5;
  const Float_t jtPtLow2 = 30;
  const Float_t jtPtHi2 = 300;
  Double_t jtPtBins2[nJtPtBins2+1];
  getLogBins(jtPtLow2, jtPtHi2, nJtPtBins2, jtPtBins2);

  const Int_t nJtEtaBins2 = 3;
  const Float_t jtEtaLow2 = -2.0;
  const Float_t jtEtaHi2 = 2.0;
  Double_t jtEtaBins2[nJtEtaBins2+1] = {-2.0, -1.0, 1.0, 2.0};
  //  getLinBins(jtEtaLow2, jtEtaHi2, nJtEtaBins2, jtEtaBins2);

  for(Int_t iter = 0; iter < nJtPtBins+1; iter++){
    std::cout << "Pt Bin, val: " << iter << ", " << jtPtBins[iter] << std::endl;
  }

  for(Int_t iter = 0; iter < nJtEtaBins+1; iter++){
    std::cout << "Eta Bin, val: " << iter << ", " << jtEtaBins[iter] << std::endl;
  }

  Int_t nCentBinsTemp = nCentBins;
  if(!isPbPb) nCentBinsTemp = 1;

  const Int_t nCentBins2 = nCentBinsTemp;

  TH1F* bkgEstimate_p[nJetAlgo][nCentBins2];
  TH1F* areaEstimate_p[nJetAlgo][nCentBins2];
  TH2F* bkgAreaEstimate_p[nJetAlgo][nCentBins2];

  TH2F* jtRecoGenDRVPt_p[nJetAlgo][nCentBins2][nQG];
  TH2F* jtRecoVGen_p[nJetAlgo][nCentBins2][nQG];
  TH2F* jtRecoOverGenVPt_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoOverGenVPt_Mean_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRecoOverGenVPt_Res_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRecoOverGenVPt_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];

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

  TH2F* jtEffVPt_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtEffVPt_Mean_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtEffVPt_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];

  TH2F* jtFakeVPt_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtFakeVPt_Mean_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtFakeVPt_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];

  TH2F* jtEffVPtEta_p[nJetAlgo][nCentBins2][nQG];
  TH2F* jtEffVPtEta_Denom_p[nJetAlgo][nCentBins2][nQG];

  TH2F* jtRecoOverGenVEta_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoOverGenVEta_Mean_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRecoOverGenVEta_Res_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRecoOverGenVEta_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];

  TH2F* jtRecoOverRawVEta_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRecoOverRawVEta_Mean_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRecoOverRawVEta_Res_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRecoOverRawVEta_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];

  TH2F* jtRawOverGenVEta_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtRawOverGenVEta_Mean_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRawOverGenVEta_Res_p[nJetAlgo][nCentBins2][nQG][nMeanFit];
  TH1F* jtRawOverGenVEta_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];

  TH2F* jtEffVEta_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtEffVEta_Mean_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtEffVEta_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];

  TH2F* jtFakeVEta_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtFakeVEta_Mean_p[nJetAlgo][nCentBins2][nQG];
  TH1F* jtFakeVEta_MeanResPts_p[nJetAlgo][nCentBins2][nQG][nJtPtBins];

  Int_t nNoLeadEle = 0;
  Int_t nNoSubleadEle = 0;

  for(Int_t iter = 0; iter < nJetAlgo; iter++){

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
	    jtRecoOverGenVCent_Mean_p[iter][ptIter][etaIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVCent_%s_%sMean_%s_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo[iter].c_str(), ptStr1.c_str(), etaStr1.c_str()), Form(";Centrality;#mu_{Reco./Gen.} (%s)", jetAlgo[iter].c_str()), nCentBins2, centBins2);
	    jtRecoOverGenVCent_Res_p[iter][ptIter][etaIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVCent_%s_%sRes_%s_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo[iter].c_str(), ptStr1.c_str(), etaStr1.c_str()), Form(";Centrality;#sigma_{Reco./Gen.} (%s)", jetAlgo[iter].c_str()), nCentBins2, centBins2);
	  }

          for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
	    std::string centStr = "PP";
	    if(isPbPb) centStr = Form("cent%dto%d", centBins[centIter+1]/2, centBins[centIter]/2);

	    jtRecoOverGenVCent_MeanResPts_p[iter][ptIter][etaIter][qgIter][centIter] = new TH1F(Form("jtRecoOverGenVCent_%s_MeanResPts_%s_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), ptStr1.c_str(), etaStr1.c_str(), centStr.c_str()), Form(";(Reco. %s Jet p_{T})/(Gen. Jet p_{T});Events (%s)", jetAlgo[iter].c_str(), centStr.c_str()), 30, 0, 3);
	  }
	}
      }
    }

    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      std::string centStr = "PP";
      if(isPbPb) centStr = Form("cent%dto%d", centBins[centIter+1]/2, centBins[centIter]/2);

      
      bkgEstimate_p[iter][centIter] = new TH1F(Form("bkgEstimate_%s_%s_h", jetAlgo[iter].c_str(), centStr.c_str()), Form(";Bkg. Estimate;Events"), 100, 0, 50);

      areaEstimate_p[iter][centIter] = new TH1F(Form("areaEstimate_%s_%s_h", jetAlgo[iter].c_str(), centStr.c_str()), Form(";Area;Events"), 100, 0, 4*2*TMath::Pi());

      bkgAreaEstimate_p[iter][centIter] = new TH2F(Form("bkgAreaEstimate_%s_%s_h", jetAlgo[iter].c_str(), centStr.c_str()), Form(";Area;Bkg."), 20, 0, 4*2*TMath::Pi(), 20, 0, 50);

      for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	jtRecoGenDRVPt_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoGenDRVPt_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; #DeltaR_{Reco. %s, Gen.}", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins, 20, 0, 0.4);

	jtRecoVGen_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoVGen_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; Reco. %s Jet p_{T}", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins, nJtPtBins, jtPtBins);

	jtRecoOverGenVPt_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoOverGenVPt_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; (Reco. %s Jet p_{T})/(Gen. Jet p_{T})", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins, 30, 0, 3);

	jtRecoOverRawVPt_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoOverRawVPt_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Reco. Jet p_{T}; (Reco %s Jet p_{T})/(Raw Jet p_{T})", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins, 30, 0, 3);

	jtRawOverGenVPt_p[iter][centIter][qgIter] = new TH2F(Form("jtRawOverGenVPt_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; (Raw %s Jet p_{T})/(Gen Jet p_{T})", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins, 30, 0, 3);

	jtEffVPt_p[iter][centIter][qgIter] = new TH2F(Form("jtEffVPt_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};Eff. (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins, 2, -0.5, 1.5);

	jtFakeVPt_p[iter][centIter][qgIter] = new TH2F(Form("jtFakeVPt_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Reco. Jet p_{T};Fake (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins, 2, -0.5, 1.5);

	jtEffVPtEta_p[iter][centIter][qgIter] = new TH2F(Form("jtEffVPtEta_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet #eta;Gen. Jet p_{T} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins, nJtPtBins, jtPtBins);
	jtEffVPtEta_Denom_p[iter][centIter][qgIter] = new TH2F(Form("jtEffVPtEta_%s_Denom_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet #eta;Gen. Jet p_{T} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins, nJtPtBins, jtPtBins);
      

	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  jtRecoOverGenVPt_Mean_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVPt_%s_%sMean_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#mu_{Reco./Gen.} (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins);

	  jtRecoOverRawVPt_Mean_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverRawVPt_%s_%sMean_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Reco. Jet p_{T};#mu_{Reco./Raw} (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins);

	  jtRawOverGenVPt_Mean_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRawOverGenVPt_%s_%sMean_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#mu_{Raw/Gen.} (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins);
	}

	jtEffVPt_Mean_p[iter][centIter][qgIter] = new TH1F(Form("jtEffVPt_%s_Mean_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};Eff. (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins);
	
	jtFakeVPt_Mean_p[iter][centIter][qgIter] = new TH1F(Form("jtFakeVPt_%s_Mean_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Reco. Jet p_{T};Fake (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins);
	
	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVPt_%s_%sRes_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#sigma_{Reco./Gen.} (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins);
	  
	  jtRecoOverRawVPt_Res_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverRawVPt_%s_%sRes_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Reco. Jet p_{T};#sigma_{Reco./Raw} (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins);
	  
	  jtRawOverGenVPt_Res_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRawOverGenVPt_%s_%sRes_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#sigma_{Raw/Gen.} (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins);
	}
	
	jtRecoOverGenVEta_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoOverGenVEta_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet #eta; (Reco. %s Jet p_{T})/(Gen. Jet p_{T})", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins, 30, 0, 3);

	jtRecoOverRawVEta_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoOverRawVEta_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Reco. Jet #eta; (Reco %s Jet p_{T})/(Raw Jet p_{T})", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins, 30, 0, 3);

	jtRawOverGenVEta_p[iter][centIter][qgIter] = new TH2F(Form("jtRawOverGenVEta_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet #eta; (Raw %s Jet p_{T})/(Gen Jet p_{T})", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins, 30, 0, 3);

	jtEffVEta_p[iter][centIter][qgIter] = new TH2F(Form("jtEffVEta_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet #eta;Eff. (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins, 2, -0.5, 1.5);
	
	jtFakeVEta_p[iter][centIter][qgIter] = new TH2F(Form("jtFakeVEta_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Reco. Jet #eta;Fake (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins, 2, -0.5, 1.5);

	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  jtRecoOverGenVEta_Mean_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVEta_%s_%sMean_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet #eta;#mu_{Reco./Gen.} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins);

	  jtRecoOverRawVEta_Mean_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverRawVEta_%s_%sMean_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Reco. Jet #eta;#mu_{Reco./Raw} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins);
	  
	  jtRawOverGenVEta_Mean_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRawOverGenVEta_%s_%sMean_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet #eta;#mu_{Raw/Gen.} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins);
	}
	
	jtEffVEta_Mean_p[iter][centIter][qgIter] = new TH1F(Form("jtEffVEta_%s_Mean_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet #eta;Eff. (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins);
	
	jtFakeVEta_Mean_p[iter][centIter][qgIter] = new TH1F(Form("jtFakeVEta_%s_Mean_%s_%s_h", qg[qgIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Reco. Jet #eta;Fake (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins);
	
	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  jtRecoOverGenVEta_Res_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVEta_%s_%sRes_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet #eta;#sigma_{Reco./Gen.} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins);
	  
	  jtRecoOverRawVEta_Res_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverRawVEta_%s_%sRes_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Reco. Jet #eta;#sigma_{Reco./Raw} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins);
	  
	  jtRawOverGenVEta_Res_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRawOverGenVEta_%s_%sRes_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo[iter].c_str(), centStr.c_str()), Form(";Gen. Jet #eta;#sigma_{Raw/Gen.} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins);
	}
	

	for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	  Int_t ptLowInt = std::trunc(jtPtBins[jtIter]);
	  Int_t ptHiInt = std::trunc(jtPtBins[jtIter+1]);
	  
	  Int_t ptLowDec = std::trunc(jtPtBins[jtIter]*10 - ptLowInt*10);
	  Int_t ptHiDec = std::trunc(jtPtBins[jtIter+1]*10 - ptHiInt*10);
	  
	  jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoOverGenVPt_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo[iter].c_str(), centStr.c_str()), Form(";(Reco. %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo[iter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 30, 0, 3);
	  
	  jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoOverRawVPt_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo[iter].c_str(), centStr.c_str()), Form(";(Reco. %s Jet p_{T})/(Raw. Jet p_{T});Events (%d.%d<p_{T,Reco.}<%d.%d)", jetAlgo[iter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 30, 0, 3);
	  
	  jtRawOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRawOverGenVPt_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo[iter].c_str(), centStr.c_str()), Form(";(Raw %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo[iter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 30, 0, 3);
	  
	  jtEffVPt_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtEffVPt_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo[iter].c_str(), centStr.c_str()), Form(";(Eff. (%s));Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo[iter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 2, -0.5, 1.5);
	  
	  jtFakeVPt_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtFakeVPt_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo[iter].c_str(), centStr.c_str()), Form(";(Fake (%s));Events (%d.%d<p_{T,Reco.}<%d.%d)", jetAlgo[iter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 2, -0.5, 1.5);
	}
	
	
	for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	  Int_t etaLowInt = std::trunc(jtEtaBins[jtIter]);
	  Int_t etaHiInt = std::trunc(jtEtaBins[jtIter+1]);
	  
	  Int_t etaLowDec = std::trunc(jtEtaBins[jtIter]*10 - etaLowInt*10);
	  Int_t etaHiDec = std::trunc(jtEtaBins[jtIter+1]*10 - etaHiInt*10);
	  
	  jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoOverGenVEta_%s_MeanResPts_Eta%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo[iter].c_str(), centStr.c_str()), Form(";(Reco. %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<#eta_{Gen.}<%d.%d)", jetAlgo[iter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec), 30, 0, 3);
	  
	  jtRecoOverRawVEta_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoOverRawVEta_%s_MeanResPts_Eta%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo[iter].c_str(), centStr.c_str()), Form(";(Reco. %s Jet p_{T})/(Raw. Jet p_{T});Events (%d.%d<#eta_{Reco.}<%d.%d)", jetAlgo[iter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec), 30, 0, 3);
	  
	  jtRawOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRawOverGenVEta_%s_MeanResPts_Eta%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo[iter].c_str(), centStr.c_str()), Form(";(Raw %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<#eta_{Gen.}<%d.%d)", jetAlgo[iter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec), 30, 0, 3);
	  
	  jtEffVEta_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtEffVEta_%s_MeanResPts_Eta%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo[iter].c_str(), centStr.c_str()), Form(";(Eff. (%s));Events (%d.%d<#eta_{Gen.}<%d.%d)", jetAlgo[iter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec), 2, -0.5, 1.5);
	  
	  jtFakeVEta_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtFakeVEta_%s_MeanResPts_Eta%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo[iter].c_str(), centStr.c_str()), Form(";(Fake (%s));Events (%d.%d<#eta_{Reco.}<%d.%d)", jetAlgo[iter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec), 2, -0.5, 1.5);
	}
      }
    }
  }

  const Int_t numberOfFiles = (Int_t)listOfFiles_p->size();
  Int_t fileDiv = ((Int_t)(numberOfFiles/10));
  if(fileDiv < 1) fileDiv = 1;

  fileDiv = 1;

  //  std::cout << "Number of files: " << numberOfFiles << std::endl;

  for(Int_t fileIter = 0; fileIter < numberOfFiles; fileIter++){
    if(fileIter%fileDiv == 0) std::cout << "File # " << fileIter << "/" << numberOfFiles << std::endl;

    TFile* inFile_p = new TFile(listOfFiles_p->at(fileIter).c_str(), "READ");

    //    std::cout << listOfFiles_p->at(fileIter).c_str() << std::endl;

    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
    TTree* genTree_p = (TTree*)inFile_p->Get("HiGenParticleAna/hi");
    TTree* jetTree_p[nJetAlgo];

    for(Int_t iter = 0; iter < nJetAlgo; iter++){
      jetTree_p[iter] = (TTree*)inFile_p->Get(Form("%sJetAnalyzer/t", jetAlgo[iter].c_str()));
    }

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
      jetTree_p[iter]->SetBranchStatus("chargedN", 1);
      jetTree_p[iter]->SetBranchStatus("refpt", 1);
      jetTree_p[iter]->SetBranchStatus("refeta", 1);
      jetTree_p[iter]->SetBranchStatus("subid", 1);
      jetTree_p[iter]->SetBranchStatus("refparton_flavor", 1);
      jetTree_p[iter]->SetBranchStatus("ngen", 1);
      jetTree_p[iter]->SetBranchStatus("genpt", 1);
      jetTree_p[iter]->SetBranchStatus("geneta", 1);
      jetTree_p[iter]->SetBranchStatus("genphi", 1);
      jetTree_p[iter]->SetBranchStatus("genmatchindex", 1);
      jetTree_p[iter]->SetBranchStatus("gensubid", 1);
      
      jetTree_p[iter]->SetBranchAddress("nref", &nJt_[iter]);
      jetTree_p[iter]->SetBranchAddress("jtpt", jtPt_[iter]);
      jetTree_p[iter]->SetBranchAddress("rawpt", jtRawPt_[iter]);
      jetTree_p[iter]->SetBranchAddress("jteta", jtEta_[iter]);
      jetTree_p[iter]->SetBranchAddress("jtphi", jtPhi_[iter]);
      jetTree_p[iter]->SetBranchAddress("hcalSum", jtHCalSum_[iter]);
      jetTree_p[iter]->SetBranchAddress("ecalSum", jtECalSum_[iter]);
      jetTree_p[iter]->SetBranchAddress("neutralSum", jtNeutralSum_[iter]);
      jetTree_p[iter]->SetBranchAddress("chargedSum", jtChargedSum_[iter]);
      jetTree_p[iter]->SetBranchAddress("chargedN", jtChargedN_[iter]);
      jetTree_p[iter]->SetBranchAddress("refpt", refPt_[iter]);
      jetTree_p[iter]->SetBranchAddress("refeta", refEta_[iter]);
      jetTree_p[iter]->SetBranchAddress("subid", refSubID_[iter]);
      jetTree_p[iter]->SetBranchAddress("refparton_flavor", refPartFlav_[iter]);
      jetTree_p[iter]->SetBranchAddress("ngen", &nGenJt_[iter]);
      jetTree_p[iter]->SetBranchAddress("genpt", genJtPt_[iter]);
      jetTree_p[iter]->SetBranchAddress("geneta", genJtEta_[iter]);
      jetTree_p[iter]->SetBranchAddress("genphi", genJtPhi_[iter]);
      jetTree_p[iter]->SetBranchAddress("genmatchindex", genJtMatchIndex_[iter]);
      jetTree_p[iter]->SetBranchAddress("gensubid", genSubId_[iter]);
    }


    //    std::cout << "Gets Here B" << std::endl;
    const Int_t nEntries = jetTree_p[0]->GetEntries();
    Int_t entryDiv = ((Int_t)(nEntries/10));

    //    std::cout << "Gets Here C" << std::endl;

    for(Int_t entry = 0; entry < nEntries; entry++){
      //      std::cout << "Gets here d" << std::endl;

      if(entry%entryDiv == 0 && nEntries >= 10000) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
      //      std::cout << "Gets here e" << std::endl;
      hiTree_p->GetEntry(entry);
      genTree_p->GetEntry(entry);

      Int_t centPos = -1;

      if(TMath::Abs(vz_) > 15) continue;

      if(isPbPb)
	for(Int_t iter = 0; iter < nCentBins2; iter++){
	  if(hiBin_ < centBins[iter] && hiBin_ >= centBins[iter+1]){
	    centPos = iter;
	    break;
	  }
	}
      else centPos = 0;

      Float_t maxElePt = -1;
      Float_t maxEleEta = -100;
      Float_t maxElePhi = -100;
      
      Float_t twoElePt = -1;
      Float_t twoEleEta = -100;
      Float_t twoElePhi = -100;
      /*
	const Int_t nMult_ = genPt_p->size();
	
	for(Int_t iter = 0; iter < nMult_; iter++){
	if(TMath::Abs(genPDG_p->at(iter)) != 11) continue;
	
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
      */
      if(maxElePt < 10) nNoLeadEle++;
      if(twoElePt < 5) nNoSubleadEle++;
      
      for(Int_t iter = 0; iter < nJetAlgo; iter++){
	jetTree_p[iter]->GetEntry(entry);
      }

      //      std::cout << "ALPHA" << std::endl;
      
      Bool_t skipEvent = false;
      
      for(Int_t iter = 0; iter < nGenJt_[0]; iter++){
	if(genJtPt_[0][iter] < 15) continue;
	if(TMath::Abs(genJtEta_[0][iter]) > 3.0){
	  skipEvent = true;
	  break;
	}
      }
      
      if(skipEvent) continue;
      
      for(Int_t algoIter = 0; algoIter < nJetAlgo; algoIter++){
	for(Int_t jtIter = 0; jtIter < nJt_[algoIter]; jtIter++){
	  if(TMath::Abs(jtEta_[algoIter][jtIter]) > 2.0) continue;
	  if(jtPt_[algoIter][jtIter] < 5.0) continue;
	  
	  if(isPbPb && !isGoodCaloJet(jetAlgoTemp[algoIter].c_str(), jtHCalSum_[algoIter][jtIter], jtECalSum_[algoIter][jtIter]) && !isGoodPFJet(jetAlgoTemp[algoIter].c_str(), jtChargedN_[algoIter][jtIter], jtChargedSum_[algoIter][jtIter], jtNeutralSum_[algoIter][jtIter], jtRawPt_[algoIter][jtIter])) continue;
	  
	  
	
	  //	if(refSubID_[algoIter][jtIter] != 0 && refSubID_[algoIter][jtIter] != -1) continue;

	  if(maxElePt > 10)
	    if(getDR(jtEta_[algoIter][jtIter], jtPhi_[algoIter][jtIter], maxEleEta, maxElePhi) < 0.4) continue;
	  
	  if(twoElePt > 5)
	    if(getDR(jtEta_[algoIter][jtIter], jtPhi_[algoIter][jtIter], twoEleEta, twoElePhi) < 0.4) continue;
	  
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
	
	if(doGetBkg) getBkgEstimate(jetAlgoR[algoIter], genPt_p, genPhi_p, genEta_p, nGenJt_[algoIter], genJtPt_[algoIter], genJtPhi_[algoIter], genJtEta_[algoIter], evtBkg, evtArea, isBkg_p);

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
	      if(isBkg_p->at(bkgIter)) std::cout << "  particle pt, phi, eta: " << genPt_p->at(bkgIter) << ", " << genPhi_p->at(bkgIter) << ", " << genEta_p->at(bkgIter) << std::endl;
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
	for(Int_t jtIter = 0; jtIter < nGenJt_[algoIter]; jtIter++){
	  if(TMath::Abs(genJtEta_[algoIter][jtIter]) > 2.0) continue;
	  
	  if(genJtPt_[algoIter][jtIter] < 5.0) break;
	  
	  if(genSubId_[algoIter][jtIter] != 0) continue;

	  if(maxElePt > 10)
	    if(getDR(genJtEta_[algoIter][jtIter], genJtPhi_[algoIter][jtIter], maxEleEta, maxElePhi) < 0.4) continue;
	  
	  if(twoElePt > 5)
	    if(getDR(genJtEta_[algoIter][jtIter], genJtPhi_[algoIter][jtIter], twoEleEta, twoElePhi) < 0.4) continue;

	  Int_t fillVal = 0;
	  Float_t tempDR = 999;
	  Int_t tempPos = -1;

	  for(Int_t jtIter2 = 0; jtIter2 < nJt_[algoIter]; jtIter2++){
	    if(isUsedRecoJet[jtIter2]) continue;
	    
	    if(isPbPb && !isGoodCaloJet(jetAlgoTemp[algoIter].c_str(), jtHCalSum_[algoIter][jtIter2], jtECalSum_[algoIter][jtIter2]) && !isGoodPFJet(jetAlgoTemp[algoIter].c_str(), jtChargedN_[algoIter][jtIter2], jtChargedSum_[algoIter][jtIter2], jtNeutralSum_[algoIter][jtIter2], jtRawPt_[algoIter][jtIter2])) continue;
	    
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
	  
	  Int_t qgPos[2] = {0, -1};
	  if(TMath::Abs(refPartFlav_[algoIter][tempPos]) < 9) qgPos[1] = 1;
	  else if(TMath::Abs(refPartFlav_[algoIter][tempPos]) == 21) qgPos[1] = 2;
	  
	  //	if(genJtMatchIndex_[algoIter][jtIter] >= 0) fillVal = 1;
	  
	  for(Int_t qgIter = 0; qgIter < 2; qgIter++){
	    if(qgPos[qgIter] == -1) continue;
	    
	    if(recoPosJet[jtIter] != -1){
	      jtRecoVGen_p[algoIter][centPos][qgPos[qgIter]]->Fill(genJtPt_[algoIter][jtIter], jtPt_[algoIter][recoPosJet[jtIter]]);

	      jtRecoGenDRVPt_p[algoIter][centPos][qgPos[qgIter]]->Fill(genJtPt_[algoIter][jtIter], getDR(jtEta_[algoIter][recoPosJet[jtIter]], jtPhi_[algoIter][recoPosJet[jtIter]], genJtEta_[algoIter][jtIter], genJtPhi_[algoIter][jtIter]));

	      jtRecoOverGenVPt_p[algoIter][centPos][qgPos[qgIter]]->Fill(genJtPt_[algoIter][jtIter], jtPt_[algoIter][recoPosJet[jtIter]]/genJtPt_[algoIter][jtIter]);
	      jtRawOverGenVPt_p[algoIter][centPos][qgPos[qgIter]]->Fill(genJtPt_[algoIter][jtIter], jtRawPt_[algoIter][recoPosJet[jtIter]]/genJtPt_[algoIter][jtIter]);
	    }
	    
	    jtEffVPt_p[algoIter][centPos][qgPos[qgIter]]->Fill(genJtPt_[algoIter][jtIter], fillVal);
	    
	    jtEffVPtEta_Denom_p[algoIter][centPos][qgPos[qgIter]]->Fill(genJtEta_[algoIter][jtIter], genJtPt_[algoIter][jtIter]);
	    if(fillVal) jtEffVPtEta_p[algoIter][centPos][qgPos[qgIter]]->Fill(genJtEta_[algoIter][jtIter], genJtPt_[algoIter][jtIter]);
	    
	    for(Int_t jtIter2 = 0; jtIter2 < nJtPtBins2; jtIter2++){
	      if(genJtPt_[algoIter][jtIter] > jtPtBins2[jtIter2] && genJtPt_[algoIter][jtIter] < jtPtBins2[jtIter2+1]){
		for(Int_t etaIter2 = 0; etaIter2 < nJtEtaBins2; etaIter2++){
		  if(genJtEta_[algoIter][jtIter] > jtEtaBins2[etaIter2] && genJtEta_[algoIter][jtIter] < jtEtaBins2[etaIter2+1]){
		    jtRecoOverGenVCent_MeanResPts_p[algoIter][jtIter2][etaIter2][qgPos[qgIter]][centPos]->Fill(jtPt_[algoIter][recoPosJet[jtIter]]/genJtPt_[algoIter][jtIter]);

		    break;
		  }
		}
		break;
	      }
	    }

	    for(Int_t jtIter2 = 0; jtIter2 < nJtPtBins; jtIter2++){
	      if(genJtPt_[algoIter][jtIter] > jtPtBins[jtIter2] && genJtPt_[algoIter][jtIter] < jtPtBins[jtIter2+1]){

		if(recoPosJet[jtIter] != -1){
		  jtRecoOverGenVPt_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(jtPt_[algoIter][recoPosJet[jtIter]]/genJtPt_[algoIter][jtIter]);
		  jtRawOverGenVPt_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(jtRawPt_[algoIter][recoPosJet[jtIter]]/genJtPt_[algoIter][jtIter]);	
		}
		
		jtEffVPt_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(fillVal);
		break;
	      }
	    }

	    if(genJtPt_[algoIter][jtIter] > 30){
	      if(recoPosJet[jtIter] != -1){
		jtRecoOverGenVEta_p[algoIter][centPos][qgPos[qgIter]]->Fill(genJtEta_[algoIter][jtIter], jtPt_[algoIter][recoPosJet[jtIter]]/genJtPt_[algoIter][jtIter]);
		jtRawOverGenVEta_p[algoIter][centPos][qgPos[qgIter]]->Fill(genJtEta_[algoIter][jtIter], jtRawPt_[algoIter][recoPosJet[jtIter]]/genJtPt_[algoIter][jtIter]);
	      }	    
	      
	      jtEffVEta_p[algoIter][centPos][qgPos[qgIter]]->Fill(genJtEta_[algoIter][jtIter], fillVal);
	      
	      for(Int_t jtIter2 = 0; jtIter2 < nJtEtaBins; jtIter2++){
		if(genJtEta_[algoIter][jtIter] > jtEtaBins[jtIter2] && genJtEta_[algoIter][jtIter] < jtEtaBins[jtIter2+1]){
		  if(recoPosJet[jtIter] != -1){
		    jtRecoOverGenVEta_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(jtPt_[algoIter][recoPosJet[jtIter]]/genJtPt_[algoIter][jtIter]);
		    jtRawOverGenVEta_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(jtRawPt_[algoIter][recoPosJet[jtIter]]/genJtPt_[algoIter][jtIter]);	
		  }
		  
		  jtEffVEta_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(fillVal);
		  break;
		}
	      }
	    }
	  }
	}
	
	for(Int_t boolIter = 0; boolIter < boolSize; boolIter++){
	  if(TMath::Abs(jtEta_[algoIter][boolIter]) > 2.0) continue;

	  if(jtPt_[algoIter][boolIter] < 5.0) continue;
	  
	  if(isPbPb && !isGoodCaloJet(jetAlgoTemp[algoIter].c_str(), jtHCalSum_[algoIter][boolIter], jtECalSum_[algoIter][boolIter]) && !isGoodPFJet(jetAlgoTemp[algoIter].c_str(), jtChargedN_[algoIter][boolIter], jtChargedSum_[algoIter][boolIter], jtNeutralSum_[algoIter][boolIter], jtRawPt_[algoIter][boolIter])) continue;
	  
	  Int_t fillVal = 1;
	  if(isUsedRecoJet[boolIter]) fillVal = 0;
	  
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
    delete inFile_p;
  }


  std::cout << "#Events with no lead electron: " << nNoLeadEle << std::endl;
  std::cout << "#Events with no sublead electron: " << nNoSubleadEle << std::endl;

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

	    FitGauss(jtRecoOverGenVCent_MeanResPts_p[iter][ptIter][etaIter][qgIter][centIter], tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);

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

    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	  Float_t tempMean[nMeanFit];
	  Float_t tempMeanErr[nMeanFit];
	  Float_t tempRes[nMeanFit];
	  Float_t tempResErr[nMeanFit];
	  
	  tempMean[0] = jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();	
	  tempRes[0] = jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();
	  
	  FitGauss(jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter], tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);
	
	  for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	    jtRecoOverGenVPt_Mean_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempMean[mIter]);
	    jtRecoOverGenVPt_Mean_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempMeanErr[mIter]);
	    jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempRes[mIter]);
	    jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempResErr[mIter]);
	  }
	  
	  tempMean[0] = jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();
	  tempRes[0] = jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();

	  FitGauss(jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter], tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);
	  
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
	  
	  FitGauss(jtRawOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter], tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);

	  for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	    jtRawOverGenVPt_Mean_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempMean[mIter]);
	    jtRawOverGenVPt_Mean_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempMeanErr[mIter]);
	    jtRawOverGenVPt_Res_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempRes[mIter]);
	    jtRawOverGenVPt_Res_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempResErr[mIter]);
	  }
	  
	  jtEffVPt_Mean_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, jtEffVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean());
	  jtEffVPt_Mean_p[iter][centIter][qgIter]->SetBinError(jtIter+1, jtEffVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError());
	  
	  jtFakeVPt_Mean_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, jtFakeVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean());
	  jtFakeVPt_Mean_p[iter][centIter][qgIter]->SetBinError(jtIter+1, jtFakeVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError());
	}

	for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	  Float_t tempMean[nMeanFit];
	  Float_t tempMeanErr[nMeanFit];
	  Float_t tempRes[nMeanFit];
	  Float_t tempResErr[nMeanFit];

	  tempMean[0] = jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();
	  tempRes[0] = jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();
	  
	  FitGauss(jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter], tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);
	  
	  
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
	  
	  FitGauss(jtRecoOverRawVEta_MeanResPts_p[iter][centIter][qgIter][jtIter], tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);
	  
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
	  
	  FitGauss(jtRawOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter], tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);

	  for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	    jtRawOverGenVEta_Mean_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempMean[mIter]);
	    jtRawOverGenVEta_Mean_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempMeanErr[mIter]);
	    jtRawOverGenVEta_Res_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempRes[mIter]);
	    jtRawOverGenVEta_Res_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempResErr[mIter]);
	  }
	  
	  jtEffVEta_Mean_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, jtEffVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean());
	  jtEffVEta_Mean_p[iter][centIter][qgIter]->SetBinError(jtIter+1, jtEffVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError());
	  
	  jtFakeVEta_Mean_p[iter][centIter][qgIter]->SetBinContent(jtIter+1, jtFakeVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean());
	  jtFakeVEta_Mean_p[iter][centIter][qgIter]->SetBinError(jtIter+1, jtFakeVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError());
	}
      }
    }
  }

  outFile_p->cd();

  TNamed pathStr("pathStr", inFileName.c_str());
  pathStr.Write("", TObject::kOverwrite);

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    TDirectory* dir_p = outFile_p->GetDirectory(Form("%s", jetAlgo[iter].c_str()));
    if(dir_p){
      dir_p->cd();
    }
    else{
      dir_p = outFile_p->mkdir(Form("%s", jetAlgo[iter].c_str()));
      dir_p->cd();
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
	jtRecoVGen_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoGenDRVPt_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoOverGenVPt_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoOverRawVPt_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRawOverGenVPt_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	
	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  jtRecoOverGenVPt_Mean_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRecoOverRawVPt_Mean_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRecoOverRawVPt_Res_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRawOverGenVPt_Mean_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRawOverGenVPt_Res_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	}
	
	jtEffVPt_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtEffVPt_Mean_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtFakeVPt_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtFakeVPt_Mean_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	
      //      jtEffVPtEta_p[iter][centIter][qgIter]->Write(Form("%s_NUM",jtEffVPtEta_p[iter][centIter][qgIter]->GetName()), TObject::kOverwrite);
      //      jtEffVPtEta_Denom_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtEffVPtEta_p[iter][centIter][qgIter]->Divide(jtEffVPtEta_Denom_p[iter][centIter][qgIter]);
	jtEffVPtEta_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);

	for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	  jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  jtRawOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  jtEffVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
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
	
	jtEffVEta_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtEffVEta_Mean_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtFakeVEta_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtFakeVEta_Mean_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	
	for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	  jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  jtRecoOverRawVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  jtRawOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  jtEffVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  jtFakeVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	}
      }
    }
  }

  for(Int_t iter = 0; iter < nJetAlgo; iter++){


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
	delete jtRecoVGen_p[iter][centIter][qgIter];
	delete jtRecoGenDRVPt_p[iter][centIter][qgIter];
	delete jtRecoOverGenVPt_p[iter][centIter][qgIter];
	delete jtRecoOverRawVPt_p[iter][centIter][qgIter];
	delete jtRawOverGenVPt_p[iter][centIter][qgIter];
	
	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  delete jtRecoOverGenVPt_Mean_p[iter][centIter][qgIter][mIter];
	  delete jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter];
	  delete jtRecoOverRawVPt_Mean_p[iter][centIter][qgIter][mIter];
	  delete jtRecoOverRawVPt_Res_p[iter][centIter][qgIter][mIter];
	  delete jtRawOverGenVPt_Mean_p[iter][centIter][qgIter][mIter];
	  delete jtRawOverGenVPt_Res_p[iter][centIter][qgIter][mIter];
	}
	
	delete jtEffVPt_p[iter][centIter][qgIter];
	delete jtEffVPt_Mean_p[iter][centIter][qgIter];
	delete jtFakeVPt_p[iter][centIter][qgIter];
	delete jtFakeVPt_Mean_p[iter][centIter][qgIter];
	
	delete jtEffVPtEta_p[iter][centIter][qgIter];
	delete jtEffVPtEta_Denom_p[iter][centIter][qgIter];
	
	for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	  delete jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtRecoOverRawVPt_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtRawOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtEffVPt_MeanResPts_p[iter][centIter][qgIter][jtIter];
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
	
	delete jtEffVEta_p[iter][centIter][qgIter];
	delete jtEffVEta_Mean_p[iter][centIter][qgIter];
	delete jtFakeVEta_p[iter][centIter][qgIter];
	delete jtFakeVEta_Mean_p[iter][centIter][qgIter];
	
	for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	  delete jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtRecoOverRawVEta_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtRawOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtEffVEta_MeanResPts_p[iter][centIter][qgIter][jtIter];
	  delete jtFakeVEta_MeanResPts_p[iter][centIter][qgIter][jtIter];
	}
      }
    }
  }
    
  outFile_p->Close();
  delete outFile_p;

  delete date;

  listOfFiles_p->clear();
  delete listOfFiles_p;

  return;
}
