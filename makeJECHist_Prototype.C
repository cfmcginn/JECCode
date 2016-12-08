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

#include "include/jecConfigParser.h"

#include "getPtEtaJetResidualCorr.h"

const std::string storeStr = "/store";
const std::string xrootdStr = "root://xrootd.unl.edu/";


const Bool_t debugMode = false;
const Bool_t doGetBkg = true;

const Int_t nJtCat = 3;
const std::string jtCat[nJtCat] = {"Eta2", "Eta1", "Dijet"};

const Int_t nMeanFit = 2;
const std::string meanFit[nMeanFit] = {"", "Fit"};

const Int_t nQG = 3;
const std::string qg[nQG] = {"Inc", "Q", "G"};

const Int_t nMaxJets = 500;

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


float findNormNcoll(int hiBin, int centPos, std::vector<unsigned int> centBins)
{
  float nCollWeight = findNcoll(hiBin);
  
  float normNColl = findNcoll(centBins.at(centPos+1)*2-1);

  return nCollWeight/normNColl;
}

void genSort(Int_t nGenJt, Float_t genJtPt[], Float_t genJtPhi[], Float_t genJtEta[], Int_t genJtMatchIndex[], Int_t genJtSubID[])
{
  Int_t genPos = 0;

  while(genPos < nGenJt){
    Bool_t doSwap = false;

    for(Int_t iter = genPos+1; iter < nGenJt; iter++){
      if(genJtPt[genPos] < genJtPt[iter]){
	Float_t tempPt = genJtPt[genPos];
	Float_t tempPhi = genJtPhi[genPos];
	Float_t tempEta = genJtEta[genPos];
	Int_t tempMatchIndex = genJtMatchIndex[genPos];
	Int_t tempSubID = genJtSubID[genPos];

	genJtPt[genPos] = genJtPt[iter];
	genJtPhi[genPos] = genJtPhi[iter];
	genJtEta[genPos] = genJtEta[iter];
	genJtMatchIndex[genPos] = genJtMatchIndex[iter];
	genJtSubID[genPos] = genJtSubID[iter];

	genJtPt[iter] = tempPt;
	genJtPhi[iter] = tempPhi;
	genJtEta[iter] = tempEta;
	genJtMatchIndex[iter] = tempMatchIndex;
	genJtSubID[iter] = tempSubID;

	doSwap = true;
      }
    }

    if(!doSwap) genPos++;
  }

  return;
}


void FitGauss(TH1F* hist_p, Float_t& mean, Float_t& meanErr, Float_t& res, Float_t& resErr, Bool_t isWeighted = true)
{
  if(hist_p->Integral() == 0) return;
  if(hist_p->GetEntries() == 0) return;

  std::string fitOpt = "Q M E";
  if(!isWeighted) fitOpt = "Q M E L";

  TF1* f1_p = new TF1("f1_p", "gaus", hist_p->GetXaxis()->GetXmin(), hist_p->GetXaxis()->GetXmax());

  hist_p->Fit("f1_p", fitOpt.c_str());

  mean = f1_p->GetParameter(1);
  res = f1_p->GetParameter(2);

  meanErr = f1_p->GetParError(1);
  resErr = f1_p->GetParError(2);

  if(f1_p->GetProb() > .01) return;

  hist_p->Rebin(2);

  hist_p->Fit("f1_p", fitOpt.c_str());

  mean = f1_p->GetParameter(1);
  res = f1_p->GetParameter(2);

  meanErr = f1_p->GetParError(1);
  resErr = f1_p->GetParError(2);

  if(f1_p->GetProb() > .01) return;

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
      hist_p->Fit("f1_p", fitOpt.c_str(), "", fitLow, fitHi);
      
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
      
      if(hist_p->Integral(tempBin, meanBin+iter) > (.95-.05*fitIter)*hist_p->Integral() || tempBin == 1){
	fitLow = hist_p->GetBinCenter(tempBin);
	fitHi = hist_p->GetBinCenter(meanBin+iter);
	break;
      }
    }
    
    
    if(nBins >= 7){
      hist_p->Fit("f1_p", fitOpt.c_str(), "", fitLow, fitHi);
      
      mean = f1_p->GetParameter(1);
      res = f1_p->GetParameter(2);
      
      meanErr = f1_p->GetParError(1);
      resErr = f1_p->GetParError(2);
    }
  }

  delete f1_p;

  return;
}


int makeJECHist_Prototype(const std::string inConfigFileName)
{
  initPtEtaJetResidualCorr();

  TFile* ratioFile_p = new TFile("outputDir/merged_dgulhan-Pythia8_Dijet30_pp_TuneCUETP8M1_Hydjet_MinBia\
s_5020GeV_RECODEBUG_758_PrivMC_forest_v28_0_20160512_QGFRACTIONHIST.root", "READ");
  TH1F* qRatio_p = (TH1F*)ratioFile_p->Get("qRatioZOverDijet_h");
  TH1F* gRatio_p = (TH1F*)ratioFile_p->Get("gRatioZOverDijet_h");


  jecConfigParser config;
  if(!config.SetConfigParser(inConfigFileName)) return 1;

  const std::string outName = config.GetOutName();
  TFile* outFile_p = new TFile(outName.c_str(), "RECREATE");
  
  TFile* tempJetInFile_p = TFile::Open(config.GetInput(0).c_str(), "READ");
  std::vector<std::string> jetAlgoInFile = config.GetJetTypesFinal();
  tempJetInFile_p->Close();

  Int_t jetAlgoCheck = 0;
  while(jetAlgoCheck < (Int_t)jetAlgoInFile.size()){
    Int_t duplicatePos = -1;
    for(Int_t iter = jetAlgoCheck+1; iter < (Int_t)jetAlgoInFile.size(); iter++){
      if(jetAlgoInFile.at(jetAlgoCheck).size() == jetAlgoInFile.at(iter).size() && jetAlgoInFile.at(jetAlgoCheck).find(jetAlgoInFile.at(iter)) != std::string::npos){
	duplicatePos = iter;
	break;
      }
    }
    
    if(duplicatePos >= 0){
      jetAlgoInFile.erase(jetAlgoInFile.begin()+jetAlgoCheck);
    }
    else jetAlgoCheck++;
  }

  const Int_t nJetAlgo = (Int_t)jetAlgoInFile.size();
  std::cout << nJetAlgo << std::endl;
  std::vector<std::string> jetAlgo;

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    jetAlgo.push_back(jetAlgoInFile.at(iter).substr(0, jetAlgoInFile.at(iter).find("JetAnalyzer")));

    if(debugMode) std::cout << __LINE__ << std::endl;
  }
  
  if(debugMode) std::cout << __LINE__ << std::endl;
  
  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    std::cout << "Algos: " << jetAlgo.at(iter) << std::endl;
  }

  Float_t vz_;
  Int_t hiBin_ = 0;

  UInt_t run_;
  ULong64_t evt_;

  std::vector<float>* genPt_p = 0;
  std::vector<float>* genEta_p = 0;
  std::vector<float>* genPhi_p = 0;
  std::vector<int>* genPDG_p = 0;

  Int_t nJt_[nJetAlgo]; 
  Float_t jtPt_[nJetAlgo][nMaxJets];
  Float_t jtEta_[nJetAlgo][nMaxJets];
  Float_t jtPhi_[nJetAlgo][nMaxJets];
  Float_t refPt_[nJetAlgo][nMaxJets];
  Float_t refEta_[nJetAlgo][nMaxJets];
  Float_t refPhi_[nJetAlgo][nMaxJets];
  Int_t refSubID_[nJetAlgo][nMaxJets];
  Int_t refPartFlav_[nJetAlgo][nMaxJets];
  Float_t ptHat_[nJetAlgo];

  Int_t nGenJt_[nJetAlgo]; 
  Float_t genJtPt_[nJetAlgo][nMaxJets];
  Float_t genJtEta_[nJetAlgo][nMaxJets];
  Float_t genJtPhi_[nJetAlgo][nMaxJets];
  Int_t genJtMatchIndex_[nJetAlgo][nMaxJets];
  Int_t genJtSubID_[nJetAlgo][nMaxJets];

  std::vector<float>* mcPt_p=0;
  std::vector<float>* mcPhi_p=0;
  std::vector<float>* mcEta_p=0;
  std::vector<int>* mcPID_p=0;
  std::vector<int>* mcMomPID_p=0;

  std::vector<float>* phoPt_p=0;
  std::vector<float>* phoPhi_p=0;
  std::vector<float>* phoEta_p=0;
  std::vector<float>* phoSeedTime_p = 0;
  std::vector<float>* phoSwissCrx_p = 0;
  std::vector<float>* phoHoverE_p = 0;
  std::vector<float>* phoSigmaIEtaIEta_2012_p=0;
  std::vector<float>* pho_ecalClusterIsoR4_p=0;
  std::vector<float>* pho_hcalRechitIsoR4_p=0;
  std::vector<float>* pho_trackIsoR4PtCut20_p=0;
  std::vector<float>* phoR9_p=0;


  //pthat30  
  const Int_t nJtPtBins = config.GetNJtPtBins();
  const Float_t jtPtLow = config.GetJtPtLow();
  const Float_t jtPtHi = config.GetJtPtHi();
  Double_t jtPtBins[nJtPtBins+1];
  if(config.GetDoJtPtLogBins()) getLogBins(jtPtLow, jtPtHi, nJtPtBins, jtPtBins);
  else getLinBins(jtPtLow, jtPtHi, nJtPtBins, jtPtBins);

  const Int_t nJtPtEtaBins = config.GetJtPtEtaBins();
  Double_t jtPtEtaBins[nJtPtEtaBins+1];
  getLinBins(0, config.GetJtEtaMax(), nJtPtEtaBins, jtPtEtaBins);

  const Int_t nJtEtaBins = config.GetJtEtaBins();
  const Float_t jtEtaLow = -(config.GetJtEtaMax());
  const Float_t jtEtaHi = config.GetJtEtaMax();
  Double_t jtEtaBins[nJtEtaBins+1];
  getLinBins(jtEtaLow, jtEtaHi, nJtEtaBins, jtEtaBins);

  const Int_t nCentBins = config.GetNCentBins();

  const Int_t nJtDPhiBins = 20;
  const Float_t jtDPhiLow = 0.;
  const Float_t jtDPhiHi = 0.4;
  Double_t jtDPhiBins[nJtDPhiBins+1];
  getLinBins(jtDPhiLow, jtDPhiHi, nJtDPhiBins, jtDPhiBins);

  const Int_t nJtDEtaBins = 20;
  const Float_t jtDEtaLow = 0.;
  const Float_t jtDEtaHi = 0.4;
  Double_t jtDEtaBins[nJtDEtaBins+1];
  getLinBins(jtDEtaLow, jtDEtaHi, nJtDEtaBins, jtDEtaBins);

  const Int_t nJtDRBins = 20;
  const Float_t jtDRLow = 0.;
  const Float_t jtDRHi = 0.4;
  Double_t jtDRBins[nJtDRBins+1];
  getLinBins(jtDRLow, jtDRHi, nJtDRBins, jtDRBins);

  config.PrintCentBins();

  TH1F* ptHat_Unweighted_h = new TH1F("ptHat_Unweighted_h", ";pthat;Events", 37, 15, 200);
  TH1F* ptHat_Weighted_h = new TH1F("ptHat_Weighted_h", ";pthat;Events", 37, 15, 200);
  ptHat_Unweighted_h->Sumw2();
  ptHat_Weighted_h->Sumw2();

  const unsigned int nPtHats = config.GetNPthats();
  TH1F* ptHat_Unweighted_PerPthat_h[nPtHats];
  TH1F* ptHat_Weighted_PerPthat_h[nPtHats];

  for(unsigned int iter = 0; iter < nPtHats; iter++){
    std::string nameUnweighted = "ptHat_Unweighted_PerPthat_" + std::to_string(config.GetPthat(iter)) + "_h";
    std::string nameWeighted = "ptHat_Weighted_PerPthat_" + std::to_string(config.GetPthat(iter)) + "_h";

    ptHat_Unweighted_PerPthat_h[iter] = new TH1F(nameUnweighted.c_str(), ";pthat;Events", 37, 15, 200);
    ptHat_Weighted_PerPthat_h[iter] = new TH1F(nameWeighted.c_str(), ";pthat;Events", 37, 15, 200);

    ptHat_Unweighted_PerPthat_h[iter]->Sumw2();
    ptHat_Weighted_PerPthat_h[iter]->Sumw2();
  }


  TH1F* genPhoPt_Unweighted_h[nCentBins];
  TH1F* genPhoPt_Weighted_h[nCentBins];
  TH1F* recoPhoPt_Unweighted_h[nCentBins];
  TH1F* recoPhoPt_Weighted_h[nCentBins];

  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    std::string centStr = "PP";
    if(config.GetIsPbPb()) centStr = "Cent" + std::to_string(config.GetCentBinFromPos(centIter)) + "to" + std::to_string(config.GetCentBinFromPos(centIter+1));

    std::string genNameUnweighted = "genPhoPt_Unweighted_" + centStr + "_h";
    std::string genNameWeighted = "genPhoPt_Weighted_" + centStr + "_h";

    genPhoPt_Unweighted_h[centIter] = new TH1F(genNameUnweighted.c_str(), ";Gen. Photon p_{T};Events", 25, config.GetMinGammaPt(), config.GetMinGammaPt()+100);
    genPhoPt_Weighted_h[centIter] = new TH1F(genNameWeighted.c_str(), ";Gen. Photon p_{T};Events", 25, config.GetMinGammaPt(), config.GetMinGammaPt()+100);

    genPhoPt_Unweighted_h[centIter]->Sumw2();
    genPhoPt_Weighted_h[centIter]->Sumw2();

    std::string recoNameUnweighted = "recoPhoPt_Unweighted_" + centStr + "_h";
    std::string recoNameWeighted = "recoPhoPt_Weighted_" + centStr + "_h";

    recoPhoPt_Unweighted_h[centIter] = new TH1F(recoNameUnweighted.c_str(), ";Reco. Photon p_{T};Events", 50, config.GetMinGammaPt(), config.GetMinGammaPt()+100);
    recoPhoPt_Weighted_h[centIter] = new TH1F(recoNameWeighted.c_str(), ";Reco. Photon p_{T};Events", 50, config.GetMinGammaPt(), config.GetMinGammaPt()+100);

    recoPhoPt_Unweighted_h[centIter]->Sumw2();
    recoPhoPt_Weighted_h[centIter]->Sumw2();
  }

  TH1F* genJtPtPerPthat_p[nJetAlgo][nCentBins][nQG][nPtHats];
  TH2F* jtRecoVGen_p[nJetAlgo][nCentBins][nQG];
  TH2F* jtRecoOverGenVPt_p[nJetAlgo][nCentBins][nQG];

  TH1F* jtRecoOverGenVPt_Mean_p[nJetAlgo][nCentBins][nJtPtEtaBins+1][nQG][nMeanFit];
  TH1F* jtRecoOverGenVPt_Res_p[nJetAlgo][nCentBins][nJtPtEtaBins+1][nQG][nMeanFit];
  TH1F* jtRecoOverGenVPt_MeanResPts_p[nJetAlgo][nCentBins][nJtPtEtaBins+1][nQG][nJtPtBins];
  Int_t jtRecoOverGenVPt_MeanResPts_COUNTS[nJetAlgo][nCentBins][nJtPtEtaBins+1][nQG][nJtPtBins];

  TH2F* jtDPhiVPt_p[nJetAlgo][nCentBins][nJtPtEtaBins+1][nQG];
  TH2F* jtDEtaVPt_p[nJetAlgo][nCentBins][nJtPtEtaBins+1][nQG];
  TH2F* jtDRVPt_p[nJetAlgo][nCentBins][nJtPtEtaBins+1][nQG];

  TH1F* jtDPhiVPt_Mean_p[nJetAlgo][nCentBins][nJtPtEtaBins+1][nQG][nMeanFit];
  TH1F* jtDPhiVPt_Res_p[nJetAlgo][nCentBins][nJtPtEtaBins+1][nQG][nMeanFit];
  TH1F* jtDPhiVPt_MeanResPts_p[nJetAlgo][nCentBins][nJtPtEtaBins+1][nQG][nJtPtBins];


  TH1F* jtEffVPt_Mean_p[nJetAlgo][nCentBins][nJtPtEtaBins+1][nJetRecoCuts];
  Double_t jtEffVPt_Num[nJetAlgo][nCentBins][nJtPtEtaBins+1][nJetRecoCuts][nJtPtBins];
  Double_t jtEffVPt_Denom[nJetAlgo][nCentBins][nJtPtEtaBins+1][nJtPtBins];

  TH1F* jtRecoOverGenVEta_Mean_p[nJetAlgo][nCentBins][nQG][nMeanFit];
  TH1F* jtRecoOverGenVEta_Res_p[nJetAlgo][nCentBins][nQG][nMeanFit];
  TH1F* jtRecoOverGenVEta_MeanResPts_p[nJetAlgo][nCentBins][nQG][nJtEtaBins];
  Int_t jtRecoOverGenVEta_MeanResPts_COUNTS[nJetAlgo][nCentBins][nQG][nJtEtaBins];



  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      std::string centStr = "PP";
      if(config.GetIsPbPb()) centStr = "Cent" + std::to_string(config.GetCentBinFromPos(centIter)) + "to" + std::to_string(config.GetCentBinFromPos(centIter+1));


      for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	jtRecoVGen_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoVGen_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; Reco. %s Jet p_{T}", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins, nJtPtBins, jtPtBins);
	
	jtRecoOverGenVPt_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoOverGenVPt_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; (Reco. %s Jet p_{T})/(Gen. Jet p_{T})", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins, 30, 0, 3);
	
	for(unsigned int pthatIter = 0; pthatIter < nPtHats; pthatIter++){
	  std::string nameGenJtPtPerPthat = "genJtPtPerPthat_" + qg[qgIter] + "_" + jetAlgo.at(iter) + "_" + centStr + "_Pthat" + std::to_string(config.GetPthat(pthatIter)) + "_h";
	  genJtPtPerPthat_p[iter][centIter][qgIter][pthatIter] = new TH1F(nameGenJtPtPerPthat.c_str(), ";Gen. Jet p_{T};Events", nJtPtBins, jtPtBins);
	}
      }

      for(Int_t ptEtaIter = 0; ptEtaIter < nJtPtEtaBins+1; ptEtaIter++){
	std::string ptEtaStr = "EtaInc";
	if(ptEtaIter > 0){
	  Int_t ptEtaLowInt = std::trunc(jtPtEtaBins[ptEtaIter-1]);
	  Int_t ptEtaHiInt = std::trunc(jtPtEtaBins[ptEtaIter]);
	  
	  Int_t ptEtaLowDec = std::trunc(jtPtEtaBins[ptEtaIter-1]*10 - ptEtaLowInt*10);
	  Int_t ptEtaHiDec = std::trunc(jtPtEtaBins[ptEtaIter]*10 - ptEtaHiInt*10);
	 
	  ptEtaStr = "Eta" + std::to_string(ptEtaLowInt) + "p" + std::to_string(ptEtaLowDec) + "to" + std::to_string(ptEtaHiInt) + "p" + std::to_string(ptEtaHiDec);
	}

	
	for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	  jtDPhiVPt_p[iter][centIter][ptEtaIter][qgIter] = new TH2F(Form("jtDPhiVPt_%s_%s_%s_%s_h", ptEtaStr.c_str(), qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#Delta#phi"), nJtPtBins, jtPtBins, nJtDPhiBins, jtDPhiBins);

	  jtDEtaVPt_p[iter][centIter][ptEtaIter][qgIter] = new TH2F(Form("jtDEtaVPt_%s_%s_%s_%s_h", ptEtaStr.c_str(), qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#Delta#eta"), nJtPtBins, jtPtBins, nJtDEtaBins, jtDEtaBins);

	  jtDRVPt_p[iter][centIter][ptEtaIter][qgIter] = new TH2F(Form("jtDRVPt_%s_%s_%s_%s_h", ptEtaStr.c_str(), qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#DeltaR"), nJtPtBins, jtPtBins, nJtDRBins, jtDRBins);
	  

	  for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	    jtRecoOverGenVPt_Mean_p[iter][centIter][ptEtaIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVPt_%s_%s_%sMean_%s_%s_h", ptEtaStr.c_str(), qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#mu_{Reco./Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);

	    jtRecoOverGenVPt_Res_p[iter][centIter][ptEtaIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVPt_%s_%s_%sRes_%s_%s_h", ptEtaStr.c_str(), qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#sigma_{Reco./Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);

	    jtDPhiVPt_Mean_p[iter][centIter][ptEtaIter][qgIter][mIter] = new TH1F(Form("jtDPhiVPt_%s_%s_%sMean_%s_%s_h", ptEtaStr.c_str(), qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#Delta#phi (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);

	    jtDPhiVPt_Res_p[iter][centIter][ptEtaIter][qgIter][mIter] = new TH1F(Form("jtDPhiVPt_%s_%s_%sRes_%s_%s_h", ptEtaStr.c_str(), qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#Delta#phi (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);
	 
	  }
	 
	  for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	    Int_t ptLowInt = std::trunc(jtPtBins[jtIter]);
	    Int_t ptHiInt = std::trunc(jtPtBins[jtIter+1]);
	    
	    Int_t ptLowDec = std::trunc(jtPtBins[jtIter]*10 - ptLowInt*10);
	    Int_t ptHiDec = std::trunc(jtPtBins[jtIter+1]*10 - ptHiInt*10);
	    
	    jtRecoOverGenVPt_MeanResPts_p[iter][centIter][ptEtaIter][qgIter][jtIter] = new TH1F(Form("jtRecoOverGenVPt_%s_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", ptEtaStr.c_str(), qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";(Reco. %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo.at(iter).c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 60, 0, 2);
	    jtRecoOverGenVPt_MeanResPts_p[iter][centIter][ptEtaIter][qgIter][jtIter]->Sumw2();
	    
	    jtRecoOverGenVPt_MeanResPts_COUNTS[iter][centIter][ptEtaIter][qgIter][jtIter] = 0;


	    jtDPhiVPt_MeanResPts_p[iter][centIter][ptEtaIter][qgIter][jtIter] = new TH1F(Form("jtDPhiVPt_%s_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", ptEtaStr.c_str(), qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";#Delta#phi (%s);Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo.at(iter).c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 120, -1, 1);
	    jtDPhiVPt_MeanResPts_p[iter][centIter][ptEtaIter][qgIter][jtIter]->Sumw2();
	  }
  

	}

	for(Int_t recoPtCutIter = 0; recoPtCutIter < nJetRecoCuts; recoPtCutIter++){
	  jtEffVPt_Mean_p[iter][centIter][ptEtaIter][recoPtCutIter] = new TH1F(Form("jtEffVPt_%s_Mean_RecoPtCut%d_%s_%s", ptEtaStr.c_str(), jetRecoCuts[recoPtCutIter], jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};Eff. (%s, p_{T} > %d)", jetAlgo.at(iter).c_str(), jetRecoCuts[recoPtCutIter]), nJtPtBins, jtPtBins);
	  
	  for(Int_t jtPtIter = 0; jtPtIter < nJtPtBins; jtPtIter++){
	    jtEffVPt_Num[iter][centIter][ptEtaIter][recoPtCutIter][jtPtIter] = 0;
	  }
	}
	
	for(Int_t jtPtIter = 0; jtPtIter < nJtPtBins; jtPtIter++){
	  jtEffVPt_Denom[iter][centIter][ptEtaIter][jtPtIter] = 0;
	}
      }
      
      for(Int_t qgIter = 0; qgIter < nQG; qgIter++){	  
	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
          jtRecoOverGenVEta_Mean_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVEta_%s_%sMean_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#mu_{Reco./Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);
        }

        for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
          jtRecoOverGenVEta_Res_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVEta_%s_%sRes_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#sigma_{Reco./Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtEtaBins, jtEtaBins);
        }



	for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	  Int_t etaLowInt = std::trunc(jtEtaBins[jtIter]);
	  Int_t etaHiInt = std::trunc(jtEtaBins[jtIter+1]);
	  
	  Int_t etaLowDec = std::trunc(jtEtaBins[jtIter]*10 - etaLowInt*10);
	  Int_t etaHiDec = std::trunc(jtEtaBins[jtIter+1]*10 - etaHiInt*10);
	  
	  
	  

	  jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoOverGenVEta_%s_MeanResPts_Eta%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";(Reco. %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo.at(iter).c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec), 60, 0, 2);
	  jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->Sumw2();

	  jtRecoOverGenVEta_MeanResPts_COUNTS[iter][centIter][qgIter][jtIter]++;
	}
      }
    }
  }
  
  const unsigned int nInputs = config.GetNInputs();
  unsigned int inputCounts[nInputs];
  unsigned int inputEvents[nInputs];
  for(unsigned int iter = 0; iter < nInputs; iter++){
    inputCounts[iter] = 0;
    inputEvents[iter] = 0;
  }

  for(unsigned int pthatIter = 0; pthatIter < nInputs; pthatIter++){
    std::cout << "Begin processing " << pthatIter+1 << "/" << nInputs << "(ptHat == " << config.GetInputPtHat(pthatIter) << "): \'" << config.GetInput(pthatIter) << "\'...." << std::endl;

    TFile* inFile_p = TFile::Open(config.GetInput(pthatIter).c_str(), "READ");
    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
    TTree* genTree_p = (TTree*)inFile_p->Get("HiGenParticleAna/hi");
    TTree* phoTree_p=0;
    if(config.GetIsPbPb()) phoTree_p = (TTree*)inFile_p->Get("ggHiNtuplizer/EventTree");
    else phoTree_p = (TTree*)inFile_p->Get("ggHiNtuplizerGED/EventTree");
    TTree* jetTree_p[nJetAlgo];
      
    for(Int_t iter = 0; iter < nJetAlgo; iter++){
      jetTree_p[iter] = (TTree*)inFile_p->Get(Form("%sJetAnalyzer/t", jetAlgo.at(iter).c_str()));
    }
      
    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("vz", 1);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("evt", 1);
    if(config.GetIsPbPb()) hiTree_p->SetBranchStatus("hiBin", 1);

    hiTree_p->SetBranchAddress("vz", &vz_);
    hiTree_p->SetBranchAddress("run", &run_);
    hiTree_p->SetBranchAddress("evt", &evt_);
    if(config.GetIsPbPb()) hiTree_p->SetBranchAddress("hiBin", &hiBin_);
      
      
    genTree_p->SetBranchStatus("*", 0);
    genTree_p->SetBranchStatus("pt", 1);
    genTree_p->SetBranchStatus("phi", 1);
    genTree_p->SetBranchStatus("eta", 1);
    genTree_p->SetBranchStatus("pdg", 1);
    
    genTree_p->SetBranchAddress("pt", &genPt_p);
    genTree_p->SetBranchAddress("phi", &genPhi_p);
    genTree_p->SetBranchAddress("eta", &genEta_p);
    genTree_p->SetBranchAddress("pdg", &genPDG_p);
      
          
    for(Int_t iter = 0; iter < nJetAlgo; iter++){
      jetTree_p[iter]->SetBranchStatus("*", 0);
      jetTree_p[iter]->SetBranchStatus("nref", 1);
      jetTree_p[iter]->SetBranchStatus("jtpt", 1);
      jetTree_p[iter]->SetBranchStatus("jteta", 1);
      jetTree_p[iter]->SetBranchStatus("jtphi", 1);
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
      
      jetTree_p[iter]->SetBranchAddress("nref", &nJt_[iter]);
      jetTree_p[iter]->SetBranchAddress("jtpt", jtPt_[iter]);
      jetTree_p[iter]->SetBranchAddress("jteta", jtEta_[iter]);
      jetTree_p[iter]->SetBranchAddress("jtphi", jtPhi_[iter]);
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
      jetTree_p[iter]->SetBranchAddress("gensubid", genJtSubID_[iter]);
    }
     
    phoTree_p->SetBranchStatus("*", 0);
    phoTree_p->SetBranchStatus("mcPt", 1);
    phoTree_p->SetBranchStatus("mcPhi", 1);
    phoTree_p->SetBranchStatus("mcEta", 1);
    phoTree_p->SetBranchStatus("mcPID", 1);
    phoTree_p->SetBranchStatus("mcMomPID", 1);
    phoTree_p->SetBranchStatus("phoEt", 1);
    phoTree_p->SetBranchStatus("phoEta", 1);
    phoTree_p->SetBranchStatus("phoPhi", 1);
    phoTree_p->SetBranchStatus("pho_seedTime", 1);
    phoTree_p->SetBranchStatus("pho_swissCrx", 1);
    phoTree_p->SetBranchStatus("phoHoverE", 1);
    phoTree_p->SetBranchStatus("phoSigmaIEtaIEta_2012", 1);
    phoTree_p->SetBranchStatus("pho_ecalClusterIsoR4", 1);
    phoTree_p->SetBranchStatus("pho_hcalRechitIsoR4", 1);
    phoTree_p->SetBranchStatus("pho_trackIsoR4PtCut20", 1);
    phoTree_p->SetBranchStatus("phoR9", 1);

    phoTree_p->SetBranchAddress("mcPt", &mcPt_p);
    phoTree_p->SetBranchAddress("mcPhi", &mcPhi_p);
    phoTree_p->SetBranchAddress("mcEta", &mcEta_p);
    phoTree_p->SetBranchAddress("mcPID", &mcPID_p);
    phoTree_p->SetBranchAddress("mcMomPID", &mcMomPID_p);
    phoTree_p->SetBranchAddress("phoEt", &phoPt_p);
    phoTree_p->SetBranchAddress("phoEta", &phoEta_p);
    phoTree_p->SetBranchAddress("phoPhi", &phoPhi_p);
    phoTree_p->SetBranchAddress("pho_seedTime", &phoSeedTime_p);
    phoTree_p->SetBranchAddress("pho_swissCrx", &phoSwissCrx_p);
    phoTree_p->SetBranchAddress("phoHoverE", &phoHoverE_p);
    phoTree_p->SetBranchAddress("phoSigmaIEtaIEta_2012", &phoSigmaIEtaIEta_2012_p);
    phoTree_p->SetBranchAddress("pho_ecalClusterIsoR4", &pho_ecalClusterIsoR4_p);
    phoTree_p->SetBranchAddress("pho_hcalRechitIsoR4", &pho_hcalRechitIsoR4_p);
    phoTree_p->SetBranchAddress("pho_trackIsoR4PtCut20", &pho_trackIsoR4PtCut20_p);
    phoTree_p->SetBranchAddress("phoR9", &phoR9_p);

 
    if(debugMode) std::cout << __LINE__ << std::endl;
    
    Int_t tempStartPos = 0;
    Int_t tempNEntries = jetTree_p[0]->GetEntries();
    
    const Int_t startPos = tempStartPos;
    const Int_t nEntries = tempNEntries;
    Int_t entryDiv = TMath::Max(1, ((Int_t)((nEntries-startPos)/20)));
    
    for(Int_t entry = startPos; entry < nEntries; entry++){
      if(entry%entryDiv == 0 && nEntries >= 10000) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
      
      hiTree_p->GetEntry(entry);
      genTree_p->GetEntry(entry);
      phoTree_p->GetEntry(entry);


      Float_t genTempLeadingPhoPt = -999;
      //      Float_t genTempLeadingPhoPhi = -999;
      //      Float_t genTempLeadingPhoEta = -999;

      Float_t recoTempLeadingPhoPt = -999;
      Float_t recoTempLeadingPhoPhi = -999;
      Float_t recoTempLeadingPhoEta = -999;

      Float_t tempLeadingMuPt_ = -999;
      Float_t tempLeadingMuPhi_ = -999;
      Float_t tempLeadingMuEta_ = -999;

      Float_t tempSubleadingMuPt_ = -999;
      Float_t tempSubleadingMuPhi_ = -999;
      Float_t tempSubleadingMuEta_ = -999;

      Float_t tempLeadingElePt_ = -999;
      Float_t tempLeadingElePhi_ = -999;
      Float_t tempLeadingEleEta_ = -999;

      Float_t tempSubleadingElePt_ = -999;
      Float_t tempSubleadingElePhi_ = -999;
      Float_t tempSubleadingEleEta_ = -999;

      Float_t tempZPhi_ = -999;

      if(config.GetIsGammaJet()){	
	const Int_t nGenPart = mcPt_p->size();

	for(Int_t genIter = 0; genIter < nGenPart; genIter++){
	  if(mcPID_p->at(genIter) != 22) continue;

	  //	  if(mcMomPID_p->at(genIter) != -999) continue;

	  if(genTempLeadingPhoPt < mcPt_p->at(genIter)){
            genTempLeadingPhoPt = mcPt_p->at(genIter);
	  }   
	}
	
	const Int_t nPho = phoPt_p->size();

	for(Int_t phoIter = 0; phoIter < nPho; phoIter++){
	  if(TMath::Abs(phoEta_p->at(phoIter)) > 1.44) continue;

	  if(phoSigmaIEtaIEta_2012_p->at(phoIter) < .002) continue;
	  if(TMath::Abs(phoSeedTime_p->at(phoIter)) > 3) continue;
	  if(phoSwissCrx_p->at(phoIter) > .9) continue;
	  if(phoSigmaIEtaIEta_2012_p->at(phoIter) > .01) continue;
	  if(pho_ecalClusterIsoR4_p->at(phoIter) + pho_hcalRechitIsoR4_p->at(phoIter) + pho_trackIsoR4PtCut20_p->at(phoIter) > 1) continue;
	  if(phoHoverE_p->at(phoIter) > .1) continue;

	  if(recoTempLeadingPhoPt < phoPt_p->at(phoIter)){
	    recoTempLeadingPhoPt = phoPt_p->at(phoIter);
	  }
	}
      }
      else if(config.GetIsZJet()){
	const Int_t nGenPart = genPt_p->size();
	
        for(Int_t genIter = 0; genIter < nGenPart; genIter++){                                                   
	  //	  if(TMath::Abs(genEta_p->at(genIter)) > 2.4) continue;
	  if(TMath::Abs(genPDG_p->at(genIter)) != 11 && TMath::Abs(genPDG_p->at(genIter)) != 13) continue; 

	  if(TMath::Abs(genPDG_p->at(genIter)) == 11){
	    if(tempLeadingElePt_ < genPt_p->at(genIter)){
	      tempSubleadingElePt_ = tempLeadingElePt_;
	      tempSubleadingElePhi_ = tempLeadingElePhi_;
	      tempSubleadingEleEta_ = tempLeadingEleEta_;

	      tempLeadingElePt_ = genPt_p->at(genIter);
	      tempLeadingElePhi_ = genPhi_p->at(genIter);
	      tempLeadingEleEta_ = genEta_p->at(genIter);
	    }
	    else if(tempSubleadingElePt_ < genPt_p->at(genIter)){
	      tempSubleadingElePt_ = genPt_p->at(genIter);
	      tempSubleadingElePhi_ = genPhi_p->at(genIter);
	      tempSubleadingEleEta_ = genEta_p->at(genIter);
	    }
	  }
	   
	  if(TMath::Abs(genPDG_p->at(genIter)) == 13){
	    if(tempLeadingMuPt_ < genPt_p->at(genIter)){
	      tempSubleadingMuPt_ = tempLeadingMuPt_;
	      tempSubleadingMuPhi_ = tempLeadingMuPhi_;
	      tempSubleadingMuEta_ = tempLeadingMuEta_;

	      tempLeadingMuPt_ = genPt_p->at(genIter);
	      tempLeadingMuPhi_ = genPhi_p->at(genIter);
	      tempLeadingMuEta_ = genEta_p->at(genIter);
	    }
	    else if(tempSubleadingMuPt_ < genPt_p->at(genIter)){
	      tempSubleadingMuPt_ = genPt_p->at(genIter);
	      tempSubleadingMuPhi_ = genPhi_p->at(genIter);
	      tempSubleadingMuEta_ = genEta_p->at(genIter);
	    }
	  }
	}                                          
      }
      
      if(config.GetIsGammaJet()){
	if(recoTempLeadingPhoPt < config.GetMinGammaPt()) continue;
	if(TMath::Abs(recoTempLeadingPhoEta) > 1.44) continue;
	if(config.GetGammaPtHatStagger() + config.GetInputPtHat(pthatIter) > recoTempLeadingPhoPt) continue;
      }
      

      if(config.GetIsZJet()){
	Bool_t isGoodMu = false;
	Bool_t isGoodEle = false;
	if(tempLeadingMuPt_ > 10 && tempSubleadingMuPt_ > 5) isGoodMu = true;
	if(tempLeadingElePt_ > 10 && tempSubleadingElePt_ > 5) isGoodEle = true;

	if(!isGoodMu && ! isGoodEle) continue;

	if(isGoodMu){
	  TLorentzVector mu1, mu2;
	  mu1.SetPtEtaPhiM(tempLeadingMuPt_, tempLeadingMuEta_, tempLeadingMuPhi_, muMass);
	  mu2.SetPtEtaPhiM(tempSubleadingMuPt_, tempSubleadingMuEta_, tempSubleadingMuPhi_, muMass);
	  
	  TLorentzVector z = mu1+mu2;

	  if(z.Pt() < 20) continue;

	  tempZPhi_ = z.Phi();
	}
	else if(isGoodEle){
	  TLorentzVector ele1, ele2;
	  ele1.SetPtEtaPhiM(tempLeadingElePt_, tempLeadingEleEta_, tempLeadingElePhi_, eleMass);
	  ele2.SetPtEtaPhiM(tempSubleadingElePt_, tempSubleadingEleEta_, tempSubleadingElePhi_, eleMass);

	  TLorentzVector z = ele1+ele2;

	  if(z.Pt() < 20) continue;

	  tempZPhi_ = z.Phi();
	}
      }

      if(TMath::Abs(vz_) > 15) continue;
      
      if(debugMode) std::cout << __LINE__ << std::endl;

      Int_t centPos = config.GetCentBinFromHiBin(hiBin_);
      if(centPos == -1) continue;

      for(Int_t iter = 0; iter < nJetAlgo; iter++){
	jetTree_p[iter]->GetEntry(entry);
      }

      if(config.GetDoPthatStagger()){
	if(config.GetInputPtHatPos(pthatIter) != config.GetNPthats()-1){
	  if(ptHat_[0] > config.GetPthat(config.GetInputPtHatPos(pthatIter)+1)) continue;
	}
      }
            
      
      if(config.GetIsPbPb()){
	for(Int_t algoIter = 0; algoIter < nJetAlgo; algoIter++){
	  for(Int_t jtIter = 0; jtIter < nJt_[algoIter]; jtIter++){
	    jtPt_[algoIter][jtIter] *= getPtEtaJetResidualCorr(jtPt_[algoIter][jtIter], jtEta_[algoIter][jtIter], hiBin_/2.);
	  }
	}
      }
      else{
	for(Int_t algoIter = 0; algoIter < nJetAlgo; algoIter++){
          for(Int_t jtIter = 0; jtIter < nJt_[algoIter]; jtIter++){
            jtPt_[algoIter][jtIter] *= .99;
          }
        }
      }
      

      Double_t ptHatWeight = config.GetPtHatWeight(ptHat_[0]);

      ptHat_Unweighted_h->Fill(ptHat_[0]);
      ptHat_Weighted_h->Fill(ptHat_[0], ptHatWeight);

      ptHat_Unweighted_PerPthat_h[config.GetInputPtHatPos(pthatIter)]->Fill(ptHat_[0]);
      ptHat_Weighted_PerPthat_h[config.GetInputPtHatPos(pthatIter)]->Fill(ptHat_[0], ptHatWeight);
      
      if(config.GetIsGammaJet()){
	if(genTempLeadingPhoPt > 0){
	  genPhoPt_Unweighted_h[centPos]->Fill(genTempLeadingPhoPt);
	  genPhoPt_Weighted_h[centPos]->Fill(genTempLeadingPhoPt, ptHatWeight);
	}
	
	if(recoTempLeadingPhoPt > 0){
	  recoPhoPt_Unweighted_h[centPos]->Fill(recoTempLeadingPhoPt);
	  recoPhoPt_Weighted_h[centPos]->Fill(recoTempLeadingPhoPt, ptHatWeight);
	}
      }

      if(debugMode) std::cout << __LINE__ << std::endl;

      Float_t hiBinWeight = 1;
      if(config.GetIsPbPb()) hiBinWeight = findNormNcoll(hiBin_, centPos, config.GetCentBins());
      //	  hiBinWeight = 1;
      
      
      inputEvents[pthatIter]++;

      for(Int_t algoIter = 0; algoIter < nJetAlgo; algoIter++){
	genSort(nGenJt_[algoIter], genJtPt_[algoIter], genJtPhi_[algoIter], genJtEta_[algoIter], genJtMatchIndex_[algoIter], genJtSubID_[algoIter]);
	
	//MODDING FOR REFPT USAGE
	for(Int_t jtIter = 0; jtIter < nJt_[algoIter]; jtIter++){
	  if(TMath::Abs(refEta_[algoIter][jtIter]) > config.GetJtEtaMax()) continue;
	  if(refPt_[algoIter][jtIter] < 5.0) continue;
	  if(refSubID_[algoIter][jtIter] != 0) continue;

	  if(refPt_[algoIter][jtIter] < config.GetJtPtLow()) continue;
	  if(refPt_[algoIter][jtIter] > config.GetJtPtHi()) continue;

	  //	  if(config.GetJtWeight(pthatIter, refPt_[algoIter][jtIter], recoTempLeadingPhoPt) < .1) continue;
	  
	  Float_t jtWeight = config.GetJtWeight(pthatIter, refPt_[algoIter][jtIter], refEta_[algoIter][jtIter], recoTempLeadingPhoPt);

	  Double_t truncPtHatWeight = config.GetTruncPtHatWeight(ptHat_[0], refPt_[algoIter][jtIter]);
	  //	  truncPtHatWeight = 1;

	  if(config.GetDoPthatStagger()){
	    if(jtWeight < .1) continue;
	  }
	  
	  if(config.GetIsGammaJet() && recoTempLeadingPhoPt > 0){
	    if(getDR(recoTempLeadingPhoEta, recoTempLeadingPhoPhi, refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter]) < 0.4) continue;
	  }
	  
	  if(config.GetIsZJet()){
	    if(tempLeadingMuPt_ > 10){
	      if(getDR(tempLeadingMuEta_, tempLeadingMuPhi_, refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter]) < 0.4) continue;
	    }
	    if(tempSubleadingMuPt_ > 5){	  
	      if(getDR(tempSubleadingMuEta_, tempSubleadingMuPhi_, refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter]) < 0.4) continue;
	    }
	    if(tempLeadingElePt_ > 10){
	      if(getDR(tempLeadingEleEta_, tempLeadingElePhi_, refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter]) < 0.4) continue;
	    }
	    if(tempSubleadingElePt_ > 5){
	      if(getDR(tempSubleadingEleEta_, tempSubleadingElePhi_, refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter]) < 0.4) continue;
	    }

	    if(TMath::Abs(getDPHI(tempZPhi_, refPhi_[algoIter][jtIter])) < 7.*TMath::Pi()/8.) continue;
	  }
	
       

	  Int_t qgPos[2] = {0, -1};
	  if(TMath::Abs(refPartFlav_[algoIter][jtIter]) < 9) qgPos[1] = 1;
	  else if(TMath::Abs(refPartFlav_[algoIter][jtIter]) == 21) qgPos[1] = 2;	
  
	  Int_t etaPos[2] = {0, -1};

	  for(Int_t jtPtEtaIter = 0; jtPtEtaIter < nJtPtEtaBins; jtPtEtaIter++){
	    if(TMath::Abs(refEta_[algoIter][jtIter]) <= jtPtEtaBins[jtPtEtaIter+1]){
	      etaPos[1] = jtPtEtaIter+1;
	      break;
	    }
	  }

	  if(etaPos[1] == -1) std::cout << "Error, no eta pos for jet w/ refeta \'" << refEta_[algoIter][jtIter] << "\'." << std::endl;

	  Float_t qgWeight = 1;

	  if(config.GetIsDijet()){
	    if(qgPos[1] == 1){
	      if(refPt_[algoIter][jtIter] < 300) qgWeight *= qRatio_p->GetBinContent(qRatio_p->FindBin(refPt_[algoIter][jtIter]));
	      else qgWeight *= qRatio_p->GetBinContent(qRatio_p->GetNbinsX());
	    }
	    else if(qgPos[1] == 2){
	      if(refPt_[algoIter][jtIter] < 300) qgWeight *= gRatio_p->GetBinContent(gRatio_p->FindBin(refPt_[algoIter][jtIter]));
	      else qgWeight *= gRatio_p->GetBinContent(gRatio_p->GetNbinsX());
	    }
	  }

	  //	  qgWeight = 1.;

      
	  for(Int_t jtPtEtaIter = 0; jtPtEtaIter < 2; jtPtEtaIter++){
	    for(Int_t qgIter = 0; qgIter < 2; qgIter++){
	      if(qgPos[qgIter] == -1) continue;
	    
	      if(jtPtEtaIter == 0){
		jtRecoVGen_p[algoIter][centPos][qgPos[qgIter]]->Fill(refPt_[algoIter][jtIter], jtPt_[algoIter][jtIter]);
		jtRecoOverGenVPt_p[algoIter][centPos][qgPos[qgIter]]->Fill(refPt_[algoIter][jtIter], jtPt_[algoIter][jtIter]/refPt_[algoIter][jtIter]);
	      }

	      jtDPhiVPt_p[algoIter][centPos][qgPos[qgIter]][etaPos[jtPtEtaIter]]->Fill(refPt_[algoIter][jtIter], TMath::Abs(getDPHI(refPhi_[algoIter][jtIter], jtPhi_[algoIter][jtIter])), qgWeight*hiBinWeight*truncPtHatWeight);

	      jtDEtaVPt_p[algoIter][centPos][qgPos[qgIter]][etaPos[jtPtEtaIter]]->Fill(refPt_[algoIter][jtIter], TMath::Abs(refEta_[algoIter][jtIter] - jtEta_[algoIter][jtIter]), qgWeight*hiBinWeight*truncPtHatWeight);

	      jtDRVPt_p[algoIter][centPos][qgPos[qgIter]][etaPos[jtPtEtaIter]]->Fill(refPt_[algoIter][jtIter], getDR(refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter], jtEta_[algoIter][jtIter], jtPhi_[algoIter][jtIter]), qgWeight*hiBinWeight*truncPtHatWeight);

	      for(Int_t jtIter2 = 0; jtIter2 < nJtPtBins; jtIter2++){
		if(refPt_[algoIter][jtIter] > jtPtBins[jtIter2] && refPt_[algoIter][jtIter] < jtPtBins[jtIter2+1]){
		  jtRecoOverGenVPt_MeanResPts_p[algoIter][centPos][etaPos[jtPtEtaIter]][qgPos[qgIter]][jtIter2]->Fill(jtPt_[algoIter][jtIter]/refPt_[algoIter][jtIter], qgWeight*hiBinWeight*truncPtHatWeight);
		  jtDPhiVPt_MeanResPts_p[algoIter][centPos][etaPos[jtPtEtaIter]][qgPos[qgIter]][jtIter2]->Fill(getDPHI(jtPhi_[algoIter][jtIter], refPhi_[algoIter][jtIter]), qgWeight*hiBinWeight*truncPtHatWeight);


		  if(qgIter == 0 && jtPtEtaIter == 0){
		    inputCounts[pthatIter]++;
		  }
		  jtRecoOverGenVPt_MeanResPts_COUNTS[algoIter][centPos][etaPos[jtPtEtaIter]][qgPos[qgIter]][jtIter2]++;
		  break;
		}
	      }
	    }
	  }
	  
	  if(refPt_[algoIter][jtIter] > config.GetJtEtaPtThresh()){
	    for(Int_t qgIter = 0; qgIter < 2; qgIter++){
	      if(qgPos[qgIter] == -1) continue;
	      
	      for(Int_t jtIter2 = 0; jtIter2 < nJtEtaBins; jtIter2++){
		if(refEta_[algoIter][jtIter] > jtEtaBins[jtIter2] && refEta_[algoIter][jtIter] < jtEtaBins[jtIter2+1]){
		  jtRecoOverGenVEta_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(jtPt_[algoIter][jtIter]/refPt_[algoIter][jtIter], qgWeight*hiBinWeight*truncPtHatWeight);
		  jtRecoOverGenVEta_MeanResPts_COUNTS[algoIter][centPos][qgPos[qgIter]][jtIter2]++;
		  break;
		}
	      }
	    }
	  } 
	}
      
	for(Int_t jtIter = 0; jtIter < nGenJt_[algoIter]; jtIter++){
	  if(TMath::Abs(genJtEta_[algoIter][jtIter]) > config.GetJtEtaMax()) continue;
	  if(genJtPt_[algoIter][jtIter] < 5.0) continue;
	  if(genJtSubID_[algoIter][jtIter] != 0) continue;
	  
	  genJtPtPerPthat_p[algoIter][centPos][0][pthatIter]->Fill(genJtPt_[algoIter][jtIter]);

	  if(genJtPt_[algoIter][jtIter] < jtPtLow || genJtPt_[algoIter][jtIter] > jtPtHi) continue;

	  if(config.GetIsZJet()){
	    if(tempLeadingMuPt_ > 10){
	      if(getDR(tempLeadingMuEta_, tempLeadingMuPhi_, genJtEta_[algoIter][jtIter], genJtPhi_[algoIter][jtIter]) < 0.4) continue;
	    }
	    if(tempSubleadingMuPt_ > 5){	  
	      if(getDR(tempSubleadingMuEta_, tempSubleadingMuPhi_, genJtEta_[algoIter][jtIter], genJtPhi_[algoIter][jtIter]) < 0.4) continue;
	    }
	    if(tempLeadingElePt_ > 10){
	      if(getDR(tempLeadingEleEta_, tempLeadingElePhi_, genJtEta_[algoIter][jtIter], genJtPhi_[algoIter][jtIter]) < 0.4) continue;
	    }
	    if(tempSubleadingElePt_ > 5){
	      if(getDR(tempSubleadingEleEta_, tempSubleadingElePhi_, genJtEta_[algoIter][jtIter], genJtPhi_[algoIter][jtIter]) < 0.4) continue;
	    }
	  }


	  Int_t etaPos[2] = {0, -1};

	  for(Int_t jtPtEtaIter = 0; jtPtEtaIter < nJtPtEtaBins; jtPtEtaIter++){
	    if(TMath::Abs(genJtEta_[algoIter][jtIter]) <= jtPtEtaBins[jtPtEtaIter+1]){
	      etaPos[1] = jtPtEtaIter+1;
	      break;
	    }
	  }

	  if(etaPos[1] == -1) std::cout << "Error, no eta pos for jet w/ geneta \'" << genJtEta_[algoIter][jtIter] << "\'." << std::endl;
	  
	  Int_t jtPtPos = -1;
	  for(Int_t jtPtIter = 0; jtPtIter < nJtPtBins; jtPtIter++){
	    if(genJtPt_[algoIter][jtIter] >= jtPtBins[jtPtIter] && genJtPt_[algoIter][jtIter] < jtPtBins[jtPtIter+1]){
	      jtPtPos = jtPtIter;
	      break;
	    }
	  }

	  if(jtPtPos == -1) std::cout << "Error, no jtpt pos for jet w/ genpt \'" << genJtPt_[algoIter][jtIter] << "\'." << std::endl;

	  for(Int_t jtPtEtaIter = 0; jtPtEtaIter < 2; jtPtEtaIter++){
	    jtEffVPt_Denom[algoIter][centPos][etaPos[jtPtEtaIter]][jtPtPos]+=hiBinWeight;
	    
	    if(genJtMatchIndex_[algoIter][jtIter] >= 0){
	      for(Int_t recoPtCutIter = 0; recoPtCutIter < nJetRecoCuts; recoPtCutIter++){
		if(jtPt_[algoIter][genJtMatchIndex_[algoIter][jtIter]] >= jetRecoCuts[recoPtCutIter]){
		  jtEffVPt_Num[algoIter][centPos][etaPos[jtPtEtaIter]][recoPtCutIter][jtPtPos]+=hiBinWeight;
		}
	      }
	    }

	  }

	}

      }
    }
    inFile_p->Close();
  }    

  if(debugMode) std::cout << __LINE__ << std::endl;

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){

      for(Int_t jtPtEtaIter = 0; jtPtEtaIter < nJtPtEtaBins+1; jtPtEtaIter++){

	for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	  for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	    Float_t tempMean[nMeanFit];
	    Float_t tempMeanErr[nMeanFit];
	    Float_t tempRes[nMeanFit];
	    Float_t tempResErr[nMeanFit];
	    
	    tempMean[0] = jtRecoOverGenVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter]->GetMean();
	    tempMeanErr[0] = jtRecoOverGenVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter]->GetMeanError();	
	    tempRes[0] = jtRecoOverGenVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter]->GetStdDev();
	    tempResErr[0] = jtRecoOverGenVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter]->GetStdDevError();
	    
	    if(jtRecoOverGenVPt_MeanResPts_COUNTS[iter][centIter][jtPtEtaIter][qgIter][jtIter] < 300. || jtRecoOverGenVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter]->GetMaximum() < 200) jtRecoOverGenVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter]->Rebin(2);
	    
	    FitGauss(jtRecoOverGenVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter], tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1], true/*config.GetDoWeights()*/);
	    
	    for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	      jtRecoOverGenVPt_Mean_p[iter][centIter][jtPtEtaIter][qgIter][mIter]->SetBinContent(jtIter+1, tempMean[mIter]);
	      jtRecoOverGenVPt_Mean_p[iter][centIter][jtPtEtaIter][qgIter][mIter]->SetBinError(jtIter+1, tempMeanErr[mIter]);
	      jtRecoOverGenVPt_Res_p[iter][centIter][jtPtEtaIter][qgIter][mIter]->SetBinContent(jtIter+1, tempRes[mIter]);
	      jtRecoOverGenVPt_Res_p[iter][centIter][jtPtEtaIter][qgIter][mIter]->SetBinError(jtIter+1, tempResErr[mIter]);
	    }
	  }

	  for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	    Float_t tempMean[nMeanFit];
	    Float_t tempMeanErr[nMeanFit];
	    Float_t tempRes[nMeanFit];
	    Float_t tempResErr[nMeanFit];
	    
	    tempMean[0] = jtDPhiVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter]->GetMean();
	    tempMeanErr[0] = jtDPhiVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter]->GetMeanError();	
	    tempRes[0] = jtDPhiVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter]->GetStdDev();
	    tempResErr[0] = jtDPhiVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter]->GetStdDevError();
	    	    
	    FitGauss(jtDPhiVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter], tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1], true/*config.GetDoWeights()*/);
	    
	    for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	      jtDPhiVPt_Mean_p[iter][centIter][jtPtEtaIter][qgIter][mIter]->SetBinContent(jtIter+1, tempMean[mIter]);
	      jtDPhiVPt_Mean_p[iter][centIter][jtPtEtaIter][qgIter][mIter]->SetBinError(jtIter+1, tempMeanErr[mIter]);
	      jtDPhiVPt_Res_p[iter][centIter][jtPtEtaIter][qgIter][mIter]->SetBinContent(jtIter+1, tempRes[mIter]);
	      jtDPhiVPt_Res_p[iter][centIter][jtPtEtaIter][qgIter][mIter]->SetBinError(jtIter+1, tempResErr[mIter]);
	    }
	  }
	}

	for(Int_t jtRecoPtCutIter = 0; jtRecoPtCutIter < nJetRecoCuts; jtRecoPtCutIter++){

	  for(Int_t jtPtIter = 0; jtPtIter < nJtPtBins; jtPtIter++){
	    jtEffVPt_Mean_p[iter][centIter][jtPtEtaIter][jtRecoPtCutIter]->SetBinContent(jtPtIter+1, jtEffVPt_Num[iter][centIter][jtPtEtaIter][jtRecoPtCutIter][jtPtIter]/jtEffVPt_Denom[iter][centIter][jtPtEtaIter][jtPtIter]);
	    jtEffVPt_Mean_p[iter][centIter][jtPtEtaIter][jtRecoPtCutIter]->SetBinError(jtPtIter+1, 0);
	  }
	}
      }

      if(debugMode) std::cout << __LINE__ << std::endl;
      for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	  Float_t tempMean[nMeanFit];
	  Float_t tempMeanErr[nMeanFit];
	  Float_t tempRes[nMeanFit];
	  Float_t tempResErr[nMeanFit];

	  if(debugMode) std::cout << __LINE__ << std::endl;

	  tempMean[0] = jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMean();
	  tempMeanErr[0] = jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMeanError();	
	  tempRes[0] = jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDev();
	  tempResErr[0] = jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetStdDevError();
	  
	  if(debugMode) std::cout << __LINE__ << std::endl;

	  if(jtRecoOverGenVEta_MeanResPts_COUNTS[iter][centIter][qgIter][jtIter] < 300. || jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->GetMaximum() < 200) jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->Rebin(2);

	  if(debugMode) std::cout << __LINE__ << std::endl;

	  FitGauss(jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter], tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1], true/*config.GetDoWeights()*/);

	  if(debugMode) std::cout << __LINE__ << std::endl;
	
	  for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	    jtRecoOverGenVEta_Mean_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempMean[mIter]);
	    jtRecoOverGenVEta_Mean_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempMeanErr[mIter]);
	    jtRecoOverGenVEta_Res_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempRes[mIter]);
	    jtRecoOverGenVEta_Res_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempResErr[mIter]);
	  }

	  if(debugMode) std::cout << __LINE__ << std::endl;
	}

	if(debugMode) std::cout << __LINE__ << std::endl;

      }
    }
  }

  if(debugMode) std::cout << __LINE__ << std::endl;    

  outFile_p->cd();

  TDirectory* dir_p = outFile_p->GetDirectory("ptHatDir");
  if(dir_p){
    dir_p->cd();
  }
  else{
    dir_p = outFile_p->mkdir("ptHatDir");
    dir_p->cd();
  }

  ptHat_Unweighted_h->Write("", TObject::kOverwrite);
  ptHat_Weighted_h->Write("", TObject::kOverwrite);
  
  for(unsigned int iter = 0; iter < nPtHats; iter++){
    ptHat_Unweighted_PerPthat_h[iter]->Write("", TObject::kOverwrite);

    unsigned int scalePos = ptHat_Weighted_PerPthat_h[iter]->GetMaximumBin();
    Double_t scaleVal = ptHat_Weighted_h->GetBinContent(scalePos);

    ptHat_Unweighted_PerPthat_h[iter]->Scale(scaleVal/ptHat_Unweighted_PerPthat_h[iter]->GetMaximum());

    std::string nameUnweighted = "ptHat_Unweighted_PerPthat_" + std::to_string(config.GetPthat(iter)) + "_RescaleToWeighted_h";

    ptHat_Unweighted_PerPthat_h[iter]->Write(nameUnweighted.c_str(), TObject::kOverwrite);

    delete ptHat_Unweighted_PerPthat_h[iter];

    ptHat_Weighted_PerPthat_h[iter]->Write("", TObject::kOverwrite);
    delete ptHat_Weighted_PerPthat_h[iter];
  }

  delete ptHat_Unweighted_h;
  delete ptHat_Weighted_h;

  outFile_p->cd();

  if(debugMode) std::cout << __LINE__ << std::endl;

  if(config.GetIsGammaJet()){
    TDirectory* dir_p = outFile_p->GetDirectory("phoDir");
    if(dir_p){
      dir_p->cd();
    }
    else{
      dir_p = outFile_p->mkdir("phoDir");
      dir_p->cd();
    }

    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      std::string centStr = "PP";
      if(config.GetIsPbPb()) centStr = "Cent" + std::to_string(config.GetCentBinFromPos(centIter)) + "to" + std::to_string(config.GetCentBinFromPos(centIter+1));

      std::string genNameUnweighted = "genPhoPt_Unweighted_Rescale_" + centStr + "_h";
      std::string genNameWeighted = "genPhoPt_Weighted_Rescale_" + centStr + "_h";

      std::string recoNameUnweighted = "recoPhoPt_Unweighted_Rescale_" + centStr + "_h";
      std::string recoNameWeighted = "recoPhoPt_Weighted_Rescale_" + centStr + "_h";
      
      genPhoPt_Unweighted_h[centIter]->Write("", TObject::kOverwrite);
      genPhoPt_Weighted_h[centIter]->Write("", TObject::kOverwrite);

      recoPhoPt_Unweighted_h[centIter]->Write("", TObject::kOverwrite);
      recoPhoPt_Weighted_h[centIter]->Write("", TObject::kOverwrite);

      genPhoPt_Unweighted_h[centIter]->Scale(1./genPhoPt_Unweighted_h[centIter]->Integral());
      genPhoPt_Weighted_h[centIter]->Scale(1./genPhoPt_Weighted_h[centIter]->Integral());

      recoPhoPt_Unweighted_h[centIter]->Scale(1./recoPhoPt_Unweighted_h[centIter]->Integral());
      recoPhoPt_Weighted_h[centIter]->Scale(1./recoPhoPt_Weighted_h[centIter]->Integral());

      genPhoPt_Unweighted_h[centIter]->Write(genNameUnweighted.c_str(), TObject::kOverwrite);
      genPhoPt_Weighted_h[centIter]->Write(genNameWeighted.c_str(), TObject::kOverwrite);

      recoPhoPt_Unweighted_h[centIter]->Write(recoNameUnweighted.c_str(), TObject::kOverwrite);
      recoPhoPt_Weighted_h[centIter]->Write(recoNameWeighted.c_str(), TObject::kOverwrite);

      delete genPhoPt_Unweighted_h[centIter];
      delete genPhoPt_Weighted_h[centIter];

      delete recoPhoPt_Unweighted_h[centIter];
      delete recoPhoPt_Weighted_h[centIter];
    }
  }

  if(debugMode) std::cout << __LINE__ << std::endl;

  outFile_p->cd();

  TF1* csn_p = new TF1("csn_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))", 30, 150);
  csn_p->FixParameter(0, 0.061);
  csn_p->FixParameter(1, 1.24);
  csn_p->FixParameter(2, 8.08);


  TF1* csn2_p = new TF1("csn2_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))", 30, 150);
  csn2_p->FixParameter(0, 0.061 - .001);
  csn2_p->FixParameter(1, 1.24 - .04);
  csn2_p->FixParameter(2, 8.08 - .15);


  TF1* csn3_p = new TF1("csn3_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))", 30, 150);
  csn3_p->FixParameter(0, 0.061 + .001);
  csn3_p->FixParameter(1, 1.24 + .04);
  csn3_p->FixParameter(2, 8.08 + .15);


  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    TDirectory* dir_p = outFile_p->GetDirectory(Form("%s", jetAlgo.at(iter).c_str()));
    if(dir_p){
      dir_p->cd();
    }
    else{
      dir_p = outFile_p->mkdir(Form("%s", jetAlgo.at(iter).c_str()));
      dir_p->cd();
    }


    if(debugMode) std::cout << __LINE__ << std::endl;

    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	jtRecoVGen_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoOverGenVPt_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	
	for(unsigned int pthatIter = 0; pthatIter < nPtHats; pthatIter++){
	  genJtPtPerPthat_p[iter][centIter][qgIter][pthatIter]->Write("", TObject::kOverwrite);
	}

	if(debugMode) std::cout << __LINE__ << std::endl;
	
	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  for(Int_t jtPtEtaIter = 0; jtPtEtaIter < nJtPtEtaBins+1; jtPtEtaIter++){

	    if(mIter == 0){
	      jtDPhiVPt_p[iter][centIter][jtPtEtaIter][qgIter]->Write("", TObject::kOverwrite);
	      jtDEtaVPt_p[iter][centIter][jtPtEtaIter][qgIter]->Write("", TObject::kOverwrite);
	      jtDRVPt_p[iter][centIter][jtPtEtaIter][qgIter]->Write("", TObject::kOverwrite);
	    }

	    jtRecoOverGenVPt_Mean_p[iter][centIter][jtPtEtaIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	    jtDPhiVPt_Mean_p[iter][centIter][jtPtEtaIter][qgIter][mIter]->Write("", TObject::kOverwrite);

	    if(mIter == 1 && qgIter == 0 && centIter == 0 && jtPtEtaIter == 0){
	      jtRecoOverGenVPt_Res_p[iter][centIter][jtPtEtaIter][qgIter][mIter]->Fit("csn_p", "Q", "", 30, 150);
	      jtRecoOverGenVPt_Res_p[iter][centIter][jtPtEtaIter][qgIter][mIter]->Fit("csn2_p", "Q", "", 30, 150);
	      jtRecoOverGenVPt_Res_p[iter][centIter][jtPtEtaIter][qgIter][mIter]->Fit("csn3_p", "Q", "", 30, 150);
	    }

	    jtRecoOverGenVPt_Res_p[iter][centIter][jtPtEtaIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	    jtDPhiVPt_Res_p[iter][centIter][jtPtEtaIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  }

	  if(debugMode) std::cout << __LINE__ << std::endl;

	  jtRecoOverGenVEta_Mean_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRecoOverGenVEta_Res_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	}

	if(debugMode) std::cout << __LINE__ << std::endl;
	
	for(Int_t jtPtEtaIter = 0; jtPtEtaIter < nJtPtEtaBins+1; jtPtEtaIter++){
	  for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	    jtRecoOverGenVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	    jtDPhiVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	  }
	}

	if(debugMode) std::cout << __LINE__ << std::endl;
	  
	for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	  jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	}	
      }

      for(Int_t jtPtEtaIter = 0; jtPtEtaIter < nJtPtEtaBins+1; jtPtEtaIter++){
	for(Int_t jtRecoPtCutIter = 0; jtRecoPtCutIter < nJetRecoCuts; jtRecoPtCutIter++){
	  jtEffVPt_Mean_p[iter][centIter][jtPtEtaIter][jtRecoPtCutIter]->Write("", TObject::kOverwrite);
	}
      }
    }
  }


  delete csn_p;
  delete csn2_p;
  delete csn3_p;

  if(debugMode) std::cout << __LINE__ << std::endl;
  
  config.WriteConfigParamsToRootFile(outFile_p);


  if(debugMode) std::cout << __LINE__ << std::endl;


  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	delete jtRecoVGen_p[iter][centIter][qgIter];
	delete jtRecoOverGenVPt_p[iter][centIter][qgIter];

        for(unsigned int pthatIter = 0; pthatIter < nPtHats; pthatIter++){
          delete genJtPtPerPthat_p[iter][centIter][qgIter][pthatIter];
        }

	if(debugMode) std::cout << __LINE__ << std::endl;

	for(Int_t jtPtEtaIter = 0; jtPtEtaIter < nJtPtEtaBins+1; jtPtEtaIter++){
	  delete jtDPhiVPt_p[iter][centIter][jtPtEtaIter][qgIter];
	  delete jtDEtaVPt_p[iter][centIter][jtPtEtaIter][qgIter];
	  delete jtDRVPt_p[iter][centIter][jtPtEtaIter][qgIter];

	  for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	    delete jtRecoOverGenVPt_Mean_p[iter][centIter][jtPtEtaIter][qgIter][mIter];
	    delete jtRecoOverGenVPt_Res_p[iter][centIter][jtPtEtaIter][qgIter][mIter];

	    delete jtDPhiVPt_Mean_p[iter][centIter][jtPtEtaIter][qgIter][mIter];
	    delete jtDPhiVPt_Res_p[iter][centIter][jtPtEtaIter][qgIter][mIter];
	  }

	  if(debugMode) std::cout << __LINE__ << std::endl;

	  for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	    delete jtRecoOverGenVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter];
	    delete jtDPhiVPt_MeanResPts_p[iter][centIter][jtPtEtaIter][qgIter][jtIter];
	  }
	}
     
	if(debugMode) std::cout << __LINE__ << std::endl;

	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  delete jtRecoOverGenVEta_Mean_p[iter][centIter][qgIter][mIter];
	  delete jtRecoOverGenVEta_Res_p[iter][centIter][qgIter][mIter];
	}
	if(debugMode) std::cout << __LINE__ << std::endl;
	
	for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	  delete jtRecoOverGenVEta_MeanResPts_p[iter][centIter][qgIter][jtIter];
	}
      }
    
      for(Int_t jtPtEtaIter = 0; jtPtEtaIter < nJtPtEtaBins+1; jtPtEtaIter++){      
	for(Int_t jtRecoPtCutIter = 0; jtRecoPtCutIter < nJetRecoCuts; jtRecoPtCutIter++){
	  delete jtEffVPt_Mean_p[iter][centIter][jtPtEtaIter][jtRecoPtCutIter];
	}
      }
    }
  }


  if(debugMode) std::cout << __LINE__ << std::endl;

  for(unsigned int iter = 0; iter < nInputs; iter++){
    std::cout << "Input " << iter << " Events, Counts: " << inputEvents[iter] << ", " << inputCounts[iter] << std::endl;
  }
    
  outFile_p->Close();
  delete outFile_p;

  ratioFile_p->Close();
  delete ratioFile_p;

  return 0;
}


int main(int argc, char *argv[])
{
  if(argc != 2){
    std::cout << "Usage: makeJECHist_Prototype.exe <inConfigFile>" << std::endl;
    std::cout << "Number of args given: " << argc << std::endl;
    for(int iter = 0; iter < argc; iter++){
      std::cout << "  argv[" << iter << "]: " << argv[iter] << std::endl;
    }
    return -1;
  }

  return makeJECHist_Prototype(argv[1]);
}
