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

const std::string storeStr = "/store";
const std::string xrootdStr = "root://xrootd.unl.edu/";


const Bool_t debugMode = false;
const Bool_t doGetBkg = false;

const Int_t nJtCat = 3;
const std::string jtCat[nJtCat] = {"Eta2", "Eta1", "Dijet"};

const Int_t nMeanFit = 2;
const std::string meanFit[nMeanFit] = {"", "Fit"};

const Int_t nQG = 3;
const std::string qg[nQG] = {"Inc", "Q", "G"};

const Int_t nMaxJets = 500;

const Int_t nCentBins = 1;

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

  //  if(f1_p->GetProb() > .01) return;
  if(!isPbPb) return;

  //  return;

  for(Int_t fitIter = 0; fitIter < 3; fitIter++){
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


int makeJECHist_Prototype(const std::string inConfigFileName)
{
  jecConfigParser config;
  if(!config.SetConfigParser(inConfigFileName)) return 1;

  const std::string outName = "fileForConfigTest.root";
  TFile* outFile_p = new TFile(outName.c_str(), "RECREATE");
  
  TFile* tempJetInFile_p = TFile::Open(config.GetInput(0).c_str(), "READ");
  std::vector<std::string> jetAlgoInFile = returnRootFileContentsList(tempJetInFile_p, "TTree", "JetAnalyzer");
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

  //pthat30  
  const Int_t nJtPtBins = config.GetNJtPtBins();
  const Float_t jtPtLow = config.GetJtPtLow();
  const Float_t jtPtHi = config.GetJtPtHi();
  Double_t jtPtBins[nJtPtBins+1];
  if(config.GetDoJtPtLogBins()) getLogBins(jtPtLow, jtPtHi, nJtPtBins, jtPtBins);
  else getLinBins(jtPtLow, jtPtHi, nJtPtBins, jtPtBins);

  TH2F* jtRecoVGen_p[nJetAlgo][nCentBins][nQG];
  TH2F* jtRecoOverGenVPt_p[nJetAlgo][nCentBins][nQG];
  TH1F* jtRecoOverGenVPt_Mean_p[nJetAlgo][nCentBins][nQG][nMeanFit];
  TH1F* jtRecoOverGenVPt_Res_p[nJetAlgo][nCentBins][nQG][nMeanFit];
  TH1F* jtRecoOverGenVPt_MeanResPts_p[nJetAlgo][nCentBins][nQG][nJtPtBins];

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      std::string centStr = "PP";

      for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	jtRecoVGen_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoVGen_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; Reco. %s Jet p_{T}", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins, nJtPtBins, jtPtBins);
	
	jtRecoOverGenVPt_p[iter][centIter][qgIter] = new TH2F(Form("jtRecoOverGenVPt_%s_%s_%s_h", qg[qgIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T}; (Reco. %s Jet p_{T})/(Gen. Jet p_{T})", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins, 30, 0, 3);
	
	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  jtRecoOverGenVPt_Mean_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVPt_%s_%sMean_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#mu_{Reco./Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);
	}
	
	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter] = new TH1F(Form("jtRecoOverGenVPt_%s_%sRes_%s_%s_h", qg[qgIter].c_str(), meanFit[mIter].c_str(), jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";Gen. Jet p_{T};#sigma_{Reco./Gen.} (%s)", jetAlgo.at(iter).c_str()), nJtPtBins, jtPtBins);	 
	}
	
	
	for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	  Int_t ptLowInt = std::trunc(jtPtBins[jtIter]);
	  Int_t ptHiInt = std::trunc(jtPtBins[jtIter+1]);
	  
	  Int_t ptLowDec = std::trunc(jtPtBins[jtIter]*10 - ptLowInt*10);
	  Int_t ptHiDec = std::trunc(jtPtBins[jtIter+1]*10 - ptHiInt*10);
	  
	  jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter] = new TH1F(Form("jtRecoOverGenVPt_%s_MeanResPts_Pt%dp%dTo%dp%d_%s_%s_h", qg[qgIter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo.at(iter).c_str(), centStr.c_str()), Form(";(Reco. %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo.at(iter).c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 45, 0, 3);
	  jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Sumw2();
	}
      }
    }
  }
  

  for(unsigned int pthatIter = 0; pthatIter < config.GetNPthats(); pthatIter++){
    TFile* inFile_p = TFile::Open(config.GetInput(pthatIter).c_str(), "READ");
    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
    TTree* genTree_p = (TTree*)inFile_p->Get("HiGenParticleAna/hi");
    TTree* jetTree_p[nJetAlgo];
      
    for(Int_t iter = 0; iter < nJetAlgo; iter++){
      jetTree_p[iter] = (TTree*)inFile_p->Get(Form("%sJetAnalyzer/t", jetAlgo.at(iter).c_str()));
    }
      
    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("vz", 1);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("evt", 1);
    
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
    }
      
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
      
      if(TMath::Abs(vz_) > 15) continue;
      
      if(debugMode) std::cout << __LINE__ << std::endl;
      
      for(Int_t iter = 0; iter < nJetAlgo; iter++){
	jetTree_p[iter]->GetEntry(entry);
      }

      Int_t centPos = 0;
      
      if(debugMode) std::cout << __LINE__ << std::endl;
      
      for(Int_t algoIter = 0; algoIter < nJetAlgo; algoIter++){
	genSort(nGenJt_[algoIter], genJtPt_[algoIter], genJtPhi_[algoIter], genJtEta_[algoIter]);
	
	//MODDING FOR REFPT USAGE
	for(Int_t jtIter = 0; jtIter < nJt_[algoIter]; jtIter++){
	  if(TMath::Abs(refEta_[algoIter][jtIter]) > 4.9) continue;
	  if(refPt_[algoIter][jtIter] < 5.0) continue;
	  if(refSubID_[algoIter][jtIter] != 0) continue;
	  
	  Int_t qgPos[2] = {0, -1};
	  if(TMath::Abs(refPartFlav_[algoIter][jtIter]) < 9) qgPos[1] = 1;
	  else if(TMath::Abs(refPartFlav_[algoIter][jtIter]) == 21) qgPos[1] = 2;
	  
	  for(Int_t qgIter = 0; qgIter < 2; qgIter++){
	    if(qgPos[qgIter] == -1) continue;
	    
	    jtRecoVGen_p[algoIter][centPos][qgPos[qgIter]]->Fill(refPt_[algoIter][jtIter], jtPt_[algoIter][jtIter]);
	    jtRecoOverGenVPt_p[algoIter][centPos][qgPos[qgIter]]->Fill(refPt_[algoIter][jtIter], jtPt_[algoIter][jtIter]/refPt_[algoIter][jtIter]);
	  
	    for(Int_t jtIter2 = 0; jtIter2 < nJtPtBins; jtIter2++){
	      if(refPt_[algoIter][jtIter] > jtPtBins[jtIter2] && refPt_[algoIter][jtIter] < jtPtBins[jtIter2+1]){
		jtRecoOverGenVPt_MeanResPts_p[algoIter][centPos][qgPos[qgIter]][jtIter2]->Fill(jtPt_[algoIter][jtIter]/refPt_[algoIter][jtIter]);
	      }
	    }
	  }
	}
      }
    }
    inFile_p->Close();
  }    

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
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
	  
	  FitGauss(jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter], config.GetIsPbPb(), tempMean[1], tempMeanErr[1], tempRes[1], tempResErr[1]);
	
	  for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	    jtRecoOverGenVPt_Mean_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempMean[mIter]);
	    jtRecoOverGenVPt_Mean_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempMeanErr[mIter]);
	    jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter]->SetBinContent(jtIter+1, tempRes[mIter]);
	    jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter]->SetBinError(jtIter+1, tempResErr[mIter]);
	  }
	}
      }
    }
  }

  outFile_p->cd();

  for(unsigned int iter = 0; iter < config.GetNPthats(); iter++){
    TNamed pathStr(Form("pathStr_%s", std::to_string(config.GetPthat(iter)).c_str()), config.GetInput(iter).c_str());
		   pathStr.Write("", TObject::kOverwrite);
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


    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	jtRecoVGen_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);
	jtRecoOverGenVPt_p[iter][centIter][qgIter]->Write("", TObject::kOverwrite);

	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  jtRecoOverGenVPt_Mean_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	  jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter]->Write("", TObject::kOverwrite);
	}
	
	for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	  jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter]->Write("", TObject::kOverwrite);
	}
      }
    }
  }


  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      for(Int_t qgIter = 0; qgIter < nQG; qgIter++){
	delete jtRecoVGen_p[iter][centIter][qgIter];
	delete jtRecoOverGenVPt_p[iter][centIter][qgIter];

	for(Int_t mIter = 0; mIter < nMeanFit; mIter++){
	  delete jtRecoOverGenVPt_Mean_p[iter][centIter][qgIter][mIter];
	  delete jtRecoOverGenVPt_Res_p[iter][centIter][qgIter][mIter];
	}

	for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	  delete jtRecoOverGenVPt_MeanResPts_p[iter][centIter][qgIter][jtIter];
	}
      }
    }
  }
    
  outFile_p->Close();
  delete outFile_p;

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
