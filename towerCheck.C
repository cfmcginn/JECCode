#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TDatime.h"
#include "TMath.h"

#include "dirent.h"
#include <vector>
#include <string>

#include "getLogBins.h"
#include "etaPhiFunc.h"

const Int_t nMaxJets = 500;
const Int_t nMaxTowers = 5000;

const Int_t nCentBins = 4;
const Int_t centBins[nCentBins+1] = {200, 100, 60, 20, 0};

void towerCheck(const std::string inFileName, Bool_t isPbPb)
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

  Int_t nJt_;
  Float_t jtPt_[nMaxJets];
  Float_t jtRawPt_[nMaxJets];
  Float_t jtEta_[nMaxJets];
  Float_t jtPhi_[nMaxJets];

  Int_t hiBin_ = -1;
  Float_t vz_;
  UInt_t run_;
  ULong64_t evt_;

  Int_t nTow_;
  Float_t towEt_[nMaxTowers];
  Float_t towVsPt_[nMaxTowers];
  Float_t towPhi_[nMaxTowers];
  Float_t towEta_[nMaxTowers];

  const Int_t nJtPtBins = 29;
  const Float_t jtPtLow = 7.999;
  const Float_t jtPtHi = 300.001;
  Double_t jtPtBins[nJtPtBins+1];
  getLogBins(jtPtLow, jtPtHi, nJtPtBins, jtPtBins);

  Int_t nCentBinsTemp = nCentBins;
  if(!isPbPb) nCentBinsTemp = 1;

  const Int_t nCentBins2 = nCentBinsTemp;

  TH2F* jtTowerOverRawVRawPt_p[nCentBins2];
  TH1F* jtTowerOverRawVRawPt_Mean_p[nCentBins2];
  TH1F* jtTowerOverRawVRawPt_MeanPts_p[nCentBins2][nJtPtBins];

  TH2F* jtVsTowerOverRawVRawPt_p[nCentBins2];
  TH1F* jtVsTowerOverRawVRawPt_Mean_p[nCentBins2];
  TH1F* jtVsTowerOverRawVRawPt_MeanPts_p[nCentBins2][nJtPtBins];

  for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
    std::string centStr = "PP";
    if(isPbPb) centStr = Form("cent%dto%d", centBins[centIter+1]/2, centBins[centIter]/2);

    jtTowerOverRawVRawPt_p[centIter] = new TH2F(Form("jtTowerOverRawVRawPt_%s_h", centStr.c_str()), Form(";Jet p_{T}^{Raw};Tower Sum/p_{T}^{Raw}"), nJtPtBins, jtPtBins, 100, 0, 10);

    jtTowerOverRawVRawPt_Mean_p[centIter] = new TH1F(Form("jtTowerOverRawVRawPt_Mean_%s_h", centStr.c_str()), Form(";Jet p_{T}^{Raw};Tower Sum/p_{T}^{Raw}"), nJtPtBins, jtPtBins);

    jtVsTowerOverRawVRawPt_p[centIter] = new TH2F(Form("jtVsTowerOverRawVRawPt_%s_h", centStr.c_str()), Form(";Jet p_{T}^{Raw};VsTower Sum/p_{T}^{Raw}"), nJtPtBins, jtPtBins, 100, 0, 10);

    jtVsTowerOverRawVRawPt_Mean_p[centIter] = new TH1F(Form("jtVsTowerOverRawVRawPt_Mean_%s_h", centStr.c_str()), Form(";Jet p_{T}^{Raw};VsTower Sum/p_{T}^{Raw}"), nJtPtBins, jtPtBins);

    for(Int_t ptIter = 0; ptIter < nJtPtBins; ptIter++){
      Int_t ptLowInt = std::trunc(jtPtBins[ptIter]);
      Int_t ptHiInt = std::trunc(jtPtBins[ptIter+1]);

      Int_t ptLowDec = std::trunc(jtPtBins[ptIter]*10 - ptLowInt*10);
      Int_t ptHiDec = std::trunc(jtPtBins[ptIter+1]*10 - ptHiInt*10);

      jtVsTowerOverRawVRawPt_MeanPts_p[centIter][ptIter] = new TH1F(Form("jtVsTowerOverRawVRawPt_MeanPts_Pt%dp%dTo%dp%d_%s_h", ptLowInt, ptLowDec, ptHiInt, ptHiDec, centStr.c_str()), Form(";VsTower Sum/p_{T}^{Raw};Events (%d.%d<p_{T,Raw}<%d.%d)", ptLowInt, ptLowDec, ptHiInt, ptHiDec), 100, 0, 10);

      jtTowerOverRawVRawPt_MeanPts_p[centIter][ptIter] = new TH1F(Form("jtTowerOverRawVRawPt_MeanPts_Pt%dp%dTo%dp%d_%s_h", ptLowInt, ptLowDec, ptHiInt, ptHiDec, centStr.c_str()), Form(";Tower Sum/p_{T}^{Raw};Events (%d.%d<p_{T,Raw}<%d.%d)", ptLowInt, ptLowDec, ptHiInt, ptHiDec), 100, 0, 10);
    }
  }

  const Int_t numberOfFiles = (Int_t)listOfFiles_p->size();
  Int_t fileDiv = ((Int_t)(numberOfFiles/10));
  if(fileDiv < 1) fileDiv = 1;

  for(Int_t fileIter = 0; fileIter < numberOfFiles; fileIter++){
    if(fileIter%fileDiv == 0) std::cout << "File # " << fileIter << "/" << numberOfFiles << std::endl;
    TFile* inFile_p = new TFile(listOfFiles_p->at(fileIter).c_str(), "READ");

    TTree* jetTree_p = (TTree*)inFile_p->Get("akVs4CaloJetAnalyzer/t");
    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
    TTree* towTree_p = (TTree*)inFile_p->Get("rechitanalyzer/tower");

    jetTree_p->SetBranchStatus("*", 0);
    jetTree_p->SetBranchStatus("nref", 1);
    jetTree_p->SetBranchStatus("jtpt", 1);
    jetTree_p->SetBranchStatus("rawpt", 1);
    jetTree_p->SetBranchStatus("jteta", 1);
    jetTree_p->SetBranchStatus("jtphi", 1);

    jetTree_p->SetBranchAddress("nref", &nJt_);
    jetTree_p->SetBranchAddress("jtpt", jtPt_);
    jetTree_p->SetBranchAddress("rawpt", jtRawPt_);
    jetTree_p->SetBranchAddress("jteta", jtEta_);
    jetTree_p->SetBranchAddress("jtphi", jtPhi_);

    hiTree_p->SetBranchStatus("*", 0);
    if(isPbPb) hiTree_p->SetBranchStatus("hiBin", 1);
    hiTree_p->SetBranchStatus("vz", 1);
    hiTree_p->SetBranchStatus("run", 1);
    hiTree_p->SetBranchStatus("evt", 1);

    if(isPbPb) hiTree_p->SetBranchAddress("hiBin", &hiBin_);
    hiTree_p->SetBranchAddress("vz", &vz_);
    hiTree_p->SetBranchAddress("run", &run_);
    hiTree_p->SetBranchAddress("evt", &evt_);

    towTree_p->SetBranchStatus("*", 0);
    towTree_p->SetBranchStatus("n", 1);
    towTree_p->SetBranchStatus("et", 1);
    towTree_p->SetBranchStatus("vsPt", 1);
    towTree_p->SetBranchStatus("phi", 1);
    towTree_p->SetBranchStatus("eta", 1);

    towTree_p->SetBranchAddress("n", &nTow_);
    towTree_p->SetBranchAddress("et", towEt_);
    towTree_p->SetBranchAddress("vsPt", towVsPt_);
    towTree_p->SetBranchAddress("phi", towPhi_);
    towTree_p->SetBranchAddress("eta", towEta_);


    const Int_t nEntries = jetTree_p->GetEntries();
    Int_t entryDiv = ((Int_t)(nEntries/10));

    for(Int_t entry = 0; entry < nEntries; entry++){
      if(entry%entryDiv == 0 && nEntries >= 10000) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
      
      jetTree_p->GetEntry(entry);
      hiTree_p->GetEntry(entry);
      towTree_p->GetEntry(entry);

      if(TMath::Abs(vz_) > 15) continue;

      Int_t centPos = -1;
      if(isPbPb)
        for(Int_t iter = 0; iter < nCentBins2; iter++){
          if(hiBin_ < centBins[iter] && hiBin_ >= centBins[iter+1]){
            centPos = iter;
            break;
          }
        }
      else centPos = 0;

      const Int_t nTowerSum = nJt_;     
      Float_t towerSums[nTowerSum];
      Float_t vsTowerSums[nTowerSum];
      for(Int_t towIter = 0; towIter < nTowerSum; towIter++){
	towerSums[towIter] = 0;
	vsTowerSums[towIter] = 0;
      }

      for(Int_t towIter = 0; towIter < nTow_; towIter++){
	if(TMath::Abs(towEta_[towIter]) > 2.0) continue;
	for(Int_t jtIter = 0; jtIter < nJt_; jtIter++){
	  if(TMath::Abs(jtEta_[jtIter]) > 1.6) continue;
	  if(getDR(towEta_[towIter], towPhi_[towIter], jtEta_[jtIter], jtPhi_[jtIter]) < 0.4){
	    towerSums[jtIter] += towEt_[towIter];
	    vsTowerSums[jtIter] += towVsPt_[towIter];
	    break;
	  }
	}
      }
     
      for(Int_t jtIter = 0; jtIter < nJt_; jtIter++){
	if(TMath::Abs(jtEta_[jtIter]) > 1.6) continue;

	jtTowerOverRawVRawPt_p[centPos]->Fill(jtRawPt_[jtIter], towerSums[jtIter]/jtRawPt_[jtIter]);
	jtVsTowerOverRawVRawPt_p[centPos]->Fill(jtRawPt_[jtIter], vsTowerSums[jtIter]/jtRawPt_[jtIter]);	

	for(Int_t jtIter2 = 0; jtIter2 < nJtPtBins; jtIter2++){
	  if(jtRawPt_[jtIter] > jtPtBins[jtIter2] && jtRawPt_[jtIter] < jtPtBins[jtIter2+1]){
	    jtTowerOverRawVRawPt_MeanPts_p[centPos][jtIter2]->Fill(towerSums[jtIter]/jtRawPt_[jtIter]);
	    jtVsTowerOverRawVRawPt_MeanPts_p[centPos][jtIter2]->Fill(vsTowerSums[jtIter]/jtRawPt_[jtIter]);
	    break;
	  }
	}
      }
      
    }

    inFile_p->Close();
    delete inFile_p;
  }
  
  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
    for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
      jtTowerOverRawVRawPt_Mean_p[centIter]->SetBinContent(jtIter+1, jtTowerOverRawVRawPt_MeanPts_p[centIter][jtIter]->GetMean());
      jtTowerOverRawVRawPt_Mean_p[centIter]->SetBinError(jtIter+1, jtTowerOverRawVRawPt_MeanPts_p[centIter][jtIter]->GetMeanError());

      jtVsTowerOverRawVRawPt_Mean_p[centIter]->SetBinContent(jtIter+1, jtVsTowerOverRawVRawPt_MeanPts_p[centIter][jtIter]->GetMean());
      jtVsTowerOverRawVRawPt_Mean_p[centIter]->SetBinError(jtIter+1, jtVsTowerOverRawVRawPt_MeanPts_p[centIter][jtIter]->GetMeanError());
    }
  }

  outFile_p->cd();

  for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
    jtTowerOverRawVRawPt_p[centIter]->Write("", TObject::kOverwrite);
    jtTowerOverRawVRawPt_Mean_p[centIter]->Write("", TObject::kOverwrite);

    jtVsTowerOverRawVRawPt_p[centIter]->Write("", TObject::kOverwrite);
    jtVsTowerOverRawVRawPt_Mean_p[centIter]->Write("", TObject::kOverwrite);

    for(Int_t ptIter = 0; ptIter < nJtPtBins; ptIter++){
      jtTowerOverRawVRawPt_MeanPts_p[centIter][ptIter]->Write("", TObject::kOverwrite);
      jtVsTowerOverRawVRawPt_MeanPts_p[centIter][ptIter]->Write("", TObject::kOverwrite);
    }
  }

  for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
    delete jtTowerOverRawVRawPt_p[centIter];
    delete jtTowerOverRawVRawPt_Mean_p[centIter];

    delete jtVsTowerOverRawVRawPt_p[centIter];
    delete jtVsTowerOverRawVRawPt_Mean_p[centIter];

    for(Int_t ptIter = 0; ptIter < nJtPtBins; ptIter++){
      delete jtTowerOverRawVRawPt_MeanPts_p[centIter][ptIter];
      delete jtVsTowerOverRawVRawPt_MeanPts_p[centIter][ptIter];
    }
  }

  outFile_p->Close();
  delete outFile_p;

  listOfFiles_p->clear();
  delete listOfFiles_p;

  return;
}
