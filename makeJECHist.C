#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TF1.h"
#include "TDirectory.h"

#include <string>
#include <iostream>
#include <math.h>

#include "getLogBins.h"
#include "getLinBins.h"
#include "etaPhiFunc.h"

const Int_t nJetAlgo = 4;
const std::string jetAlgo[nJetAlgo] = {"akVs4Calo", "akPu4Calo", "akCs4PF", "akPu4PF"};

const Int_t nMaxJets = 500;

const Int_t nCentBins = 8;
const Int_t centBins[nCentBins+1] = {200, 140, 100, 80, 60, 40, 20, 10, 0};

const Int_t nMaxGen = 100000;

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

  if(TMath::Abs(mean - 1.0) < 0.01) return;
  if(f1_p->GetProb() > .01) return;

  Float_t max = -1;
  Int_t maxBin = -1;

  for(Int_t iter = 0; iter < hist_p->GetNbinsX(); iter++){
    if(hist_p->GetBinContent(iter+1) > max){
      max = hist_p->GetBinContent(iter+1);
      maxBin = iter+1;
    }
  }

  Float_t fitLow = -1;
  Float_t fitHi = -1;

  for(Int_t iter = 1; iter < hist_p->GetNbinsX(); iter++){
    Int_t tempBin = maxBin - iter;
    if(tempBin < 1) tempBin = 1;

    if(hist_p->Integral(tempBin, maxBin+iter) > .9*hist_p->Integral() || tempBin == 1){
      fitLow = hist_p->GetBinCenter(tempBin);
      fitHi = hist_p->GetBinCenter(maxBin+iter);
      break;
    }
  }

  hist_p->Fit("f1_p", "LL Q M", "", fitLow, fitHi);

  mean = f1_p->GetParameter(1);
  res = f1_p->GetParameter(2);

  meanErr = f1_p->GetParError(1);
  resErr = f1_p->GetParError(2);

  if(TMath::Abs(mean - 1.0) < 0.01) return;
  if(f1_p->GetProb() > .01) return;

  Int_t meanBin = hist_p->FindBin(hist_p->GetMean());
  fitLow = -1;
  fitHi = -1;

  for(Int_t iter = 1; iter < hist_p->GetNbinsX(); iter++){
    Int_t tempBin = meanBin - iter;
    if(tempBin < 1) tempBin = 1;
    
    if(hist_p->Integral(tempBin, meanBin+iter) > .9*hist_p->Integral() || tempBin == 1){
      fitLow = hist_p->GetBinCenter(tempBin);
      fitHi = hist_p->GetBinCenter(meanBin+iter);
      break;
    }
  }

  hist_p->Fit("f1_p", "LL Q M", "", fitLow, fitHi);

  mean = f1_p->GetParameter(1);
  res = f1_p->GetParameter(2);

  meanErr = f1_p->GetParError(1);
  resErr = f1_p->GetParError(2);

  delete f1_p;

  return;
}


void makeJECHist(const std::string inFileName)
{
  /*
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");

  TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
  TTree* genTree_p = (TTree*)inFile_p->Get("HiGenParticleAna/hi");
  TTree* jetTree_p[nJetAlgo];

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    jetTree_p[iter] = (TTree*)inFile_p->Get(Form("%sJetAnalyzer/t", jetAlgo[iter].c_str()));
  }
  */


  TChain* hiTree_p = new TChain("hiEvtAnalyzer/HiTree");
  TChain* genTree_p = new TChain("HiGenParticleAna/hi");
  TChain* jetTree_p[nJetAlgo];

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    jetTree_p[iter] = new TChain(Form("%sJetAnalyzer/t", jetAlgo[iter].c_str()));
  }

  hiTree_p->Add(inFileName.c_str());
  genTree_p->Add(inFileName.c_str());

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    jetTree_p[iter]->Add(inFileName.c_str());
  }

  Int_t hiBin_;

  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("hiBin", 1);

  hiTree_p->SetBranchAddress("hiBin", &hiBin_);
  /*
  std::vector<float>* genPt_p = 0;
  std::vector<float>* genEta_p = 0;
  std::vector<float>* genPhi_p = 0;
  std::vector<int>* genPDG_p = 0;

  genTree_p->SetBranchStatus("*", 0);
  genTree_p->SetBranchStatus("pt", 1);
  genTree_p->SetBranchStatus("phi", 1);
  genTree_p->SetBranchStatus("eta", 1);
  genTree_p->SetBranchStatus("pdg", 1);

  genTree_p->SetBranchAddress("pt", &genPt_p);
  genTree_p->SetBranchAddress("phi", &genPhi_p);
  genTree_p->SetBranchAddress("eta", &genEta_p);
  genTree_p->SetBranchAddress("pdg", &genPDG_p);
  */
  Int_t nJt_[nJetAlgo]; 
  Float_t jtPt_[nJetAlgo][nMaxJets];
  Float_t jtRawPt_[nJetAlgo][nMaxJets];
  Float_t jtEta_[nJetAlgo][nMaxJets];
  Float_t jtPhi_[nJetAlgo][nMaxJets];
  Float_t refPt_[nJetAlgo][nMaxJets];
  Float_t refEta_[nJetAlgo][nMaxJets];
  Int_t refSubID_[nJetAlgo][nMaxJets];

  Int_t nGenJt_[nJetAlgo]; 
  Float_t genJtPt_[nJetAlgo][nMaxJets];
  Float_t genJtEta_[nJetAlgo][nMaxJets];
  Float_t genJtPhi_[nJetAlgo][nMaxJets];
  Int_t genJtMatchIndex_[nJetAlgo][nMaxJets];

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    jetTree_p[iter]->SetBranchStatus("*", 0);
    jetTree_p[iter]->SetBranchStatus("nref", 1);
    jetTree_p[iter]->SetBranchStatus("jtpt", 1);
    jetTree_p[iter]->SetBranchStatus("rawpt", 1);
    jetTree_p[iter]->SetBranchStatus("jteta", 1);
    jetTree_p[iter]->SetBranchStatus("jtphi", 1);
    jetTree_p[iter]->SetBranchStatus("refpt", 1);
    jetTree_p[iter]->SetBranchStatus("refeta", 1);
    jetTree_p[iter]->SetBranchStatus("subid", 1);
    jetTree_p[iter]->SetBranchStatus("ngen", 1);
    jetTree_p[iter]->SetBranchStatus("genpt", 1);
    jetTree_p[iter]->SetBranchStatus("geneta", 1);
    jetTree_p[iter]->SetBranchStatus("genphi", 1);
    jetTree_p[iter]->SetBranchStatus("genmatchindex", 1);

    jetTree_p[iter]->SetBranchAddress("nref", &nJt_[iter]);
    jetTree_p[iter]->SetBranchAddress("jtpt", jtPt_[iter]);
    jetTree_p[iter]->SetBranchAddress("rawpt", jtRawPt_[iter]);
    jetTree_p[iter]->SetBranchAddress("jteta", jtEta_[iter]);
    jetTree_p[iter]->SetBranchAddress("jtphi", jtPhi_[iter]);
    jetTree_p[iter]->SetBranchAddress("refpt", refPt_[iter]);
    jetTree_p[iter]->SetBranchAddress("refeta", refEta_[iter]);
    jetTree_p[iter]->SetBranchAddress("subid", refSubID_[iter]);
    jetTree_p[iter]->SetBranchAddress("ngen", &nGenJt_[iter]);
    jetTree_p[iter]->SetBranchAddress("genpt", genJtPt_[iter]);
    jetTree_p[iter]->SetBranchAddress("geneta", genJtEta_[iter]);
    jetTree_p[iter]->SetBranchAddress("genphi", genJtPhi_[iter]);
    jetTree_p[iter]->SetBranchAddress("genmatchindex", genJtMatchIndex_[iter]);
  }
  
  const Int_t nJtPtBins = 29;
  const Float_t jtPtLow = 9.999;
  const Float_t jtPtHi = 300.001;
  Double_t jtPtBins[nJtPtBins+1];
  getLogBins(jtPtLow, jtPtHi, nJtPtBins, jtPtBins);

  const Int_t nJtEtaBins = 20;
  const Float_t jtEtaLow = -2.0;
  const Float_t jtEtaHi = 2.0;
  Double_t jtEtaBins[nJtEtaBins+1];
  getLinBins(jtEtaLow, jtEtaHi, nJtEtaBins, jtEtaBins);

  for(Int_t iter = 0; iter < nJtPtBins+1; iter++){
    std::cout << "Pt Bin, val: " << iter << ", " << jtPtBins[iter] << std::endl;
  }

  for(Int_t iter = 0; iter < nJtEtaBins+1; iter++){
    std::cout << "Eta Bin, val: " << iter << ", " << jtEtaBins[iter] << std::endl;
  }

  TH2F* jtRecoOverGenVPt_p[nJetAlgo][nCentBins];
  TH1F* jtRecoOverGenVPt_Mean_p[nJetAlgo][nCentBins];
  TH1F* jtRecoOverGenVPt_Res_p[nJetAlgo][nCentBins];
  TH1F* jtRecoOverGenVPt_MeanResPts_p[nJetAlgo][nCentBins][nJtPtBins];

  TH2F* jtRecoOverRawVPt_p[nJetAlgo][nCentBins];
  TH1F* jtRecoOverRawVPt_Mean_p[nJetAlgo][nCentBins];
  TH1F* jtRecoOverRawVPt_Res_p[nJetAlgo][nCentBins];
  TH1F* jtRecoOverRawVPt_MeanResPts_p[nJetAlgo][nCentBins][nJtPtBins];

  TH2F* jtRawOverGenVPt_p[nJetAlgo][nCentBins];
  TH1F* jtRawOverGenVPt_Mean_p[nJetAlgo][nCentBins];
  TH1F* jtRawOverGenVPt_Res_p[nJetAlgo][nCentBins];
  TH1F* jtRawOverGenVPt_MeanResPts_p[nJetAlgo][nCentBins][nJtPtBins];

  TH2F* jtEffVPt_p[nJetAlgo][nCentBins];
  TH1F* jtEffVPt_Mean_p[nJetAlgo][nCentBins];
  TH1F* jtEffVPt_MeanResPts_p[nJetAlgo][nCentBins][nJtPtBins];

  TH2F* jtEffVPtEta_p[nJetAlgo][nCentBins];
  TH2F* jtEffVPtEta_Denom_p[nJetAlgo][nCentBins];

  TH2F* jtRecoOverGenVEta_p[nJetAlgo][nCentBins];
  TH1F* jtRecoOverGenVEta_Mean_p[nJetAlgo][nCentBins];
  TH1F* jtRecoOverGenVEta_Res_p[nJetAlgo][nCentBins];
  TH1F* jtRecoOverGenVEta_MeanResPts_p[nJetAlgo][nCentBins][nJtPtBins];

  TH2F* jtRecoOverRawVEta_p[nJetAlgo][nCentBins];
  TH1F* jtRecoOverRawVEta_Mean_p[nJetAlgo][nCentBins];
  TH1F* jtRecoOverRawVEta_Res_p[nJetAlgo][nCentBins];
  TH1F* jtRecoOverRawVEta_MeanResPts_p[nJetAlgo][nCentBins][nJtPtBins];

  TH2F* jtRawOverGenVEta_p[nJetAlgo][nCentBins];
  TH1F* jtRawOverGenVEta_Mean_p[nJetAlgo][nCentBins];
  TH1F* jtRawOverGenVEta_Res_p[nJetAlgo][nCentBins];
  TH1F* jtRawOverGenVEta_MeanResPts_p[nJetAlgo][nCentBins][nJtPtBins];

  TH2F* jtEffVEta_p[nJetAlgo][nCentBins];
  TH1F* jtEffVEta_Mean_p[nJetAlgo][nCentBins];
  TH1F* jtEffVEta_MeanResPts_p[nJetAlgo][nCentBins][nJtPtBins];

  Int_t nNoLeadEle = 0;
  Int_t nNoSubleadEle = 0;

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      jtRecoOverGenVPt_p[iter][centIter] = new TH2F(Form("jtRecoOverGenVPt_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet p_{T}; (Reco. %s Jet p_{T})/(Gen. Jet p_{T})", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins, 30, 0, 3);

      jtRecoOverRawVPt_p[iter][centIter] = new TH2F(Form("jtRecoOverRawVPt_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Reco. Jet p_{T}; (Reco %s Jet p_{T})/(Raw Jet p_{T})", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins, 30, 0, 3);

      jtRawOverGenVPt_p[iter][centIter] = new TH2F(Form("jtRawOverGenVPt_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet p_{T}; (Raw %s Jet p_{T})/(Gen Jet p_{T})", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins, 30, 0, 3);

      jtEffVPt_p[iter][centIter] = new TH2F(Form("jtEffVPt_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet p_{T};Eff. (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins, 2, -0.5, 1.5);

      jtEffVPtEta_p[iter][centIter] = new TH2F(Form("jtEffVPtEta_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet #eta;Gen. Jet p_{T} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins, nJtPtBins, jtPtBins);
      jtEffVPtEta_Denom_p[iter][centIter] = new TH2F(Form("jtEffVPtEta_Denom_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet #eta;Gen. Jet p_{T} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins, nJtPtBins, jtPtBins);

      jtRecoOverGenVPt_Mean_p[iter][centIter] = new TH1F(Form("jtRecoOverGenVPt_Mean_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet p_{T};#mu_{Reco./Gen.} (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins);

      jtRecoOverRawVPt_Mean_p[iter][centIter] = new TH1F(Form("jtRecoOverRawVPt_Mean_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Reco. Jet p_{T};#mu_{Reco./Raw} (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins);

      jtRawOverGenVPt_Mean_p[iter][centIter] = new TH1F(Form("jtRawOverGenVPt_Mean_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet p_{T};#mu_{Raw/Gen.} (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins);

      jtEffVPt_Mean_p[iter][centIter] = new TH1F(Form("jtEffVPt_Mean_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet p_{T};Eff. (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins);

      jtRecoOverGenVPt_Res_p[iter][centIter] = new TH1F(Form("jtRecoOverGenVPt_Res_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet p_{T};#sigma_{Reco./Gen.} (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins);

      jtRecoOverRawVPt_Res_p[iter][centIter] = new TH1F(Form("jtRecoOverRawVPt_Res_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Reco. Jet p_{T};#sigma_{Reco./Raw} (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins);

      jtRawOverGenVPt_Res_p[iter][centIter] = new TH1F(Form("jtRawOverGenVPt_Res_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet p_{T};#sigma_{Raw/Gen.} (%s)", jetAlgo[iter].c_str()), nJtPtBins, jtPtBins);

      jtRecoOverGenVEta_p[iter][centIter] = new TH2F(Form("jtRecoOverGenVEta_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet #eta; (Reco. %s Jet p_{T})/(Gen. Jet p_{T})", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins, 30, 0, 3);

      jtRecoOverRawVEta_p[iter][centIter] = new TH2F(Form("jtRecoOverRawVEta_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Reco. Jet #eta; (Reco %s Jet p_{T})/(Raw Jet p_{T})", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins, 30, 0, 3);

      jtRawOverGenVEta_p[iter][centIter] = new TH2F(Form("jtRawOverGenVEta_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet #eta; (Raw %s Jet p_{T})/(Gen Jet p_{T})", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins, 30, 0, 3);

      jtEffVEta_p[iter][centIter] = new TH2F(Form("jtEffVEta_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet #eta;Eff. (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins, 2, -0.5, 1.5);

      jtRecoOverGenVEta_Mean_p[iter][centIter] = new TH1F(Form("jtRecoOverGenVEta_Mean_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet #eta;#mu_{Reco./Gen.} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins);

      jtRecoOverRawVEta_Mean_p[iter][centIter] = new TH1F(Form("jtRecoOverRawVEta_Mean_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Reco. Jet #eta;#mu_{Reco./Raw} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins);

      jtRawOverGenVEta_Mean_p[iter][centIter] = new TH1F(Form("jtRawOverGenVEta_Mean_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet #eta;#mu_{Raw/Gen.} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins);

      jtEffVEta_Mean_p[iter][centIter] = new TH1F(Form("jtEffVEta_Mean_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet #eta;Eff. (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins);

      jtRecoOverGenVEta_Res_p[iter][centIter] = new TH1F(Form("jtRecoOverGenVEta_Res_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet #eta;#sigma_{Reco./Gen.} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins);

      jtRecoOverRawVEta_Res_p[iter][centIter] = new TH1F(Form("jtRecoOverRawVEta_Res_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Reco. Jet #eta;#sigma_{Reco./Raw} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins);

      jtRawOverGenVEta_Res_p[iter][centIter] = new TH1F(Form("jtRawOverGenVEta_Res_%s_cent%dto%d_h", jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";Gen. Jet #eta;#sigma_{Raw/Gen.} (%s)", jetAlgo[iter].c_str()), nJtEtaBins, jtEtaBins);

      for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	Int_t ptLowInt = std::trunc(jtPtBins[jtIter]);
	Int_t ptHiInt = std::trunc(jtPtBins[jtIter+1]);
	
	Int_t ptLowDec = std::trunc(jtPtBins[jtIter]*10 - ptLowInt*10);
	Int_t ptHiDec = std::trunc(jtPtBins[jtIter+1]*10 - ptHiInt*10);
	
	jtRecoOverGenVPt_MeanResPts_p[iter][centIter][jtIter] = new TH1F(Form("jtRecoOverGenVPt_MeanResPts_Pt%dp%dTo%dp%d_%s_cent%dto%d_h", ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";(Reco. %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo[iter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 30, 0, 3);

	jtRecoOverRawVPt_MeanResPts_p[iter][centIter][jtIter] = new TH1F(Form("jtRecoOverRawVPt_MeanResPts_Pt%dp%dTo%dp%d_%s_cent%dto%d_h", ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";(Reco. %s Jet p_{T})/(Raw. Jet p_{T});Events (%d.%d<p_{T,Reco.}<%d.%d)", jetAlgo[iter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 30, 0, 3);

	jtRawOverGenVPt_MeanResPts_p[iter][centIter][jtIter] = new TH1F(Form("jtRawOverGenVPt_MeanResPts_Pt%dp%dTo%dp%d_%s_cent%dto%d_h", ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";(Raw %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo[iter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 30, 0, 3);

	jtEffVPt_MeanResPts_p[iter][centIter][jtIter] = new TH1F(Form("jtEffVPt_MeanResPts_Pt%dp%dTo%dp%d_%s_cent%dto%d_h", ptLowInt, ptLowDec, ptHiInt, ptHiDec, jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";(Eff. (%s));Events (%d.%d<p_{T,Gen.}<%d.%d)", jetAlgo[iter].c_str(), ptLowInt, ptLowDec, ptHiInt, ptHiDec), 2, -0.5, 1.5);
      }


      for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	Int_t etaLowInt = std::trunc(jtEtaBins[jtIter]);
	Int_t etaHiInt = std::trunc(jtEtaBins[jtIter+1]);

	Int_t etaLowDec = std::trunc(jtEtaBins[jtIter]*10 - etaLowInt*10);
	Int_t etaHiDec = std::trunc(jtEtaBins[jtIter+1]*10 - etaHiInt*10);
	
	jtRecoOverGenVEta_MeanResPts_p[iter][centIter][jtIter] = new TH1F(Form("jtRecoOverGenVEta_MeanResPts_Eta%dp%dTo%dp%d_%s_cent%dto%d_h", etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";(Reco. %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<#eta_{Gen.}<%d.%d)", jetAlgo[iter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec), 30, 0, 3);

	jtRecoOverRawVEta_MeanResPts_p[iter][centIter][jtIter] = new TH1F(Form("jtRecoOverRawVEta_MeanResPts_Eta%dp%dTo%dp%d_%s_cent%dto%d_h", etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";(Reco. %s Jet p_{T})/(Raw. Jet p_{T});Events (%d.%d<#eta_{Reco.}<%d.%d)", jetAlgo[iter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec), 30, 0, 3);

	jtRawOverGenVEta_MeanResPts_p[iter][centIter][jtIter] = new TH1F(Form("jtRawOverGenVEta_MeanResPts_Eta%dp%dTo%dp%d_%s_cent%dto%d_h", etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";(Raw %s Jet p_{T})/(Gen. Jet p_{T});Events (%d.%d<#eta_{Gen.}<%d.%d)", jetAlgo[iter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec), 30, 0, 3);

	jtEffVEta_MeanResPts_p[iter][centIter][jtIter] = new TH1F(Form("jtEffVEta_MeanResPts_Eta%dp%dTo%dp%d_%s_cent%dto%d_h", etaLowInt, etaLowDec, etaHiInt, etaHiDec, jetAlgo[iter].c_str(), centBins[centIter+1]/2, centBins[centIter]/2), Form(";(Eff. (%s));Events (%d.%d<#eta_{Gen.}<%d.%d)", jetAlgo[iter].c_str(), etaLowInt, etaLowDec, etaHiInt, etaHiDec), 2, -0.5, 1.5);
      }
    }
  }

  std::cout << "Getting entries..." << std::endl;

  const Int_t nEntries = jetTree_p[0]->GetEntries();

  std::cout << "Number of Events: " << nEntries << std::endl;

  for(Int_t entry = 0; entry < nEntries; entry++){
    if(entry%10000 == 0) std::cout << "Entry: " << entry << std::endl;

    hiTree_p->GetEntry(entry);
    //    genTree_p->GetEntry(entry);

    Int_t centPos = -1;

    for(Int_t iter = 0; iter < nCentBins; iter++){
      if(hiBin_ < centBins[iter] && hiBin_ >= centBins[iter+1]){
	centPos = iter;
	break;
      }
    }

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

    for(Int_t algoIter = 0; algoIter < nJetAlgo; algoIter++){
      for(Int_t jtIter = 0; jtIter < nJt_[algoIter]; jtIter++){
	if(TMath::Abs(jtEta_[algoIter][jtIter]) > 2.0) continue;
       
	if(jtPt_[algoIter][jtIter] < 5.0) continue;

	//	if(refSubID_[algoIter][jtIter] != 0 && refSubID_[algoIter][jtIter] != -1) continue;

	if(maxElePt > 10)
	  if(getDR(jtEta_[algoIter][jtIter], jtPhi_[algoIter][jtIter], maxEleEta, maxElePhi) < 0.4) continue;

	if(twoElePt > 5)
	  if(getDR(jtEta_[algoIter][jtIter], jtPhi_[algoIter][jtIter], twoEleEta, twoElePhi) < 0.4) continue;

	jtRecoOverGenVPt_p[algoIter][centPos]->Fill(refPt_[algoIter][jtIter], jtPt_[algoIter][jtIter]/refPt_[algoIter][jtIter]);

	jtRecoOverRawVPt_p[algoIter][centPos]->Fill(jtPt_[algoIter][jtIter], jtPt_[algoIter][jtIter]/jtRawPt_[algoIter][jtIter]);

	jtRawOverGenVPt_p[algoIter][centPos]->Fill(refPt_[algoIter][jtIter], jtRawPt_[algoIter][jtIter]/refPt_[algoIter][jtIter]);
	
	for(Int_t jtIter2 = 0; jtIter2 < nJtPtBins; jtIter2++){
	  if(refPt_[algoIter][jtIter] > jtPtBins[jtIter2] && refPt_[algoIter][jtIter] < jtPtBins[jtIter2+1]){
	    jtRecoOverGenVPt_MeanResPts_p[algoIter][centPos][jtIter2]->Fill(jtPt_[algoIter][jtIter]/refPt_[algoIter][jtIter]);
	    jtRawOverGenVPt_MeanResPts_p[algoIter][centPos][jtIter2]->Fill(jtRawPt_[algoIter][jtIter]/refPt_[algoIter][jtIter]);	
	    break;
	  }
	}

	for(Int_t jtIter2 = 0; jtIter2 < nJtPtBins; jtIter2++){
	  if(jtPt_[algoIter][jtIter] > jtPtBins[jtIter2] && jtPt_[algoIter][jtIter] < jtPtBins[jtIter2+1]){
	    jtRecoOverRawVPt_MeanResPts_p[algoIter][centPos][jtIter2]->Fill(jtPt_[algoIter][jtIter]/jtRawPt_[algoIter][jtIter]);
	    break;
	  }
	}


	if(refPt_[algoIter][jtIter] > 30){
	  jtRecoOverGenVEta_p[algoIter][centPos]->Fill(refEta_[algoIter][jtIter], jtPt_[algoIter][jtIter]/refPt_[algoIter][jtIter]);
	  jtRawOverGenVEta_p[algoIter][centPos]->Fill(refEta_[algoIter][jtIter], jtRawPt_[algoIter][jtIter]/refPt_[algoIter][jtIter]);
	
	  for(Int_t jtIter2 = 0; jtIter2 < nJtEtaBins; jtIter2++){
	    if(refEta_[algoIter][jtIter] > jtEtaBins[jtIter2] && refEta_[algoIter][jtIter] < jtEtaBins[jtIter2+1]){
	      jtRecoOverGenVEta_MeanResPts_p[algoIter][centPos][jtIter2]->Fill(jtPt_[algoIter][jtIter]/refPt_[algoIter][jtIter]);
	      jtRawOverGenVEta_MeanResPts_p[algoIter][centPos][jtIter2]->Fill(jtRawPt_[algoIter][jtIter]/refPt_[algoIter][jtIter]);	
	      break;
	    }
	  }
	}

	if(jtPt_[algoIter][jtIter] > 30){
	  jtRecoOverRawVEta_p[algoIter][centPos]->Fill(jtEta_[algoIter][jtIter], jtPt_[algoIter][jtIter]/jtRawPt_[algoIter][jtIter]);

	  for(Int_t jtIter2 = 0; jtIter2 < nJtEtaBins; jtIter2++){
	    if(jtEta_[algoIter][jtIter] > jtEtaBins[jtIter2] && jtEta_[algoIter][jtIter] < jtEtaBins[jtIter2+1]){
	      jtRecoOverRawVEta_MeanResPts_p[algoIter][centPos][jtIter2]->Fill(jtPt_[algoIter][jtIter]/jtRawPt_[algoIter][jtIter]);
	      break;
	    }
	  }
	}
      }

      genSort(nGenJt_[algoIter], genJtPt_[algoIter], genJtPhi_[algoIter], genJtEta_[algoIter]);
      const Int_t boolSize = nJt_[algoIter];
      Bool_t isUsedRecoJet[boolSize];
      for(Int_t boolIter = 0; boolIter < boolSize; boolIter++){
	isUsedRecoJet[boolIter] = false;
      }

      //      if(nGen_[algoIter] != 0) std::cout << "NGENJT: " << nGen_[algoIter] << std::endl;
      for(Int_t jtIter = 0; jtIter < nGenJt_[algoIter]; jtIter++){
	if(TMath::Abs(genJtEta_[algoIter][jtIter]) > 2.0) continue;

	if(genJtPt_[algoIter][jtIter] < 5.0) break;

        if(maxElePt > 10)
          if(getDR(genJtEta_[algoIter][jtIter], genJtPhi_[algoIter][jtIter], maxEleEta, maxElePhi) < 0.4) continue;

        if(twoElePt > 5)
          if(getDR(genJtEta_[algoIter][jtIter], genJtPhi_[algoIter][jtIter], twoEleEta, twoElePhi) < 0.4) continue;

	Int_t fillVal = 0;
	Float_t tempDR = 999;
	Int_t tempPos = -1;

	for(Int_t jtIter2 = 0; jtIter2 < nJt_[algoIter]; jtIter2++){
	  if(isUsedRecoJet[jtIter2]) continue;

	  if(getDR(genJtEta_[algoIter][jtIter], genJtPhi_[algoIter][jtIter], jtEta_[algoIter][jtIter2], jtPhi_[algoIter][jtIter2]) < tempDR){
	    tempPos = jtIter2;
	    tempDR = getDR(genJtEta_[algoIter][jtIter], genJtPhi_[algoIter][jtIter], jtEta_[algoIter][jtIter2], jtPhi_[algoIter][jtIter2]);
	  }
	}

	if(tempDR < 0.4){
	  fillVal = 1;
	  isUsedRecoJet[tempPos] = true;
	}
       
	//	if(genJtMatchIndex_[algoIter][jtIter] >= 0) fillVal = 1;

	jtEffVPt_p[algoIter][centPos]->Fill(genJtPt_[algoIter][jtIter], fillVal);

	jtEffVPtEta_Denom_p[algoIter][centPos]->Fill(genJtEta_[algoIter][jtIter], genJtPt_[algoIter][jtIter]);
	if(fillVal) jtEffVPtEta_p[algoIter][centPos]->Fill(genJtEta_[algoIter][jtIter], genJtPt_[algoIter][jtIter]);
	
	for(Int_t jtIter2 = 0; jtIter2 < nJtPtBins; jtIter2++){
	  if(genJtPt_[algoIter][jtIter] > jtPtBins[jtIter2] && genJtPt_[algoIter][jtIter] < jtPtBins[jtIter2+1]){
	    jtEffVPt_MeanResPts_p[algoIter][centPos][jtIter2]->Fill(fillVal);
	    break;
	  }
	}

	if(genJtPt_[algoIter][jtIter] > 30){
	  jtEffVEta_p[algoIter][centPos]->Fill(genJtEta_[algoIter][jtIter], fillVal);

	  for(Int_t jtIter2 = 0; jtIter2 < nJtEtaBins; jtIter2++){
	    if(genJtEta_[algoIter][jtIter] > jtEtaBins[jtIter2] && genJtEta_[algoIter][jtIter] < jtEtaBins[jtIter2+1]){
	      jtEffVEta_MeanResPts_p[algoIter][centPos][jtIter2]->Fill(fillVal);
	      break;
	    }
	  }
	}
      }
    }
  }

  std::cout << "#Events with no lead electron: " << nNoLeadEle << std::endl;
  std::cout << "#Events with no sublead electron: " << nNoSubleadEle << std::endl;

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	Float_t tempMean = -999;
	Float_t tempMeanErr = -999;
	Float_t tempRes = -999;
	Float_t tempResErr = -999;

	FitGauss(jtRecoOverGenVPt_MeanResPts_p[iter][centIter][jtIter], tempMean, tempMeanErr, tempRes, tempResErr);
	
	jtRecoOverGenVPt_Mean_p[iter][centIter]->SetBinContent(jtIter+1, tempMean);
	jtRecoOverGenVPt_Mean_p[iter][centIter]->SetBinError(jtIter+1, tempMeanErr);
	jtRecoOverGenVPt_Res_p[iter][centIter]->SetBinContent(jtIter+1, tempRes);
	jtRecoOverGenVPt_Res_p[iter][centIter]->SetBinError(jtIter+1, tempResErr);

	FitGauss(jtRecoOverRawVPt_MeanResPts_p[iter][centIter][jtIter], tempMean, tempMeanErr, tempRes, tempResErr);

	jtRecoOverRawVPt_Mean_p[iter][centIter]->SetBinContent(jtIter+1, tempMean);
	jtRecoOverRawVPt_Mean_p[iter][centIter]->SetBinError(jtIter+1, tempMeanErr);
	jtRecoOverRawVPt_Res_p[iter][centIter]->SetBinContent(jtIter+1, tempRes);
	jtRecoOverRawVPt_Res_p[iter][centIter]->SetBinError(jtIter+1, tempResErr);

	FitGauss(jtRawOverGenVPt_MeanResPts_p[iter][centIter][jtIter], tempMean, tempMeanErr, tempRes, tempResErr);

	jtRawOverGenVPt_Mean_p[iter][centIter]->SetBinContent(jtIter+1, tempMean);
	jtRawOverGenVPt_Mean_p[iter][centIter]->SetBinError(jtIter+1, tempMeanErr);
	jtRawOverGenVPt_Res_p[iter][centIter]->SetBinContent(jtIter+1, tempRes);
	jtRawOverGenVPt_Res_p[iter][centIter]->SetBinError(jtIter+1, tempResErr);

	jtEffVPt_Mean_p[iter][centIter]->SetBinContent(jtIter+1, jtEffVPt_MeanResPts_p[iter][centIter][jtIter]->GetMean());
	jtEffVPt_Mean_p[iter][centIter]->SetBinError(jtIter+1, jtEffVPt_MeanResPts_p[iter][centIter][jtIter]->GetMeanError());
      }

      for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	Float_t tempMean = -999;
	Float_t tempMeanErr = -999;
	Float_t tempRes = -999;
	Float_t tempResErr = -999;

        FitGauss(jtRecoOverGenVEta_MeanResPts_p[iter][centIter][jtIter], tempMean, tempMeanErr, tempRes, tempResErr);

	jtRecoOverGenVEta_Mean_p[iter][centIter]->SetBinContent(jtIter+1, tempMean);
	jtRecoOverGenVEta_Mean_p[iter][centIter]->SetBinError(jtIter+1, tempMeanErr);
	jtRecoOverGenVEta_Res_p[iter][centIter]->SetBinContent(jtIter+1, tempRes);
	jtRecoOverGenVEta_Res_p[iter][centIter]->SetBinError(jtIter+1, tempResErr);

        FitGauss(jtRecoOverRawVEta_MeanResPts_p[iter][centIter][jtIter], tempMean, tempMeanErr, tempRes, tempResErr);

	jtRecoOverRawVEta_Mean_p[iter][centIter]->SetBinContent(jtIter+1, tempMean);
	jtRecoOverRawVEta_Mean_p[iter][centIter]->SetBinError(jtIter+1, tempMeanErr);
	jtRecoOverRawVEta_Res_p[iter][centIter]->SetBinContent(jtIter+1, tempRes);
	jtRecoOverRawVEta_Res_p[iter][centIter]->SetBinError(jtIter+1, tempResErr);

        FitGauss(jtRawOverGenVEta_MeanResPts_p[iter][centIter][jtIter], tempMean, tempMeanErr, tempRes, tempResErr);

	jtRawOverGenVEta_Mean_p[iter][centIter]->SetBinContent(jtIter+1, tempMean);
	jtRawOverGenVEta_Mean_p[iter][centIter]->SetBinError(jtIter+1, tempMeanErr);
	jtRawOverGenVEta_Res_p[iter][centIter]->SetBinContent(jtIter+1, tempRes);
	jtRawOverGenVEta_Res_p[iter][centIter]->SetBinError(jtIter+1, tempResErr);

	jtEffVEta_Mean_p[iter][centIter]->SetBinContent(jtIter+1, jtEffVEta_MeanResPts_p[iter][centIter][jtIter]->GetMean());
	jtEffVEta_Mean_p[iter][centIter]->SetBinError(jtIter+1, jtEffVEta_MeanResPts_p[iter][centIter][jtIter]->GetMeanError());
      }
    }
  }

  std::string outName = inFileName;
  const std::string inString = ".root";
  const std::string inString2 = "*";
  const std::string outString = "_HIST.root";
  std::size_t strIndex = 0;

  while(outName.find("/") != std::string::npos){
    outName.replace(0, outName.find("/")+1, "");
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

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    TDirectory* dir_p = outFile_p->GetDirectory(Form("%s", jetAlgo[iter].c_str()));
    if(dir_p){
      dir_p->cd();
    }
    else{
      dir_p = outFile_p->mkdir(Form("%s", jetAlgo[iter].c_str()));
      dir_p->cd();
    }

    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      jtRecoOverGenVPt_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRecoOverGenVPt_Mean_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRecoOverGenVPt_Res_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRecoOverRawVPt_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRecoOverRawVPt_Mean_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRecoOverRawVPt_Res_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRawOverGenVPt_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRawOverGenVPt_Mean_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRawOverGenVPt_Res_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtEffVPt_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtEffVPt_Mean_p[iter][centIter]->Write("", TObject::kOverwrite);
      

      //      jtEffVPtEta_p[iter][centIter]->Write(Form("%s_NUM",jtEffVPtEta_p[iter][centIter]->GetName()), TObject::kOverwrite);
      //      jtEffVPtEta_Denom_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtEffVPtEta_p[iter][centIter]->Divide(jtEffVPtEta_Denom_p[iter][centIter]);
      jtEffVPtEta_p[iter][centIter]->Write("", TObject::kOverwrite);

      for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	jtRecoOverGenVPt_MeanResPts_p[iter][centIter][jtIter]->Write("", TObject::kOverwrite);
	jtRecoOverRawVPt_MeanResPts_p[iter][centIter][jtIter]->Write("", TObject::kOverwrite);
	jtRawOverGenVPt_MeanResPts_p[iter][centIter][jtIter]->Write("", TObject::kOverwrite);
	jtEffVPt_MeanResPts_p[iter][centIter][jtIter]->Write("", TObject::kOverwrite);
      }

      jtRecoOverGenVEta_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRecoOverGenVEta_Mean_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRecoOverGenVEta_Res_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRecoOverRawVEta_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRecoOverRawVEta_Mean_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRecoOverRawVEta_Res_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRawOverGenVEta_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRawOverGenVEta_Mean_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtRawOverGenVEta_Res_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtEffVEta_p[iter][centIter]->Write("", TObject::kOverwrite);
      jtEffVEta_Mean_p[iter][centIter]->Write("", TObject::kOverwrite);
      
      for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	jtRecoOverGenVEta_MeanResPts_p[iter][centIter][jtIter]->Write("", TObject::kOverwrite);
	jtRecoOverRawVEta_MeanResPts_p[iter][centIter][jtIter]->Write("", TObject::kOverwrite);
	jtRawOverGenVEta_MeanResPts_p[iter][centIter][jtIter]->Write("", TObject::kOverwrite);
	jtEffVEta_MeanResPts_p[iter][centIter][jtIter]->Write("", TObject::kOverwrite);
      }
    }
  }

  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      delete jtRecoOverGenVPt_p[iter][centIter];
      delete jtRecoOverGenVPt_Mean_p[iter][centIter];
      delete jtRecoOverGenVPt_Res_p[iter][centIter];
      delete jtRecoOverRawVPt_p[iter][centIter];
      delete jtRecoOverRawVPt_Mean_p[iter][centIter];
      delete jtRecoOverRawVPt_Res_p[iter][centIter];
      delete jtRawOverGenVPt_p[iter][centIter];
      delete jtRawOverGenVPt_Mean_p[iter][centIter];
      delete jtRawOverGenVPt_Res_p[iter][centIter];
      delete jtEffVPt_p[iter][centIter];
      delete jtEffVPt_Mean_p[iter][centIter];

      delete jtEffVPtEta_p[iter][centIter];
      delete jtEffVPtEta_Denom_p[iter][centIter];
      
      for(Int_t jtIter = 0; jtIter < nJtPtBins; jtIter++){
	delete jtRecoOverGenVPt_MeanResPts_p[iter][centIter][jtIter];
	delete jtRecoOverRawVPt_MeanResPts_p[iter][centIter][jtIter];
	delete jtRawOverGenVPt_MeanResPts_p[iter][centIter][jtIter];
	delete jtEffVPt_MeanResPts_p[iter][centIter][jtIter];
      }

      delete jtRecoOverGenVEta_p[iter][centIter];
      delete jtRecoOverGenVEta_Mean_p[iter][centIter];
      delete jtRecoOverGenVEta_Res_p[iter][centIter];
      delete jtRecoOverRawVEta_p[iter][centIter];
      delete jtRecoOverRawVEta_Mean_p[iter][centIter];
      delete jtRecoOverRawVEta_Res_p[iter][centIter];
      delete jtRawOverGenVEta_p[iter][centIter];
      delete jtRawOverGenVEta_Mean_p[iter][centIter];
      delete jtRawOverGenVEta_Res_p[iter][centIter];
      delete jtEffVEta_p[iter][centIter];
      delete jtEffVEta_Mean_p[iter][centIter];

      for(Int_t jtIter = 0; jtIter < nJtEtaBins; jtIter++){
	delete jtRecoOverGenVEta_MeanResPts_p[iter][centIter][jtIter];
	delete jtRecoOverRawVEta_MeanResPts_p[iter][centIter][jtIter];
	delete jtRawOverGenVEta_MeanResPts_p[iter][centIter][jtIter];
	delete jtEffVEta_MeanResPts_p[iter][centIter][jtIter];
      }
    }
  }
    
  //  inFile_p->Close();
  //  delete inFile_p;

  return;
}
