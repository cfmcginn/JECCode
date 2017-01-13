#include "TFile.h"
#include "TDirectoryFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TKey.h"
#include "TLine.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TStyle.h"
#include "TF1.h"

#include <string>
#include <iostream>
#include <vector>

#include "include/doGlobalDebug.h"
#include "include/checkMakeDir.h"
#include "include/jecConfigParser.h"

#include "SuperCanvas/include/superCanvas.h"

Bool_t debugMode = false;

Bool_t plotTrue = false;

Bool_t addResCorrStr = true;

const Int_t nHistName = 3;
const std::string inHistName[nHistName] = {"RecoOverGen", "Eff", "DPhi"};
const std::string xAxisLabel[nHistName] = {"Reco./Gen.", "Eff.", "<#Delta#phi>"};

const Int_t nMeanRes = 3;
const std::string meanResStr[nMeanRes] = {"Mean", "Res", "ResOverMean"};
const std::string meanResStr2[nMeanRes] = {"#mu", "#sigma", "#sigma/#mu"};

const Int_t nQG = 4;
const std::string qgStr[nQG] = {"Inc", "Q", "G", "Untagged"};
const std::string qgStr2[nQG] = {"Inc.", "Quarks", "Gluons", "Untagged"};
const Int_t qgCol[nQG] = {kGray+1, kBlue, kRed, kYellow+2};

const std::string meanResPtsStr = "MeanResPts";

void claverCanvasSaving(TCanvas* c, TString s,TString format="gif"){
  TDatime* date = new TDatime();
  c->SaveAs(Form("%s_%d.%s",s.Data(),date->GetDate(), format.Data()));
  return;
}

Float_t setMaxMinNice(Float_t inMaxMin, Bool_t isMax){
  for(Int_t iter = 0; iter < 5; iter++){
    if(inMaxMin*TMath::Power(10, iter+1) > 10){
      inMaxMin = inMaxMin*TMath::Power(10, iter+1);
      if(isMax) inMaxMin = std::ceil(inMaxMin);
      else inMaxMin = std::floor(inMaxMin);
      inMaxMin = inMaxMin/TMath::Power(10, iter+1);
      if(isMax) inMaxMin -= .01/TMath::Power(10, iter+1);
      else inMaxMin += .01/TMath::Power(10, iter+1);

      break;
    }
  }
  
  return inMaxMin;
}


int makeJECPlotBasicQGDistrib(const std::string inFileName, jecConfigParser config, std::string distribStartString, const bool doBinWidthNorm, const bool doIncNorm, std::vector<std::string>* pdfList_p)
{
  const Bool_t isPbPb = config.GetIsPbPb();
  const Int_t nCentBins = config.GetNCentBins();
  std::vector<unsigned int> centBins = config.GetCentBins();

  if(isPbPb){
    if((unsigned int)(nCentBins+1) != centBins.size()){
      std::cout << "Config nCentBins \'" << nCentBins << "\' != centBins.size \'" << centBins.size() << "\'. return 1" << std::endl;
      return 1;
    }
  }
  
  std::string centStrings[nCentBins];
  std::string centStrings2[nCentBins];
  
  for(Int_t iter = 0; iter < nCentBins; iter++){
    if(isPbPb){
      centStrings[iter] = "Cent" + std::to_string(centBins.at(iter)) + "to" + std::to_string(centBins.at(iter+1));
      centStrings2[iter] = std::to_string(centBins.at(iter)) + "-" + std::to_string(centBins.at(iter+1)) + "%";
    }
    else{
      centStrings[iter] = "PP";
      centStrings2[iter] = "PP";
    }
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  const Int_t nContents = inFile_p->GetListOfKeys()->GetEntries();

  Int_t nTH1Temp = 0;
  std::vector<std::string>* th1Names_p = new std::vector<std::string>;

  Int_t nDirTemp = 0;
  std::vector<std::string>* dirNames_p = new std::vector<std::string>;

  for(Int_t iter = 0; iter < nContents; iter++){
    TString className = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetClassName();
    TString name = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetName();

    if(className.Index("TDirectoryFile") >= 0 && name.Index("ak") >= 0){
      nDirTemp++;
      dirNames_p->push_back(name.Data());

      TDirectoryFile* tempDir_p = (TDirectoryFile*)inFile_p->Get(name);

      const Int_t nDirContents = tempDir_p->GetListOfKeys()->GetEntries();
      
      for(Int_t dirIter = 0; dirIter < nDirContents; dirIter++){
	TString name2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetName();
	TString className2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetClassName();

	if(className2.Index("TH1") >= 0){
	  if(name2.Index(distribStartString.c_str()) != 0) continue;

	  if(name2.Index("jtPtVPt") >= 0) continue;

	  nTH1Temp++;
	  th1Names_p->push_back(Form("%s/%s", name.Data(), name2.Data()));
	}
      }
    }
  }
  
  const Int_t nTH1 = nTH1Temp;
  const Int_t nDir = nDirTemp;
  TH1F* th1_p[nTH1];  
  
  if(debugMode) std::cout << __LINE__ << std::endl;

  std::vector<superCanvas*> th1Canv_p(nDir);
  //  TCanvas* th1Canv_p[nDir];

  Int_t nXPanel = nCentBins;
  Int_t nYPanel = 1;
  if(!isPbPb){
    nXPanel = 1;
    nYPanel = 1;
  }

  Int_t th1CanvTot[nDir][nXPanel][nYPanel];
  Int_t th1CanvCount[nDir][nXPanel][nYPanel];
  
  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t xIter = 0; xIter < nXPanel; xIter++){
      for(Int_t yIter = 0; yIter < nYPanel; yIter++){
        th1CanvTot[iter][xIter][yIter] = 0;
        th1CanvCount[iter][xIter][yIter] = 0;
      }
    }
  }


  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    if(doIncNorm && th1Name.find("_Inc_") != std::string::npos) continue;

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
      if(th1Name.substr(0, th1Name.find("/")).size() != dirNames_p->at(dirIter).size()) continue;
      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    Int_t centPos = -1;
    if(isPbPb){
      for(Int_t centIter = 0; centIter < nCentBins; centIter++){
	std::size_t pos = th1Name.find(centStrings[centIter]);
        if(pos != std::string::npos){
          centPos = centIter;
          break;
        }
      }
    }
    else centPos = 0;

    th1CanvTot[dirPos][centPos][0]++;
  }  

  if(debugMode) std::cout << __LINE__ << std::endl;

  for(Int_t iter = 0; iter < nDir; iter++){
    th1Canv_p[iter] = new superCanvas(nXPanel, nYPanel, th1CanvTot[iter][0][0], 1000, dirNames_p->at(iter));

    th1Canv_p[iter]->SetGlobalMaxMin();
    
    th1Canv_p[iter]->canv_p->cd();
    th1Canv_p[iter]->canv_p->Draw();

    for(Int_t xIter = 0; xIter < nXPanel; xIter++){
      for(Int_t yIter = 0; yIter < nYPanel; yIter++){
	th1Canv_p[iter]->canv_p->cd();
	th1Canv_p[iter]->pads_p[xIter][yIter]->Draw();
	th1Canv_p[iter]->pads_p[xIter][yIter]->cd();
      }
    }
  }

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);
    th1_p[iter] = (TH1F*)inFile_p->Get(th1Name.c_str());
  }

  if(doIncNorm){
    for(Int_t iter = 0; iter < nTH1; iter++){
      std::string th1Name = th1Names_p->at(iter);
      if(th1Name.find("_Inc_") == std::string::npos) continue;
      
      Int_t dirPos = -1;
      for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
	std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
	
	if(debugMode) std::cout << __LINE__ << std::endl;
	
	if(th1Name.substr(0, th1Name.find("/")).size() != dirNames_p->at(dirIter).size()) continue;
	
	if(pos != std::string::npos){
	  dirPos = dirIter;
	  break;
	}
      }
      
      Int_t centPos = -1;
      if(isPbPb){
	for(Int_t centIter = 0; centIter < nCentBins; centIter++){
	  std::size_t pos = th1Name.find(centStrings[centIter]);
	  if(pos != std::string::npos){
	    centPos = centIter;
	    break;
	  }
	}
      }
      else centPos = 0;
      
      for(Int_t iter2 = 0; iter2 < nTH1; iter2++){
	if(iter == iter2) continue;
	
	Int_t dirPos2 = -1;
	for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
	  std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
	  
	  if(debugMode) std::cout << __LINE__ << std::endl;
	  
	  if(th1Name.substr(0, th1Name.find("/")).size() != dirNames_p->at(dirIter).size()) continue;
	  
	  if(pos != std::string::npos){
	    dirPos2 = dirIter;
	    break;
	  }
	}
	
	if(dirPos != dirPos2) continue;
	
	Int_t centPos2 = -1;
	if(isPbPb){
	  for(Int_t centIter = 0; centIter < nCentBins; centIter++){
	    std::size_t pos = th1Name.find(centStrings[centIter]);
	    if(pos != std::string::npos){
	      centPos2 = centIter;
	      break;
	    }
	  }
	}
	else centPos2 = 0;
	
	if(centPos != centPos2) continue;
	
	th1_p[iter2]->Divide(th1_p[iter]);
      }
    }
  }

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    if(doIncNorm && th1Name.find("_Inc_") != std::string::npos) continue;

    if(doBinWidthNorm){
      std::string newTitle = th1_p[iter]->GetYaxis()->GetTitle();
      newTitle = "#frac{1.}{Bin Width} " + newTitle;
      th1_p[iter]->GetYaxis()->SetTitle(newTitle.c_str());

      for(Int_t binIter = 0; binIter < th1_p[iter]->GetNbinsX(); binIter++){
	Double_t binWidth = th1_p[iter]->GetBinWidth(binIter+1);

	th1_p[iter]->SetBinContent(binIter+1, th1_p[iter]->GetBinContent(binIter+1)/binWidth);
	th1_p[iter]->SetBinError(binIter+1, th1_p[iter]->GetBinError(binIter+1)/binWidth);
      }
    }
    else if(doIncNorm){
      th1_p[iter]->GetYaxis()->SetTitle("Ratio #frac{Flavor}{Inclusive}");
    }

    th1_p[iter]->SetMarkerSize(0.5);

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
      
      if(debugMode) std::cout << __LINE__ << std::endl;

      if(th1Name.substr(0, th1Name.find("/")).size() != dirNames_p->at(dirIter).size()) continue;

      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    th1Canv_p[dirPos]->canv_p->cd();
  
    Int_t centPos = -1;
    if(isPbPb){
      for(Int_t centIter = 0; centIter < nCentBins; centIter++){
	std::size_t pos = th1Name.find(centStrings[centIter]);
	if(pos != std::string::npos){
	  centPos = centIter;
	  break;
	}
      }
    }
    else centPos = 0;

    std::string legStr = "Inc.";
    if(th1Name.find("_Q_") != std::string::npos) legStr = "Quark";
    else if(th1Name.find("_G_") != std::string::npos) legStr = "Gluon";
    else if(th1Name.find("_Untagged_") != std::string::npos) legStr = "Untagged";

    th1Canv_p[dirPos]->SetHist(th1_p[iter], centPos, 0, th1CanvCount[dirPos][centPos][0], legStr);
    th1CanvCount[dirPos][centPos][0]++;
  }


  for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
    th1Canv_p[dirIter]->canv_p->cd();
    th1Canv_p[dirIter]->MakeHistMaxMinNice();
    th1Canv_p[dirIter]->SetHistMaxMin();

    gStyle->SetOptStat(0);

    th1Canv_p[dirIter]->SetPanelWhiteSpace();

    if(debugMode) std::cout << __LINE__ << ", " << dirIter << std::endl;
	  
    for(Int_t iter = 0; iter < nXPanel; iter++){
      for(Int_t iter2 = 0; iter2 < nYPanel; iter2++){
	if(debugMode) std::cout << __LINE__ << ", " << dirIter << ", " << iter << ", " << iter2 << std::endl;

	th1Canv_p[dirIter]->canv_p->cd();
	th1Canv_p[dirIter]->pads_p[iter][iter2]->cd();

	bool isDrawn = false;
	
	for(Int_t histIter = 0; histIter < th1CanvCount[dirIter][iter][iter2]; histIter++){
	  if(debugMode) std::cout << __LINE__ << ", " << dirIter << ", " << iter << ", " << iter2 << ", " << histIter << std::endl;

	  th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->SetMarkerSize(2);
	  th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->SetMarkerStyle(20);	  

	  std::string th1Name = th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->GetName();
	  if(doIncNorm && th1Name.find("_Inc_") != std::string::npos) continue;

	  Int_t col = kGray+1;
	  if(th1Name.find("_Q_") != std::string::npos) col = kBlue;
	  else if(th1Name.find("_G_") != std::string::npos) col = kRed;
	  else if(th1Name.find("_Untagged_") != std::string::npos) col = kYellow+2;

	  th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->SetMarkerColor(col);	  
	  th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->SetFillColor(col);	  
	  th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->SetLineColor(1);

	 

	  if(!isDrawn){
	    isDrawn = true;
	    th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->GetXaxis()->SetRange(1, th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->GetNbinsX());
	    th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->DrawCopy("E1 P");
	  }
	  else th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->DrawCopy("SAME E1 P");
	}

	if(debugMode) std::cout << __LINE__ << std::endl;

	if(iter == 0) th1Canv_p[dirIter]->DrawLabel1(iter, 0, dirNames_p->at(dirIter));
	if(isPbPb) th1Canv_p[dirIter]->DrawLabel2(iter, 0, centStrings2[iter]);
	else th1Canv_p[dirIter]->DrawLabel2(iter, 0, "PP");
      }
    }
  }

  std::string outName = inFileName;
  const std::string inString = "_HIST.root";
  const std::string outString = "_PLOT.root";
  std::size_t strIndex = 0;

  strIndex = outName.find(inString);
  if(!(strIndex == std::string::npos)){
    outName.replace(strIndex, inString.length(), outString);
  }

  std::string binWidthNorm = "";
  if(doBinWidthNorm) binWidthNorm = "_BinWidthNorm";

  std::string incNorm = "";
  if(doIncNorm) incNorm = "_IncNorm";

  TFile* outFile_p = new TFile(outName.c_str(), "UPDATE");

  std::string tempFileName = config.GetConfigFileNameNoExt().c_str();
  while(tempFileName.find("/") != std::string::npos){
    tempFileName.replace(0, tempFileName.find("/")+1, "");
  }

  for(Int_t iter = 0; iter < nDir; iter++){

    std::string pdfName =  Form("pdfDir/%s_%sQGDistrib%s%s_%s", dirNames_p->at(iter).c_str(), distribStartString.c_str(), binWidthNorm.c_str(), incNorm.c_str(), tempFileName.c_str());
    
    gStyle->SetOptStat(0);
   
    claverCanvasSaving(th1Canv_p[iter]->canv_p, pdfName.c_str(), "pdf");
    TDatime* date = new TDatime();
    pdfName = pdfName + "_" + std::to_string(date->GetDate()) + ".pdf";
    delete date;
    
    pdfList_p->push_back(pdfName);
  }

  if(debugMode) std::cout << __LINE__ << std::endl;

  outFile_p->Close();
  delete outFile_p;

  for(Int_t iter = 0; iter < nDir; iter++){
    delete th1Canv_p[iter]->canv_p;
    delete th1Canv_p[iter];
  }

  th1Names_p->clear();
  delete th1Names_p;

  dirNames_p->clear();
  delete dirNames_p;

  inFile_p->Close();
  delete inFile_p;

  return 0;
}


int makeJECPlotMeanRes(const std::string inFileName, jecConfigParser config, const std::string fitMeanStr, const Int_t inHistNum, const Int_t meanResNum, const std::string ptEtaStr, std::vector<std::string>* pdfList_p)
{
  std::string ptEtaStr2 = "Pt";
  std::string ptEtaStr3 = "p_{T}";
  if(ptEtaStr.find("VEta") != std::string::npos){
    ptEtaStr2 = "Eta";
    ptEtaStr3 = "#eta";
  }

  const Bool_t isPbPb = config.GetIsPbPb();

  const Int_t nCentBins = config.GetNCentBins();
  std::vector<unsigned int> centBins = config.GetCentBins();

  if(isPbPb){
    if((unsigned int)(nCentBins+1) != centBins.size()){
      std::cout << "Config nCentBins \'" << nCentBins << "\' != centBins.size \'" << centBins.size() << "\'. return 1" << std::endl;
      return 1;
    }
  }

  std::string centStrings[nCentBins];
  std::string centStrings2[nCentBins];

  for(Int_t iter = 0; iter < nCentBins; iter++){
    if(isPbPb){
      centStrings[iter] = "Cent" + std::to_string(centBins.at(iter)) + "to" + std::to_string(centBins.at(iter+1));
      centStrings2[iter] = std::to_string(centBins.at(iter)) + "-" + std::to_string(centBins.at(iter+1)) + "%";
    }
    else{
      centStrings[iter] = "PP";
      centStrings2[iter] = "PP";
    }
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  const Int_t nContents = inFile_p->GetListOfKeys()->GetEntries();

  Int_t nTH1Temp = 0;
  std::vector<std::string>* th1Names_p = new std::vector<std::string>;

  Int_t nDirTemp = 0;
  std::vector<std::string>* dirNames_p = new std::vector<std::string>;

  Bool_t isMean = false;
  if(meanResStr[meanResNum].size() == 4 && meanResStr[meanResNum].find("Mean") != std::string::npos) isMean = true;

  Bool_t isRes = false;
  if(meanResStr[meanResNum].size() == 3 && meanResStr[meanResNum].find("Res") != std::string::npos) isRes = true;

  Bool_t isResOverMean = false;
  if(meanResStr[meanResNum].size() == 11 && meanResStr[meanResNum].find("ResOverMean") != std::string::npos) isResOverMean = true;

  for(Int_t iter = 0; iter < nContents; iter++){
    TString className = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetClassName();
    TString name = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetName();

    if(className.Index("TDirectoryFile") >= 0 && name.Index("ak") >= 0){
      nDirTemp++;
      dirNames_p->push_back(name.Data());

      TDirectoryFile* tempDir_p = (TDirectoryFile*)inFile_p->Get(name);

      const Int_t nDirContents = tempDir_p->GetListOfKeys()->GetEntries();
      
      for(Int_t dirIter = 0; dirIter < nDirContents; dirIter++){
	TString name2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetName();
	TString className2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetClassName();

	if(className2.Index("TH1") >= 0){
	  if(name2.Index(Form("%s_", meanResStr[meanResNum].c_str())) >= 0){
	
	    if((isMean || isRes) && name2.Index("ResOverMean_") >= 0) continue;

	    if(name2.Index(ptEtaStr.c_str()) < 0) continue;

	    if(isMean && name2.Index("Mean") >= 0 && name2.Index("_Q_") >= 0) continue;
	    if(isMean && name2.Index("Mean") >= 0 && name2.Index("_G_") >= 0) continue;
	    if(isMean && name2.Index("Mean") >= 0 && name2.Index("_Untagged_") >= 0) continue;

	    if(isRes && name2.Index("Res") >= 0 && name2.Index("_Q_") >= 0) continue;
	    if(isRes && name2.Index("Res") >= 0 && name2.Index("_G_") >= 0) continue;
	    if(isRes && name2.Index("Res") >= 0 && name2.Index("_Untagged_") >= 0) continue;

	    if(isResOverMean && name2.Index("ResOverMean") >= 0 && name2.Index("_Q_") >= 0) continue;
	    if(isResOverMean && name2.Index("ResOverMean") >= 0 && name2.Index("_G_") >= 0) continue;
	    if(isResOverMean && name2.Index("ResOverMean") >= 0 && name2.Index("_Untagged_") >= 0) continue;

	    //	    if(isMean && name2.Index("Mean") >= 0 && name2.Index("_G_") >= 0) continue;

	    if(!strcmp("Fake", inHistName[inHistNum].c_str()) && (name2.Index("_Q_") >= 0 || name2.Index("_G_") >= 0)) continue;
	    if(!strcmp("Eff", inHistName[inHistNum].c_str()) && (name2.Index("_Q_") >= 0 || name2.Index("_G_") >= 0)) continue;

	    if(fitMeanStr.size() != 0){
	      if(name2.Index("RecoOverGen") >= 0 && name2.Index("_Mean_") >= 0) continue;
	      if(name2.Index("RecoOverGen") >= 0 && name2.Index("_Res_") >= 0) continue;
	      if(name2.Index("RecoOverGen") >= 0 && name2.Index("_ResOverMean_") >= 0) continue;
	    }
	    else{
	      if(name2.Index("RecoOverGen") >= 0 && name2.Index("_FitMean_") >= 0) continue;
	      if(name2.Index("RecoOverGen") >= 0 && name2.Index("_FitRes_") >= 0) continue;
	      if(name2.Index("RecoOverGen") >= 0 && name2.Index("_FitResOverMean_") >= 0) continue;
	    }

	    if(name2.Index("RecoPtCut5") >= 0) continue;
	    if(name2.Index("RecoPtCut30") >= 0) continue;
	    //	    if(name2.Index("RecoPtCut10") >= 0) continue;

	    if(inHistName[inHistNum].find("DPhi") != std::string::npos && name2.Index("_Q_") >= 0) continue;
	    if(inHistName[inHistNum].find("DPhi") != std::string::npos && name2.Index("_G_") >= 0) continue;

	    if(name2.Index(inHistName[inHistNum].c_str()) >= 0){
	      nTH1Temp++;
	      th1Names_p->push_back(Form("%s/%s", name.Data(), name2.Data()));
	    }
	  }
	}
      }
    }
  }

  const Int_t nTH1 = nTH1Temp;
  const Int_t nDir = nDirTemp;
  TH1F* th1_p[nTH1];  
  
  if(debugMode) std::cout << __LINE__ << std::endl;

  std::vector<superCanvas*> th1Canv_p(nDir);
  //  TCanvas* th1Canv_p[nDir];

  Int_t nXPanel = nCentBins;
  Int_t nYPanel = 1;
  if(!isPbPb){
    nXPanel = 1;
    nYPanel = 1;
  }

  Int_t th1CanvTot[nDir][nXPanel][nYPanel];
  Int_t th1CanvCount[nDir][nXPanel][nYPanel];
  
  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t xIter = 0; xIter < nXPanel; xIter++){
      for(Int_t yIter = 0; yIter < nYPanel; yIter++){
        th1CanvTot[iter][xIter][yIter] = 0;
        th1CanvCount[iter][xIter][yIter] = 0;
      }
    }
  }

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
      if(th1Name.substr(0, th1Name.find("/")).size() != dirNames_p->at(dirIter).size()) continue;
      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    Int_t centPos = -1;
    if(isPbPb){
      for(Int_t centIter = 0; centIter < nCentBins; centIter++){
	std::size_t pos = th1Name.find(centStrings[centIter]);
        if(pos != std::string::npos){
          centPos = centIter;
          break;
        }
      }
    }
    else centPos = 0;

    th1CanvTot[dirPos][centPos][0]++;
  }  

  if(debugMode) std::cout << __LINE__ << std::endl;

  for(Int_t iter = 0; iter < nDir; iter++){
    th1Canv_p[iter] = new superCanvas(nXPanel, nYPanel, th1CanvTot[iter][0][0], 1000, dirNames_p->at(iter));

    th1Canv_p[iter]->SetGlobalMaxMin();
    
    if(isMean){
      if(inHistName[inHistNum].find("RecoOverGen") != std::string::npos){
	th1Canv_p[iter]->CapGlobalMaxMin(1.3, 0.8);
	th1Canv_p[iter]->UnderCapGlobalMaxMin(1.3, 0.8);
      }
      if(inHistName[inHistNum].find("Eff") != std::string::npos){
	th1Canv_p[iter]->CapGlobalMaxMin(1.05, 0.0);
	th1Canv_p[iter]->UnderCapGlobalMaxMin(1.05, 0.0);
      }
    }
    if(isRes || isResOverMean){
      if(inHistName[inHistNum].find("DPhi") != std::string::npos){
	th1Canv_p[iter]->CapGlobalMaxMin(0.2, 0.0);
	th1Canv_p[iter]->UnderCapGlobalMaxMin(0.2, 0.0);
      }
      else th1Canv_p[iter]->CapGlobalMaxMin(0.7, 0.0);
    }

    th1Canv_p[iter]->canv_p->cd();
    th1Canv_p[iter]->canv_p->Draw();

    for(Int_t xIter = 0; xIter < nXPanel; xIter++){
      for(Int_t yIter = 0; yIter < nYPanel; yIter++){
	th1Canv_p[iter]->canv_p->cd();
	th1Canv_p[iter]->pads_p[xIter][yIter]->Draw();
	th1Canv_p[iter]->pads_p[xIter][yIter]->cd();

      }
    }
  }

  if(debugMode) std::cout << __LINE__ << std::endl;

  //  Bool_t legAdded[2][nQG] = {{false}, {false}};

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  //  Float_t dirMax[nDir];
  //  Float_t dirMin[nDir];
  /*
  for(Int_t iter = 0; iter < nDir; iter++){
  dirMax[iter] = -1;
    dirMin[iter] = 100;
    }*/

  Float_t th1XMin = 30;

  if(debugMode) std::cout << __LINE__ << std::endl;

  TLegend* meanLeg_p;
  if(isPbPb){
    if(!strcmp("Eff", inHistName[inHistNum].c_str())) meanLeg_p = new TLegend(.70, .08, .90, .55);
    if(isRes || isResOverMean) meanLeg_p = new TLegend(.70, .08, .90, .35);
    else meanLeg_p = new TLegend(.50, .68, .90, .95);
  }
  else if(isRes || isResOverMean) meanLeg_p = new TLegend(.3, .38, .90, .5);
  else if(!strcmp("Eff", inHistName[inHistNum].c_str())) meanLeg_p = new TLegend(.6, .25, .95, .75);
  else meanLeg_p = new TLegend(.50, .58, .90, .85);

  meanLeg_p->SetBorderSize(0);
  meanLeg_p->SetFillColor(0);
  meanLeg_p->SetFillStyle(0);
  meanLeg_p->SetTextFont(43);
  meanLeg_p->SetTextSize(18);


  TLegend* meanLeg2_p;
  if(!strcmp("Eff", inHistName[inHistNum].c_str())){
    if(isPbPb) meanLeg2_p = new TLegend(.60, .25, .95, .75);
    else  meanLeg2_p = new TLegend(.50, .35, .95, .85);
  }
  else if(isRes || isResOverMean) meanLeg2_p = new TLegend(.50, .38, .90, .65);
  else meanLeg2_p = new TLegend(.50, .68, .90, .95);
  meanLeg2_p->SetBorderSize(0);
  meanLeg2_p->SetFillColor(0);
  meanLeg2_p->SetFillStyle(0);
  meanLeg2_p->SetTextFont(43);
  meanLeg2_p->SetTextSize(18);

  if(debugMode) std::cout << __LINE__ << std::endl;

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    th1_p[iter] = (TH1F*)inFile_p->Get(th1Name.c_str());

    th1_p[iter]->SetMarkerSize(0.5);

    if(iter == 0 && strcmp("Pt", ptEtaStr2.c_str()) != 0) th1XMin = th1_p[iter]->GetBinCenter(1);

    Int_t xMinBin = th1_p[iter]->FindBin(th1XMin);

    if(iter == 0) th1XMin = th1_p[iter]->GetBinLowEdge(xMinBin);

    th1XMin = th1_p[iter]->GetXaxis()->GetXmin();
    //    if(!strcmp("Eta", ptEtaStr2[ptEtaNum].c_str())) th1XMin += .1;
    th1_p[iter]->SetAxisRange(th1XMin, th1_p[iter]->GetXaxis()->GetXmax(), "X");

    if(debugMode) std::cout << __LINE__ << std::endl;

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
      
      if(debugMode) std::cout << __LINE__ << std::endl;

      if(th1Name.substr(0, th1Name.find("/")).size() != dirNames_p->at(dirIter).size()) continue;

      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    if(debugMode) std::cout << __LINE__ << std::endl;


    th1Canv_p[dirPos]->canv_p->cd();
  
    Int_t centPos = -1;
    if(isPbPb){
      for(Int_t centIter = 0; centIter < nCentBins; centIter++){
	std::size_t pos = th1Name.find(centStrings[centIter]);
	if(pos != std::string::npos){
	  centPos = centIter;
	  break;
	}
      }
    }
    else centPos = 0;


    if(debugMode) std::cout << __LINE__ << std::endl;
    

    std::string legStr = "Incl.";
    if(inHistName[inHistNum].find("RecoOverGen") != std::string::npos){
      if(th1Name.find("_Q_") != std::string::npos) legStr = "Quark";
      else if(th1Name.find("_G_") != std::string::npos) legStr = "Gluon";
    }
    else if(inHistName[inHistNum].find("Eff") != std::string::npos){
      if(th1Name.find("RecoPtCut5") != std::string::npos) legStr = "Reco. p_{T} > 5";
      else if(th1Name.find("RecoPtCut10") != std::string::npos) legStr = "Reco. p_{T} > 10";
      else if(th1Name.find("RecoPtCut15") != std::string::npos) legStr = "Reco. p_{T} > 15";
      else if(th1Name.find("RecoPtCut20") != std::string::npos) legStr = "Reco. p_{T} > 20";
      else if(th1Name.find("RecoPtCut25") != std::string::npos) legStr = "Reco. p_{T} > 25";
      else if(th1Name.find("RecoPtCut30") != std::string::npos) legStr = "Reco. p_{T} > 30";
    }
    else if(inHistName[inHistNum].find("DPhi") != std::string::npos){
      if(th1Name.find("_Mean_") != std::string::npos) legStr = "Hist Mean";
      else if(th1Name.find("_FitMean_") != std::string::npos) legStr = "Fit Mean";
      else if(th1Name.find("_Res_") != std::string::npos) legStr = "Hist Res.";
      else if(th1Name.find("_FitRes_") != std::string::npos) legStr = "Fit Res.";
      else if(th1Name.find("_ResOverMean_") != std::string::npos) legStr = "Hist #frac{Res.}{Mean}";
      else if(th1Name.find("_FitResOverMean_") != std::string::npos) legStr = "Fit #frac{Res.}{Mean}";
    }

    if(debugMode) std::cout << __LINE__ << std::endl;

    th1Canv_p[dirPos]->SetHist(th1_p[iter], centPos, 0, th1CanvCount[dirPos][centPos][0], legStr);
    th1CanvCount[dirPos][centPos][0]++;

    if(debugMode) std::cout << __LINE__ << std::endl;
  }

  if(debugMode) std::cout << __LINE__ << std::endl;

  for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
    th1Canv_p[dirIter]->canv_p->cd();
    th1Canv_p[dirIter]->MakeHistMaxMinNice();
    th1Canv_p[dirIter]->SetHistMaxMin();

    gStyle->SetOptStat(0);

    th1Canv_p[dirIter]->SetPanelWhiteSpace();

    if(debugMode) std::cout << __LINE__ << ", " << dirIter << std::endl;
	  
    for(Int_t iter = 0; iter < nXPanel; iter++){
      for(Int_t iter2 = 0; iter2 < nYPanel; iter2++){
	if(debugMode) std::cout << __LINE__ << ", " << dirIter << ", " << iter << ", " << iter2 << std::endl;

	th1Canv_p[dirIter]->canv_p->cd();
	th1Canv_p[dirIter]->pads_p[iter][iter2]->cd();
	
	for(Int_t histIter = 0; histIter < th1CanvCount[dirIter][iter][iter2]; histIter++){
	  if(debugMode) std::cout << __LINE__ << ", " << dirIter << ", " << iter << ", " << iter2 << ", " << histIter << std::endl;

	  th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->SetMarkerSize(2);
	  th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->SetMarkerStyle(20);

	  

	  if(histIter == 0){

	    if(isRes || isResOverMean){
	      th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->SetMaximum(0.5);
	      th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->SetMinimum(0.0);
	    }

	    th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->GetXaxis()->SetRange(1, th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->GetNbinsX());
	    th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->DrawCopy("E1 P");

	    if(isRes || isResOverMean){
	      if(ptEtaStr.find("Pt") != std::string::npos){
		std::string th1Name = th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->GetName();
	       
		if(th1Name.find("VPt") != std::string::npos){
		  if(th1Name.find("RecoOverGen") != std::string::npos){
		    th1Name.replace(th1Name.find("RecoOverGen"), 11, "Pt");
		    if(isRes){
		      if(th1Name.find("FitRes") != std::string::npos) th1Name.replace(th1Name.find("FitRes"), 6, "Mean");
		      if(th1Name.find("Res") != std::string::npos) th1Name.replace(th1Name.find("Res"), 3, "Mean");
		      
		    }
		    else if(isResOverMean){
		      if(th1Name.find("FitResOverMean") != std::string::npos) th1Name.replace(th1Name.find("FitResOverMean"), 14, "Mean");
		      if(th1Name.find("ResOverMean") != std::string::npos) th1Name.replace(th1Name.find("ResOverMean"), 11, "Mean");
		    }

		    th1Name = dirNames_p->at(dirIter) + "/" + th1Name;

		    TH1F* meanLineHist_p = (TH1F*)inFile_p->Get(th1Name.c_str());

		    for(Int_t binIter = 0; binIter < meanLineHist_p->GetNbinsX(); binIter++){
		      TLine* line_p = new TLine();
		      line_p->SetLineColor(1);
		      line_p->SetLineStyle(2);

		      line_p->DrawLine(meanLineHist_p->GetBinContent(binIter+1), th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->GetBinContent(binIter+1) - 0.05, meanLineHist_p->GetBinContent(binIter+1), th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->GetBinContent(binIter+1) + 0.05);
		    }
		  }
		}
	      }

	    }
	  }
	  else th1Canv_p[dirIter]->hists_p[iter][iter2][histIter]->DrawCopy("SAME E1 P");
	}

	if(debugMode) std::cout << __LINE__ << std::endl;

	if(iter == 0) th1Canv_p[dirIter]->DrawLabel1(iter, 0, dirNames_p->at(dirIter));
	if(isPbPb) th1Canv_p[dirIter]->DrawLabel2(iter, 0, centStrings2[iter]);
	else th1Canv_p[dirIter]->DrawLabel2(iter, 0, "PP");
      }
    }
    
    if(debugMode) std::cout << __LINE__ << std::endl;

    //    th1Canv_p[dirIter]->DrawLegend();

    if(debugMode) std::cout << __LINE__ << std::endl;

    if(ptEtaStr.find("VPt") != std::string::npos){
      for(Int_t iter = 0; iter < nXPanel; iter++){
       	th1Canv_p[dirIter]->SetColumnLogX(iter); 
      }
    }
    
  }

  if(debugMode) std::cout << __LINE__ << std::endl;


  std::string outName = inFileName;
  const std::string inString = "_HIST.root";
  const std::string outString = "_PLOT.root";
  std::size_t strIndex = 0;

  strIndex = outName.find(inString);
  if(!(strIndex == std::string::npos)){
    outName.replace(strIndex, inString.length(), outString);
  }

  TF1* csn_p = 0;
  TF1* csn2_p = 0;
  TF1* csn3_p = 0;

  TF1* csn_PP_p = 0;
  TF1* csn2_PP_p = 0;
  TF1* csn3_PP_p = 0;

  if(debugMode) std::cout << __LINE__ << std::endl;

  if(isPbPb){
    csn_p = new TF1("csn_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))", 30, 150);
    csn_p->FixParameter(0, 0.061);
    csn_p->FixParameter(1, 1.24);
    csn_p->FixParameter(2, 8.08);
    
    csn2_p= new TF1("csn2_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))", 30, 150);
    csn2_p->FixParameter(0, 0.061 - .001);
    csn2_p->FixParameter(1, 1.24 - .04);
    csn2_p->FixParameter(2, 8.08 - .15);
    csn2_p->SetLineStyle(2);
    
    csn3_p = new TF1("csn3_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))", 30, 150);
    csn3_p->FixParameter(0, 0.061 + .001);
    csn3_p->FixParameter(1, 1.24 + .04);
    csn3_p->FixParameter(2, 8.08 + .15);
    csn3_p->SetLineStyle(2);
  }
  else{
    csn_p = new TF1("csn_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x)", 30, 150);
    csn_p->FixParameter(0, 0.061);
    csn_p->FixParameter(1, 1.24);
    
    csn2_p= new TF1("csn2_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x)", 30, 150);
    csn2_p->FixParameter(0, 0.061 - .001);
    csn2_p->FixParameter(1, 1.24 - .04);
    csn2_p->SetLineStyle(2);
    
    csn3_p = new TF1("csn3_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x)", 30, 150);
    csn3_p->FixParameter(0, 0.061 + .001);
    csn3_p->FixParameter(1, 1.24 + .04);
    csn3_p->SetLineStyle(2);

    csn_PP_p = new TF1("csn_PP_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x)", 30, 150);
    csn_PP_p->FixParameter(0, 0.061);
    csn_PP_p->FixParameter(1, .95);
    csn_PP_p->SetLineColor(kBlue);
    csn_PP_p->SetMarkerColor(kBlue);

    csn2_PP_p= new TF1("csn2_PP_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x)", 30, 150);
    csn2_PP_p->FixParameter(0, 0.061 - .001);
    csn2_PP_p->FixParameter(1, .95 - .04);
    csn2_PP_p->SetLineStyle(2);
    csn2_PP_p->SetLineColor(kBlue);
    csn2_PP_p->SetMarkerColor(kBlue);
    
    csn3_PP_p = new TF1("csn3_PP_p", "TMath::Sqrt([0]*[0] + [1]*[1]/x)", 30, 150);
    csn3_PP_p->FixParameter(0, 0.061 + .001);
    csn3_PP_p->FixParameter(1, .95 + .04);
    csn3_PP_p->SetLineStyle(2);
    csn3_PP_p->SetLineColor(kBlue);
    csn3_PP_p->SetMarkerColor(kBlue);
  }

  TFile* outFile_p = new TFile(outName.c_str(), "UPDATE");
  if(debugMode) std::cout << __LINE__ << std::endl;

  for(Int_t iter = 0; iter < nDir; iter++){
    if(inHistName[inHistNum].find("RecoOverGen") != std::string::npos){
      if(isMean){
	th1Canv_p[iter]->DrawGlobalHorizontalLine(1);
	th1Canv_p[iter]->DrawGlobalHorizontalLine(.95);
	th1Canv_p[iter]->DrawGlobalHorizontalLine(1.05);
      }
      else if(isRes || isResOverMean){
	th1Canv_p[iter]->pads_p[0][0]->cd();
	csn_p->Draw("SAME C");
	csn2_p->Draw("SAME C");
	csn3_p->Draw("SAME C");
	
	if(!isPbPb){
	  csn_PP_p->Draw("SAME C");
	  csn2_PP_p->Draw("SAME C");
	  csn3_PP_p->Draw("SAME C");
	}
      }

      th1Canv_p[iter]->DrawGlobalVerticalLine(30);
    }
    if(debugMode) std::cout << __LINE__ << std::endl;

    if(inHistName[inHistNum].find("Eff") != std::string::npos && isMean){
      th1Canv_p[iter]->DrawGlobalHorizontalLine(1);
      th1Canv_p[iter]->DrawGlobalVerticalLine(30);
    }

    if(inHistName[inHistNum].find("DPhi") != std::string::npos){
      th1Canv_p[iter]->DrawGlobalVerticalLine(30);
    }

    if(debugMode) std::cout << __LINE__ << std::endl;

    std::string resCorrStr = "";
    if(addResCorrStr) resCorrStr = "_RESCORR";

    if(debugMode) std::cout << __LINE__ << std::endl;

    //    th1Canv_p[iter]->canv_p->Write(Form("%s_%s%s%sc", dirNames_p->at(iter).c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr.c_str()), TObject::kOverwrite);
    
    if(debugMode) std::cout << __LINE__ << std::endl;

    std::string tempFileName = config.GetConfigFileNameNoExt().c_str();
    while(tempFileName.find("/") != std::string::npos){
      tempFileName.replace(0, tempFileName.find("/")+1, "");
    }

    std::string pdfName =  Form("pdfDir/%s_%s%s%s%sc_%s", dirNames_p->at(iter).c_str(), fitMeanStr.c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr.c_str(), tempFileName.c_str());

    gStyle->SetOptStat(0);

    claverCanvasSaving(th1Canv_p[iter]->canv_p, pdfName.c_str(), "pdf");
    TDatime* date = new TDatime();
    pdfName = pdfName + "_" + std::to_string(date->GetDate()) + ".pdf";
    delete date;
    pdfList_p->push_back(pdfName);
  }

  if(debugMode) std::cout << __LINE__ << std::endl;

  outFile_p->Close();
  delete outFile_p;

  delete meanLeg_p;

  delete label_p;

  

  for(Int_t iter = 0; iter < nDir; iter++){
    delete th1Canv_p[iter]->canv_p;
    delete th1Canv_p[iter];
  }

  th1Names_p->clear();
  delete th1Names_p;

  dirNames_p->clear();
  delete dirNames_p;

  inFile_p->Close();
  delete inFile_p;

  return 0;
}


int makeJECPlotMeanPts(const std::string inFileName, jecConfigParser config, const Int_t inHistNum, const std::string ptEtaStr, const Int_t qgNum, std::vector<std::string>* pdfList_p)
{
  std::string ptEtaStr2 = "Pt";
  std::string ptEtaStr3 = "p_{T}";
  if(ptEtaStr.find("VEta") != std::string::npos){
    ptEtaStr2 = "Eta";
    ptEtaStr3 = "#eta";
  }

  const Bool_t isPbPb = config.GetIsPbPb();

  const Int_t nCentBins = config.GetNCentBins();
  std::vector<unsigned int> centBins = config.GetCentBins();

  if(isPbPb){
    if((unsigned int)(nCentBins+1) != centBins.size()){
      std::cout << "Config nCentBins \'" << nCentBins << "\' != centBins.size \'" << centBins.size() << "\'. return 1" << std::endl;
      return 1;
    }
  }

  std::string centStrings[nCentBins];
  std::string centStrings2[nCentBins];

  for(Int_t iter = 0; iter < nCentBins; iter++){
    if(isPbPb){
      centStrings[iter] = "Cent" + std::to_string(centBins.at(iter)) + "to" + std::to_string(centBins.at(iter+1));
      centStrings2[iter] = std::to_string(centBins.at(iter)) + "-" + std::to_string(centBins.at(iter+1)) + "%";
    }
    else{
      centStrings[iter] = "PP";
      centStrings2[iter] = "PP";
    }
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  const Int_t nContents = inFile_p->GetListOfKeys()->GetEntries();

  Int_t nTH1Temp = 0;
  std::vector<std::string>* th1Names_p = new std::vector<std::string>;

  Int_t nDirTemp = 0;
  std::vector<std::string>* dirNames_p = new std::vector<std::string>;

  Int_t nPtBinsTemp = 0;
  std::vector<std::string>* ptBins_p = new std::vector<std::string>;
  std::vector<Int_t>* ptInts_p = new std::vector<Int_t>;
  std::vector<Int_t>* ptDec_p = new std::vector<Int_t>;

  for(Int_t iter = 0; iter < nContents; iter++){
    TString className = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetClassName();
    TString name = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetName();

    if(className.Index("TDirectoryFile") >= 0 && name.Index("ak") >= 0){
      nDirTemp++;
      dirNames_p->push_back(name.Data());

      TDirectoryFile* tempDir_p = (TDirectoryFile*)inFile_p->Get(name);

      const Int_t nDirContents = tempDir_p->GetListOfKeys()->GetEntries();
      
      for(Int_t dirIter = 0; dirIter < nDirContents; dirIter++){
	TString name2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetName();
	TString className2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetClassName();

	if(className2.Index("TH1") >= 0){
	  if(name2.Index(ptEtaStr.c_str()) < 0) continue;
	  if(name2.Index(Form("_%s_", qgStr[qgNum].c_str())) < 0) continue;
	  if(name2.Index(Form("%s", meanResPtsStr.c_str())) < 0) continue;
	  if(name2.Index(inHistName[inHistNum].c_str()) < 0) continue;
	      

	  TString name3 = name2(name2.Index(Form("%s", meanResPtsStr.c_str()))+meanResPtsStr.size(), name2.Length());
	  name3 = name3(name3.Index(Form("_%s", ptEtaStr2.c_str()))+1, name3.Length());
	  TString name4 = name3(0, name3.Index("_"));
	  std::string ptRange = name4.Data();
	  
	  Bool_t addPtRange = true;
	  
	  for(Int_t ptIter = 0; ptIter < (Int_t)ptBins_p->size(); ptIter++){
	    if(!strcmp(ptRange.c_str(), ptBins_p->at(ptIter).c_str())){
	      addPtRange = false;
	      break;
	    }
	  }
	  
	  if(addPtRange){
	    ptBins_p->push_back(ptRange);
	    ptInts_p->push_back(0);
	    ptDec_p->push_back(0);
	    nPtBinsTemp++;
	  }

	  nTH1Temp++;
	  th1Names_p->push_back(Form("%s/%s", name.Data(), name2.Data()));
	}
      }
    }
  }


  ptInts_p->push_back(0);
  ptDec_p->push_back(0);

  Int_t sortIter = 0; 

  while(sortIter < (Int_t)(ptBins_p->size()-1)){
    std::string ptBins1Temp = ptBins_p->at(sortIter);
    std::string ptBins2Temp = ptBins_p->at(sortIter+1);

    while(ptBins1Temp.find("Neg") != std::string::npos){
      ptBins1Temp.replace(ptBins1Temp.find("Neg"), 3, "-");
    }

    while(ptBins2Temp.find("Neg") != std::string::npos){
      ptBins2Temp.replace(ptBins2Temp.find("Neg"), 3, "-");
    }


    TString ptBins1 = ptBins1Temp;
    TString ptBins2 = ptBins2Temp;   

    TString pt1Part1 = ptBins1(ptBins1.Index(ptEtaStr2.c_str()) + ptEtaStr2.size(), TString(ptBins1(ptBins1.Index(ptEtaStr2.c_str())+ptEtaStr2.size(), ptBins1.Length() - ptEtaStr2.size())).Index("p"));
    TString pt1Part2 = ptBins1(ptBins1.Index("p")+1, 1);
    if(!strcmp(pt1Part2.Data(), "-")) pt1Part2 = ptBins1(ptBins1.Index("p")+1, 2);

    TString pt2Part1 = ptBins2(ptBins2.Index(ptEtaStr2.c_str()) + ptEtaStr2.size(), TString(ptBins2(ptBins2.Index(ptEtaStr2.c_str()) + ptEtaStr2.size(), ptBins2.Length() - ptEtaStr2.size())).Index("p"));
    TString pt2Part2 = ptBins2(ptBins2.Index("p")+1, 1);
    if(!strcmp(pt2Part2.Data(), "-")) pt2Part2 = ptBins2(ptBins2.Index("p")+1, 2);

    Float_t pt1 = 0;
    Float_t pt2 = 0;
    Bool_t pt1IsNeg = false;
    Bool_t pt2IsNeg = false;

    if(pt1Part1.Index("-") >= 0){
      pt1Part1 = pt1Part1(1, pt1Part1.Length() - 1);
      pt1IsNeg = true;
    }
    if(pt1Part2.Index("-") >= 0){
      pt1Part2 = pt1Part2(1, pt1Part2.Length() - 1);
      pt1IsNeg = true;
    }
    if(pt2Part1.Index("-") >= 0){
      pt2Part1 = pt2Part1(1, pt2Part1.Length() - 1);
      pt2IsNeg = true;
    }
    if(pt2Part2.Index("-") >= 0){
      pt2Part2 = pt2Part2(1, pt2Part2.Length() - 1);
      pt2IsNeg = true;
    }
  
    for(Int_t iter = 0; iter < 10; iter++){
      for(Int_t iter2 = 0; iter2 < Int_t(pt1Part1.Length()); iter2++){
	if(((TString)(pt1Part1(iter2, 1))).Index(Form("%d", iter)) >= 0) pt1 += TMath::Power(10, pt1Part1.Length() - iter2 - 1)*iter;
      }

      for(Int_t iter2 = 0; iter2 < Int_t(pt2Part1.Length()); iter2++){
	if(((TString)(pt2Part1(iter2, 1))).Index(Form("%d", iter)) >= 0) pt2 += TMath::Power(10, pt2Part1.Length() - iter2 - 1)*iter;
      }


      if(((TString)(pt1Part2(0, 1))).Index(Form("%d", iter)) >= 0) pt1 += iter/10.;
      if(((TString)(pt2Part2(0, 1))).Index(Form("%d", iter)) >= 0) pt2 += iter/10.;
    }



    Float_t tempPt1 = pt1;
    Float_t tempPt2 = pt2;

    if(pt1IsNeg) pt1 = -pt1;
    if(pt2IsNeg) pt2 = -pt2;
   
    if(pt1 > pt2){
      std::string tempStr = ptBins_p->at(sortIter); 
      ptBins_p->at(sortIter) = ptBins_p->at(sortIter+1);
      ptBins_p->at(sortIter+1) = tempStr;
      
      ptInts_p->at(sortIter) = std::floor(tempPt1);
      ptDec_p->at(sortIter) = tempPt1*10 - std::floor(tempPt1)*10;

      if(pt1IsNeg) ptInts_p->at(sortIter) = -ptInts_p->at(sortIter);

      ptInts_p->at(sortIter+1) = std::floor(tempPt2);
      ptDec_p->at(sortIter+1) = tempPt2*10 - std::floor(tempPt2)*10;

      if(pt2IsNeg) ptInts_p->at(sortIter+1) = -ptInts_p->at(sortIter+1);

      sortIter = 0;
    }
    else{
      ptInts_p->at(sortIter) = std::floor(tempPt1);
      ptDec_p->at(sortIter) = tempPt1*10 - std::floor(tempPt1)*10;

      if(pt1IsNeg) ptInts_p->at(sortIter) = -ptInts_p->at(sortIter);

      ptInts_p->at(sortIter+1) = std::floor(tempPt2);
      ptDec_p->at(sortIter+1) = tempPt2*10 - std::floor(tempPt2)*10;

      if(pt2IsNeg) ptInts_p->at(sortIter+1) = -ptInts_p->at(sortIter+1);

      sortIter++;
    }
  }

  Int_t nCentBinsTemp = 1;
  if(isPbPb) nCentBinsTemp = nCentBins;
  const Int_t nCentBins2 = nCentBinsTemp;

  const Int_t nTH1 = nTH1Temp;
  const Int_t nDir = nDirTemp;
  const Int_t nPtBins = nPtBinsTemp;
  TH1F* th1_p[nTH1];
  TCanvas* th1Canv_p[nDir][nCentBins2];

  Int_t xBins = 0;
  Int_t yBins = 0;
  Float_t breakBool = false;

  for(Int_t iter = 0; iter < nPtBins; iter++){
    if(breakBool) break;
    for(Int_t iter2 = 0; iter2 < 4; iter2++){
   
      xBins = iter + 4 - iter2;
      yBins = iter;

      //      if(xBins == nPtBins-1) continue;

      if(xBins >= 6){
	if(xBins*yBins > nPtBins && xBins*yBins - nPtBins < 4){
	  breakBool = true;
	  break;
	}
      }
      else if(xBins*yBins > nPtBins && xBins*yBins - nPtBins < ((Float_t)xBins)/2.){
	breakBool = true;
	break;
      }
    }
  }

  if(xBins == nPtBins && nPtBins >= 8){
    xBins = nPtBins/2;
    yBins = 2;

    if(xBins*yBins < nPtBins) yBins = 3;
  }

  std::string tempFileName = config.GetConfigFileNameNoExt().c_str();
  while(tempFileName.find("/") != std::string::npos){
    tempFileName.replace(0, tempFileName.find("/")+1, "");
  }

  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      std::string centStr = centStrings[centIter];
      if(!isPbPb) centStr = "PP";

      th1Canv_p[iter][centIter] = new TCanvas(Form("%s_%s_MeanPts%s%s%s_c_%s", dirNames_p->at(iter).c_str(), qgStr[qgNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr.c_str(), centStr.c_str(), tempFileName.c_str()), Form("%s_%s_MeanPts%s%s%s_c_%s", dirNames_p->at(iter).c_str(), qgStr[qgNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr.c_str(), centStr.c_str(), tempFileName.c_str()), xBins*300, yBins*325);
      //      th1Canv_p[iter][centIter]->Divide(xBins, yBins, 0.0, 0.0);
      th1Canv_p[iter][centIter]->Divide(xBins, yBins, .000005, 0.000005);
    }
  }

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  Float_t dirMax[nDir][nCentBins2];
  Float_t dirMin[nDir][nCentBins2];

  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      dirMax[iter][centIter] = -1;
      dirMin[iter][centIter] = 100;
    }
  }

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    th1_p[iter] = (TH1F*)inFile_p->Get(th1Name.c_str());

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));

      if(th1Name.substr(0, th1Name.find("/")).size() != dirNames_p->at(dirIter).size()) continue;

      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    Int_t centPos = -1;
    if(isPbPb){
      for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
	std::size_t pos = th1Name.find(centStrings[centIter]);
	if(pos != std::string::npos){
	  centPos = centIter;
	  break;
	}
      }
    }
    else centPos = 0;



    for(Int_t binIter = 0; binIter < th1_p[iter]->GetNbinsX(); binIter++){
      if(th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1) > dirMax[dirPos][centPos]) dirMax[dirPos][centPos] = th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1);

      if(th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1) < dirMin[dirPos][centPos]) dirMin[dirPos][centPos] = th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1);
    }
  }


  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      dirMax[iter][centIter] = setMaxMinNice(dirMax[iter][centIter], true);
      dirMin[iter][centIter] = setMaxMinNice(dirMin[iter][centIter], false);
    }
  }

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));

      if(th1Name.substr(0, th1Name.find("/")).size() != dirNames_p->at(dirIter).size()) continue;

      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    Int_t centPos = -1;
    if(isPbPb){
      for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
	std::size_t pos = th1Name.find(centStrings[centIter]);
	if(pos != std::string::npos){
	  centPos = centIter;
	  break;
	}
      }
    }
    else centPos = 0;

    Int_t ptPos = -1;
    for(Int_t ptIter = 0; ptIter < (Int_t)ptBins_p->size(); ptIter++){
      std::size_t pos = th1Name.find(ptBins_p->at(ptIter).c_str());
      if(pos != std::string::npos){
        ptPos = ptIter;
        break;
      }
    }

    th1Canv_p[dirPos][centPos]->cd(ptPos+1);

    gPad->SetLeftMargin(gPad->GetLeftMargin()/2);
    gPad->SetBottomMargin(gPad->GetBottomMargin()/2);

    gStyle->SetOptStat(0);

    //    th1_p[iter]->SetMaximum(dirMax[dirPos][centPos]);
    //    th1_p[iter]->SetMinimum(dirMin[dirPos][centPos]);

    th1_p[iter]->GetYaxis()->SetNdivisions(404);

    th1_p[iter]->GetXaxis()->CenterTitle();
    th1_p[iter]->GetXaxis()->SetTitleOffset(th1_p[iter]->GetXaxis()->GetTitleOffset() + 2.0);
    th1_p[iter]->GetXaxis()->SetTitleFont(43);
    th1_p[iter]->GetXaxis()->SetTitleSize(26);
    th1_p[iter]->GetYaxis()->CenterTitle();
    th1_p[iter]->GetYaxis()->SetTitleOffset(th1_p[iter]->GetYaxis()->GetTitleOffset() + 2.5);
    
    th1_p[iter]->GetXaxis()->SetTitle(xAxisLabel[inHistNum].c_str());
    th1_p[iter]->GetYaxis()->SetTitle("Events");

    gPad->SetTicks(gPad->GetTickx(), 1);

    th1_p[iter]->DrawCopy("E1 P Y+");
    
    TLine* meanLine_p = new TLine(th1_p[iter]->GetMean(), th1_p[iter]->GetMinimum(), th1_p[iter]->GetMean(), th1_p[iter]->GetMaximum()*.8);
    meanLine_p->SetLineStyle(2);
    meanLine_p->DrawClone();
    delete meanLine_p;
    
    if(th1_p[iter]->GetListOfFunctions()->GetSize() != 0){
      TLine* meanLine2_p = new TLine(th1_p[iter]->GetFunction("f1_p")->GetParameter(1), th1_p[iter]->GetMinimum(), th1_p[iter]->GetFunction("f1_p")->GetParameter(1), th1_p[iter]->GetMaximum()*.8);
      meanLine2_p->SetLineStyle(2);
      meanLine2_p->SetLineColor(kRed);
      meanLine2_p->DrawClone();
      delete meanLine2_p;
    }

    //    if(centPos == 0) label_p->DrawLatex(.30, .9, Form("#bf{%s}", dirNames_p->at(dirPos).c_str()));
    
    if(ptPos%xBins == 0) label_p->DrawLatex(.30, .88, Form("%d.%d<%s<%d.%d", ptInts_p->at(ptPos), ptDec_p->at(ptPos), ptEtaStr3.c_str(), ptInts_p->at(ptPos+1), ptDec_p->at(ptPos+1)));
    else label_p->DrawLatex(.20, .88, Form("%d.%d<%s<%d.%d", ptInts_p->at(ptPos), ptDec_p->at(ptPos), ptEtaStr3.c_str(), ptInts_p->at(ptPos+1), ptDec_p->at(ptPos+1)));


    if(ptPos == 0){
      if(isPbPb) label_p->DrawLatex(.30, .78, Form("#bf{#color[2]{%s (%s)}}", dirNames_p->at(dirPos).c_str(), centStrings2[centPos].c_str()));
      else label_p->DrawLatex(.30, .78, Form("#bf{#color[2]{%s (%s)}}", dirNames_p->at(dirPos).c_str(), "PP"));
    }
    if(ptPos == 1) label_p->DrawLatex(.15, .78, Form("#bf{#color[2]{|#eta_{jet}|<1.6}}"));
    if(ptPos == 2 && !strcmp(ptEtaStr2.c_str(), "Eta" )) label_p->DrawLatex(.15, .78, Form("#bf{#color[2]{Gen. p_{T}>30}}"));
    if(ptPos == 3) label_p->DrawLatex(.15, .78, Form("#bf{#color[2]{%s}}", qgStr2[qgNum].c_str()));

    //    if(centPos%xBins == 0) label_p->DrawLatex(.68, .9, centStrings2[centPos].c_str());
    //    else label_p->DrawLatex(.60, .9, centStrings2[centPos].c_str());
  }

  std::string outName = inFileName;
  const std::string inString = "_HIST.root";
  const std::string outString = "_PLOT.root";
  std::size_t strIndex = 0;

  strIndex = outName.find(inString);
  if(!(strIndex == std::string::npos)){
    outName.replace(strIndex, inString.length(), outString);
  }


  TFile* outFile_p = new TFile(outName.c_str(), "UPDATE");

  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      //      th1Canv_p[iter][centIter]->Write("", TObject::kOverwrite);

      std::string pdfName = Form("pdfDir/%s", th1Canv_p[iter][centIter]->GetName());
      claverCanvasSaving(th1Canv_p[iter][centIter], pdfName.c_str(), "pdf");
      TDatime* date = new TDatime();
      pdfName = pdfName + "_" + std::to_string(date->GetDate()) + ".pdf";
      pdfList_p->push_back(pdfName);
    }
  }

  outFile_p->Close();
  delete outFile_p;

  delete label_p;

  delete ptBins_p;
  delete ptInts_p;
  delete ptDec_p;

  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t centIter = 0; centIter < nCentBins2; centIter++){
      delete th1Canv_p[iter][centIter];
    }
  }

  th1Names_p->clear();
  delete th1Names_p;

  dirNames_p->clear();
  delete dirNames_p;

  inFile_p->Close();
  delete inFile_p;

  return 0;
}


int makeTEXPlots(const std::string inFileName, std::string texFileName, std::vector<std::string>* pdfList_p)
{
  std::string texDirName = texFileName;
  while(texDirName.find(".") != std::string::npos){
    texDirName.replace(texDirName.find("."), texDirName.size() - texDirName.find("."), "");
  }

  checkMakeDir("texDir");

  texDirName = "texDir/" + texDirName + "_TeXDir";
  checkMakeDir(texDirName);
  std::string texDirNamePDF = texDirName + "/pdfDir";
  checkMakeDir(texDirNamePDF);

  texFileName = texDirName + "/" + texFileName;

  std::ofstream texFile(texFileName.c_str());
  texFile.close();

  texFile.open(texFileName.c_str(), std::ios::app);

  texFile << "\\RequirePackage{xspace}" << std::endl;
  texFile << "\\RequirePackage{amsmath}" << std::endl;
  texFile << std::endl;

  texFile << "\\documentclass[xcolor=dvipsnames]{beamer}" << std::endl;
  texFile << "\\usetheme{Warsaw}" << std::endl;
  texFile << "\\setbeamercolor{structure}{fg=NavyBlue!90!NavyBlue}" << std::endl;
  texFile << "\\setbeamercolor{footlinecolor}{fg=white,bg=lightgray}" << std::endl;
  texFile << std::endl;

  texFile << "\\newcommand{\\pt}{\\ensuremath{p_{\\mathrm{T}}}\\xspace}" << std::endl;
  texFile << std::endl;

  texFile << "\\setbeamersize{text margin left=5pt,text margin right=5pt}" << std::endl;
  texFile << std::endl;

  texFile << "\\setbeamertemplate{frametitle}" << std::endl;
  texFile << "{" << std::endl;
  texFile << "    \\nointerlineskip" << std::endl;
  texFile << "    \\begin{beamercolorbox}[sep=0.3cm, ht=1.8em, wd=\\paperwidth]{frametitle}" << std::endl;
  texFile << "        \\vbox{}\\vskip-2ex%" << std::endl;
  texFile << "        \\strut\\insertframetitle\\strut" << std::endl;
  texFile << "        \\vskip-0.8ex%" << std::endl;
  texFile << "    \\end{beamercolorbox}" << std::endl;
  texFile << "}" << std::endl;
  texFile << std::endl;

  texFile << "\\setbeamertemplate{footline}{%" << std::endl;
  texFile << "  \\begin{beamercolorbox}[sep=.8em,wd=\\paperwidth,leftskip=0.5cm,rightskip=0.5cm]{footlinecolor}" << std::endl;
  texFile << "    \\hspace{0.3cm}%" << std::endl;
  texFile << "    \\hfill\\insertauthor \\hfill\\insertpagenumber" << std::endl;
  texFile << "  \\end{beamercolorbox}%" << std::endl;
  texFile << "}" << std::endl;
  texFile << "\\setbeamertemplate{navigation symbols}{}" << std::endl;
  texFile << std::endl;

  texFile << "\\setbeamertemplate{itemize item}[circle]" << std::endl;
  texFile << "\\setbeamertemplate{itemize subitem}[circle]" << std::endl;
  texFile << "\\setbeamertemplate{itemize subsubitem}[circle]" << std::endl;
  texFile << "\\setbeamercolor{itemize item}{fg=black}" << std::endl;
  texFile << "\\setbeamercolor{itemize subitem}{fg=black}" << std::endl;
  texFile << "\\setbeamercolor{itemize subsubitem}{fg=black}" << std::endl;
  texFile << std::endl;

  texFile << "\\definecolor{links}{HTML}{00BFFF}" << std::endl;
  texFile << "\\hypersetup{colorlinks,linkcolor=,urlcolor=links}" << std::endl;
  texFile << std::endl;

  texFile << "\\author[CM]{Placeholder}" << std::endl;
  texFile << std::endl;

  texFile << "\\begin{document}" << std::endl;
  texFile << "\\begin{frame}" << std::endl;
  TDatime* date = new TDatime();
  texFile << "\\frametitle{\\centerline{JEC Validation (" << date->GetYear() << "." << date->GetMonth() << "." << date->GetDay() << ")}}" << std::endl;
  delete date;
  texFile << " \\begin{itemize}" << std::endl;
  texFile << "  \\fontsize{10}{10}\\selectfont" << std::endl;
  texFile << "  \\item{Placeholder}" << std::endl;
  texFile << "  \\begin{itemize}" << std::endl;
  texFile << "   \\fontsize{10}{10}\\selectfont" << std::endl;
  texFile << "   \\item{Placeholder}" << std::endl;
  texFile << "  \\end{itemize}" << std::endl;
  texFile << " \\end{itemize}" << std::endl;
  texFile << "\\end{frame}" << std::endl;

  texFile << std::endl;
  
  for(unsigned int iter = 0; iter < pdfList_p->size(); iter++){
    std::ifstream src(pdfList_p->at(iter).c_str(), std::ios::binary);

    std::string destName = pdfList_p->at(iter);
    while(destName.find("/") != std::string::npos){
      destName.replace(0, destName.find("/")+1, "");
    }
    destName = texDirNamePDF + "/" + destName;

    std::ofstream dest(destName.c_str(), std::ios::binary);
    dest << src.rdbuf();

    dest.close();
    src.close();
  }

  for(unsigned int iter = 0; iter < pdfList_p->size(); iter++){
    //  if(pdfList_p->at(iter).find("EtaInc") == std::string::npos && pdfList_p->at(iter).find("PtInc") == std::string::npos && pdfList_p->at(iter).find("QGDistrib") == std::string::npos) continue;

    if(pdfList_p->at(iter).find("_Eta") == std::string::npos && pdfList_p->at(iter).find("_Pt") == std::string::npos) continue;

    if(pdfList_p->at(iter).find("_Q_") != std::string::npos) continue;
    if(pdfList_p->at(iter).find("_G_") != std::string::npos) continue;
    if(pdfList_p->at(iter).find("_Untagged_") != std::string::npos) continue;

    std::string placeHolderName = pdfList_p->at(iter);
    while(placeHolderName.find("/") != std::string::npos){
      placeHolderName.replace(0, placeHolderName.find("/")+1, "");
    }
    while(placeHolderName.find(".") != std::string::npos){
      placeHolderName.replace(placeHolderName.find("."), placeHolderName.size() - placeHolderName.find("."), "");
    }

    int charIter = 0;
    while(charIter < (int)placeHolderName.size()){
      if(placeHolderName.substr(charIter, 1).find("_") != std::string::npos){
	if(charIter == 0) placeHolderName.replace(charIter, 1, "\\_");
	else if(placeHolderName.substr(charIter-1, 1).find("\\") == std::string::npos) placeHolderName.replace(charIter, 1, "\\_");
	else charIter++;
      }
      else charIter++;
    }

    charIter = placeHolderName.size()-1;
    while(charIter >= 0){
      if(placeHolderName.substr(charIter, 1).find("_") != std::string::npos){
	placeHolderName.replace(charIter-1, placeHolderName.size()-(charIter-1), "");
	break;
      }
      else charIter--;
    }

    texFile << "\\begin{frame}" << std::endl;
    texFile << "\\frametitle{\\centerline{Placeholder}}" << std::endl;
    texFile << "\\begin{center}" << std::endl;
    texFile << "\\includegraphics[width=.6\\textwidth]{" << pdfList_p->at(iter) << "}" << std::endl;
    texFile << "\\end{center}" << std::endl;
    texFile << "\\begin{itemize}" << std::endl;
    texFile << "\\fontsize{8}{8}\\selectfont" << std::endl;
    texFile << "\\item{" << placeHolderName << "}" << std::endl;
    texFile << "\\end{itemize}" << std::endl;
    texFile << "\\end{frame}" << std::endl;

    texFile << std::endl;
  }


  texFile << "\\begin{frame}" << std::endl;
  texFile << "\\frametitle{\\centerline{Conclusions}}" << std::endl;
  texFile << " \\begin{itemize}" << std::endl;
  texFile << "  \\fontsize{10}{10}\\selectfont" << std::endl;
  texFile << "  \\item{Placeholder}" << std::endl;
  texFile << "  \\begin{itemize}" << std::endl;
  texFile << "   \\fontsize{10}{10}\\selectfont" << std::endl;
  texFile << "   \\item{Placeholder}" << std::endl;
  texFile << "  \\end{itemize}" << std::endl;
  texFile << " \\end{itemize}" << std::endl;
  texFile << "\\end{frame}" << std::endl;
  
  texFile << std::endl;

  texFile << "\\begin{frame}" << std::endl;
  texFile << "\\frametitle{\\centerline{Backup}}" << std::endl;
  texFile << "\\end{frame}" << std::endl;

  texFile << std::endl;


  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  std::vector<std::string> tempConfigParams = returnRootFileContentsList(inFile_p, "TNamed", "configParamsDir");
  const Int_t nConfigParams = tempConfigParams.size();
  const Int_t nParamPerSlide = 16;

  for(Int_t iter = 0; iter < nConfigParams/nParamPerSlide + 1; iter++){
    if(iter == nConfigParams/nParamPerSlide && nConfigParams%nParamPerSlide == 0) break;

    texFile << "\\begin{frame}" << std::endl;
    texFile << "\\frametitle{\\centerline{Config. Params (" << iter+1  << "/" << nConfigParams/nParamPerSlide + 1 << ")}}" << std::endl;
    texFile << "\\begin{itemize}" << std::endl;
    texFile << "\\fontsize{6}{6}\\selectfont" << std::endl;
    
    for(Int_t iter2 = nParamPerSlide*iter; iter2 < TMath::Min(nParamPerSlide*(iter+1), nConfigParams); iter2++){
      TNamed* tempName_p = (TNamed*)inFile_p->Get(tempConfigParams.at(iter2).c_str());
      std::string tempConfigParam = tempConfigParams.at(iter2);
      std::string tempParamName = tempName_p->GetTitle();
      while(tempConfigParam.find("/") != std::string::npos){
	tempConfigParam.replace(0, tempConfigParam.find("/")+1, "");
      }
      
      int charIter = 0;
      while(charIter < (int)tempConfigParam.size()){
	if(tempConfigParam.substr(charIter, 1).find("_") != std::string::npos){
	  if(charIter == 0) tempConfigParam.replace(charIter, 1, "\\_");
	  else if(tempConfigParam.substr(charIter-1, 1).find("\\") == std::string::npos) tempConfigParam.replace(charIter, 1, "\\_");
	  else charIter++;
	}
	else charIter++;
      }

      charIter = 0;
      while(charIter < (int)tempParamName.size()){
	if(tempParamName.substr(charIter, 1).find("_") != std::string::npos){
	  if(charIter == 0) tempParamName.replace(charIter, 1, "\\_");
	  else if(tempParamName.substr(charIter-1, 1).find("\\") == std::string::npos) tempParamName.replace(charIter, 1, "\\_");
	  else charIter++;
	}
	else charIter++;
      }
      
      texFile << "\\item{" << tempConfigParam << "=" << tempParamName << "}" << std::endl;
    }

    texFile << "\\end{itemize}" << std::endl;
    texFile << "\\end{frame}" << std::endl;
    texFile << std::endl;
  }

  texFile << "\\end{document}" << std::endl;

  inFile_p->Close();
  delete inFile_p;

  texFile.close();

  return 0;
}


int makeJECPlot(const std::string inFileName, const bool produceTeX)
{
  jecConfigParser config;
  if(!config.SetConfigParser(inFileName)) return 1;

  const Int_t nJtPtEtaBins = config.GetNJtPtEtaBins();
  Double_t jtPtEtaBins[nJtPtEtaBins+1];
  if(config.GetDoJtPtEtaCustomBins()) config.FillJtPtEtaCustomBins(jtPtEtaBins);
  else getLinBins(config.GetJtPtEtaMin(), config.GetJtPtEtaMax(), nJtPtEtaBins, jtPtEtaBins);

  std::string jtPtEtaBinStrings[nJtPtEtaBins+1];
  jtPtEtaBinStrings[0] = "VPt_EtaInc_";
  for(Int_t ptEtaIter = 0; ptEtaIter < nJtPtEtaBins; ptEtaIter++){
    jtPtEtaBinStrings[ptEtaIter+1] = "VPt_Eta" + config.FloatRangeToTitleString(jtPtEtaBins[ptEtaIter], jtPtEtaBins[ptEtaIter+1]) + "_";
  }

  const Int_t nJtEtaPtBins = config.GetNJtEtaPtBins();
  const Float_t jtEtaPtMin = config.GetJtEtaPtMin();
  const Float_t jtEtaPtMax = config.GetJtEtaPtMax();  
  Double_t jtEtaPtBins[nJtEtaPtBins+1];
  if(config.GetDoJtEtaPtLogBins()) getLogBins(jtEtaPtMin, jtEtaPtMax, nJtEtaPtBins, jtEtaPtBins);
  else if(config.GetDoJtEtaPtCustomBins()) config.FillJtEtaPtCustomBins(jtEtaPtBins);
  else getLinBins(jtEtaPtMin, jtEtaPtMax, nJtEtaPtBins, jtEtaPtBins);

  std::string jtEtaPtBinStrings[nJtEtaPtBins+1];
  jtEtaPtBinStrings[0] = "VEta_PtInc_";
  for(Int_t etaPtIter = 0; etaPtIter < nJtEtaPtBins; etaPtIter++){
    jtEtaPtBinStrings[etaPtIter+1] = "VEta_Pt" + config.FloatRangeToTitleString(jtEtaPtBins[etaPtIter], jtEtaPtBins[etaPtIter+1]) + "_";
  }

  std::cout << "Eta partial strings: ";
  for(Int_t ptEtaIter = 0; ptEtaIter < nJtPtEtaBins+1; ptEtaIter++){
    std::cout << jtPtEtaBinStrings[ptEtaIter] << ", ";
  }
  std::cout << std::endl;


  std::cout << "Pt partial strings: ";
  for(Int_t etaPtIter = 0; etaPtIter < nJtEtaPtBins+1; etaPtIter++){
    std::cout << jtEtaPtBinStrings[etaPtIter] << ", ";
  }
  std::cout << std::endl;

  checkMakeDir("pdfDir");
  
  std::vector<std::string>* pdfList_p = new std::vector<std::string>;

  int retVal = 0;

  retVal += makeJECPlotBasicQGDistrib(inFileName, config, "jtPt_All_", false, false, pdfList_p);
  retVal += makeJECPlotBasicQGDistrib(inFileName, config, "jtPt_All_", false, true, pdfList_p);
  retVal += makeJECPlotBasicQGDistrib(inFileName, config, "jtPt_All_", true, false, pdfList_p);

  if(config.GetDoGenGammaCutOverride()){
    retVal += makeJECPlotBasicQGDistrib(inFileName, config, "jtPt_GammaOverride_", false, false, pdfList_p);
    retVal += makeJECPlotBasicQGDistrib(inFileName, config, "jtPt_GammaOverride_", false, true, pdfList_p);
    retVal += makeJECPlotBasicQGDistrib(inFileName, config, "jtPt_GammaOverride_", true, false, pdfList_p);
  }

  for(Int_t iter = 0; iter < nHistName; iter++){
    for(Int_t ptEtaIter = 0; ptEtaIter < nJtPtEtaBins+1; ptEtaIter++){
      for(Int_t iter2 = 0; iter2 < nMeanRes; iter2++){

	if(iter == 0){
	  if(debugMode) std::cout << __LINE__ << std::endl;

	  retVal += makeJECPlotMeanRes(inFileName, config, "", iter, iter2, jtPtEtaBinStrings[ptEtaIter], pdfList_p);
	  retVal += makeJECPlotMeanRes(inFileName, config, "Fit", iter, iter2, jtPtEtaBinStrings[ptEtaIter], pdfList_p);
	}
      }
      
      if(iter == 0){
	for(Int_t qgIter = 0; qgIter < 1; qgIter++){
	  if(debugMode) std::cout << __LINE__ << std::endl;
	  makeJECPlotMeanPts(inFileName, config, iter, jtPtEtaBinStrings[ptEtaIter], qgIter, pdfList_p);
	}
      }
    }

    for(Int_t etaPtIter = 0; etaPtIter < nJtEtaPtBins+1; etaPtIter++){
      for(Int_t iter2 = 0; iter2 < nMeanRes; iter2++){
	if(iter == 0){
	  retVal += makeJECPlotMeanRes(inFileName, config, "", iter, iter2, jtEtaPtBinStrings[etaPtIter], pdfList_p);
	  retVal += makeJECPlotMeanRes(inFileName, config, "Fit", iter, iter2, jtEtaPtBinStrings[etaPtIter], pdfList_p);
	}
      }
      
      if(iter == 0){
	for(Int_t qgIter = 0; qgIter < 1; qgIter++){
	  if(qgIter == 0) retVal += makeJECPlotMeanPts(inFileName, config, iter, jtEtaPtBinStrings[etaPtIter], qgIter, pdfList_p);
	}
      }
    }
  }

  std::string tempFileName = config.GetConfigFileNameNoExt().c_str();
  while(tempFileName.find("/") != std::string::npos){
    tempFileName.replace(0, tempFileName.find("/")+1, "");
  }

  if(produceTeX){
    TDatime* date = new TDatime();
    const std::string texName = tempFileName + "_" + std::to_string(date->GetDate()) + ".tex";
    std::cout << "TexName: " << texName << std::endl;
    delete date;
    makeTEXPlots(inFileName, texName, pdfList_p);
  }
  
  pdfList_p->clear();
  delete pdfList_p;
  
  return retVal;
}

int main(int argc, char *argv[])
{
  if(argc != 3){
    std::cout << "Usage: makeJECPlot_Prototype.exe <inHistFile> <produceTeX>" << std::endl;
    std::cout << "Number of args given: " << argc << std::endl;
    for(int iter = 0; iter < argc; iter++){
      std::cout << "  argv[" << iter << "]: " << argv[iter] << std::endl;
    }
    std::cout << "return -1" << std::endl;
    return -1;
  }

  std::string texBool = argv[2];

  if(texBool.size() != 1){
    std::cout << "argument <produceTeX> is not 0 or 1. return -1" << std::endl;
    return -1;
  }

  if(texBool.find("0") == std::string::npos && texBool.find("1") == std::string::npos){
    std::cout << "argument <produceTeX> is not 0 or 1. return -1" << std::endl;
    return -1;
  }

  int retVal = makeJECPlot(argv[1], std::stoi(argv[2]));
  return retVal;
}
