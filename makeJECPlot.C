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

#include <string>
#include <iostream>
#include <vector>

const Int_t nJetAlgo = 4;
const std::string jetAlgo[nJetAlgo] = {"akVs4Calo", "akPu4Calo", "akVs4PF", "akPu4PF"};


const Int_t nCentBins = 8;
const std::string centStrings[nCentBins] = {"cent70to100", "cent50to70", "cent40to50", "cent30to40", "cent20to30", "cent10to20", "cent5to10", "cent0to5"};
const std::string centStrings2[nCentBins] = {"70-100%", "50-70%", "40-50%", "30-40%", "20-30%", "10-20%", "5-10%", "0-5%"};

const Int_t nHistName = 4;
const std::string inHistName[nHistName] = {"RecoOverGen", "RecoOverRaw", "RawOverGen", "Eff"};
const std::string xAxisLabel[nHistName] = {"Reco./Gen.", "Reco./Raw", "Raw/Gen.", "Eff."};

const Int_t nMeanRes = 2;
const std::string meanResStr[nMeanRes] = {"Mean", "Res"};

const Int_t nPtEta = 3;
const std::string ptEtaStr[nPtEta] = {"VPt_", "VEta_", "VPtEta_"};
const std::string ptEtaStr2[nPtEta] = {"Pt", "Eta", "PtEta"};
const std::string ptEtaStr3[nPtEta] = {"p_{T}", "#eta", ""};


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


void makeJECPlotMeanRes(const std::string inFileName, const Int_t inHistNum, const Int_t meanResNum, const Int_t ptEtaNum)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  const Int_t nContents = inFile_p->GetListOfKeys()->GetEntries();

  Int_t nTH1Temp = 0;
  std::vector<std::string>* th1Names_p = new std::vector<std::string>;

  Int_t nDirTemp = 0;
  std::vector<std::string>* dirNames_p = new std::vector<std::string>;

  for(Int_t iter = 0; iter < nContents; iter++){
    TString className = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetClassName();
    TString name = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetName();

    if(className.Index("TDirectoryFile") >= 0){
      nDirTemp++;
      dirNames_p->push_back(name.Data());

      TDirectoryFile* tempDir_p = (TDirectoryFile*)inFile_p->Get(name);

      const Int_t nDirContents = tempDir_p->GetListOfKeys()->GetEntries();
      
      for(Int_t dirIter = 0; dirIter < nDirContents; dirIter++){
	TString name2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetName();
	TString className2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetClassName();

	if(className2.Index("TH1") >= 0){
	  if(name2.Index(Form("_%s_", meanResStr[meanResNum].c_str())) >= 0){

	    if(name2.Index(ptEtaStr[ptEtaNum].c_str()) < 0) continue;

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
  
  Float_t th1XMin = 30;

  TCanvas* th1Canv_p[nDir];

  for(Int_t iter = 0; iter < nDir; iter++){
    th1Canv_p[iter] = new TCanvas(Form("%s_%s%s%sc", dirNames_p->at(iter).c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str()), Form("%s_%s%s%sc", dirNames_p->at(iter).c_str(), meanResStr[meanResNum].c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str()), 4*300, 2*325);
    th1Canv_p[iter]->Divide(4, 2, 0.0, 0.0);
  }

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  Float_t dirMax[nDir];
  Float_t dirMin[nDir];

  for(Int_t iter = 0; iter < nDir; iter++){
    dirMax[iter] = -1;
    dirMin[iter] = 100;
  }


  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    th1_p[iter] = (TH1F*)inFile_p->Get(th1Name.c_str());

    Int_t xMinBin = th1_p[iter]->FindBin(th1XMin);

    if(iter == 0) th1XMin = th1_p[iter]->GetBinLowEdge(xMinBin);

    th1_p[iter]->SetAxisRange(th1XMin, th1_p[iter]->GetXaxis()->GetXmax(), "X");

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    for(Int_t binIter = xMinBin-1; binIter < th1_p[iter]->GetNbinsX(); binIter++){
      if(th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1) > dirMax[dirPos]) dirMax[dirPos] = th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1);
      
      if(th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1) < dirMin[dirPos] && th1_p[iter]->GetBinContent(binIter+1) != 0) dirMin[dirPos] = th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1);
    }
  }

  for(Int_t iter = 0; iter < nDir; iter++){
    dirMax[iter] = setMaxMinNice(dirMax[iter], true);
    dirMin[iter] = setMaxMinNice(dirMin[iter], false);
  }

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    Int_t centPos = -1;
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      std::size_t pos = th1Name.find(centStrings[centIter]);
      if(pos != std::string::npos){
        centPos = centIter;
        break;
      }
    }

    th1Canv_p[dirPos]->cd(centPos+1);

    th1_p[iter]->SetMaximum(dirMax[dirPos]);
    th1_p[iter]->SetMinimum(dirMin[dirPos]);

    th1_p[iter]->GetXaxis()->CenterTitle();
    th1_p[iter]->GetXaxis()->SetTitleOffset(th1_p[iter]->GetXaxis()->GetTitleOffset() + .9);
    th1_p[iter]->GetXaxis()->SetTitleFont(43);
    th1_p[iter]->GetXaxis()->SetTitleSize(26);
    th1_p[iter]->GetYaxis()->CenterTitle();
    th1_p[iter]->GetYaxis()->SetTitleOffset(th1_p[iter]->GetYaxis()->GetTitleOffset() + .5);
    
    if(!strcmp(ptEtaStr2[ptEtaNum].c_str(), "Pt")) gPad->SetLogx();
    
    th1_p[iter]->DrawCopy("E1 P");
    
    if(dirMin[dirPos] < 1 && dirMax[dirPos] > 1){
      TLine* oneLine_p = new TLine(th1XMin, 1, th1_p[iter]->GetXaxis()->GetXmax(), 1);
      oneLine_p->SetLineStyle(2);
      oneLine_p->DrawClone();
      delete oneLine_p;
    }    

    if(centPos == 0) label_p->DrawLatex(.30, .9, Form("#bf{#color[2]{%s}}", dirNames_p->at(dirPos).c_str()));

    if(centPos == 1) label_p->DrawLatex(.10, .9, Form("#bf{#color[2]{|#eta_{jet}|<2.0}}"));
    if(centPos == 2 && !strcmp(ptEtaStr2[ptEtaNum].c_str(), "Eta" )) label_p->DrawLatex(.10, .9, Form("#bf{#color[2]{Gen. p_{T}>30}}"));

    if(centPos%4 == 0) label_p->DrawLatex(.68, .9, centStrings2[centPos].c_str());
    else label_p->DrawLatex(.60, .9, centStrings2[centPos].c_str());
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
    th1Canv_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(th1Canv_p[iter], Form("pdfDir/%s", th1Canv_p[iter]->GetName()), "pdf");
  }

  outFile_p->Close();
  delete outFile_p;

  delete label_p;

  for(Int_t iter = 0; iter < nDir; iter++){
    delete th1Canv_p[iter];
  }

  th1Names_p->clear();
  delete th1Names_p;

  dirNames_p->clear();
  delete dirNames_p;

  inFile_p->Close();
  delete inFile_p;

  return;
}



void makeJECPlotMeanPts(const std::string inFileName, const Int_t inHistNum, const Int_t ptEtaNum)
{
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

    if(className.Index("TDirectoryFile") >= 0){
      nDirTemp++;
      dirNames_p->push_back(name.Data());

      TDirectoryFile* tempDir_p = (TDirectoryFile*)inFile_p->Get(name);

      const Int_t nDirContents = tempDir_p->GetListOfKeys()->GetEntries();
      
      for(Int_t dirIter = 0; dirIter < nDirContents; dirIter++){
	TString name2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetName();
	TString className2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetClassName();

	if(className2.Index("TH1") >= 0){
	  if(name2.Index(ptEtaStr[ptEtaNum].c_str()) < 0) continue;

	  if(name2.Index(Form("_%s", ptEtaStr2[ptEtaNum].c_str())) >= 0){
	    if(name2.Index(inHistName[inHistNum].c_str()) >= 0){
	      
	      TString name3 = name2(name2.Index(Form("_%s", ptEtaStr2[ptEtaNum].c_str()))+1, name2.Length()); 
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
    }
  }

  ptInts_p->push_back(0);
  ptDec_p->push_back(0);

  Int_t sortIter = 0; 

  while(sortIter < (Int_t)(ptBins_p->size() - 1)){
    TString ptBins1 = ptBins_p->at(sortIter).c_str();
    TString ptBins2 = ptBins_p->at(sortIter+1).c_str();

    TString pt1Part1 = ptBins1(ptBins1.Index(ptEtaStr2[ptEtaNum].c_str()) + ptEtaStr2[ptEtaNum].size(), TString(ptBins1(ptBins1.Index(ptEtaStr2[ptEtaNum].c_str())+ptEtaStr2[ptEtaNum].size(), ptBins1.Length() - ptEtaStr2[ptEtaNum].size())).Index("p"));
    TString pt1Part2 = ptBins1(ptBins1.Index("p")+1, 1);
    if(!strcmp(pt1Part2.Data(), "-")) pt1Part2 = ptBins1(ptBins1.Index("p")+1, 2);

    TString pt2Part1 = ptBins2(ptBins2.Index(ptEtaStr2[ptEtaNum].c_str()) + ptEtaStr2[ptEtaNum].size(), TString(ptBins2(ptBins2.Index(ptEtaStr2[ptEtaNum].c_str()) + ptEtaStr2[ptEtaNum].size(), ptBins2.Length() - ptEtaStr2[ptEtaNum].size())).Index("p"));
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

  const Int_t nTH1 = nTH1Temp;
  const Int_t nDir = nDirTemp;
  const Int_t nPtBins = nPtBinsTemp;
  TH1F* th1_p[nTH1];
  TCanvas* th1Canv_p[nDir][nCentBins];

  Int_t xBins = 0;
  Int_t yBins = 0;
  Float_t breakBool = false;

  for(Int_t iter = 0; iter < nPtBins; iter++){
    if(breakBool) break;
    for(Int_t iter2 = 0; iter2 < 4; iter2++){
      xBins = iter + 4 - iter2;
      yBins = iter;

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

  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      th1Canv_p[iter][centIter] = new TCanvas(Form("%s_MeanPts%s%s%s_c", dirNames_p->at(iter).c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str(), centStrings[centIter].c_str()), Form("%s_MeanPts%s%s%s_c", dirNames_p->at(iter).c_str(), inHistName[inHistNum].c_str(), ptEtaStr[ptEtaNum].c_str(), centStrings[centIter].c_str()), xBins*300, yBins*325);
      th1Canv_p[iter][centIter]->Divide(xBins, yBins, 0.0, 0.0);
    }
  }

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  Float_t dirMax[nDir][nCentBins];
  Float_t dirMin[nDir][nCentBins];

  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
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
      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    Int_t centPos = -1;
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      std::size_t pos = th1Name.find(centStrings[centIter]);
      if(pos != std::string::npos){
        centPos = centIter;
        break;
      }
    }



    for(Int_t binIter = 0; binIter < th1_p[iter]->GetNbinsX(); binIter++){
      if(th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1) > dirMax[dirPos][centPos]) dirMax[dirPos][centPos] = th1_p[iter]->GetBinContent(binIter+1) + th1_p[iter]->GetBinError(binIter+1);

      if(th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1) < dirMin[dirPos][centPos]) dirMin[dirPos][centPos] = th1_p[iter]->GetBinContent(binIter+1) - th1_p[iter]->GetBinError(binIter+1);
    }
  }


  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      dirMax[iter][centIter] = setMaxMinNice(dirMax[iter][centIter], true);
      dirMin[iter][centIter] = setMaxMinNice(dirMin[iter][centIter], false);
    }
  }

  for(Int_t iter = 0; iter < nTH1; iter++){
    std::string th1Name = th1Names_p->at(iter);

    Int_t dirPos = -1;
    for(Int_t dirIter = 0; dirIter < nDir; dirIter++){
      std::size_t pos = th1Name.find(dirNames_p->at(dirIter));
      if(pos != std::string::npos){
        dirPos = dirIter;
        break;
      }
    }

    Int_t centPos = -1;
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      std::size_t pos = th1Name.find(centStrings[centIter]);
      if(pos != std::string::npos){
        centPos = centIter;
        break;
      }
    }

    Int_t ptPos = -1;
    for(Int_t ptIter = 0; ptIter < (Int_t)ptBins_p->size(); ptIter++){
      std::size_t pos = th1Name.find(ptBins_p->at(ptIter).c_str());
      if(pos != std::string::npos){
        ptPos = ptIter;
        break;
      }
    }

    th1Canv_p[dirPos][centPos]->cd(ptPos+1);

    gStyle->SetOptStat(0);

    th1_p[iter]->SetMaximum(dirMax[dirPos][centPos]);
    th1_p[iter]->SetMinimum(dirMin[dirPos][centPos]);

    th1_p[iter]->GetXaxis()->CenterTitle();
    th1_p[iter]->GetXaxis()->SetTitleOffset(th1_p[iter]->GetXaxis()->GetTitleOffset() + 2.0);
    th1_p[iter]->GetXaxis()->SetTitleFont(43);
    th1_p[iter]->GetXaxis()->SetTitleSize(26);
    th1_p[iter]->GetYaxis()->CenterTitle();
    th1_p[iter]->GetYaxis()->SetTitleOffset(th1_p[iter]->GetYaxis()->GetTitleOffset() + 2.5);
    
    th1_p[iter]->GetXaxis()->SetTitle(xAxisLabel[inHistNum].c_str());
    th1_p[iter]->GetYaxis()->SetTitle("Events");

    th1_p[iter]->DrawCopy("E1 P");
    
    TLine* meanLine_p = new TLine(th1_p[iter]->GetMean(), th1_p[iter]->GetMinimum(), th1_p[iter]->GetMean(), th1_p[iter]->GetMaximum()*.8);
    meanLine_p->SetLineStyle(2);
    meanLine_p->DrawClone();
    delete meanLine_p;
    
    //    if(centPos == 0) label_p->DrawLatex(.30, .9, Form("#bf{%s}", dirNames_p->at(dirPos).c_str()));
    
    if(ptPos%xBins == 0) label_p->DrawLatex(.30, .88, Form("%d.%d<%s<%d.%d", ptInts_p->at(ptPos), ptDec_p->at(ptPos), ptEtaStr3[ptEtaNum].c_str(), ptInts_p->at(ptPos+1), ptDec_p->at(ptPos+1)));
    else label_p->DrawLatex(.20, .88, Form("%d.%d<%s<%d.%d", ptInts_p->at(ptPos), ptDec_p->at(ptPos), ptEtaStr3[ptEtaNum].c_str(), ptInts_p->at(ptPos+1), ptDec_p->at(ptPos+1)));

    if(ptPos == 0) label_p->DrawLatex(.30, .78, Form("#bf{#color[2]{%s (%s)}}", dirNames_p->at(dirPos).c_str(), centStrings2[centPos].c_str()));
    if(ptPos == 1) label_p->DrawLatex(.15, .78, Form("#bf{#color[2]{|#eta_{jet}|<2.0}}"));
    if(ptPos == 2 && !strcmp(ptEtaStr2[ptEtaNum].c_str(), "Eta" )) label_p->DrawLatex(.15, .78, Form("#bf{#color[2]{Gen. p_{T}>30}}"));

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
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      th1Canv_p[iter][centIter]->Write("", TObject::kOverwrite);
      claverCanvasSaving(th1Canv_p[iter][centIter], Form("pdfDir/%s", th1Canv_p[iter][centIter]->GetName()), "pdf");
    }
  }

  outFile_p->Close();
  delete outFile_p;

  delete label_p;

  delete ptBins_p;
  delete ptInts_p;
  delete ptDec_p;


  for(Int_t iter = 0; iter < nDir; iter++){
    for(Int_t centIter = 0; centIter < nCentBins; centIter++){
      delete th1Canv_p[iter][centIter];
    }
  }

  th1Names_p->clear();
  delete th1Names_p;

  dirNames_p->clear();
  delete dirNames_p;

  inFile_p->Close();
  delete inFile_p;

  return;
}


void makeJECPlotScatter(const std::string inFileName, const Int_t inHistNum, const Int_t ptEtaNum)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  const Int_t nContents = inFile_p->GetListOfKeys()->GetEntries();

  Int_t nTH2Temp = 0;
  std::vector<std::string>* th2Names_p = new std::vector<std::string>;

  Int_t nDirTemp = 0;
  std::vector<std::string>* dirNames_p = new std::vector<std::string>;

  for(Int_t iter = 0; iter < nContents; iter++){
    TString className = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetClassName();
    TString name = ((TKey*)inFile_p->GetListOfKeys()->At(iter))->GetName();

    if(className.Index("TDirectoryFile") >= 0){
      nDirTemp++;
      dirNames_p->push_back(name.Data());

      TDirectoryFile* tempDir_p = (TDirectoryFile*)inFile_p->Get(name);

      const Int_t nDirContents = tempDir_p->GetListOfKeys()->GetEntries();
      for(Int_t dirIter = 0; dirIter < nDirContents; dirIter++){
        TString name2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetName();
        TString className2 = ((TKey*)tempDir_p->GetListOfKeys()->At(dirIter))->GetClassName();

        if(className2.Index("TH2") >= 0){
	  if(name2.Index(ptEtaStr[ptEtaNum].c_str()) < 0) continue;
          if(name2.Index(inHistName[inHistNum].c_str()) >= 0){
            nTH2Temp++;
            th2Names_p->push_back(Form("%s/%s", name.Data(), name2.Data()));
          }
        }
      }
    }    
  }
  
  const Int_t nTH2 = nTH2Temp;
  TH2F* th2_p[nTH2];
  TCanvas* th2Canv_p[nTH2];

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSizePixels(28);

  for(Int_t iter = 0; iter < nTH2; iter++){
    std::string th2Name = th2Names_p->at(iter);

    th2_p[iter] = (TH2F*)inFile_p->Get(th2Name.c_str());

    std::size_t pos = th2Name.find("_h");
    if(pos != std::string::npos) th2Name.replace(pos, pos+2, "_c");

    pos = th2Name.find("/");
    if(pos != std::string::npos) th2Name.replace(0, pos+1, "");

    th2Canv_p[iter] = new TCanvas(th2Name.c_str(), th2Name.c_str(), 700, 700);

    gPad->SetLeftMargin(gPad->GetLeftMargin() - .06);
    gPad->SetBottomMargin(gPad->GetBottomMargin() - .06);
    gPad->SetRightMargin(gPad->GetRightMargin() + .11);

    if(!strcmp(ptEtaStr2[ptEtaNum].c_str(), "Pt")) gPad->SetLogx();
    if(strcmp(ptEtaStr2[ptEtaNum].c_str(), "PtEta") != 0) gPad->SetLogz();
    if(!strcmp(ptEtaStr2[ptEtaNum].c_str(), "PtEta")) gPad->SetLogy();

    th2_p[iter]->GetXaxis()->CenterTitle();
    th2_p[iter]->GetXaxis()->SetTitleOffset(th2_p[iter]->GetXaxis()->GetTitleOffset() - .20);
    th2_p[iter]->GetYaxis()->CenterTitle();
    th2_p[iter]->GetYaxis()->SetTitleOffset(th2_p[iter]->GetYaxis()->GetTitleOffset() - .20);

    th2_p[iter]->DrawCopy("COLZ");

    pos = th2Name.find("cent");
    if(pos != std::string::npos){
      std::string centStr = th2Name.substr(pos);
      pos = centStr.find("_");
      if(pos != std::string::npos) centStr = centStr.substr(4, pos-4);

      pos = centStr.find("to");

      if(pos != std::string::npos){
	if(strcmp(ptEtaStr2[ptEtaNum].c_str(), "PtEta") != 0) label_p->DrawLatex(.65, .85, Form("#bf{#color[2]{%s-%s%%}}", centStr.substr(0, pos).c_str(), centStr.substr(pos+2).c_str()));
	else label_p->DrawLatex(.65, .85, Form("#bf{%s-%s%%}", centStr.substr(0, pos).c_str(), centStr.substr(pos+2).c_str()));
      }
    }

    pos = th2Name.find("ak");
    if(pos != std::string::npos){
      std::string centStr = th2Name.substr(pos);
      pos = centStr.find("_");
      if(pos != std::string::npos) centStr = centStr.substr(0, pos);
      if(pos != std::string::npos) label_p->DrawLatex(.65, .92, Form("#bf{#color[2]{%s}}", centStr.c_str()));
    }

    label_p->DrawLatex(.65, .78, "#bf{#color[2]{|#eta_{jet}|<2.0}}");
    if(!strcmp(ptEtaStr2[ptEtaNum].c_str(), "Eta" )) label_p->DrawLatex(.65, .71, Form("#bf{#color[2]{Gen. p_{T}>30}}"));
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

  for(Int_t iter = 0; iter < nTH2; iter++){
    th2Canv_p[iter]->Write("", TObject::kOverwrite);
    claverCanvasSaving(th2Canv_p[iter], Form("pdfDir/%s", th2Canv_p[iter]->GetName()), "pdf");
  }

  outFile_p->Close();
  delete outFile_p;

  delete label_p;

  for(Int_t iter = 0; iter < nTH2; iter++){
    delete th2Canv_p[iter];
  }

  dirNames_p->clear();
  delete dirNames_p;

  th2Names_p->clear();
  delete th2Names_p;

  inFile_p->Close();
  delete inFile_p;

  return;
}




void makeJECPlot(const std::string inFileName)
{
  for(Int_t iter = 0; iter < nHistName; iter++){
    for(Int_t ptEtaIter = 0; ptEtaIter < nPtEta; ptEtaIter++){
      for(Int_t iter2 = 0; iter2 < nMeanRes; iter2++){
	if(!strcmp(inHistName[iter].c_str(), "Eff") && !strcmp(meanResStr[iter2].c_str(), "Res")) continue;

	if(ptEtaIter < 2) makeJECPlotMeanRes(inFileName, iter, iter2, ptEtaIter);
      }

      if(ptEtaIter < 2) makeJECPlotMeanPts(inFileName, iter, ptEtaIter);
      makeJECPlotScatter(inFileName, iter, ptEtaIter);
    }
  }

  return;
}
