#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TDatime.h"
#include "TLorentzVector.h"
#include "TNamed.h"

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


int makeJECTree_Prototype(const std::string inConfigFileName)
{
  jecConfigParser config;
  if(!config.SetConfigParser(inConfigFileName)) return 1;

  if(config.GetDoCorrections()) initPtEtaJetResidualCorr(config.GetCorrFileName(), config.GetCorrForm());
  //if(config.GetDoCorrections()) initGetResidualJetCorr(config.GetCorrFileName());

  std::string outName = config.GetOutName();

  if(outName.find(".root") != std::string::npos) outName.replace(outName.find(".root"), 5, "_TREE.root");
  else outName = outName + "_TREE.root";

  TFile* outFile_p = new TFile(outName.c_str(), "RECREATE");
  
  std::vector<std::string> jetAlgoInFile = config.GetJetTypesFinal();

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


  std::ofstream outFileRunLumiEvt("runLumiEvtFileForXCheck.txt");
  outFileRunLumiEvt << "Run,Lumi,Evt" << std::endl;
  outFileRunLumiEvt.close();
    

  const Int_t nJetAlgo = (Int_t)jetAlgoInFile.size();
  std::cout << nJetAlgo << std::endl;
  std::vector<std::string> jetAlgo;

  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    jetAlgo.push_back(jetAlgoInFile.at(iter).substr(0, jetAlgoInFile.at(iter).find("JetAnalyzer")));
  }
  
  for(Int_t iter = 0; iter < nJetAlgo; iter++){
    std::cout << "Algos: " << jetAlgo.at(iter) << std::endl;
  }

  outFile_p->cd();

  Int_t RunOut_[nJetAlgo];
  Int_t LumiOut_[nJetAlgo];
  ULong64_t EventOut_[nJetAlgo];

  Double_t leptPt1Out_[nJetAlgo];
  Double_t leptEta1Out_[nJetAlgo];
  Double_t leptPhi1Out_[nJetAlgo];
  Double_t leptPt2Out_[nJetAlgo];
  Double_t leptEta2Out_[nJetAlgo];
  Double_t leptPhi2Out_[nJetAlgo];
  Double_t zPtOut_[nJetAlgo];
  Double_t zEtaOut_[nJetAlgo];
  Double_t zPhiOut_[nJetAlgo];
  Double_t zMOut_[nJetAlgo];

  Double_t phoPtOut_[nJetAlgo];
  Double_t phoEtaOut_[nJetAlgo];
  Double_t phoPhiOut_[nJetAlgo];

  Double_t pho2PtOut_[nJetAlgo];
  Double_t pho2EtaOut_[nJetAlgo];
  Double_t pho2PhiOut_[nJetAlgo];

  Double_t genPhoPtOut_[nJetAlgo];
  Double_t genPhoEtaOut_[nJetAlgo];
  Double_t genPhoPhiOut_[nJetAlgo];

  Int_t hiBinOut_[nJetAlgo];
  Double_t hiBinWeightOut_[nJetAlgo];
  Double_t hiBinWeightNormOut_[nJetAlgo];
  Double_t ptHatOut_[nJetAlgo];
  Double_t ptHatWeightOut_[nJetAlgo];
  Double_t fullWeightOut_[nJetAlgo];
  Double_t jtPtOut_[nJetAlgo];
  Double_t rawPtOut_[nJetAlgo];
  Double_t jtPtOrigOut_[nJetAlgo];
  Double_t jtEtaOut_[nJetAlgo];
  Double_t jtPhiOut_[nJetAlgo];
  Double_t jtDelZPhiOut_[nJetAlgo];
  Double_t jtDelPhoPhiOut_[nJetAlgo];
  Double_t jtDelPho2PhiOut_[nJetAlgo];
  Double_t jtOverRefPtOut_[nJetAlgo];
  Double_t jtOrigOverRefPtOut_[nJetAlgo];
  Double_t jtDelRefPhiOut_[nJetAlgo];
  Double_t jtDelRefEtaOut_[nJetAlgo];
  Double_t jtChargedSumOut_[nJetAlgo];
  Double_t jtChargedMaxOut_[nJetAlgo];
  Double_t jtTrackSumOut_[nJetAlgo];
  Double_t jtTrackMaxOut_[nJetAlgo];
  Double_t jtPhotonSumOut_[nJetAlgo];
  Double_t jtNeutralSumOut_[nJetAlgo];
  Double_t jtMuSumOut_[nJetAlgo];
  Double_t jtESumOut_[nJetAlgo];
  Double_t jtMaxTrackInConePtOut_[nJetAlgo];
  Bool_t jtMaxTrackInConeHPOut_[nJetAlgo];
  Double_t refPtOut_[nJetAlgo];
  Double_t refEtaOut_[nJetAlgo];
  Double_t refPhiOut_[nJetAlgo];
  Int_t refSubIDOut_[nJetAlgo];
  Int_t refPartFlavOut_[nJetAlgo];

  Int_t RunOutBad_[nJetAlgo];
  Int_t LumiOutBad_[nJetAlgo];
  ULong64_t EventOutBad_[nJetAlgo];

  Double_t leptPt1OutBad_[nJetAlgo];
  Double_t leptEta1OutBad_[nJetAlgo];
  Double_t leptPhi1OutBad_[nJetAlgo];
  Double_t leptPt2OutBad_[nJetAlgo];
  Double_t leptEta2OutBad_[nJetAlgo];
  Double_t leptPhi2OutBad_[nJetAlgo];
  Double_t zPtOutBad_[nJetAlgo];
  Double_t zEtaOutBad_[nJetAlgo];
  Double_t zPhiOutBad_[nJetAlgo];
  Double_t zMOutBad_[nJetAlgo];

  Double_t phoPtOutBad_[nJetAlgo];
  Double_t phoEtaOutBad_[nJetAlgo];
  Double_t phoPhiOutBad_[nJetAlgo];

  Double_t pho2PtOutBad_[nJetAlgo];
  Double_t pho2EtaOutBad_[nJetAlgo];
  Double_t pho2PhiOutBad_[nJetAlgo];

  Double_t genPhoPtOutBad_[nJetAlgo];
  Double_t genPhoEtaOutBad_[nJetAlgo];
  Double_t genPhoPhiOutBad_[nJetAlgo];

  Int_t hiBinOutBad_[nJetAlgo];
  Double_t hiBinWeightOutBad_[nJetAlgo];
  Double_t hiBinWeightNormOutBad_[nJetAlgo];
  Double_t ptHatOutBad_[nJetAlgo];
  Double_t ptHatWeightOutBad_[nJetAlgo];
  Double_t fullWeightOutBad_[nJetAlgo];
  Double_t jtPtOutBad_[nJetAlgo];
  Double_t rawPtOutBad_[nJetAlgo];
  Double_t jtPtOrigOutBad_[nJetAlgo];
  Double_t jtEtaOutBad_[nJetAlgo];
  Double_t jtPhiOutBad_[nJetAlgo];
  Double_t jtDelZPhiOutBad_[nJetAlgo];
  Double_t jtDelPhoPhiOutBad_[nJetAlgo];
  Double_t jtDelPho2PhiOutBad_[nJetAlgo];
  Double_t jtOverRefPtOutBad_[nJetAlgo];
  Double_t jtOrigOverRefPtOutBad_[nJetAlgo];
  Double_t jtDelRefPhiOutBad_[nJetAlgo];
  Double_t jtDelRefEtaOutBad_[nJetAlgo];
  Double_t jtChargedSumOutBad_[nJetAlgo];
  Double_t jtTrackSumOutBad_[nJetAlgo];
  Double_t jtTrackMaxOutBad_[nJetAlgo];
  Double_t jtPhotonSumOutBad_[nJetAlgo];
  Double_t jtNeutralSumOutBad_[nJetAlgo];
  Double_t jtMuSumOutBad_[nJetAlgo];
  Double_t jtESumOutBad_[nJetAlgo];
  Double_t refPtOutBad_[nJetAlgo];
  Double_t refEtaOutBad_[nJetAlgo];
  Double_t refPhiOutBad_[nJetAlgo];
  Int_t refSubIDOutBad_[nJetAlgo];
  Int_t refPartFlavOutBad_[nJetAlgo];

  Int_t RunOutGen_[nJetAlgo];
  Int_t LumiOutGen_[nJetAlgo];
  ULong64_t EventOutGen_[nJetAlgo];

  Int_t hiBinOutGen_[nJetAlgo];
  Double_t hiBinWeightOutGen_[nJetAlgo];
  Double_t hiBinWeightNormOutGen_[nJetAlgo];
  Double_t ptHatOutGen_[nJetAlgo];
  Double_t ptHatWeightOutGen_[nJetAlgo];
  Double_t fullWeightOutGen_[nJetAlgo];  
  Double_t genJtPtOut_[nJetAlgo];
  Double_t genJtEtaOut_[nJetAlgo];
  Double_t genJtPhiOut_[nJetAlgo];
  Double_t genJtDelZPhiOut_[nJetAlgo];
  Int_t genJtMatchIndexOut_[nJetAlgo];
  Double_t genJtMatchRecoPtOut_[nJetAlgo];
  Double_t genJtMatchRecoPhiOut_[nJetAlgo];
  Double_t genJtMatchRecoEtaOut_[nJetAlgo];
  Int_t genJtManualMatchIndexOut_[nJetAlgo];
  Double_t genJtManualMatchRecoPtOut_[nJetAlgo];
  Double_t genJtManualMatchRecoPhiOut_[nJetAlgo];
  Double_t genJtManualMatchRecoEtaOut_[nJetAlgo];
  Double_t genJtMaxPartPtOut_[nJetAlgo];
  Double_t genJtMaxPartPhiOut_[nJetAlgo];
  Double_t genJtMaxPartEtaOut_[nJetAlgo];
  Int_t genJtMaxPartPDGIDOut_[nJetAlgo];
  Double_t genJtSumPartPtOut_[nJetAlgo];
  Int_t genJtNPartPtOut_[nJetAlgo];
  Int_t genJtSubIDOut_[nJetAlgo];
  Int_t genJtPartFlavOut_[nJetAlgo];
  Int_t genJtManualMatchPartFlavOut_[nJetAlgo];


  TTree* outTree_p[nJetAlgo];
  TTree* outTreeBad_p[nJetAlgo];
  TTree* outTreeGen_p[nJetAlgo];

  for(int iter = 0; iter < nJetAlgo; iter++){
    std::string treeName = "jecTree_" + jetAlgo.at(iter);
    outTree_p[iter] = new TTree(treeName.c_str(), treeName.c_str());

    std::string leptPt1Str1 = "leptPt1";
    std::string leptPt1Str2 = "leptPt1/D";
    std::string leptPhi1Str1 = "leptPhi1";
    std::string leptPhi1Str2 = "leptPhi1/D";
    std::string leptEta1Str1 = "leptEta1";
    std::string leptEta1Str2 = "leptEta1/D";

    std::string leptPt2Str1 = "leptPt2";
    std::string leptPt2Str2 = "leptPt2/D";
    std::string leptPhi2Str1 = "leptPhi2";
    std::string leptPhi2Str2 = "leptPhi2/D";
    std::string leptEta2Str1 = "leptEta2";
    std::string leptEta2Str2 = "leptEta2/D";

    std::string zPtStr1 = "zPt";
    std::string zPtStr2 = "zPt/D";
    std::string zPhiStr1 = "zPhi";
    std::string zPhiStr2 = "zPhi/D";
    std::string zEtaStr1 = "zEta";
    std::string zEtaStr2 = "zEta/D";
    std::string zMStr1 = "zM";
    std::string zMStr2 = "zM/D";

    std::string phoPtStr1 = "phoPt";
    std::string phoPtStr2 = "phoPt/D";
    std::string phoPhiStr1 = "phoPhi";
    std::string phoPhiStr2 = "phoPhi/D";
    std::string phoEtaStr1 = "phoEta";
    std::string phoEtaStr2 = "phoEta/D";

    std::string pho2PtStr1 = "pho2Pt";
    std::string pho2PtStr2 = "pho2Pt/D";
    std::string pho2PhiStr1 = "pho2Phi";
    std::string pho2PhiStr2 = "pho2Phi/D";
    std::string pho2EtaStr1 = "pho2Eta";
    std::string pho2EtaStr2 = "pho2Eta/D";


    std::string genPhoPtStr1 = "genPhoPt";
    std::string genPhoPtStr2 = "genPhoPt/D";
    std::string genPhoPhiStr1 = "genPhoPhi";
    std::string genPhoPhiStr2 = "genPhoPhi/D";
    std::string genPhoEtaStr1 = "genPhoEta";
    std::string genPhoEtaStr2 = "genPhoEta/D";

    outTree_p[iter]->Branch(leptPt1Str1.c_str(), &leptPt1Out_[iter], leptPt1Str2.c_str());
    outTree_p[iter]->Branch(leptPhi1Str1.c_str(), &leptPhi1Out_[iter], leptPhi1Str2.c_str());
    outTree_p[iter]->Branch(leptEta1Str1.c_str(), &leptEta1Out_[iter], leptEta1Str2.c_str());

    outTree_p[iter]->Branch(leptPt2Str1.c_str(), &leptPt2Out_[iter], leptPt2Str2.c_str());
    outTree_p[iter]->Branch(leptPhi2Str1.c_str(), &leptPhi2Out_[iter], leptPhi2Str2.c_str());
    outTree_p[iter]->Branch(leptEta2Str1.c_str(), &leptEta2Out_[iter], leptEta2Str2.c_str());

    outTree_p[iter]->Branch(zPtStr1.c_str(), &zPtOut_[iter], zPtStr2.c_str());
    outTree_p[iter]->Branch(zPhiStr1.c_str(), &zPhiOut_[iter], zPhiStr2.c_str());
    outTree_p[iter]->Branch(zEtaStr1.c_str(), &zEtaOut_[iter], zEtaStr2.c_str());
    outTree_p[iter]->Branch(zMStr1.c_str(), &zMOut_[iter], zMStr2.c_str());

    outTree_p[iter]->Branch(phoPtStr1.c_str(), &phoPtOut_[iter], phoPtStr2.c_str());
    outTree_p[iter]->Branch(phoPhiStr1.c_str(), &phoPhiOut_[iter], phoPhiStr2.c_str());
    outTree_p[iter]->Branch(phoEtaStr1.c_str(), &phoEtaOut_[iter], phoEtaStr2.c_str());

    outTree_p[iter]->Branch(pho2PtStr1.c_str(), &pho2PtOut_[iter], pho2PtStr2.c_str());
    outTree_p[iter]->Branch(pho2PhiStr1.c_str(), &pho2PhiOut_[iter], pho2PhiStr2.c_str());
    outTree_p[iter]->Branch(pho2EtaStr1.c_str(), &pho2EtaOut_[iter], pho2EtaStr2.c_str());

    outTree_p[iter]->Branch(genPhoPtStr1.c_str(), &genPhoPtOut_[iter], genPhoPtStr2.c_str());
    outTree_p[iter]->Branch(genPhoPhiStr1.c_str(), &genPhoPhiOut_[iter], genPhoPhiStr2.c_str());
    outTree_p[iter]->Branch(genPhoEtaStr1.c_str(), &genPhoEtaOut_[iter], genPhoEtaStr2.c_str());

    std::string runStr1 = "Run";
    std::string runStr2 = "Run/I";
    std::string lumiStr1 = "Lumi";
    std::string lumiStr2 = "Lumi/I";
    std::string eventStr1 = "Event";
    std::string eventStr2 = "Event/l";

    outTree_p[iter]->Branch(runStr1.c_str(), &RunOut_[iter], runStr2.c_str());
    outTree_p[iter]->Branch(lumiStr1.c_str(), &LumiOut_[iter], lumiStr2.c_str());
    outTree_p[iter]->Branch(eventStr1.c_str(), &EventOut_[iter], eventStr2.c_str());

    std::string hiBinStr1 = "hiBin";
    std::string hiBinStr2 = "hiBin/I";
    std::string hiBinWeightStr1 = "hiBinWeight";
    std::string hiBinWeightStr2 = "hiBinWeight/D";
    std::string hiBinWeightNormStr1 = "hiBinWeightNorm";
    std::string hiBinWeightNormStr2 = "hiBinWeightNorm/D";
    std::string ptHatStr1 = "ptHat";
    std::string ptHatStr2 = "ptHat/D";
    std::string ptHatWeightStr1 = "ptHatWeight";
    std::string ptHatWeightStr2 = "ptHatWeight/D";
    std::string fullWeightStr1 = "fullWeight";
    std::string fullWeightStr2 = "fullWeight/D";

    outTree_p[iter]->Branch(hiBinStr1.c_str(), &hiBinOut_[iter], hiBinStr2.c_str());
    outTree_p[iter]->Branch(hiBinWeightStr1.c_str(), &hiBinWeightOut_[iter], hiBinWeightStr2.c_str());
    outTree_p[iter]->Branch(hiBinWeightNormStr1.c_str(), &hiBinWeightNormOut_[iter], hiBinWeightNormStr2.c_str());
    outTree_p[iter]->Branch(ptHatStr1.c_str(), &ptHatOut_[iter], ptHatStr2.c_str());
    outTree_p[iter]->Branch(ptHatWeightStr1.c_str(), &ptHatWeightOut_[iter], ptHatWeightStr2.c_str());
    outTree_p[iter]->Branch(fullWeightStr1.c_str(), &fullWeightOut_[iter], fullWeightStr2.c_str());

    std::string jtPtStr1 = "jtPt";
    std::string jtPtStr2 = "jtPt/D";
    std::string rawPtStr1 = "rawPt";
    std::string rawPtStr2 = "rawPt/D";
    std::string jtPtOrigStr1 = "jtPtOrig";
    std::string jtPtOrigStr2 = "jtPtOrig/D";
    std::string jtEtaStr1 = "jtEta";
    std::string jtEtaStr2 = "jtEta/D";
    std::string jtPhiStr1 = "jtPhi";
    std::string jtPhiStr2 = "jtPhi/D";
    std::string jtDelZPhiStr1 = "jtDelZPhi";
    std::string jtDelZPhiStr2 = "jtDelZPhi/D";
    std::string jtDelPhoPhiStr1 = "jtDelPhoPhi";
    std::string jtDelPhoPhiStr2 = "jtDelPhoPhi/D";
    std::string jtDelPho2PhiStr1 = "jtDelPho2Phi";
    std::string jtDelPho2PhiStr2 = "jtDelPho2Phi/D";

    std::string jtOverRefPtStr1 = "jtOverRefPt";
    std::string jtOverRefPtStr2 = "jtOverRefPt/D";
    std::string jtOrigOverRefPtStr1 = "jtOrigOverRefPt";
    std::string jtOrigOverRefPtStr2 = "jtOrigOverRefPt/D";

    std::string jtDelRefPhiStr1 = "jtDelRefPhi";
    std::string jtDelRefPhiStr2 = "jtDelRefPhi/D";
    std::string jtDelRefEtaStr1 = "jtDelRefEta";
    std::string jtDelRefEtaStr2 = "jtDelRefEta/D";
    std::string jtChargedSumStr1 = "jtChargedSum";
    std::string jtChargedSumStr2 = "jtChargedSum/D";
    std::string jtChargedMaxStr1 = "jtChargedMax";
    std::string jtChargedMaxStr2 = "jtChargedMax/D";
    std::string jtTrackSumStr1 = "jtTrackSum";
    std::string jtTrackSumStr2 = "jtTrackSum/D";
    std::string jtTrackMaxStr1 = "jtTrackMax";
    std::string jtTrackMaxStr2 = "jtTrackMax/D";
    std::string jtPhotonSumStr1 = "jtPhotonSum";
    std::string jtPhotonSumStr2 = "jtPhotonSum/D";
    std::string jtNeutralSumStr1 = "jtNeutralSum";
    std::string jtNeutralSumStr2 = "jtNeutralSum/D";
    std::string jtMuSumStr1 = "jtMuSum";
    std::string jtMuSumStr2 = "jtMuSum/D";
    std::string jtESumStr1 = "jtESum";
    std::string jtESumStr2 = "jtESum/D";
    std::string jtMaxTrackInConePtStr1 = "jtMaxTrackInConePt";
    std::string jtMaxTrackInConePtStr2 = "jtMaxTrackInConePt/D";
    std::string jtMaxTrackInConeHPStr1 = "jtMaxTrackInConeHP";
    std::string jtMaxTrackInConeHPStr2 = "jtMaxTrackInConeHP/O";

    outTree_p[iter]->Branch(jtPtStr1.c_str(), &jtPtOut_[iter], jtPtStr2.c_str());
    outTree_p[iter]->Branch(rawPtStr1.c_str(), &rawPtOut_[iter], rawPtStr2.c_str());
    outTree_p[iter]->Branch(jtPtOrigStr1.c_str(), &jtPtOrigOut_[iter], jtPtOrigStr2.c_str());
    outTree_p[iter]->Branch(jtEtaStr1.c_str(), &jtEtaOut_[iter], jtEtaStr2.c_str());
    outTree_p[iter]->Branch(jtPhiStr1.c_str(), &jtPhiOut_[iter], jtPhiStr2.c_str());
    outTree_p[iter]->Branch(jtDelZPhiStr1.c_str(), &jtDelZPhiOut_[iter], jtDelZPhiStr2.c_str());
    outTree_p[iter]->Branch(jtDelPhoPhiStr1.c_str(), &jtDelPhoPhiOut_[iter], jtDelPhoPhiStr2.c_str());
    outTree_p[iter]->Branch(jtDelPho2PhiStr1.c_str(), &jtDelPho2PhiOut_[iter], jtDelPho2PhiStr2.c_str());

    outTree_p[iter]->Branch(jtOverRefPtStr1.c_str(), &jtOverRefPtOut_[iter], jtOverRefPtStr2.c_str());
    outTree_p[iter]->Branch(jtOrigOverRefPtStr1.c_str(), &jtOrigOverRefPtOut_[iter], jtOrigOverRefPtStr2.c_str());
    outTree_p[iter]->Branch(jtDelRefPhiStr1.c_str(), &jtDelRefPhiOut_[iter], jtDelRefPhiStr2.c_str());
    outTree_p[iter]->Branch(jtDelRefEtaStr1.c_str(), &jtDelRefEtaOut_[iter], jtDelRefEtaStr2.c_str());
    outTree_p[iter]->Branch(jtChargedSumStr1.c_str(), &jtChargedSumOut_[iter], jtChargedSumStr2.c_str());
    outTree_p[iter]->Branch(jtChargedMaxStr1.c_str(), &jtChargedMaxOut_[iter], jtChargedMaxStr2.c_str());
    outTree_p[iter]->Branch(jtTrackSumStr1.c_str(), &jtTrackSumOut_[iter], jtTrackSumStr2.c_str());
    outTree_p[iter]->Branch(jtTrackMaxStr1.c_str(), &jtTrackMaxOut_[iter], jtTrackMaxStr2.c_str());
    outTree_p[iter]->Branch(jtPhotonSumStr1.c_str(), &jtPhotonSumOut_[iter], jtPhotonSumStr2.c_str());
    outTree_p[iter]->Branch(jtNeutralSumStr1.c_str(), &jtNeutralSumOut_[iter], jtNeutralSumStr2.c_str());
    outTree_p[iter]->Branch(jtMuSumStr1.c_str(), &jtMuSumOut_[iter], jtMuSumStr2.c_str());
    outTree_p[iter]->Branch(jtESumStr1.c_str(), &jtESumOut_[iter], jtESumStr2.c_str());
    outTree_p[iter]->Branch(jtMaxTrackInConePtStr1.c_str(), &jtMaxTrackInConePtOut_[iter], jtMaxTrackInConePtStr2.c_str());
    outTree_p[iter]->Branch(jtMaxTrackInConeHPStr1.c_str(), &jtMaxTrackInConeHPOut_[iter], jtMaxTrackInConeHPStr2.c_str());

    std::string refPtStr1 = "refPt";
    std::string refPtStr2 = "refPt/D";
    std::string refEtaStr1 = "refEta";
    std::string refEtaStr2 = "refEta/D";
    std::string refPhiStr1 = "refPhi";
    std::string refPhiStr2 = "refPhi/D";
    std::string refSubIDStr1 = "refSubID";
    std::string refSubIDStr2 = "refSubID/I";
    std::string refPartFlavStr1 = "refPartFlav";
    std::string refPartFlavStr2 = "refPartFlav/I";

    outTree_p[iter]->Branch(refPtStr1.c_str(), &refPtOut_[iter], refPtStr2.c_str());
    outTree_p[iter]->Branch(refEtaStr1.c_str(), &refEtaOut_[iter], refEtaStr2.c_str());
    outTree_p[iter]->Branch(refPhiStr1.c_str(), &refPhiOut_[iter], refPhiStr2.c_str());
    outTree_p[iter]->Branch(refSubIDStr1.c_str(), &refSubIDOut_[iter], refSubIDStr2.c_str());
    outTree_p[iter]->Branch(refPartFlavStr1.c_str(), &refPartFlavOut_[iter], refPartFlavStr2.c_str());

    std::string treeNameBad = "jecTree_Bad_" + jetAlgo.at(iter);
    outTreeBad_p[iter] = new TTree(treeNameBad.c_str(), treeNameBad.c_str());

    outTreeBad_p[iter]->Branch(leptPt1Str1.c_str(), &leptPt1OutBad_[iter], leptPt1Str2.c_str());
    outTreeBad_p[iter]->Branch(leptPhi1Str1.c_str(), &leptPhi1OutBad_[iter], leptPhi1Str2.c_str());
    outTreeBad_p[iter]->Branch(leptEta1Str1.c_str(), &leptEta1OutBad_[iter], leptEta1Str2.c_str());

    outTreeBad_p[iter]->Branch(leptPt2Str1.c_str(), &leptPt2OutBad_[iter], leptPt2Str2.c_str());
    outTreeBad_p[iter]->Branch(leptPhi2Str1.c_str(), &leptPhi2OutBad_[iter], leptPhi2Str2.c_str());
    outTreeBad_p[iter]->Branch(leptEta2Str1.c_str(), &leptEta2OutBad_[iter], leptEta2Str2.c_str());

    outTreeBad_p[iter]->Branch(zPtStr1.c_str(), &zPtOutBad_[iter], zPtStr2.c_str());
    outTreeBad_p[iter]->Branch(zPhiStr1.c_str(), &zPhiOutBad_[iter], zPhiStr2.c_str());
    outTreeBad_p[iter]->Branch(zEtaStr1.c_str(), &zEtaOutBad_[iter], zEtaStr2.c_str());
    outTreeBad_p[iter]->Branch(zMStr1.c_str(), &zMOutBad_[iter], zMStr2.c_str());

    outTreeBad_p[iter]->Branch(phoPtStr1.c_str(), &phoPtOutBad_[iter], phoPtStr2.c_str());
    outTreeBad_p[iter]->Branch(phoPhiStr1.c_str(), &phoPhiOutBad_[iter], phoPhiStr2.c_str());
    outTreeBad_p[iter]->Branch(phoEtaStr1.c_str(), &phoEtaOutBad_[iter], phoEtaStr2.c_str());

    outTreeBad_p[iter]->Branch(pho2PtStr1.c_str(), &pho2PtOutBad_[iter], pho2PtStr2.c_str());
    outTreeBad_p[iter]->Branch(pho2PhiStr1.c_str(), &pho2PhiOutBad_[iter], pho2PhiStr2.c_str());
    outTreeBad_p[iter]->Branch(pho2EtaStr1.c_str(), &pho2EtaOutBad_[iter], pho2EtaStr2.c_str());

    outTreeBad_p[iter]->Branch(genPhoPtStr1.c_str(), &genPhoPtOutBad_[iter], genPhoPtStr2.c_str());
    outTreeBad_p[iter]->Branch(genPhoPhiStr1.c_str(), &genPhoPhiOutBad_[iter], genPhoPhiStr2.c_str());
    outTreeBad_p[iter]->Branch(genPhoEtaStr1.c_str(), &genPhoEtaOutBad_[iter], genPhoEtaStr2.c_str());

    outTree_p[iter]->Branch(runStr1.c_str(), &RunOutBad_[iter], runStr2.c_str());
    outTree_p[iter]->Branch(lumiStr1.c_str(), &LumiOutBad_[iter], lumiStr2.c_str());
    outTree_p[iter]->Branch(eventStr1.c_str(), &EventOutBad_[iter], eventStr2.c_str());


    outTreeBad_p[iter]->Branch(hiBinStr1.c_str(), &hiBinOutBad_[iter], hiBinStr2.c_str());
    outTreeBad_p[iter]->Branch(hiBinWeightStr1.c_str(), &hiBinWeightOutBad_[iter], hiBinWeightStr2.c_str());
    outTreeBad_p[iter]->Branch(hiBinWeightNormStr1.c_str(), &hiBinWeightNormOutBad_[iter], hiBinWeightNormStr2.c_str());
    outTreeBad_p[iter]->Branch(ptHatStr1.c_str(), &ptHatOutBad_[iter], ptHatStr2.c_str());
    outTreeBad_p[iter]->Branch(ptHatWeightStr1.c_str(), &ptHatWeightOutBad_[iter], ptHatWeightStr2.c_str());
    outTreeBad_p[iter]->Branch(fullWeightStr1.c_str(), &fullWeightOutBad_[iter], fullWeightStr2.c_str());

    outTreeBad_p[iter]->Branch(jtPtStr1.c_str(), &jtPtOutBad_[iter], jtPtStr2.c_str());
    outTreeBad_p[iter]->Branch(rawPtStr1.c_str(), &rawPtOutBad_[iter], rawPtStr2.c_str());
    outTreeBad_p[iter]->Branch(jtPtOrigStr1.c_str(), &jtPtOrigOutBad_[iter], jtPtOrigStr2.c_str());
    outTreeBad_p[iter]->Branch(jtEtaStr1.c_str(), &jtEtaOutBad_[iter], jtEtaStr2.c_str());
    outTreeBad_p[iter]->Branch(jtPhiStr1.c_str(), &jtPhiOutBad_[iter], jtPhiStr2.c_str());
    outTreeBad_p[iter]->Branch(jtDelZPhiStr1.c_str(), &jtDelZPhiOutBad_[iter], jtDelZPhiStr2.c_str());
    outTreeBad_p[iter]->Branch(jtDelPhoPhiStr1.c_str(), &jtDelPhoPhiOutBad_[iter], jtDelPhoPhiStr2.c_str());
    outTreeBad_p[iter]->Branch(jtDelPho2PhiStr1.c_str(), &jtDelPho2PhiOutBad_[iter], jtDelPho2PhiStr2.c_str());

    outTreeBad_p[iter]->Branch(jtOverRefPtStr1.c_str(), &jtOverRefPtOutBad_[iter], jtOverRefPtStr2.c_str());
    outTreeBad_p[iter]->Branch(jtOrigOverRefPtStr1.c_str(), &jtOrigOverRefPtOutBad_[iter], jtOrigOverRefPtStr2.c_str());
    outTreeBad_p[iter]->Branch(jtDelRefPhiStr1.c_str(), &jtDelRefPhiOutBad_[iter], jtDelRefPhiStr2.c_str());
    outTreeBad_p[iter]->Branch(jtDelRefEtaStr1.c_str(), &jtDelRefEtaOutBad_[iter], jtDelRefEtaStr2.c_str());
    outTreeBad_p[iter]->Branch(jtChargedSumStr1.c_str(), &jtChargedSumOutBad_[iter], jtChargedSumStr2.c_str());
    outTreeBad_p[iter]->Branch(jtTrackSumStr1.c_str(), &jtTrackSumOutBad_[iter], jtTrackSumStr2.c_str());
    outTreeBad_p[iter]->Branch(jtTrackMaxStr1.c_str(), &jtTrackMaxOutBad_[iter], jtTrackMaxStr2.c_str());
    outTreeBad_p[iter]->Branch(jtPhotonSumStr1.c_str(), &jtPhotonSumOutBad_[iter], jtPhotonSumStr2.c_str());
    outTreeBad_p[iter]->Branch(jtNeutralSumStr1.c_str(), &jtNeutralSumOutBad_[iter], jtNeutralSumStr2.c_str());
    outTreeBad_p[iter]->Branch(jtMuSumStr1.c_str(), &jtMuSumOutBad_[iter], jtMuSumStr2.c_str());
    outTreeBad_p[iter]->Branch(jtESumStr1.c_str(), &jtESumOutBad_[iter], jtESumStr2.c_str());

    outTreeBad_p[iter]->Branch(refPtStr1.c_str(), &refPtOutBad_[iter], refPtStr2.c_str());
    outTreeBad_p[iter]->Branch(refEtaStr1.c_str(), &refEtaOutBad_[iter], refEtaStr2.c_str());
    outTreeBad_p[iter]->Branch(refPhiStr1.c_str(), &refPhiOutBad_[iter], refPhiStr2.c_str());
    outTreeBad_p[iter]->Branch(refSubIDStr1.c_str(), &refSubIDOutBad_[iter], refSubIDStr2.c_str());
    outTreeBad_p[iter]->Branch(refPartFlavStr1.c_str(), &refPartFlavOutBad_[iter], refPartFlavStr2.c_str());

    std::string genJtPtStr1 = "genJtPt";
    std::string genJtPtStr2 = "genJtPt/D";
    std::string genJtPhiStr1 = "genJtPhi";
    std::string genJtPhiStr2 = "genJtPhi/D";
    std::string genJtDelZPhiStr1 = "genJtDelZPhi";
    std::string genJtDelZPhiStr2 = "genJtDelZPhi/D";
    std::string genJtEtaStr1 = "genJtEta";
    std::string genJtEtaStr2 = "genJtEta/D";
    std::string genJtMatchIndexStr1 = "genJtMatchIndex";
    std::string genJtMatchIndexStr2 = "genJtMatchIndex/I";
    std::string genJtMatchRecoPtStr1 = "genJtMatchRecoPt";
    std::string genJtMatchRecoPtStr2 = "genJtMatchRecoPt/D";
    std::string genJtMatchRecoPhiStr1 = "genJtMatchRecoPhi";
    std::string genJtMatchRecoPhiStr2 = "genJtMatchRecoPhi/D";
    std::string genJtMatchRecoEtaStr1 = "genJtMatchRecoEta";
    std::string genJtMatchRecoEtaStr2 = "genJtMatchRecoEta/D";

    std::string genJtManualMatchIndexStr1 = "genJtManualMatchIndex";
    std::string genJtManualMatchIndexStr2 = "genJtManualMatchIndex/I";
    std::string genJtManualMatchRecoPtStr1 = "genJtManualMatchRecoPt";
    std::string genJtManualMatchRecoPtStr2 = "genJtManualMatchRecoPt/D";
    std::string genJtManualMatchRecoPhiStr1 = "genJtManualMatchRecoPhi";
    std::string genJtManualMatchRecoPhiStr2 = "genJtManualMatchRecoPhi/D";
    std::string genJtManualMatchRecoEtaStr1 = "genJtManualMatchRecoEta";
    std::string genJtManualMatchRecoEtaStr2 = "genJtManualMatchRecoEta/D";

    std::string genJtMaxPartPtStr1 = "genJtMaxPartPt";
    std::string genJtMaxPartPtStr2 = "genJtMaxPartPt/D";
    std::string genJtMaxPartPhiStr1 = "genJtMaxPartPhi";
    std::string genJtMaxPartPhiStr2 = "genJtMaxPartPhi/D";
    std::string genJtMaxPartEtaStr1 = "genJtMaxPartEta";
    std::string genJtMaxPartEtaStr2 = "genJtMaxPartEta/D";
    std::string genJtMaxPartPDGIDStr1 = "genJtMaxPartPDGID";
    std::string genJtMaxPartPDGIDStr2 = "genJtMaxPartPDGID/I";
    std::string genJtSumPartPtStr1 = "genJtSumPartPt";
    std::string genJtSumPartPtStr2 = "genJtSumPartPt/D";
    std::string genJtNPartPtStr1 = "genJtNPartPt";
    std::string genJtNPartPtStr2 = "genJtNPartPt/I";

    std::string genJtSubIDStr1 = "genJtSubID";
    std::string genJtSubIDStr2 = "genJtSubID/I";

    std::string genJtPartFlavStr1 = "genJtPartFlav";
    std::string genJtPartFlavStr2 = "genJtPartFlav/I";
    std::string genJtManualMatchPartFlavStr1 = "genJtManualMatchPartFlav";
    std::string genJtManualMatchPartFlavStr2 = "genJtManualMatchPartFlav/I";

    std::string treeNameGen = "jecTree_Gen_" + jetAlgo.at(iter);
    outTreeGen_p[iter] = new TTree(treeNameGen.c_str(), treeNameGen.c_str());

    outTree_p[iter]->Branch(runStr1.c_str(), &RunOutGen_[iter], runStr2.c_str());
    outTree_p[iter]->Branch(lumiStr1.c_str(), &LumiOutGen_[iter], lumiStr2.c_str());
    outTree_p[iter]->Branch(eventStr1.c_str(), &EventOutGen_[iter], eventStr2.c_str());    

    outTreeGen_p[iter]->Branch(hiBinStr1.c_str(), &hiBinOutGen_[iter], hiBinStr2.c_str());
    outTreeGen_p[iter]->Branch(hiBinWeightStr1.c_str(), &hiBinWeightOutGen_[iter], hiBinWeightStr2.c_str());
    outTreeGen_p[iter]->Branch(hiBinWeightNormStr1.c_str(), &hiBinWeightNormOutGen_[iter], hiBinWeightNormStr2.c_str());
    outTreeGen_p[iter]->Branch(ptHatStr1.c_str(), &ptHatOutGen_[iter], ptHatStr2.c_str());
    outTreeGen_p[iter]->Branch(ptHatWeightStr1.c_str(), &ptHatWeightOutGen_[iter], ptHatWeightStr2.c_str());
    outTreeGen_p[iter]->Branch(fullWeightStr1.c_str(), &fullWeightOutGen_[iter], fullWeightStr2.c_str());
			       
    outTreeGen_p[iter]->Branch(genJtPtStr1.c_str(), &genJtPtOut_[iter], genJtPtStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtPhiStr1.c_str(), &genJtPhiOut_[iter], genJtPhiStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtDelZPhiStr1.c_str(), &genJtDelZPhiOut_[iter], genJtDelZPhiStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtEtaStr1.c_str(), &genJtEtaOut_[iter], genJtEtaStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtMatchIndexStr1.c_str(), &genJtMatchIndexOut_[iter], genJtMatchIndexStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtMatchRecoPtStr1.c_str(), &genJtMatchRecoPtOut_[iter], genJtMatchRecoPtStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtMatchRecoPhiStr1.c_str(), &genJtMatchRecoPhiOut_[iter], genJtMatchRecoPhiStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtMatchRecoEtaStr1.c_str(), &genJtMatchRecoEtaOut_[iter], genJtMatchRecoEtaStr2.c_str());

    outTreeGen_p[iter]->Branch(genJtManualMatchIndexStr1.c_str(), &genJtManualMatchIndexOut_[iter], genJtManualMatchIndexStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtManualMatchRecoPtStr1.c_str(), &genJtManualMatchRecoPtOut_[iter], genJtManualMatchRecoPtStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtManualMatchRecoPhiStr1.c_str(), &genJtManualMatchRecoPhiOut_[iter], genJtManualMatchRecoPhiStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtManualMatchRecoEtaStr1.c_str(), &genJtManualMatchRecoEtaOut_[iter], genJtManualMatchRecoEtaStr2.c_str());

    outTreeGen_p[iter]->Branch(genJtMaxPartPtStr1.c_str(), &genJtMaxPartPtOut_[iter], genJtMaxPartPtStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtMaxPartPhiStr1.c_str(), &genJtMaxPartPhiOut_[iter], genJtMaxPartPhiStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtMaxPartEtaStr1.c_str(), &genJtMaxPartEtaOut_[iter], genJtMaxPartEtaStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtMaxPartPDGIDStr1.c_str(), &genJtMaxPartPDGIDOut_[iter], genJtMaxPartPDGIDStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtSumPartPtStr1.c_str(), &genJtSumPartPtOut_[iter], genJtSumPartPtStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtNPartPtStr1.c_str(), &genJtNPartPtOut_[iter], genJtNPartPtStr2.c_str());

    outTreeGen_p[iter]->Branch(genJtSubIDStr1.c_str(), &genJtSubIDOut_[iter], genJtSubIDStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtPartFlavStr1.c_str(), &genJtPartFlavOut_[iter], genJtPartFlavStr2.c_str());
    outTreeGen_p[iter]->Branch(genJtManualMatchPartFlavStr1.c_str(), &genJtManualMatchPartFlavOut_[iter], genJtManualMatchPartFlavStr2.c_str());

  }

  Int_t nZMoreThanOne = 0;

  Int_t Run_;
  Int_t Lumi_;
  ULong64_t Event_;

  Float_t vz_;
  Int_t hiBin_ = 0;

  UInt_t run_;
  UInt_t lumi_;
  ULong64_t evt_;

  Int_t pcollisionEventSelection_;

  std::vector<float>* genPt_p = 0;
  std::vector<float>* genEta_p = 0;
  std::vector<float>* genPhi_p = 0;
  std::vector<int>* genPDG_p = 0;

  Int_t nJt_[nJetAlgo]; 
  Float_t jtPt_[nJetAlgo][nMaxJets];
  Float_t jtCorrFact_[nJetAlgo][nMaxJets];
  Float_t rawPt_[nJetAlgo][nMaxJets];
  Float_t jtEta_[nJetAlgo][nMaxJets];
  Float_t jtPhi_[nJetAlgo][nMaxJets];
  Float_t jtTrackSum_[nJetAlgo][nMaxJets];
  Float_t jtTrackMax_[nJetAlgo][nMaxJets];
  Float_t jtChargedSum_[nJetAlgo][nMaxJets];
  Float_t jtChargedMax_[nJetAlgo][nMaxJets];
  Float_t jtPhotonSum_[nJetAlgo][nMaxJets];
  Float_t jtNeutralSum_[nJetAlgo][nMaxJets];
  Float_t jtMuSum_[nJetAlgo][nMaxJets];
  Float_t jtESum_[nJetAlgo][nMaxJets];

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

  Int_t nMaxTrk = 10000;
  Int_t nTrk_;
  Float_t trkPt_[nMaxTrk];
  Float_t trkPhi_[nMaxTrk];
  Float_t trkEta_[nMaxTrk];
  Bool_t highPurity_[nMaxTrk];
  Float_t trkPtError_[nMaxTrk];
  Float_t trkDz1_[nMaxTrk];
  Float_t trkDzError1_[nMaxTrk];
  Float_t trkDxy1_[nMaxTrk];
  Float_t trkDxyError1_[nMaxTrk];

  int goodZEvents = 0;
  int goodVzEvents = 0;
  int goodPCollEvents = 0;
  int goodHiBinEvents = 0;

  int totalJets = 0;
  int jets5PtCut = 0;
  int jetsSubidCut = 0;
  int jetsLowCut = 0;
  int jetsHiCut = 0;
  int jetsEtaCut = 0;
  int jetsZCut = 0;
  int jetsMuCut = 0;
  int finalJets = 0;

  const unsigned int nInputs = config.GetNInputs();
  unsigned int inputCounts[nInputs];
  unsigned int inputEvents[nInputs];
  for(unsigned int iter = 0; iter < nInputs; iter++){
    inputCounts[iter] = 0;
    inputEvents[iter] = 0;
  }

  int nPho = 0;
  int nGenPho = 0;
  int nPhoAndGenPho = 0;
  int nPhoOrGenPho = 0;

  int nTruePho = 0;

  for(unsigned int pthatIter = 0; pthatIter < nInputs; pthatIter++){
    std::cout << "Begin processing " << pthatIter+1 << "/" << nInputs << "(ptHat == " << config.GetInputPtHat(pthatIter) << "): \'" << config.GetInput(pthatIter) << "\'...." << std::endl;

    std::vector<std::string> fileList;

    outFileRunLumiEvt.open("runLumiEvtFileForXCheck.txt", std::ios_base::app);
    outFileRunLumiEvt << "From " << pthatIter+1 << "/" << nInputs << "(ptHat == " << config.GetInputPtHat(pthatIter) << ")" << std::endl;
    outFileRunLumiEvt.close();


    if(checkDir(config.GetInput(pthatIter))){
      fileList = returnFileList(config.GetInput(pthatIter), "HiForest");
      
      int filePos = 0;
      while((int)fileList.size() > filePos){
	if(fileList.at(filePos).find("/failed/") != std::string::npos) fileList.erase(fileList.begin() + filePos);
	else filePos++;
      }
    }
    else if(checkFile(config.GetInput(pthatIter))) fileList.push_back(config.GetInput(pthatIter));
  
    for(unsigned int fileIter = 0; fileIter < fileList.size(); fileIter++){
        
      if(fileIter%10 == 0) std::cout << "File " << fileIter << "/" << fileList.size() << ": " << fileList.at(fileIter) << std::endl;

      

      TFile* inFile_p = TFile::Open(fileList.at(fileIter).c_str(), "READ");
      TTree* hltTree_p = (TTree*)inFile_p->Get("hltanalysis/HltTree");
      TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
      TTree* skimTree_p = (TTree*)inFile_p->Get("skimanalysis/HltTree");
      TTree* genTree_p = (TTree*)inFile_p->Get("HiGenParticleAna/hi");
      TTree* phoTree_p=0;
      if(config.GetIsPbPb()) phoTree_p = (TTree*)inFile_p->Get("ggHiNtuplizer/EventTree");
      else phoTree_p = (TTree*)inFile_p->Get("ggHiNtuplizerGED/EventTree");
      TTree* trkTree_p=0;
      if(config.GetIsPbPb()) trkTree_p = (TTree*)inFile_p->Get("anaTrack/trackTree");
      else trkTree_p = (TTree*)inFile_p->Get("ppTrack/trackTree");
      TTree* jetTree_p[nJetAlgo];
      
      if(debugMode) std::cout << __LINE__ << ", nJetAlgo: " << nJetAlgo << std::endl;

      for(Int_t iter = 0; iter < nJetAlgo; iter++){
	jetTree_p[iter] = (TTree*)inFile_p->Get(Form("%sJetAnalyzer/t", jetAlgo.at(iter).c_str()));
      }
      
      hltTree_p->SetBranchStatus("*", 0);
      hltTree_p->SetBranchStatus("Run", 1);
      hltTree_p->SetBranchStatus("LumiBlock", 1);
      hltTree_p->SetBranchStatus("Event", 1);

      hltTree_p->SetBranchAddress("Run", &Run_);
      hltTree_p->SetBranchAddress("LumiBlock", &Lumi_);
      hltTree_p->SetBranchAddress("Event", &Event_);

      hiTree_p->SetBranchStatus("*", 0);
      hiTree_p->SetBranchStatus("vz", 1);
      hiTree_p->SetBranchStatus("run", 1);
      hiTree_p->SetBranchStatus("lumi", 1);
      hiTree_p->SetBranchStatus("evt", 1);
      if(config.GetIsPbPb()) hiTree_p->SetBranchStatus("hiBin", 1);
      
      hiTree_p->SetBranchAddress("vz", &vz_);
      hiTree_p->SetBranchAddress("run", &run_);
      hiTree_p->SetBranchAddress("lumi", &lumi_);
      hiTree_p->SetBranchAddress("evt", &evt_);
      if(config.GetIsPbPb()) hiTree_p->SetBranchAddress("hiBin", &hiBin_);


      skimTree_p->SetBranchStatus("*", 0);

      if(config.GetIsPbPb()){
	skimTree_p->SetBranchStatus("pcollisionEventSelection", 1);
	
	skimTree_p->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection_);
      }
  
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
	jetTree_p[iter]->SetBranchStatus("rawpt", 1);
	jetTree_p[iter]->SetBranchStatus("jteta", 1);
	jetTree_p[iter]->SetBranchStatus("jtphi", 1);
	jetTree_p[iter]->SetBranchStatus("trackSum", 1);
	jetTree_p[iter]->SetBranchStatus("trackMax", 1);
	jetTree_p[iter]->SetBranchStatus("chargedSum", 1);
	jetTree_p[iter]->SetBranchStatus("chargedMax", 1);
	jetTree_p[iter]->SetBranchStatus("photonSum", 1);
	jetTree_p[iter]->SetBranchStatus("neutralSum", 1);
	jetTree_p[iter]->SetBranchStatus("eSum", 1);
	jetTree_p[iter]->SetBranchStatus("muSum", 1);

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
	jetTree_p[iter]->SetBranchAddress("rawpt", rawPt_[iter]);
	jetTree_p[iter]->SetBranchAddress("jteta", jtEta_[iter]);
	jetTree_p[iter]->SetBranchAddress("jtphi", jtPhi_[iter]);
	jetTree_p[iter]->SetBranchAddress("trackSum", jtTrackSum_[iter]);
	jetTree_p[iter]->SetBranchAddress("trackMax", jtTrackMax_[iter]);
	jetTree_p[iter]->SetBranchAddress("chargedSum", jtChargedSum_[iter]);
	jetTree_p[iter]->SetBranchAddress("chargedMax", jtChargedMax_[iter]);
	jetTree_p[iter]->SetBranchAddress("neutralSum", jtNeutralSum_[iter]);
	jetTree_p[iter]->SetBranchAddress("photonSum", jtPhotonSum_[iter]);
	jetTree_p[iter]->SetBranchAddress("eSum", jtESum_[iter]);
	jetTree_p[iter]->SetBranchAddress("muSum", jtMuSum_[iter]);

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

      trkTree_p->SetBranchStatus("*", 0);
      trkTree_p->SetBranchStatus("nTrk", 1);
      trkTree_p->SetBranchStatus("trkPt", 1);
      trkTree_p->SetBranchStatus("trkPhi", 1);
      trkTree_p->SetBranchStatus("trkEta", 1);
      trkTree_p->SetBranchStatus("highPurity", 1);
      trkTree_p->SetBranchStatus("trkPtError", 1);
      trkTree_p->SetBranchStatus("trkDz1", 1);
      trkTree_p->SetBranchStatus("trkDzError1", 1);
      trkTree_p->SetBranchStatus("trkDxy1", 1);
      trkTree_p->SetBranchStatus("trkDxyError1", 1);

      trkTree_p->SetBranchAddress("nTrk", &nTrk_);
      trkTree_p->SetBranchAddress("trkPt", trkPt_);
      trkTree_p->SetBranchAddress("trkPhi", trkPhi_);
      trkTree_p->SetBranchAddress("trkEta", trkEta_);
      trkTree_p->SetBranchAddress("highPurity", highPurity_);
      trkTree_p->SetBranchAddress("trkPtError", trkPtError_);
      trkTree_p->SetBranchAddress("trkDz1", trkDz1_);
      trkTree_p->SetBranchAddress("trkDzError1", trkDzError1_);
      trkTree_p->SetBranchAddress("trkDxy1", trkDxy1_);
      trkTree_p->SetBranchAddress("trkDxyError1", trkDxyError1_);

      if(debugMode) std::cout << __LINE__ << std::endl;

      Int_t tempStartPos = 0;
      Int_t tempNEntries = jetTree_p[0]->GetEntries();

      if(debugMode) std::cout << __LINE__ << std::endl;
      
      const Int_t startPos = tempStartPos;
      const Int_t nEntries = tempNEntries;
      Int_t entryDiv = TMath::Max(1, ((Int_t)((nEntries-startPos)/20)));
      
      if(entryDiv < 25000) entryDiv = 25000;

      if(debugMode) std::cout << __LINE__ << std::endl;
    
      for(Int_t entry = startPos; entry < nEntries; entry++){
	if(entry%entryDiv == 0 && nEntries >= 10000){
	  std::cout << "Entry # " << entry << "/" << nEntries << std::endl;
	}

	if(debugMode) std::cout << __LINE__ << std::endl;
	
	if(debugMode) std::cout << __LINE__ << std::endl;

	hltTree_p->GetEntry(entry);
	skimTree_p->GetEntry(entry);
	hiTree_p->GetEntry(entry);
	genTree_p->GetEntry(entry);
	phoTree_p->GetEntry(entry);

	if(TMath::Abs(vz_) > 15) continue;

	goodVzEvents++;

	if(!pcollisionEventSelection_ && config.GetIsPbPb()) continue;

	goodPCollEvents++;

	Int_t centPos = config.GetCentBinFromHiBin(hiBin_);
	if(centPos == -1){
	  continue;
	}

	goodHiBinEvents++;


	if(debugMode) std::cout << __LINE__ << std::endl;
	
	
	Float_t recoTempLeadingPhoPt = -999;
	Float_t recoTempLeadingPhoPhi = -999;
	Float_t recoTempLeadingPhoEta = -999;

	Float_t recoTempLeadingPho2Pt = -999;
	Float_t recoTempLeadingPho2Phi = -999;
	Float_t recoTempLeadingPho2Eta = -999;


	Float_t genTempLeadingPhoPt = -999;
	Float_t genTempLeadingPhoPhi = -999;
	Float_t genTempLeadingPhoEta = -999;
	
	std::vector<Float_t> passingMuPt_p;
	std::vector<Float_t> passingMuPhi_p;
	std::vector<Float_t> passingMuEta_p;
	std::vector<Int_t> passingMuPDG_p;
	
	std::vector<Float_t> passingElePt_p;
	std::vector<Float_t> passingElePhi_p;
	std::vector<Float_t> passingEleEta_p;
	std::vector<Int_t> passingElePDG_p;
	
	std::vector<Float_t> zPt_p;
	std::vector<Float_t> zPhi_p;
	std::vector<Float_t> zEta_p;
	std::vector<Float_t> zM_p;

	std::vector<Float_t> zLeadLeptPt_p;
	std::vector<Float_t> zLeadLeptPhi_p;
	std::vector<Float_t> zLeadLeptEta_p;
	
	std::vector<Float_t> zSubleadLeptPt_p;
	std::vector<Float_t> zSubleadLeptPhi_p;
	std::vector<Float_t> zSubleadLeptEta_p;
	
	bool isLead = true;

	if(config.GetIsGammaJet() || config.GetDoGenGammaCutOverride()){		  
	  const Int_t nPho = phoPt_p->size();
	  
	  double tempLeadPhoPt = -999;
	  int tempLeadPhoPos = -1;
	  
	  for(Int_t phoIter = 0; phoIter < nPho; ++phoIter){
            if(TMath::Abs(phoEta_p->at(phoIter)) > 1.44) continue;

            if(phoSigmaIEtaIEta_2012_p->at(phoIter) < .002) continue;
            if(TMath::Abs(phoSeedTime_p->at(phoIter)) > 3) continue;
            if(phoSwissCrx_p->at(phoIter) > .9) continue;

	    if(phoPt_p->at(phoIter) > tempLeadPhoPt){
	      tempLeadPhoPt = phoPt_p->at(phoIter);
	      tempLeadPhoPos = phoIter;
	    }

	  }

	  for(Int_t phoIter = 0; phoIter < nPho; ++phoIter){
	    if(TMath::Abs(phoEta_p->at(phoIter)) > 1.44) continue;
	    
	    if(phoSigmaIEtaIEta_2012_p->at(phoIter) < .002) continue;
	    if(TMath::Abs(phoSeedTime_p->at(phoIter)) > 3) continue;
	    if(phoSwissCrx_p->at(phoIter) > .9) continue;

	    

	    if(phoSigmaIEtaIEta_2012_p->at(phoIter) > .01 && phoIter == tempLeadPhoPos){
	      isLead = false;
	      continue;
	    }
	    if(pho_ecalClusterIsoR4_p->at(phoIter) + pho_hcalRechitIsoR4_p->at(phoIter) + pho_trackIsoR4PtCut20_p->at(phoIter) > 1 && phoIter == tempLeadPhoPos){
	      isLead = false;
	      continue;
	    }
	    if(phoHoverE_p->at(phoIter) > .1 && phoIter == tempLeadPhoPos){
	      isLead = false;
	      continue;
	    }
	    
	    if(recoTempLeadingPhoPt < phoPt_p->at(phoIter)){
	      if(isLead && phoPt_p->at(phoIter) > 40 && phoIter == tempLeadPhoPos){
		nTruePho++;

		outFileRunLumiEvt.open("runLumiEvtFileForXCheck.txt", std::ios_base::app);
		outFileRunLumiEvt << run_ << "," << lumi_ << "," << evt_ << std::endl;
		outFileRunLumiEvt.close();

		recoTempLeadingPho2Pt = phoPt_p->at(phoIter);
		recoTempLeadingPho2Eta = phoEta_p->at(phoIter);
		recoTempLeadingPho2Phi = phoPhi_p->at(phoIter);

		isLead = false;
	      }

	      recoTempLeadingPhoPt = phoPt_p->at(phoIter);
	      recoTempLeadingPhoEta = phoEta_p->at(phoIter);
	      recoTempLeadingPhoPhi = phoPhi_p->at(phoIter);
	    }
	  }

	  const Int_t nGenPho = mcPt_p->size();

          for(Int_t phoIter = 0; phoIter < nGenPho; phoIter++){
            if(TMath::Abs(mcEta_p->at(phoIter)) > 1.44) continue;

	    if(TMath::Abs(mcPID_p->at(phoIter)) != 22) continue;

            if(genTempLeadingPhoPt < mcPt_p->at(phoIter)){
              genTempLeadingPhoPt = mcPt_p->at(phoIter);
              genTempLeadingPhoEta = mcEta_p->at(phoIter);
              genTempLeadingPhoPhi = mcPhi_p->at(phoIter);
            }
          }
	}
	if(config.GetIsZJet()){
	  const Int_t nGenPart = genPt_p->size();
	  
	  for(Int_t genIter = 0; genIter < nGenPart; genIter++){                                                   
	    if(TMath::Abs(genPDG_p->at(genIter)) != 11 && TMath::Abs(genPDG_p->at(genIter)) != 13) continue; 
	    
	    if(TMath::Abs(genPDG_p->at(genIter)) == 11){
	      if(config.KeepLepton(genPt_p->at(genIter), genEta_p->at(genIter), 11)){
		passingElePt_p.push_back(genPt_p->at(genIter)); 
		passingElePhi_p.push_back(genPhi_p->at(genIter)); 
		passingEleEta_p.push_back(genEta_p->at(genIter)); 
		passingElePDG_p.push_back(genPDG_p->at(genIter)); 
	      }
	    }
	    else if(TMath::Abs(genPDG_p->at(genIter)) == 13){
	      if(config.KeepLepton(genPt_p->at(genIter), genEta_p->at(genIter), 13)){
		passingMuPt_p.push_back(genPt_p->at(genIter));
		passingMuPhi_p.push_back(genPhi_p->at(genIter));
		passingMuEta_p.push_back(genEta_p->at(genIter));
		passingMuPDG_p.push_back(genPDG_p->at(genIter));
	      }
	    }
	  }                                          
	}
            	
	if(config.GetIsGammaJet()){
	  if(!config.KeepEventGamma(recoTempLeadingPhoPt, recoTempLeadingPhoEta) && !config.KeepEventGamma(genTempLeadingPhoPt, genTempLeadingPhoEta)) continue;
	}

	if(config.GetIsGammaJet()){
	  if(config.KeepEventGamma(recoTempLeadingPhoPt, recoTempLeadingPhoEta)) ++nPho;
	  if(config.KeepEventGamma(genTempLeadingPhoPt, genTempLeadingPhoEta)) ++nGenPho;
	  if(config.KeepEventGamma(recoTempLeadingPhoPt, recoTempLeadingPhoEta) && config.KeepEventGamma(genTempLeadingPhoPt, genTempLeadingPhoEta)) ++nPhoAndGenPho;
	  if(config.KeepEventGamma(recoTempLeadingPhoPt, recoTempLeadingPhoEta) || config.KeepEventGamma(genTempLeadingPhoPt, genTempLeadingPhoEta)) ++nPhoOrGenPho;
	}
	
	if(config.GetIsZJet()){
	  unsigned int eleSortIter = 0;
	  while(eleSortIter < passingElePt_p.size()){
	    bool isLeading = true;
	    
	    for(unsigned int eleIter = eleSortIter+1; eleIter < passingElePt_p.size(); eleIter++){
	      if(passingElePt_p.at(eleSortIter) < passingElePt_p.at(eleIter)){
		Float_t tempElePt = passingElePt_p.at(eleSortIter);
		Float_t tempElePhi = passingElePhi_p.at(eleSortIter);
		Float_t tempEleEta = passingEleEta_p.at(eleSortIter);
		Int_t tempElePDG = passingElePDG_p.at(eleSortIter);
		
		passingElePt_p.at(eleSortIter) = passingElePt_p.at(eleIter);
		passingElePhi_p.at(eleSortIter) = passingElePhi_p.at(eleIter);
		passingEleEta_p.at(eleSortIter) = passingEleEta_p.at(eleIter);
		passingElePDG_p.at(eleSortIter) = passingElePDG_p.at(eleIter);
		
		passingElePt_p.at(eleIter) = tempElePt;
		passingElePhi_p.at(eleIter) = tempElePhi;
		passingEleEta_p.at(eleIter) = tempEleEta;
		passingElePDG_p.at(eleIter) = tempElePDG;

		isLeading = false;
	      }
	    }
	    
	    if(isLeading) eleSortIter++;
	  }
		  
	  unsigned int muSortIter = 0;
	  while(muSortIter < passingMuPt_p.size()){
	    bool isLeading = true;
	    
	    for(unsigned int muIter = muSortIter+1; muIter < passingMuPt_p.size(); muIter++){
	      if(passingMuPt_p.at(muSortIter) < passingMuPt_p.at(muIter)){
		Float_t tempMuPt = passingMuPt_p.at(muSortIter);
		Float_t tempMuPhi = passingMuPhi_p.at(muSortIter);
		Float_t tempMuEta = passingMuEta_p.at(muSortIter);
		Int_t tempMuPDG = passingMuPDG_p.at(muSortIter);
		
		passingMuPt_p.at(muSortIter) = passingMuPt_p.at(muIter);
		passingMuPhi_p.at(muSortIter) = passingMuPhi_p.at(muIter);
		passingMuEta_p.at(muSortIter) = passingMuEta_p.at(muIter);
		passingMuPDG_p.at(muSortIter) = passingMuPDG_p.at(muIter);
		
		passingMuPt_p.at(muIter) = tempMuPt;
		passingMuPhi_p.at(muIter) = tempMuPhi;
		passingMuEta_p.at(muIter) = tempMuEta;
		passingMuPDG_p.at(muIter) = tempMuPDG;
		
		isLeading = false;
	      }
	    }
	    
	    if(isLeading) muSortIter++;
	  }
	  
	  Bool_t isGoodZ = false;
	  
	  if(passingMuPt_p.size() > 1){
	    TLorentzVector mu1, mu2;
	    
	    for(unsigned int muIter = 0; muIter < passingMuPt_p.size(); muIter++){
	      mu1.SetPtEtaPhiM(passingMuPt_p.at(muIter), passingMuEta_p.at(muIter), passingMuPhi_p.at(muIter), muMass);
	      
	      for(unsigned int muIter2 = muIter+1; muIter2 < passingMuPt_p.size(); muIter2++){
		mu2.SetPtEtaPhiM(passingMuPt_p.at(muIter2), passingMuEta_p.at(muIter2), passingMuPhi_p.at(muIter2), muMass);

		if(passingMuPDG_p.at(muIter) == passingMuPDG_p.at(muIter2)) continue;
		
		TLorentzVector z = mu1+mu2;
		if(config.KeepEventZ(z.Pt(), z.M())){
		  isGoodZ = true;
		  zPt_p.push_back(z.Pt());
		  zPhi_p.push_back(z.Phi());
		  zEta_p.push_back(z.Eta());
		  zM_p.push_back(z.M());
		  
		  zLeadLeptPt_p.push_back(passingMuPt_p.at(muIter));
		  zLeadLeptPhi_p.push_back(passingMuPhi_p.at(muIter));
		  zLeadLeptEta_p.push_back(passingMuEta_p.at(muIter));
		  
		  zSubleadLeptPt_p.push_back(passingMuPt_p.at(muIter2));
		  zSubleadLeptPhi_p.push_back(passingMuPhi_p.at(muIter2));
		  zSubleadLeptEta_p.push_back(passingMuEta_p.at(muIter2));		
		}	 
	      }
	    }
	  }
		  
	  if(passingElePt_p.size() > 1){
	    TLorentzVector ele1, ele2;
	    
	    for(unsigned int eleIter = 0; eleIter < passingElePt_p.size(); eleIter++){
	      ele1.SetPtEtaPhiM(passingElePt_p.at(eleIter), passingEleEta_p.at(eleIter), passingElePhi_p.at(eleIter), eleMass);
	      
	      for(unsigned int eleIter2 = eleIter+1; eleIter2 < passingElePt_p.size(); eleIter2++){
		ele2.SetPtEtaPhiM(passingElePt_p.at(eleIter2), passingEleEta_p.at(eleIter2), passingElePhi_p.at(eleIter2), eleMass);

		if(passingElePDG_p.at(eleIter) == passingElePDG_p.at(eleIter2)) continue;
		
		TLorentzVector z = ele1+ele2;
		if(config.KeepEventZ(z.Pt(), z.M())){
		  isGoodZ = true;
		  zPt_p.push_back(z.Pt());
		  zPhi_p.push_back(z.Phi());
		  zEta_p.push_back(z.Eta());
		  zM_p.push_back(z.M());
		  
		  zLeadLeptPt_p.push_back(passingElePt_p.at(eleIter));
		  zLeadLeptPhi_p.push_back(passingElePhi_p.at(eleIter));
		  zLeadLeptEta_p.push_back(passingEleEta_p.at(eleIter));
		  
		  zSubleadLeptPt_p.push_back(passingElePt_p.at(eleIter2));
		  zSubleadLeptPhi_p.push_back(passingElePhi_p.at(eleIter2));
		  zSubleadLeptEta_p.push_back(passingEleEta_p.at(eleIter2));		
		}	 
	      }
	    }
	  }
	    
	  if(!isGoodZ) continue;
	  
	  goodZEvents++;

	  if(zPt_p.size() > 1) nZMoreThanOne++;
	}
      

	if(debugMode) std::cout << __LINE__ << std::endl;
	
	trkTree_p->GetEntry(entry);

	for(Int_t iter = 0; iter < nJetAlgo; iter++){
	  jetTree_p[iter]->GetEntry(entry);
	}		

	if(debugMode) std::cout << __LINE__ << std::endl;

	for(Int_t algoIter = 0; algoIter < nJetAlgo; algoIter++){
	  for(Int_t jtIter = 0; jtIter < nJt_[algoIter]; jtIter++){
	    jtCorrFact_[algoIter][jtIter] = 1.;
	  }
	}
      	
	if(config.GetDoCorrections()){
	  if(config.GetIsPbPb()){
	    for(Int_t algoIter = 0; algoIter < nJetAlgo; algoIter++){
	      for(Int_t jtIter = 0; jtIter < nJt_[algoIter]; jtIter++){
		//		jtPt_[algoIter][jtIter] *= getResCorrJetPt(jtPt_[algoIter][jtIter], hiBin_);                      
		jtCorrFact_[algoIter][jtIter] = config.GetConstCorrFactor()*getPtEtaJetResidualCorr(jtPt_[algoIter][jtIter], jtEta_[algoIter][jtIter], hiBin_/2., config.GetCorrForm());

		jtPt_[algoIter][jtIter] *= config.GetConstCorrFactor()*getPtEtaJetResidualCorr(jtPt_[algoIter][jtIter], jtEta_[algoIter][jtIter], hiBin_/2., config.GetCorrForm());
	      }
	    }
	  }
	  else{
	    for(Int_t algoIter = 0; algoIter < nJetAlgo; algoIter++){
	      for(Int_t jtIter = 0; jtIter < nJt_[algoIter]; jtIter++){
		jtCorrFact_[algoIter][jtIter] = config.GetConstCorrFactor();
		jtPt_[algoIter][jtIter] *= config.GetConstCorrFactor();
	      }
	    }
	  }
	}

	if(debugMode) std::cout << __LINE__ << std::endl;

	if(debugMode) std::cout << __LINE__ << std::endl;

	inputEvents[pthatIter]++;
      
	for(Int_t algoIter = 0; algoIter < nJetAlgo; algoIter++){
	  RunOut_[algoIter] = Run_;
	  LumiOut_[algoIter] = Lumi_;
	  EventOut_[algoIter] = Event_;

	  hiBinOut_[algoIter] = hiBin_;
	  ptHatOut_[algoIter] = ptHat_[algoIter];
	  ptHatWeightOut_[algoIter] = config.GetPtHatWeight(ptHat_[algoIter]);

	  if(config.GetIsPbPb()){
	    hiBinWeightOut_[algoIter] = findNcoll(hiBin_);	
	    hiBinWeightNormOut_[algoIter] = findNormNcoll(hiBin_, centPos, config.GetCentBins());
	  }
	  else{
	    hiBinWeightOut_[algoIter] = 1;
            hiBinWeightNormOut_[algoIter] = 1;
	  }

	  fullWeightOut_[algoIter] = ptHatWeightOut_[algoIter]*hiBinWeightOut_[algoIter];

	  RunOutBad_[algoIter] = Run_;
	  LumiOutBad_[algoIter] = Lumi_;
	  EventOutBad_[algoIter] = Event_;

	  hiBinOutBad_[algoIter] = hiBin_;
	  ptHatOutBad_[algoIter] = ptHat_[algoIter];
	  ptHatWeightOutBad_[algoIter] = config.GetPtHatWeight(ptHat_[algoIter]);

	  if(config.GetIsPbPb()){
	    hiBinWeightOutBad_[algoIter] = findNcoll(hiBin_);	
	    hiBinWeightNormOutBad_[algoIter] = findNormNcoll(hiBin_, centPos, config.GetCentBins());
	  }
	  else{
	    hiBinWeightOutBad_[algoIter] = 1;
            hiBinWeightNormOutBad_[algoIter] = 1;
	  }

	  fullWeightOutBad_[algoIter] = ptHatWeightOutBad_[algoIter]*hiBinWeightOutBad_[algoIter];

	  if(debugMode) std::cout << __LINE__ << std::endl;

	  if(config.GetIsZJet()){
	    leptPt1Out_[algoIter] = zLeadLeptPt_p.at(0);
	    leptPhi1Out_[algoIter] = zLeadLeptPhi_p.at(0);
	    leptEta1Out_[algoIter] = zLeadLeptEta_p.at(0);

	    leptPt2Out_[algoIter] = zSubleadLeptPt_p.at(0);
	    leptPhi2Out_[algoIter] = zSubleadLeptPhi_p.at(0);
	    leptEta2Out_[algoIter] = zSubleadLeptEta_p.at(0);

	    zPtOut_[algoIter] = zPt_p.at(0);
	    zPhiOut_[algoIter] = zPhi_p.at(0);
	    zEtaOut_[algoIter] = zEta_p.at(0);
	    zMOut_[algoIter] = zM_p.at(0);	    

	    leptPt1OutBad_[algoIter] = zLeadLeptPt_p.at(0);
	    leptPhi1OutBad_[algoIter] = zLeadLeptPhi_p.at(0);
	    leptEta1OutBad_[algoIter] = zLeadLeptEta_p.at(0);

	    leptPt2OutBad_[algoIter] = zSubleadLeptPt_p.at(0);
	    leptPhi2OutBad_[algoIter] = zSubleadLeptPhi_p.at(0);
	    leptEta2OutBad_[algoIter] = zSubleadLeptEta_p.at(0);

	    zPtOutBad_[algoIter] = zPt_p.at(0);
	    zPhiOutBad_[algoIter] = zPhi_p.at(0);
	    zEtaOutBad_[algoIter] = zEta_p.at(0);
	    zMOutBad_[algoIter] = zM_p.at(0);	    
	  }
	  else if(config.GetIsGammaJet()){
	    phoPtOut_[algoIter] = recoTempLeadingPhoPt;
	    phoPhiOut_[algoIter] = recoTempLeadingPhoPhi;
	    phoEtaOut_[algoIter] = recoTempLeadingPhoEta;

	    phoPtOutBad_[algoIter] = recoTempLeadingPhoPt;
	    phoPhiOutBad_[algoIter] = recoTempLeadingPhoPhi;
	    phoEtaOutBad_[algoIter] = recoTempLeadingPhoEta;

	    pho2PtOut_[algoIter] = recoTempLeadingPho2Pt;
	    pho2PhiOut_[algoIter] = recoTempLeadingPho2Phi;
	    pho2EtaOut_[algoIter] = recoTempLeadingPho2Eta;

	    pho2PtOutBad_[algoIter] = recoTempLeadingPho2Pt;
	    pho2PhiOutBad_[algoIter] = recoTempLeadingPho2Phi;
	    pho2EtaOutBad_[algoIter] = recoTempLeadingPho2Eta;

	    genPhoPtOut_[algoIter] = genTempLeadingPhoPt;
	    genPhoPhiOut_[algoIter] = genTempLeadingPhoPhi;
	    genPhoEtaOut_[algoIter] = genTempLeadingPhoEta;

	    genPhoPtOutBad_[algoIter] = genTempLeadingPhoPt;
	    genPhoPhiOutBad_[algoIter] = genTempLeadingPhoPhi;
	    genPhoEtaOutBad_[algoIter] = genTempLeadingPhoEta;
	  }


	  if(debugMode) std::cout << __LINE__ << std::endl;

	  genSort(nGenJt_[algoIter], genJtPt_[algoIter], genJtPhi_[algoIter], genJtEta_[algoIter], genJtMatchIndex_[algoIter], genJtSubID_[algoIter]);

	  if(debugMode) std::cout << __LINE__ << std::endl;
	 
	  //MODDING FOR REFPT USAGE
	  for(Int_t jtIter = 0; jtIter < nJt_[algoIter]; jtIter++){
	    jtPtOutBad_[algoIter] = jtPt_[algoIter][jtIter];
	    rawPtOutBad_[algoIter] = rawPt_[algoIter][jtIter];
	    if(config.GetDoCorrections()) jtPtOrigOutBad_[algoIter] = jtPt_[algoIter][jtIter]/jtCorrFact_[algoIter][jtIter];
	    else jtPtOrigOutBad_[algoIter] = jtPt_[algoIter][jtIter];
	    jtPhiOutBad_[algoIter] = jtPhi_[algoIter][jtIter];
	    if(config.GetIsZJet()) jtDelZPhiOutBad_[algoIter] = TMath::Abs(getDPHI(jtPhi_[algoIter][jtIter], zPhi_p.at(0)));
	    else jtDelZPhiOutBad_[algoIter] = -999;

	    if(config.GetIsGammaJet()){
	      if(recoTempLeadingPhoPt > 0) jtDelPhoPhiOutBad_[algoIter] = TMath::Abs(getDPHI(jtPhi_[algoIter][jtIter], recoTempLeadingPhoPhi));
	      if(recoTempLeadingPho2Pt > 0) jtDelPho2PhiOutBad_[algoIter] = TMath::Abs(getDPHI(jtPhi_[algoIter][jtIter], recoTempLeadingPho2Phi));
	    }
	    else{
	      jtDelPhoPhiOutBad_[algoIter] = -999;
	      jtDelPho2PhiOutBad_[algoIter] = -999;
	    }

	    jtEtaOutBad_[algoIter] = jtEta_[algoIter][jtIter];
	    jtTrackSumOutBad_[algoIter] = jtTrackSum_[algoIter][jtIter];
	    jtTrackMaxOutBad_[algoIter] = jtTrackMax_[algoIter][jtIter];
	    jtChargedSumOutBad_[algoIter] = jtChargedSum_[algoIter][jtIter];
	    jtNeutralSumOutBad_[algoIter] = jtNeutralSum_[algoIter][jtIter];
	    jtPhotonSumOutBad_[algoIter] = jtPhotonSum_[algoIter][jtIter];
	    jtESumOutBad_[algoIter] = jtESum_[algoIter][jtIter];
	    jtMuSumOutBad_[algoIter] = jtMuSum_[algoIter][jtIter];

	    if(debugMode) std::cout << __LINE__ << std::endl;

	    jtOverRefPtOutBad_[algoIter] = jtPt_[algoIter][jtIter]/refPt_[algoIter][jtIter];
	    jtOrigOverRefPtOutBad_[algoIter] = (jtPt_[algoIter][jtIter]/jtCorrFact_[algoIter][jtIter])/refPt_[algoIter][jtIter];
	    if(refPt_[algoIter][jtIter] > 0) jtDelRefPhiOutBad_[algoIter] = getDPHI(jtPhi_[algoIter][jtIter], refPhi_[algoIter][jtIter]);
	    else jtDelRefPhiOutBad_[algoIter] = -999;
	    if(refPt_[algoIter][jtIter] > 0) jtDelRefEtaOutBad_[algoIter] = jtEta_[algoIter][jtIter] -refEta_[algoIter][jtIter];
	    else jtDelRefEtaOutBad_[algoIter] = -999;

	    refPtOutBad_[algoIter] = refPt_[algoIter][jtIter];
	    refPhiOutBad_[algoIter] = refPhi_[algoIter][jtIter];
	    refEtaOutBad_[algoIter] = refEta_[algoIter][jtIter];
	    refSubIDOutBad_[algoIter] = refSubID_[algoIter][jtIter];
	    refPartFlavOutBad_[algoIter] = refPartFlav_[algoIter][jtIter];

	    if(debugMode) std::cout << __LINE__ << std::endl;

	    totalJets++;
	    if(refEta_[algoIter][jtIter] > config.GetJtEtaMax()){
	      outTreeBad_p[algoIter]->Fill();
	      continue;
	    }
	    if(refEta_[algoIter][jtIter] < config.GetJtEtaMin()){
	      outTreeBad_p[algoIter]->Fill();
	      continue;
	    }
	    jetsEtaCut++;

	    if(debugMode) std::cout << __LINE__ << std::endl;

	    if(refPt_[algoIter][jtIter] < 5.0){
	      outTreeBad_p[algoIter]->Fill();
	      continue;
	    }
	    jets5PtCut++;
	    if(refSubID_[algoIter][jtIter] != 0){
	      outTreeBad_p[algoIter]->Fill();
	      continue;
	    }
	    jetsSubidCut++;
	    
	    if(refPt_[algoIter][jtIter] < config.GetJtPtMin()){
	      outTreeBad_p[algoIter]->Fill();
	      continue;
	    }
	    jetsLowCut++;
	    if(refPt_[algoIter][jtIter] > config.GetJtPtMax()){
	      outTreeBad_p[algoIter]->Fill();
	      continue;
	    }

	    jetsHiCut++;
	  	
	    if(config.GetIsGammaJet() && recoTempLeadingPhoPt > 0){
	      if(getDR(recoTempLeadingPhoEta, recoTempLeadingPhoPhi, refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter]) < 0.4){
		outTreeBad_p[algoIter]->Fill();
		continue;
	      }
	    }	  
	    else if(config.GetIsZJet()){
	      if(getDR(zLeadLeptEta_p.at(0), zLeadLeptPhi_p.at(0), refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter]) < 0.4){
		outTreeBad_p[algoIter]->Fill();
		continue;
	      }
	      if(getDR(zSubleadLeptEta_p.at(0), zSubleadLeptPhi_p.at(0), refEta_[algoIter][jtIter], refPhi_[algoIter][jtIter]) < 0.4){
		outTreeBad_p[algoIter]->Fill();
		continue;
	      }   
	      
	      jetsMuCut++;
	      
	      if(config.GetDoZJtDPhiCut()){
		if(!config.PassesZJetDPhiCut(zPhi_p.at(0), refPhi_[algoIter][jtIter])){
		  outTreeBad_p[algoIter]->Fill();
		  continue;
		} 
	      }
	      jetsZCut++;

	      finalJets++;
	    }

	    if(debugMode) std::cout << __LINE__ << std::endl;

	    jtPtOut_[algoIter] = jtPt_[algoIter][jtIter];
	    rawPtOut_[algoIter] = rawPt_[algoIter][jtIter];
	    if(config.GetDoCorrections()) jtPtOrigOut_[algoIter] = jtPt_[algoIter][jtIter]/jtCorrFact_[algoIter][jtIter];
	    else jtPtOrigOut_[algoIter] = jtPt_[algoIter][jtIter];
	    jtPhiOut_[algoIter] = jtPhi_[algoIter][jtIter];
	    if(config.GetIsZJet()) jtDelZPhiOut_[algoIter] = TMath::Abs(getDPHI(jtPhi_[algoIter][jtIter], zPhi_p.at(0)));
	    else jtDelZPhiOut_[algoIter] = -999;

	    if(config.GetIsGammaJet()){
	      if(recoTempLeadingPhoPt > 0) jtDelPhoPhiOut_[algoIter] = TMath::Abs(getDPHI(jtPhi_[algoIter][jtIter], recoTempLeadingPhoPhi));
	      if(recoTempLeadingPho2Pt > 0) jtDelPho2PhiOut_[algoIter] = TMath::Abs(getDPHI(jtPhi_[algoIter][jtIter], recoTempLeadingPho2Phi));

	    }
	    else{
	      jtDelPhoPhiOut_[algoIter] = -999;
	      jtDelPho2PhiOut_[algoIter] = -999;
	    }

	    jtEtaOut_[algoIter] = jtEta_[algoIter][jtIter];
	    jtTrackSumOut_[algoIter] = jtTrackSum_[algoIter][jtIter];
	    jtTrackMaxOut_[algoIter] = jtTrackMax_[algoIter][jtIter];
	    jtChargedSumOut_[algoIter] = jtChargedSum_[algoIter][jtIter];
	    jtChargedMaxOut_[algoIter] = jtChargedMax_[algoIter][jtIter];
	    jtNeutralSumOut_[algoIter] = jtNeutralSum_[algoIter][jtIter];
	    jtPhotonSumOut_[algoIter] = jtPhotonSum_[algoIter][jtIter];
	    jtESumOut_[algoIter] = jtESum_[algoIter][jtIter];
	    jtMuSumOut_[algoIter] = jtMuSum_[algoIter][jtIter];

	    jtOverRefPtOut_[algoIter] = jtPt_[algoIter][jtIter]/refPt_[algoIter][jtIter];
	    jtOrigOverRefPtOut_[algoIter] = (jtPt_[algoIter][jtIter]/jtCorrFact_[algoIter][jtIter])/refPt_[algoIter][jtIter];
	    jtDelRefPhiOut_[algoIter] = getDPHI(jtPhi_[algoIter][jtIter], refPhi_[algoIter][jtIter]);
	    jtDelRefEtaOut_[algoIter] = jtEta_[algoIter][jtIter] - refEta_[algoIter][jtIter];

	    refPtOut_[algoIter] = refPt_[algoIter][jtIter];
	    refPhiOut_[algoIter] = refPhi_[algoIter][jtIter];
	    refEtaOut_[algoIter] = refEta_[algoIter][jtIter];
	    refSubIDOut_[algoIter] = refSubID_[algoIter][jtIter];
	    refPartFlavOut_[algoIter] = refPartFlav_[algoIter][jtIter];

	    jtMaxTrackInConePtOut_[algoIter] = -999;
	    jtMaxTrackInConeHPOut_[algoIter] = false;

	    for(int trkIter = 0; trkIter < nTrk_; trkIter++){
	      if(getDR(trkEta_[trkIter], trkPhi_[trkIter], jtEta_[algoIter][jtIter], jtPhi_[algoIter][jtIter]) > .3) continue;
	      if(!highPurity_[trkIter]) continue;
	      if(trkPtError_[trkIter]/trkPt_[trkIter] > 0.3) continue;
	      if(TMath::Abs(trkDz1_[trkIter]/trkDzError1_[trkIter]) > 3) continue;
	      if(TMath::Abs(trkDxy1_[trkIter]/trkDxyError1_[trkIter]) > 3) continue;

	      if(TMath::Abs(trkEta_[trkIter]) > 2.4) continue;

	      if(trkPt_[trkIter] > jtMaxTrackInConePtOut_[algoIter]){
		jtMaxTrackInConePtOut_[algoIter] = trkPt_[trkIter];
		jtMaxTrackInConeHPOut_[algoIter] = highPurity_[trkIter];
	      }

	    } 

	    outTree_p[algoIter]->Fill();
	    
	    inputCounts[pthatIter]++;
	  }


	  if(debugMode) std::cout << __LINE__ << std::endl;

	  const int nRecoJets = nJt_[algoIter];
          Bool_t isJetMatched[nRecoJets];

          for(int recoIter = 0; recoIter < nRecoJets; recoIter++){
            isJetMatched[recoIter] = false;
          }

          for(Int_t genJtIter = 0; genJtIter < nGenJt_[algoIter]; genJtIter++){
	    if(genJtEta_[algoIter][genJtIter] > config.GetJtEtaMax()) continue;
            if(genJtEta_[algoIter][genJtIter] < config.GetJtEtaMin()) continue;
            if(genJtPt_[algoIter][genJtIter] < 5.0) continue;
            if(genJtSubID_[algoIter][genJtIter] != 0) continue;
            if(genJtPt_[algoIter][genJtIter] < config.GetJtPtMin() || genJtPt_[algoIter][genJtIter] > config.GetJtPtMax()) continue;

	    if(debugMode) std::cout << __LINE__ << std::endl;

            if(config.GetIsZJet()){
              if(getDR(zLeadLeptEta_p.at(0), zLeadLeptPhi_p.at(0), genJtEta_[algoIter][genJtIter], genJtPhi_[algoIter][genJtIter]) < 0.4) continue;
              if(getDR(zSubleadLeptEta_p.at(0), zSubleadLeptPhi_p.at(0), genJtEta_[algoIter][genJtIter], genJtPhi_[algoIter][genJtIter]) < 0.4) continue;
	      
	      if(config.GetDoZJtDPhiCut()){
                if(!config.PassesZJetDPhiCut(zPhi_p.at(0), genJtPhi_[algoIter][genJtIter])) continue;
              }
            }

            int recoJetPos = -1;
            for(int recoIter = 0; recoIter < nRecoJets; recoIter++){
              if(isJetMatched[recoIter]) continue;

              if(getDR(jtEta_[algoIter][recoIter], jtPhi_[algoIter][recoIter], genJtEta_[algoIter][genJtIter], genJtPhi_[algoIter][genJtIter]) < 0.15){
                recoJetPos = recoIter;
                isJetMatched[recoIter] = true;
                break;
              }
            }

	    

	    RunOutGen_[algoIter] = Run_;
	    LumiOutGen_[algoIter] = Lumi_;
	    EventOutGen_[algoIter] = Event_;

	    hiBinOutGen_[algoIter] = hiBin_;
	    ptHatOutGen_[algoIter] = ptHat_[algoIter];
	    ptHatWeightOutGen_[algoIter] = config.GetPtHatWeight(ptHat_[algoIter]);
	    if(config.GetIsPbPb()){
	      hiBinWeightOutGen_[algoIter] = findNcoll(hiBin_);
	      hiBinWeightNormOutGen_[algoIter] = findNormNcoll(hiBin_, centPos, config.GetCentBins());
	    }
	    else{
	      hiBinWeightOutGen_[algoIter] = 1;
	      hiBinWeightNormOutGen_[algoIter] = 1;
	    }

	    fullWeightOutGen_[algoIter] = ptHatWeightOutGen_[algoIter]*hiBinWeightOutGen_[algoIter];


	    genJtPtOut_[algoIter] = genJtPt_[algoIter][genJtIter];
	    genJtPhiOut_[algoIter] = genJtPhi_[algoIter][genJtIter];
	    if(config.GetIsZJet()) genJtDelZPhiOut_[algoIter] = TMath::Abs(getDPHI(genJtPhi_[algoIter][genJtIter], zPhi_p.at(0)));
	    else genJtDelZPhiOut_[algoIter] = -999;
	    genJtEtaOut_[algoIter] = genJtEta_[algoIter][genJtIter];
	    genJtMatchIndexOut_[algoIter] = genJtMatchIndex_[algoIter][genJtIter];
	    if(genJtMatchIndex_[algoIter][genJtIter] >= 0){
	      genJtMatchRecoPtOut_[algoIter] = jtPt_[algoIter][genJtMatchIndex_[algoIter][genJtIter]];
	      genJtMatchRecoPhiOut_[algoIter] = jtPhi_[algoIter][genJtMatchIndex_[algoIter][genJtIter]];
	      genJtMatchRecoEtaOut_[algoIter] = jtEta_[algoIter][genJtMatchIndex_[algoIter][genJtIter]];
	    }
	    else{
	      genJtMatchRecoPtOut_[algoIter] = -999;
              genJtMatchRecoPhiOut_[algoIter] = -999;
              genJtMatchRecoEtaOut_[algoIter] = -999;
	    }
	    genJtManualMatchIndexOut_[algoIter] = recoJetPos;
	    if(recoJetPos >= 0){
	      genJtManualMatchRecoPtOut_[algoIter] = jtPt_[algoIter][recoJetPos];
	      genJtManualMatchRecoPhiOut_[algoIter] = jtPhi_[algoIter][recoJetPos];
	      genJtManualMatchRecoEtaOut_[algoIter] = jtEta_[algoIter][recoJetPos];
	    }
	    else{
	      genJtManualMatchRecoPtOut_[algoIter] = -999;
	      genJtManualMatchRecoPhiOut_[algoIter] = -999;
	      genJtManualMatchRecoEtaOut_[algoIter] = -999;
	    }

	    genJtMaxPartPtOut_[algoIter] = -999;
	    genJtMaxPartPhiOut_[algoIter] = -999;
	    genJtMaxPartEtaOut_[algoIter] = -999;
	    genJtMaxPartPDGIDOut_[algoIter] = -999;
	    genJtSumPartPtOut_[algoIter] = 0;
	    genJtNPartPtOut_[algoIter] = 0;

	    const int nGenPart_ = genPt_p->size();
	    for(int partIter = 0; partIter < nGenPart_; partIter++){
	      if(getDR(genEta_p->at(partIter), genPhi_p->at(partIter), genJtEta_[algoIter][genJtIter], genJtPhi_[algoIter][genJtIter]) > .3) continue;

	      if(TMath::Abs(genPDG_p->at(partIter)) == 12) continue;
	      if(TMath::Abs(genPDG_p->at(partIter)) == 14) continue;
	      if(TMath::Abs(genPDG_p->at(partIter)) == 16) continue;

	      genJtSumPartPtOut_[algoIter] += genPt_p->at(partIter);
	      genJtNPartPtOut_[algoIter]++;

	      if(genPt_p->at(partIter) > genJtMaxPartPtOut_[algoIter]){
		genJtMaxPartPtOut_[algoIter] = genPt_p->at(partIter);
		genJtMaxPartPhiOut_[algoIter] = genPhi_p->at(partIter);
		genJtMaxPartEtaOut_[algoIter] = genEta_p->at(partIter);
		genJtMaxPartPDGIDOut_[algoIter] = genPDG_p->at(partIter);

	      }
	    }

	    genJtSubIDOut_[algoIter] = genJtSubID_[algoIter][genJtIter];
	    genJtPartFlavOut_[algoIter] = -999;
	    if(genJtMatchIndex_[algoIter][genJtIter] >= 0) genJtPartFlavOut_[algoIter] = refPartFlav_[algoIter][genJtMatchIndex_[algoIter][genJtIter]];
	    genJtManualMatchPartFlavOut_[algoIter] = -999;
	    if(recoJetPos >= 0) genJtManualMatchPartFlavOut_[algoIter] = refPartFlav_[algoIter][recoJetPos];

	    outTreeGen_p[algoIter]->Fill();
	  }

	}
      	
	passingMuPt_p.clear();
	passingMuPhi_p.clear();
	passingMuEta_p.clear();
	
	passingElePt_p.clear();
	passingElePhi_p.clear();
	passingEleEta_p.clear();
	
	zPt_p.clear();
	zPhi_p.clear();
	zEta_p.clear();
	
	zLeadLeptPt_p.clear();
	zLeadLeptPhi_p.clear();
	zLeadLeptEta_p.clear();
	
	zSubleadLeptPt_p.clear();
	zSubleadLeptPhi_p.clear();
	zSubleadLeptEta_p.clear();	
      }
    
      hiTree_p->ResetBranchAddresses();
      genTree_p->ResetBranchAddresses();
      phoTree_p->ResetBranchAddresses();
      for(Int_t iter = 0; iter < nJetAlgo; iter++){
        jetTree_p[iter]->ResetBranchAddresses();
      }

      inFile_p->Close();
    }

    if(config.GetIsZJet()) std::cout << "Finished file, nEvents w/ more than 1 Z: " << nZMoreThanOne << std::endl;

    std::cout << "End of pthat pho stats: " << std::endl;
    std::cout << " nPho: " << nPho << std::endl;
    std::cout << " nTruePho: " << nTruePho << std::endl;
    std::cout << " nGenPho: " << nGenPho << std::endl;
    std::cout << " nPhoAndGenPho: " << nPhoAndGenPho << std::endl;
    std::cout << " nPhoOrGenPho: " << nPhoOrGenPho << std::endl;
    
  }

  std::cout << "Events with good vz: " << goodVzEvents << std::endl;
  std::cout << "Events with good pcoll: " << goodPCollEvents << std::endl;
  std::cout << "Events with good hibin: " << goodHiBinEvents << std::endl;
  std::cout << "Events with good z: " << goodZEvents << std::endl;

  std::cout << "Jets total: " << totalJets << std::endl;
  std::cout << "Jets eta cut: " << jetsEtaCut << std::endl;
  std::cout << "Jets 5pt cut: " << jets5PtCut << std::endl;
  std::cout << "Jets subid cut: " << jetsSubidCut << std::endl;
  std::cout << "Jets low cut: " << jetsLowCut << std::endl;
  std::cout << "Jets hi cut: " << jetsHiCut << std::endl;
  std::cout << "Jets mu cut: " << jetsMuCut << std::endl;
  std::cout << "Jets z cut: " << jetsZCut << std::endl;
  std::cout << "Jets finale: " << finalJets << std::endl;

    
  outFile_p->cd();

  TNamed nPhoName("nPhoName", std::to_string(nPho).c_str());
  nPhoName.Write("", TObject::kOverwrite);

  TNamed nGenPhoName("nGenPhoName", std::to_string(nGenPho).c_str());
  nGenPhoName.Write("", TObject::kOverwrite);

  TNamed nPhoAndGenPhoName("nPhoAndGenPhoName", std::to_string(nPhoAndGenPho).c_str());
  nPhoAndGenPhoName.Write("", TObject::kOverwrite);

  TNamed nPhoOrGenPhoName("nPhoOrGenPhoName", std::to_string(nPhoOrGenPho).c_str());
  nPhoOrGenPhoName.Write("", TObject::kOverwrite);

  for(int iter = 0; iter < nJetAlgo; iter++){
    outTree_p[iter]->Write("", TObject::kOverwrite);
    outTreeBad_p[iter]->Write("", TObject::kOverwrite);
    outTreeGen_p[iter]->Write("", TObject::kOverwrite);

    delete outTree_p[iter];
    delete outTreeBad_p[iter];
    delete outTreeGen_p[iter];
  }
  config.WriteConfigParamsToRootFile(outFile_p);

  outFile_p->Close();
  delete outFile_p;

  return 0;
}


int main(int argc, char *argv[])
{
  if(argc != 2){
    std::cout << "Usage: makeJECTree_Prototype.exe <inConfigFile>" << std::endl;
    std::cout << "Number of args given: " << argc << std::endl;
    for(int iter = 0; iter < argc; iter++){
      std::cout << "  argv[" << iter << "]: " << argv[iter] << std::endl;
    }
    return -1;
  }

  return makeJECTree_Prototype(argv[1]);
}
