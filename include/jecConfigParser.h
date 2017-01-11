#ifndef JECCONFIGPARSER_H
#define JECCONFIGPARSER_H

//c standard
#include <iostream> 
#include <string>
#include <vector>
#include <fstream>

//root headers
#include "TFile.h"
#include "TNamed.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TMath.h"

//include headers
#include "include/doGlobalDebug.h"
#include "include/checkMakeDir.h"
#include "include/getLogBins.h"
#include "include/getLinBins.h"
#include "include/returnRootFileContentsList.h"
#include "include/etaPhiFunc.h"


class jecConfigParser{
 private:
  const static unsigned int nEventTypes = 3;
  const std::string validEventTypes[nEventTypes] = {"DIJET", "GAMMAJET", "ZJET"};
  const static unsigned int nMaxPtHat = 10;
  const std::string txtStr = ".txt";
  const std::string rootStr = ".root";
  const std::string numStr = "0123456789";
  const std::string commaStr = ",";
  const std::string dotStr = ".";
  const std::string minusStr = "-";
  const std::string alphaUpperStr = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  const std::string alphaLowerStr = "abcdefghijklmnopqrstuvwxyz";
  const static unsigned int nValidTrueFalse = 2;
  const std::string validTrue[nValidTrueFalse] = {"true", "1"};
  const std::string validFalse[nValidTrueFalse] = {"false", "0"};


  const static unsigned int nValidConfigVals = 69;
  enum configIter {EVENTTYPE, //0
		   OUTNAME, //1
		   NPTHAT, //2
		   PTHAT, //3
		   INPUT, //4
		   ISPBPB, //5
		   JETTYPES, //6
		   NJTPTBINS, //7
		   JTPTMIN, //8
		   JTPTMAX, //9
		   DOJTPTLOGBINS, //10
		   DOJTPTCUSTOMBINS, //11
		   JTPTCUSTOMBINS, //12
		   NJTETABINS, //13
		   JTETAMIN, //14
		   JTETAMAX, //15
		   DOJTETACUSTOMBINS, //16
		   JTETACUSTOMBINS, //17
		   NJTETAPTBINS, //18
		   JTETAPTMIN, //19
		   JTETAPTMAX, //20
		   DOJTETAPTLOGBINS, //21
		   DOJTETAPTCUSTOMBINS, //22
		   JTETAPTCUSTOMBINS, //23
		   NJTPTETABINS, //24
		   JTPTETAMIN, //25
		   JTPTETAMAX, //26
		   DOJTPTETAABS, //27
		   DOJTPTETACUSTOMBINS, //28
		   JTPTETACUSTOMBINS, //29
		   FITACCEPTPROBABILITY, //30
		   DOITERATIVEFIT, //31
		   FITITERATIONS, //32
		   FITITERATIONINTERVAL, //33
		   DOWEIGHTS, //34
		   DOWEIGHTTRUNC, //35
		   PTHATWEIGHTS, //36
		   DOPTHATSTAGGER, //37
		   STAGGEROFFSET, //38
		   KEEPFLAVORLESSJETS, //39
		   SETFLAVORLESSJETSASGLUON, //40
		   DOQGREWEIGHT, //41
		   QGREWEIGHTFILENAME, //42
		   DOCORRECTIONS, //43
		   NCENTBINS, //44
		   CENTBINS, //45
		   MINGAMMAPT, //46
		   MAXGAMMAPT, //47
		   MINGAMMAETA, //48
		   MAXGAMMAETA, //49
		   GAMMAPTHATSTAGGER, //50
		   DOGAMMAJTDPHICUT, // 51
		   GAMMAJTDPHICUT, //52
		   MINLEPTONPT, //53
		   MINLEPTONETA, //54
		   MAXLEPTONETA, //55
		   DOLEPTONFLAVORCUTS, //56
		   MINELECTRONPT, //57
		   MINELECTRONETA, //58
		   MAXELECTRONETA, //59
		   MINMUONPT, //60
		   MINMUONETA, //61
		   MAXMUONETA, //62
		   MINZPT, //63
		   MAXZPT, //64
		   MINZM, //65
		   MAXZM, //66
		   DOZJTDPHICUT, //67
		   ZJTDPHICUT}; //68
 
  const std::string validConfigVals[nValidConfigVals] = {"EVENTTYPE", //0
							 "OUTNAME", //1
							 "NPTHAT", //2
							 "PTHAT", //3
							 "INPUT", //4
							 "ISPBPB", //5
							 "JETTYPES", //6
							 "NJTPTBINS", //7
							 "JTPTMIN", //8
							 "JTPTMAX", //9
							 "DOJTPTLOGBINS", //10
							 "DOJTPTCUSTOMBINS", //11
							 "JTPTCUSTOMBINS", //12
							 "NJTETABINS", //13
							 "JTETAMIN", //14
							 "JTETAMAX", //15
							 "DOJTETACUSTOMBINS", //16
							 "JTETACUSTOMBINS", //17
							 "NJTETAPTBINS", //18
							 "JTETAPTMIN", //19
							 "JTETAPTMAX", //20
							 "DOJTETAPTLOGBINS", //21
							 "DOJTETAPTCUSTOMBINS", //22
							 "JTETAPTCUSTOMBINS", //23
							 "NJTPTETABINS", //24
							 "JTPTETAMIN", //25
							 "JTPTETAMAX", //26
							 "DOJTPTETAABS", //27
							 "DOJTPTETACUSTOMBINS", //28
							 "JTPTETACUSTOMBINS", //29
							 "FITACCEPTPROBABILITY", //30
							 "DOITERATIVEFIT", //31
							 "FITITERATIONS", //32
							 "FITITERATIONINTERVAL", //33
							 "DOWEIGHTS", //34
							 "DOWEIGHTTRUNC", //35
							 "PTHATWEIGHTS", //36
							 "DOPTHATSTAGGER", //37
							 "STAGGEROFFSET", //38
							 "KEEPFLAVORLESSJETS", //39
							 "SETFLAVORLESSJETSASGLUON", //40
							 "DOQGREWEIGHT", //41
							 "QGREWEIGHTFILENAME", //42
							 "DOCORRECTIONS", //43
							 "NCENTBINS", //44
							 "CENTBINS", //45
							 "MINGAMMAPT", //46
							 "MAXGAMMAPT", //47
							 "MINGAMMAETA", //48
							 "MAXGAMMAETA", //49
							 "GAMMAPTHATSTAGGER", //50
							 "DOGAMMAJTDPHICUT", //51
							 "GAMMAJTDPHICUT", //52
							 "MINLEPTONPT", //53
							 "MINLEPTONETA", //54
							 "MAXLEPTONETA", //55
							 "DOLEPTONFLAVORCUTS", //56
							 "MINELECTRONPT", //57
							 "MINELECTRONETA", //58
							 "MAXELECTRONETA", //59
							 "MINMUONPT", //60
							 "MINMUONETA", //61
							 "MAXMUONETA", //62
							 "MINZPT", //63
							 "MAXZPT", //64
							 "MINZM", //65
							 "MAXZM", //66
		       					 "DOZJTDPHICUT", //67
							 "ZJTDPHICUT"}; //68


  const std::string configTypes[nValidConfigVals] = {"std::string", //0
						     "std::string", //1
						     "unsigned int", //2
						     "std::vector<unsigned int>", //3
						     "std::vector<std::string>", //4
						     "bool", //5
						     "std::string", //6
						     "unsigned int", //7
						     "unsigned float", //8
						     "unsigned float", //9
						     "bool", //10
						     "bool", //11
						     "std::vector<float>", //12
						     "unsigned int", //13
						     "float", //14
						     "float", //15
						     "bool", //16
						     "std::vector<float>" //17
						     "unsigned int", //18
						     "unsigned float", //19
						     "unsigned float", //20
						     "bool", //21
						     "bool", //22
						     "std::vector<float>", //23
						     "unsigned int", //24
						     "float", //25
						     "float", //26
						     "bool", //27
						     "bool", //28
						     "std::vector<float>", //29
						     "unsigned float", //30
						     "bool", //31
						     "unsigned int", //32
						     "unsigned float", //33
						     "bool", //34
						     "bool", //35
						     "std::vector<float>", //36
						     "bool", //37
						     "unsigned float", //38
						     "bool", //39
						     "bool", //40
						     "bool", //41
						     "std::string", //42
						     "bool", //43
						     "unsgined int", //44
						     "std::vector<unsigned int>", //45
						     "unsigned float", //46
						     "unsigned float", //47
						     "float", //48
						     "float", //49
						     "unsigned float", //50
						     "bool", //51
						     "unsigned float", //52
						     "unsigned float", //53
						     "float", //54
						     "float", //55
						     "bool", //56
						     "unsigned float", //57
						     "float", //58
						     "float", //59
						     "unsigned float", //60
						     "float", //61
						     "float", //62
						     "unsigned float", //63
						     "unsigned float", //64
						     "unsigned float", //65
						     "unsigned float", //66
						     "bool", //67
						     "unsigned float"}; //68

  const std::string defaultConfigInputs[nValidConfigVals] = {"", //0
							     "", //1
							     "0", //2
							     "", //3
							     "", //4
							     "FALSE", //5
							     "", //6
							     "0", //7
							     "30", //8
							     "100", //9
							     "FALSE", //10
							     "FALSE", //11
							     "", //12
							     "16", //13
							     "-1.6", //14
							     "1.6", //15
							     "FALSE", //16
							     "", //17
							     "3", //18
							     "35", //19
							     "10000", //20
							     "FALSE", //21
							     "FALSE", //22
							     "", //23
							     "3", //24
							     "0.0", //25
							     "1.6", //26
							     "TRUE", //27
							     "FALSE", //28
							     "", //29
							     "0.05", //30
							     "TRUE", //31
							     "3", //32
							     "0.05", //33
							     "FALSE", //34
							     "FALSE", //35
							     "", //36
							     "TRUE", //37
							     "20", //38
							     "TRUE", //39
							     "FALSE", //40
							     "FALSE", //41
							     "", //42
							     "FALSE", //43
							     "2", //44
							     "100,30,0", //45
							     "40.", //46
							     "10000.", //47
							     "-1.44", //48
							     "1.44", //49
							     "0.", //50
							     "FALSE", //51
							     "7./8.*PI", //52
							     "15", //53
							     "-1.44", //54
							     "1.44", //55
							     "FALSE", //56
							     "15", //57
							     "-1.44", //58
							     "1.44", //59
							     "10", //60
							     "-2.5", //61
							     "2.5", //62
							     "20", //63
							     "10000", //64
							     "50", //65
							     "130", //66
							     "FALSE", //67
							     "7./8.*PI"}; //68

  unsigned int nConfigInputs[nValidConfigVals] = {0, //0
						  0, //1
						  0, //2
						  0, //3
						  0, //4
						  0, //5
						  0, //6
						  0, //7
						  0, //8
						  0, //9
						  0, //10
						  0, //11
						  0, //12
						  0, //13
						  0, //14
						  0, //15
						  0, //16
						  0, //17
						  0, //18
						  0, //19
						  0, //20
						  0, //21
						  0, //22
						  0, //23
						  0, //24
						  0, //25
						  0, //26
						  0, //27
						  0, //28
						  0, //29
						  0, //30
						  0, //31
						  0, //32
						  0, //33
						  0, //34
						  0, //35
						  0, //36
						  0, //37
						  0, //38
						  0, //39
						  0, //40
						  0, //41
						  0, //42
						  0, //43
						  0, //44
						  0, //45
						  0, //46
						  0, //47
						  0, //48
						  0, //49
						  0, //50
						  0, //51
						  0, //52
						  0, //53
						  0, //54
						  0, //55
						  0, //56
						  0, //57
						  0, //58
						  0, //59
						  0, //60
						  0, //61
						  0, //62
						  0, //63
						  0, //64
						  0, //65
						  0, //66
						  0, //67
						  0}; //68
  

  std::string configInputs[nValidConfigVals] = {"", //0
						"", //1
						"", //2
						"", //3
						"", //4
						"", //5
						"", //6
						"", //7
						"", //8
						"", //9
						"", //10
						"", //11
						"", //12
						"", //13
						"", //14
						"", //15
						"", //16
						"", //17
						"", //18
						"", //19
						"", //20
						"", //21
						"", //22
						"", //23
						"", //24
						"", //25
						"", //26
						"", //27
						"", //28
						"", //29
						"", //30
						"", //31
						"", //32
						"", //33
						"", //34
						"", //35
						"", //36
						"", //37
						"", //38
						"", //39
						"", //40
						"", //41
						"", //42
						"", //43
						"", //44
						"", //45
						"", //46
						"", //47
						"", //48
						"", //49
						"", //50
						"", //51
						"", //52
						"", //53
						"", //54
						"", //55
						"", //56
						"", //57
						"", //58
						"", //59
						"", //60
						"", //61
						"", //62
						"", //63
						"", //64
						"", //65
						"", //66
						"", //67
						""}; //68

  std::string configFileName = "";
  std::string eventTypeStr = "";
  std::string outNameStr = "";
  bool isDijet = false;
  bool isGammaJet = false;
  bool isZJet = false;

  unsigned int nPthats = 0;
  std::vector<unsigned int> pthats;
  std::vector<std::string> inputStrings;
  std::vector<unsigned int> inputFilePtHats;
  std::vector<std::string> inputFileList;
  bool isPbPb = false;

  std::string jetTypes = "";
  std::vector<std::string> jetTypesKeep;
  std::vector<std::string> jetTypesRemove;
  std::vector<std::string> jetTypesFinal;

  unsigned int nJtPtBins = 0;
  float jtPtMin = 30.;
  float jtPtMax = 100.;
  bool doJtPtLogBins = false;
  bool doJtPtCustomBins = false;
  std::vector<float> jtPtCustomBins;

  unsigned int nJtEtaBins = 16;
  float jtEtaMin = -1.6;
  float jtEtaMax = 1.6;
  bool doJtEtaCustomBins = false;
  std::vector<float> jtEtaCustomBins;

  unsigned int nJtEtaPtBins = 3;
  float jtEtaPtMin = 35;
  float jtEtaPtMax = 10000;
  bool doJtEtaPtLogBins = false;
  bool doJtEtaPtCustomBins = false;
  std::vector<float> jtEtaPtCustomBins;

  unsigned int nJtPtEtaBins = 3;
  float jtPtEtaMin = 0.0;
  float jtPtEtaMax = 1.6;
  bool doJtPtEtaAbs = true;
  bool doJtPtEtaCustomBins = true;
  std::vector<float> jtPtEtaCustomBins;

  float fitAcceptProbability = 0.05;
  bool doIterativeFit = true;
  unsigned int fitIterations = 3;
  float fitIterationInterval= 0.05;

  bool doWeights = false;
  bool doWeightTrunc = false;
  std::vector<double> pthatWeights;
  bool doPthatStagger = true;
  float staggerOffset = 20.;

  bool keepFlavorlessJets = true;
  bool setFlavorlessJetsAsGluon = false;

  bool doQGReweight = false;
  std::string qgReweightFileName = "";

  bool doCorrections = false;

  unsigned int nCentBins = 2;
  std::vector<unsigned int> centBins = {100, 30, 0};

  float minGammaPt = 40.;
  float maxGammaPt = 10000.;
  float minGammaEta = -1.44;
  float maxGammaEta = 1.44;
  float gammaPtHatStagger = 0.;

  bool doGammaJtDPhiCut = false;
  float gammaJtDPhiCut = 7.*TMath::Pi()/8.;

  float minLeptonPt = 15.;
  float minLeptonEta = -1.44;
  float maxLeptonEta = 1.44;

  bool doLeptonFlavorCuts = false;

  float minElectronPt = 15.;
  float minElectronEta = -1.44;
  float maxElectronEta = 1.44;

  float minMuonPt = 10.;
  float minMuonEta = -2.5;
  float maxMuonEta = 2.5;

  float minZPt = 40.;
  float maxZPt = 10000.;
  float minZM = 50;
  float maxZM = 130.;
  bool doZJtDPhiCut = false;
  float zJtDPhiCut = 7.*TMath::Pi()/8.;

  std::string returnLowerStr(std::string);
  bool isTrueFalseStr(std::string);
  bool parseTrueFalseStr(std::string);

 public:
  jecConfigParser();
  jecConfigParser(const std::string);
  bool StringIsGoodFloat(const std::string);
  bool StringIsGoodUFloat(const std::string);
  bool StringIsGoodInt(const std::string);
  bool StringIsGoodUInt(const std::string);
  bool ProcessUInt(const std::string, unsigned int &, jecConfigParser::configIter);
  bool ProcessUFloat(const std::string, float &, jecConfigParser::configIter);
  bool ProcessFloat(const std::string, float &, jecConfigParser::configIter);
  bool ProcessBool(const std::string, bool &, jecConfigParser::configIter);
  std::string FloatRangeToTitleString(const float, const float);
  std::string FloatRangeToLabelString(const float, const float, const std::string);
  bool SetConfigParser(const std::string);
  void ResetConfigParser();
  std::string GetConfigFileName();
  std::string GetConfigFileNameNoExt();
  std::string GetEventType();
  std::string GetOutName();
  bool GetIsDijet();
  bool GetIsGammaJet();
  bool GetIsZJet();
  unsigned int GetNPthats();
  unsigned int GetNInputs();
  unsigned int GetPthat(const unsigned int);
  void PrintPthats();
  void PrintInputs();
  std::string GetInput(const unsigned int);
  unsigned int GetInputPtHat(const unsigned int);
  unsigned int GetInputPtHatPos(const unsigned int);
  bool GetIsPbPb();
  std::string GetJetTypes();
  void PrintJetTypesKeep();
  void PrintJetTypesRemove();
  void PrintJetTypesFinal();
  std::vector<std::string> GetJetTypesFinal();
  unsigned int GetNJtPtBins();
  float GetJtPtMin();
  float GetJtPtMax();
  bool GetDoJtPtLogBins();
  bool GetDoJtPtCustomBins();
  void FillJtPtCustomBins(Float_t jtPtBinArr[]);
  void FillJtPtCustomBins(Double_t jtPtBinArr[]);
  unsigned int GetNJtEtaBins();
  float GetJtEtaMin();
  float GetJtEtaMax();
  bool GetDoJtEtaCustomBins();
  void FillJtEtaCustomBins(Float_t jtEtaBinArr[]);
  void FillJtEtaCustomBins(Double_t jtEtaBinArr[]);
  unsigned int GetNJtEtaPtBins();
  float GetJtEtaPtMin();
  float GetJtEtaPtMax();
  bool GetDoJtEtaPtLogBins();
  bool GetDoJtEtaPtCustomBins();
  void FillJtEtaPtCustomBins(Float_t jtEtaPtBinArr[]);
  void FillJtEtaPtCustomBins(Double_t jtEtaPtBinArr[]);
  int GetJtEtaPtBinPos(const float);
  unsigned int GetNJtPtEtaBins();
  float GetJtPtEtaMin();
  float GetJtPtEtaMax();
  bool GetDoJtPtEtaAbs();
  bool GetDoJtPtEtaCustomBins();
  void FillJtPtEtaCustomBins(Float_t jtPtEtaBinArr[]);
  void FillJtPtEtaCustomBins(Double_t jtPtEtaBinArr[]);
  int GetJtPtEtaBinPos(const float);
  float GetFitAcceptProbability();
  bool GetDoIterativeFit();
  unsigned int GetFitIterations();
  float GetFitIterationInterval();
  bool GetDoWeights();
  bool GetDoWeightTrunc();
  bool GetDoPthatStagger();
  bool GetKeepFlavorlessJets();
  bool GetSetFlavorlessJetsAsGluon();
  bool GetDoQGReweight();
  std::string GetQGReweightFileName();
  bool GetDoCorrections();
  bool KeepEventGamma(const float, const float);
  bool KeepEventLeptons(const float, const float, const float, const float, const int);
  bool KeepEventZ(const float, const float);
  bool PassesZJetDPhiCut(const float, const float);
  float GetJtWeight(const unsigned int, const float, const float);
  double GetPtHatWeight(const float);
  double GetTruncPtHatWeight(const float, const float);
  unsigned int GetNCentBins();
  std::vector<unsigned int> GetCentBins();
  unsigned int GetCentBinFromPos(const unsigned int);
  unsigned int GetCentBinFromCent(const unsigned int);
  int GetCentBinFromHiBin(const unsigned int);
  void PrintCentBins();
  float GetMinGammaPt();
  float GetMaxGammaPt();
  float GetMinGammaEta();
  float GetMaxGammaEta();
  float GetGammaPtHatStagger();
  bool GetDoGammaJtDPhiCut();
  float GetGammaJtDPhiCut();
  float GetMinLeptonPt();
  float GetMinLeptonEta();
  float GetMaxLeptonEta();
  bool GetDoLeptonFlavorCuts();
  float GetMinElectronPt();
  float GetMinElectronEta();
  float GetMaxElectronEta();
  float GetMinMuonPt();
  float GetMinMuonEta();
  float GetMaxMuonEta();
  float GetMinZPt();
  float GetMaxZPt();
  float GetMinZM();
  float GetMaxZM();
  bool GetDoZJtDPhiCut();
  float GetZJtDPhiCut();


  void WriteConfigParamsToRootFile(TFile*);
};

jecConfigParser::jecConfigParser()
{
  nPthats = 0;
  return;
}

jecConfigParser::jecConfigParser(const std::string inConfigFile)
{
  if(SetConfigParser(inConfigFile)) return;
  else{
    std::cout << "Initializing jecConfigParser failed. nPthats set to 0." << std::endl;
    nPthats = 0;
  }

  return;
}


bool jecConfigParser::StringIsGoodFloat(std::string floatString)
{
  while(floatString.find(" ") != std::string::npos){
    floatString.replace(floatString.find(" "), 1, "");
  }

  if(floatString.size() == 0) return false;

  bool isNum = true; // check if number                                                                            
  for(unsigned int iter = 0; iter < floatString.size(); iter++){
    if(numStr.find(floatString.at(iter)) == std::string::npos && dotStr.find(floatString.at(iter)) == std::string::npos && minusStr.find(floatString.at(iter)) == std::string::npos){
      isNum = false;
      break;
    }
    else if(minusStr.find(floatString.at(iter)) != std::string::npos && iter != 0){
      isNum = false;
      break;
    }
  }
  if(!isNum){
    std::cout << "\'" << floatString << "\' is invalid, must be number from \'" << numStr << "\' or \'" << dotStr << "\' or \'" << minusStr << "\'. If \'" << minusStr << "\', must be first char. return false" << std::endl;
  }

  return isNum;
}

bool jecConfigParser::StringIsGoodUFloat(const std::string floatString)
{
  bool isUFloat = StringIsGoodFloat(floatString);
  if(!isUFloat) return isUFloat;

  if(floatString.find(minusStr) != std::string::npos){
    isUFloat = false;
    std::cout << "\'" << floatString << "\' is invalid, must be number from \'" << numStr << "\' or \'" << dotStr << "\', non-negative. return false" << std::endl;
  }

  return isUFloat;
}

bool jecConfigParser::StringIsGoodInt(const std::string intString)
{
  bool isInt = StringIsGoodFloat(intString);
  if(!isInt) return isInt;

  if(intString.find(dotStr) != std::string::npos){
    std::cout << "\'" << intString << "\' is invalid, must be number from \'" << numStr << "\' or \'" << minusStr << "\'. If \'" << minusStr << "\', must be first char. return false" << std::endl;
    isInt = false;
  }
  return isInt;
}

bool jecConfigParser::StringIsGoodUInt(std::string intString)
{
  bool isUInt = StringIsGoodUFloat(intString);
  if(!isUInt) return isUInt;

  if(intString.find(dotStr) != std::string::npos && intString.find(dotStr) != intString.size()-1){
    std::cout << "\'" << intString << "\' is invalid, must be number from \'" << numStr << "\', non-negative and non-decimal. return false" << std::endl;
    isUInt = false;
  }
  return isUInt;
}


bool jecConfigParser::ProcessUInt(std::string intString, unsigned int& returnInt, jecConfigParser::configIter enumVal)
{
  if(!StringIsGoodUInt(intString)){
    std::cout << validConfigVals[enumVal] << " value \'" << intString << "\' is invalid. Setting to default." << std::endl;
    returnInt = (unsigned int)std::stoi(defaultConfigInputs[enumVal]);
    return false;
  }
  // setting                                                                                                                    
  configInputs[enumVal] = intString;
  returnInt = (unsigned int)std::stoi(intString);
  nConfigInputs[enumVal]++;

  return true;
}


bool jecConfigParser::ProcessUFloat(std::string floatString, float& returnFloat, jecConfigParser::configIter enumVal)
{
  std::string originalFloatString = floatString;
  bool doPiMult = false;
  if(floatString.find("*PI") != std::string::npos){
    floatString.replace(floatString.find("*PI"), 3, "");
    doPiMult = true;
  }

  std::vector<std::string> floatStrings;
  std::vector<std::string> divOrMult;
  int splitIter = 0;

  while(splitIter < (int)floatString.size()){
    if(floatString.substr(splitIter, 1).find("/") != std::string::npos){
      floatStrings.push_back(floatString.substr(0, splitIter));
      divOrMult.push_back("DIV");
      splitIter = 0;
      floatString.replace(0, floatString.find("/")+1, "");
    }
    else if(floatString.substr(splitIter, 1).find("*") != std::string::npos){
      floatStrings.push_back(floatString.substr(0, splitIter));
      divOrMult.push_back("MULT");
      splitIter = 0;
      floatString.replace(0, floatString.find("*")+1, "");
    }
    else splitIter++;
  }
  if(floatString.size() != 0) floatStrings.push_back(floatString);


  for(int iter = 0; iter < (int)floatStrings.size(); iter++){
    if(!StringIsGoodUFloat(floatStrings.at(iter))){
      std::cout << validConfigVals[enumVal] << " value \'" << originalFloatString << "\', specifically \'" << floatStrings.at(iter) << "\', is invalid. Setting to default." << std::endl;
      returnFloat = std::stof(defaultConfigInputs[enumVal]);
      return false;
    }
  }
  // setting                                                                                                                    
  configInputs[enumVal] = originalFloatString;
  returnFloat = std::stof(floatStrings.at(0));
  
  if(floatStrings.size() > 1){
    for(int iter = 1; iter < (int)floatStrings.size(); iter++){
      if(divOrMult.at(iter-1).find("MULT") != std::string::npos) returnFloat *= std::stof(floatStrings.at(iter));
      else if(divOrMult.at(iter-1).find("DIV") != std::string::npos) returnFloat /= std::stof(floatStrings.at(iter));
    }
  }

  if(doPiMult) returnFloat *= TMath::Pi();
  nConfigInputs[enumVal]++;

  return true;
}

bool jecConfigParser::ProcessFloat(std::string floatString, float& returnFloat, jecConfigParser::configIter enumVal)
{
  std::string originalFloatString = floatString;
  bool doPiMult = false;
  if(floatString.find("*PI") != std::string::npos){
    floatString.replace(floatString.find("*PI"), 3, "");
    doPiMult = true;
  }  

  std::vector<std::string> floatStrings;
  std::vector<std::string> divOrMult;
  int splitIter = 0;

  while(splitIter < (int)floatString.size()){
    if(floatString.substr(splitIter, 1).find("/") != std::string::npos){
      floatStrings.push_back(floatString.substr(0, splitIter));
      divOrMult.push_back("DIV");
      splitIter = 0;
      floatString.replace(0, floatString.find("/")+1, "");
    }
    else if(floatString.substr(splitIter, 1).find("*") != std::string::npos){
      floatStrings.push_back(floatString.substr(0, splitIter));
      divOrMult.push_back("MULT");
      splitIter = 0;
      floatString.replace(0, floatString.find("*")+1, "");
    }
    else splitIter++;
  }
  if(floatString.size() != 0) floatStrings.push_back(floatString);

  for(int iter = 0; iter < (int)floatStrings.size(); iter++){
    if(!StringIsGoodFloat(floatStrings.at(iter))){
      std::cout << validConfigVals[enumVal] << " value \'" << originalFloatString << "\', specifically \'" << floatStrings.at(iter) << "\', is invalid. Setting to default." << std::endl;
      returnFloat = std::stof(defaultConfigInputs[enumVal]);
      return false;
    }
  }
  // setting                                                                                                                    
  configInputs[enumVal] = originalFloatString;
  returnFloat = std::stof(floatStrings.at(0));
  
  if(floatStrings.size() > 1){
    for(int iter = 1; iter < (int)floatStrings.size(); iter++){
      if(divOrMult.at(iter-1).find("MULT") != std::string::npos) returnFloat *= std::stof(floatStrings.at(iter));
      else if(divOrMult.at(iter-1).find("DIV") != std::string::npos) returnFloat /= std::stof(floatStrings.at(iter));
    }
  }

  if(doPiMult) returnFloat *= TMath::Pi();
  nConfigInputs[enumVal]++;

  return true;
}

bool jecConfigParser::ProcessBool(std::string boolString, bool& returnBool, jecConfigParser::configIter enumVal)
{
  if(!isTrueFalseStr(boolString)){
    std::cout << validConfigVals[enumVal] << " value \'" << boolString << "\' is invalid. Please give \'TRUE\' or \'FALSE\'. Return false" << std::endl;
    return false;
  }

  configInputs[enumVal] = boolString;
  nConfigInputs[enumVal]++;

  if(parseTrueFalseStr(boolString)) returnBool = true;
  else returnBool = false;

  return true;
}

std::string jecConfigParser::FloatRangeToTitleString(const float rangeMin, const float rangeMax)
{
  std::string negMinString = "";
  if(rangeMin < 0) negMinString = "Neg";

  std::string negMaxString = "";
  if(rangeMax < 0) negMaxString = "Neg";

  Int_t rangeLowInt = std::trunc(TMath::Abs(rangeMin));
  Int_t rangeHiInt = std::trunc(TMath::Abs(rangeMax));

  Int_t rangeLowDec = std::trunc(TMath::Abs(rangeMin*10) - rangeLowInt*10);
  Int_t rangeHiDec = std::trunc(TMath::Abs(rangeMax*10) - rangeHiInt*10);

  std::string rangeString = negMinString + std::to_string(rangeLowInt) + "p" + std::to_string(rangeLowDec) + "to" + negMaxString + std::to_string(rangeHiInt) + "p" + std::to_string(rangeHiDec);

  return rangeString;
}

std::string jecConfigParser::FloatRangeToLabelString(const float rangeMin, const float rangeMax, const std::string inBetweenString)
{
  std::string negMinString = "";
  if(rangeMin < 0) negMinString = "-";

  std::string negMaxString = "";
  if(rangeMax < 0) negMaxString = "-";

  Int_t rangeLowInt = std::trunc(TMath::Abs(rangeMin));
  Int_t rangeHiInt = std::trunc(TMath::Abs(rangeMax));

  Int_t rangeLowDec = std::trunc(TMath::Abs(rangeMin*10) - rangeLowInt*10);
  Int_t rangeHiDec = std::trunc(TMath::Abs(rangeMax*10) - rangeHiInt*10);

  std::string rangeString = negMinString + std::to_string(rangeLowInt) + "." + std::to_string(rangeLowDec) + "<" + inBetweenString + "<" + negMaxString + std::to_string(rangeHiInt) + "." + std::to_string(rangeHiDec);

  return rangeString;
}

bool jecConfigParser::SetConfigParser(const std::string inConfigFile)
{
  std::vector<std::string> lines;
  if(!checkFile(inConfigFile)){
    std::cout << "Input jecConfig txt file, \'" << inConfigFile << "\', is not valid file. return false" << std::endl;
    return false;
  }

  bool fileIsTxt = false;  
  bool fileIsRoot = false;
  
  if(inConfigFile.size() > txtStr.size()){
    if(inConfigFile.substr(inConfigFile.size()-txtStr.size(), txtStr.size()).find(txtStr) != std::string::npos) fileIsTxt = true;
  }

  if(!fileIsTxt){
    if(inConfigFile.size() > rootStr.size()){
      if(inConfigFile.substr(inConfigFile.size()-rootStr.size(), rootStr.size()).find(rootStr) != std::string::npos) fileIsRoot = true;
    }
  }


  if(!fileIsTxt && !fileIsRoot){
    std::cout << "Input jecConfig txt file, \'" << inConfigFile << "\', doesn't end in \'" << txtStr << "\' and not in \'" << rootStr << "\'. return false" << std::endl;
    return false;
  }

  std::string readFile = "";

  if(fileIsRoot){
    std::cout << "Creating dummy text file \'tempOutFile.txt\' for .root reading..." << std::endl;
    std::ofstream tempOutFile("tempOutFile.txt");
    tempOutFile.close();
    tempOutFile.open("tempOutFile.txt", std::ios::app);

    readFile = "tempOutFile.txt";

    TFile* tempInFile_p = new TFile(inConfigFile.c_str(), "READ");
    std::vector<std::string> tempConfigParams = returnRootFileContentsList(tempInFile_p, "TNamed", "configParamsDir");
    
    const Int_t nConfigParams = tempConfigParams.size();

    for(Int_t iter = 0; iter < nConfigParams; iter++){
      TNamed* tempName_p = (TNamed*)tempInFile_p->Get(tempConfigParams.at(iter).c_str());
      std::string tempConfigParam = tempConfigParams.at(iter);
      while(tempConfigParam.find("/") != std::string::npos){
	tempConfigParam.replace(0, tempConfigParam.find("/")+1, "");
      }

      //      std::cout << "Iter " << iter << "/" << nConfigParams << ": " << tempConfigParam << ", " << tempName_p->GetTitle() << std::endl;
      
      if(tempConfigParam.find("ConfigFileName") != std::string::npos) configFileName = tempName_p->GetTitle();
      else if(tempConfigParam.find("INPUT_PTHAT") != std::string::npos){
	tempConfigParam = tempConfigParam.substr(11, tempConfigParam.size()-11);
	tempOutFile << "INPUT=" << tempConfigParam << "," << tempName_p->GetTitle() << std::endl;
      }
      else{
	tempOutFile << tempConfigParam << "=" << tempName_p->GetTitle() << std::endl;
      }
    }

    tempOutFile.close();
  }
  else{
    configFileName = inConfigFile;
    readFile = inConfigFile;
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::ifstream inFile(readFile.c_str());
  std::string tempStr;
  while(std::getline(inFile, tempStr)){
    if(tempStr.size() == 0) continue; //skip the empties
    if(tempStr.substr(0, 1).find("#") != std::string::npos) continue; //comment handling

    while(tempStr.find(" ") != std::string::npos){
      tempStr.replace(tempStr.find(" "), 1, "");
    }

    //skip invalid inputs and print skip
    bool isValid = false;
    if(tempStr.find("=") != std::string::npos){
      std::string tempStr2 = tempStr.substr(0, tempStr.find("="));
      while(tempStr2.find(" ") != std::string::npos){
	tempStr2.replace(tempStr2.find(" "), 1, "");
      }

      for(unsigned int iter = 0; iter < tempStr2.size(); iter++){
	if(alphaLowerStr.find(tempStr2.at(iter)) != std::string::npos) tempStr2.replace(iter, 1, alphaUpperStr.substr(alphaLowerStr.find(tempStr2.at(iter)), 1));
      }

      for(unsigned int iter = 0; iter < nValidConfigVals; iter++){
	if(tempStr2.find(validConfigVals[iter]) != std::string::npos && tempStr2.size() == validConfigVals[iter].size()){
	  isValid = true;
	  break;
	}
      }
    }
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(!isValid){
      std::cout << "Config line \'" << tempStr << "\' is invalid, will be skipped. Please use valid inputs: ";
      for(unsigned int iter = 0; iter < nValidConfigVals; iter++){
	if(iter < nValidConfigVals-1) std::cout << validConfigVals[iter] << ", ";
	else std::cout << validConfigVals[iter] << ", and \'=\'." << std::endl;
      }
      continue;
    }

    lines.push_back(tempStr);
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  //begin setting class params
  const unsigned int nLines = lines.size();  
  unsigned int replaceIter = 0;

  while(replaceIter < nLines){
    tempStr = lines.at(replaceIter);
    std::string valStr = tempStr.substr(tempStr.find("=")+1, tempStr.size() - tempStr.find("=")+1);
    tempStr.replace(tempStr.find("="), tempStr.size() - tempStr.find("="), "");

    while(tempStr.find(" ") != std::string::npos){
      tempStr.replace(tempStr.find(" "), 1, "");
    }

    while(valStr.find(" ") != std::string::npos){
      valStr.replace(valStr.find(" "), 1, "");
    }

    for(unsigned int iter2 = 0; iter2 < tempStr.size(); iter2++){
      if(alphaLowerStr.find(tempStr.at(iter2)) != std::string::npos) tempStr.replace(iter2, 1, alphaUpperStr.substr(alphaLowerStr.find(tempStr.at(iter2)), 1));
    }

    for(unsigned int iter2 = 0; iter2 < valStr.size(); iter2++){
      if(alphaLowerStr.find(valStr.at(iter2)) != std::string::npos) valStr.replace(iter2, 1, alphaUpperStr.substr(alphaLowerStr.find(valStr.at(iter2)), 1));
    }

    int valPos = -1;
    for(unsigned iter2 = 0; iter2 < nValidConfigVals; iter2++){
      if(valStr.find(validConfigVals[iter2]) != std::string::npos){
	valPos = iter2;
	break;
      }
    }

    if(valPos == -1){
      replaceIter++;
      continue;
    }

    int tempPos = -1;
    for(unsigned iter2 = 0; iter2 < nValidConfigVals; iter2++){
      if(tempStr.size() == validConfigVals[iter2].size() && tempStr.find(validConfigVals[iter2]) != std::string::npos){
        tempPos = iter2;
        break;
      }
    }

    if(valPos == tempPos){
      std::cout << "Request to input value at \'" << validConfigVals[valPos] << "\' for line \'" << lines.at(replaceIter) << "\' is self-referential. Fix config, return empty." << std::endl;
      ResetConfigParser();
      return false;
    }
    else{
      std::string newValStr = "";
      for(unsigned iter2 = 0; iter2 < nLines; iter2++){
	std::string tempStr2 = lines.at(iter2);
	std::string valStr2 = tempStr2.substr(tempStr2.find("=")+1, tempStr2.size() - tempStr2.find("=")+1);
	tempStr2.replace(tempStr2.find("="), tempStr2.size() - tempStr2.find("="), "");

	while(tempStr2.find(" ") != std::string::npos){
	  tempStr2.replace(tempStr2.find(" "), 1, "");
	}

	while(valStr2.find(" ") != std::string::npos){
	  valStr2.replace(valStr2.find(" "), 1, "");
	}

	for(unsigned int iter3 = 0; iter3 < tempStr2.size(); iter3++){
	  if(alphaLowerStr.find(tempStr2.at(iter3)) != std::string::npos) tempStr2.replace(iter3, 1, alphaUpperStr.substr(alphaLowerStr.find(tempStr2.at(iter3)), 1));
	}

	if(tempStr2.size() == validConfigVals[valPos].size() && tempStr2.find(validConfigVals[valPos]) != std::string::npos){
	  newValStr = valStr2;
	  break;
	}
      }

      if(newValStr.size() == 0){
	std::cout << "Request to input value at \'" << validConfigVals[valPos] << "\' for line \'" << lines.at(replaceIter) << "\' is invalid. \'" << validConfigVals[valPos] << "\' not specified in config file. Fix config, return empty." << std::endl;
	ResetConfigParser();
	return false;
      }

      std::string replaceStr = lines.at(replaceIter);
      replaceStr.replace(replaceStr.find(validConfigVals[valPos]), validConfigVals[valPos].size(), newValStr);

      std::cout << "Replacing \'" << lines.at(replaceIter) << "\' with \'" << replaceStr << "\'." << std::endl;
      lines.at(replaceIter) = replaceStr;
    }
  }


  for(unsigned int iter = 0; iter < nLines; iter++){
    tempStr = lines.at(iter);
    std::string valStr = tempStr.substr(tempStr.find("=")+1, tempStr.size() - tempStr.find("=")+1);
    tempStr.replace(tempStr.find("="), tempStr.size() - tempStr.find("="), "");

    for(unsigned int iter = 0; iter < tempStr.size(); iter++){
      if(alphaLowerStr.find(tempStr.at(iter)) != std::string::npos) tempStr.replace(iter, 1, alphaUpperStr.substr(alphaLowerStr.find(tempStr.at(iter)), 1));
    }

    if(tempStr.substr(0, validConfigVals[EVENTTYPE].size()).find(validConfigVals[EVENTTYPE]) != std::string::npos){
      eventTypeStr = valStr;
      configInputs[EVENTTYPE] = eventTypeStr;
    }

    if(tempStr.substr(0, validConfigVals[OUTNAME].size()).find(validConfigVals[OUTNAME]) != std::string::npos){
      outNameStr = valStr;
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(tempStr.substr(0, validConfigVals[NPTHAT].size()).find(validConfigVals[NPTHAT]) != std::string::npos){
      if(!ProcessUInt(valStr, nPthats, NPTHAT)) continue;

      // finally check its below class cap
      if(nPthats > nMaxPtHat){
	std::cout << tempStr << " value \'" << valStr << "\' is greater than class cap \'" << nMaxPtHat << "\'. Consider raising cap or using fewer pthats. Setting to 0." << std::endl;
	nPthats = 0;
	valStr = "";

	configInputs[NPTHAT] = "";
	nConfigInputs[NPTHAT]--;
      }
    }
    if(tempStr.size() == validConfigVals[PTHAT].size() && tempStr.substr(0, validConfigVals[PTHAT].size()).find(validConfigVals[PTHAT]) != std::string::npos){
      valStr = valStr + ",";
      while(valStr.find(",,") != std::string::npos){
	valStr.replace(valStr.find(",,"), 2, ",");
      }

      std::string tempValStr = valStr;
      bool isGoodString = true;
      while(tempValStr.find(",") != std::string::npos){
	if(!StringIsGoodUInt(tempValStr.substr(0, tempValStr.find(",")))){
	  std::cout << tempStr << " value \'" << valStr << "\', specifically \'" << tempValStr.substr(0, tempValStr.find(",")) << "\' is invalid. Return empty." << std::endl;
	  isGoodString = false;
	  break;
	}
	else tempValStr.replace(0, tempValStr.find(",")+1, "");
      }
      if(!isGoodString) continue;

      // setting
      while(valStr.find(",") != std::string::npos){
	pthats.push_back((unsigned int)std::stoi(valStr.substr(0, valStr.find(","))));
	valStr.replace(0, valStr.find(",")+1, "");
      }
      if(valStr.size() != 0) pthats.push_back((unsigned int)std::stoi(valStr));

      nConfigInputs[PTHAT]++;      
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(tempStr.substr(0, validConfigVals[INPUT].size()).find(validConfigVals[INPUT]) != std::string::npos){
      inputStrings.push_back(valStr);
      nConfigInputs[INPUT]++;      
    }

    if(tempStr.substr(0, validConfigVals[ISPBPB].size()).find(validConfigVals[ISPBPB]) != std::string::npos){
      if(!ProcessBool(valStr, isPbPb, ISPBPB)){
	ResetConfigParser();
	return false;
      }
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(tempStr.substr(0, validConfigVals[JETTYPES].size()).find(validConfigVals[JETTYPES]) != std::string::npos){
      jetTypes = valStr;
      configInputs[JETTYPES] = jetTypes;
      nConfigInputs[JETTYPES]++;

      while(valStr.find("*") != std::string::npos){
	if(valStr.substr(0, 1).find("!") != std::string::npos){
	  std::string tempValStr = valStr.substr(1, valStr.find("*")-1);
	  if(tempValStr.size() != 0) jetTypesRemove.push_back(tempValStr);
	  valStr.replace(0, valStr.find("*")+1, "");
	}
	else{
	  std::string tempValStr = valStr.substr(0, valStr.find("*"));
	  if(tempValStr.size() != 0) jetTypesKeep.push_back(tempValStr);
	  valStr.replace(0, valStr.find("*")+1, "");
	}
      }

      if(valStr.size() != 0){
	if(valStr.substr(0, 1).find("!") != std::string::npos){
	  std::string tempValStr = valStr.substr(1, valStr.find("*")-1);
	  if(tempValStr.size() != 0) jetTypesRemove.push_back(tempValStr);
	}
	else{
	  std::string tempValStr = valStr.substr(0, valStr.find("*"));
          if(tempValStr.size() != 0) jetTypesKeep.push_back(tempValStr);
	}
      }      
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(tempStr.substr(0, validConfigVals[NJTPTBINS].size()).find(validConfigVals[NJTPTBINS]) != std::string::npos){
      if(!ProcessUInt(valStr, nJtPtBins, NJTPTBINS)) continue;
    }

    if(tempStr.substr(0, validConfigVals[JTPTMIN].size()).find(validConfigVals[JTPTMIN]) != std::string::npos){
      if(!ProcessUFloat(valStr, jtPtMin, JTPTMIN)) continue;
    }

    if(tempStr.substr(0, validConfigVals[JTPTMAX].size()).find(validConfigVals[JTPTMAX]) != std::string::npos){
      if(!ProcessUFloat(valStr, jtPtMax, JTPTMAX)) continue;
    }

    if(tempStr.substr(0, validConfigVals[DOJTPTLOGBINS].size()).find(validConfigVals[DOJTPTLOGBINS]) != std::string::npos){
      if(!ProcessBool(valStr, doJtPtLogBins, DOJTPTLOGBINS)){
	doJtPtLogBins = false;
	continue;
      }
    }  

    if(tempStr.substr(0, validConfigVals[DOJTPTCUSTOMBINS].size()).find(validConfigVals[DOJTPTCUSTOMBINS]) != std::string::npos){
      if(!ProcessBool(valStr, doJtPtCustomBins, DOJTPTCUSTOMBINS)){
	doJtPtCustomBins = false;
	continue;
      }
    }  

    if(tempStr.substr(0, validConfigVals[JTPTCUSTOMBINS].size()).find(validConfigVals[JTPTCUSTOMBINS]) != std::string::npos){
      valStr = valStr + ",";
      while(valStr.find(",,") != std::string::npos){
	valStr.replace(valStr.find(",,"), 2, ",");
      }

      std::string tempValStr = valStr;
      bool isGoodString = true;
      while(tempValStr.find(",") != std::string::npos){
        if(!StringIsGoodUFloat(tempValStr.substr(0, tempValStr.find(",")))){
	  std::cout << tempStr << " value \'" << valStr << "\', specifically \'" << tempValStr.substr(0, tempValStr.find(",")) << "\' is invalid. Return empty." << std::endl;
          isGoodString = false;
          break;
        }
        else tempValStr.replace(0, tempValStr.find(",")+1, "");
      }
      if(!isGoodString) continue;

      // setting
      jtPtCustomBins.clear();
      configInputs[JTPTCUSTOMBINS] = valStr;
      while(valStr.find(",") != std::string::npos){
	jtPtCustomBins.push_back(std::stof(valStr.substr(0, valStr.find(","))));
	valStr.replace(0, valStr.find(",")+1, "");
      }
      if(valStr.size() != 0) jtPtCustomBins.push_back(std::stof(valStr));

      nConfigInputs[JTPTCUSTOMBINS]++;      
    }




    if(tempStr.substr(0, validConfigVals[NJTETABINS].size()).find(validConfigVals[NJTETABINS]) != std::string::npos){
      if(!ProcessUInt(valStr, nJtEtaBins, NJTETABINS)) continue;
    }

    if(tempStr.substr(0, validConfigVals[JTETAMIN].size()).find(validConfigVals[JTETAMIN]) != std::string::npos){
      if(!ProcessFloat(valStr, jtEtaMin, JTETAMIN)) continue;
    }

    if(tempStr.substr(0, validConfigVals[JTETAMAX].size()).find(validConfigVals[JTETAMAX]) != std::string::npos){
      if(!ProcessFloat(valStr, jtEtaMax, JTETAMAX)) continue;
    }

    if(tempStr.substr(0, validConfigVals[DOJTETACUSTOMBINS].size()).find(validConfigVals[DOJTETACUSTOMBINS]) != std::string::npos){
      if(!ProcessBool(valStr, doJtEtaCustomBins, DOJTETACUSTOMBINS)){
	doJtEtaCustomBins = false;
	continue;
      }
    }  

    if(tempStr.substr(0, validConfigVals[JTETACUSTOMBINS].size()).find(validConfigVals[JTETACUSTOMBINS]) != std::string::npos){
      valStr = valStr + ",";
      while(valStr.find(",,") != std::string::npos){
	valStr.replace(valStr.find(",,"), 2, ",");
      }

      std::string tempValStr = valStr;
      bool isGoodString = true;
      while(tempValStr.find(",") != std::string::npos){
        if(!StringIsGoodFloat(tempValStr.substr(0, tempValStr.find(",")))){
	  std::cout << tempStr << " value \'" << valStr << "\', specifically \'" << tempValStr.substr(0, tempValStr.find(",")) << "\' is invalid. Return empty." << std::endl;
          isGoodString = false;
          break;
        }
        else tempValStr.replace(0, tempValStr.find(",")+1, "");
      }
      if(!isGoodString) continue;

      // setting
      jtEtaCustomBins.clear();
      configInputs[JTETACUSTOMBINS] = valStr;
      while(valStr.find(",") != std::string::npos){
	jtEtaCustomBins.push_back(std::stof(valStr.substr(0, valStr.find(","))));
	valStr.replace(0, valStr.find(",")+1, "");
      }
      if(valStr.size() != 0) jtEtaCustomBins.push_back(std::stof(valStr));

      nConfigInputs[JTETACUSTOMBINS]++;      
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(tempStr.substr(0, validConfigVals[NJTETAPTBINS].size()).find(validConfigVals[NJTETAPTBINS]) != std::string::npos){
      if(!ProcessUInt(valStr, nJtEtaPtBins, NJTETAPTBINS)) continue;
    }

    if(tempStr.substr(0, validConfigVals[JTETAPTMIN].size()).find(validConfigVals[JTETAPTMIN]) != std::string::npos){
      if(!ProcessUFloat(valStr, jtEtaPtMin, JTETAPTMIN)) continue;
    }

    if(tempStr.substr(0, validConfigVals[JTETAPTMAX].size()).find(validConfigVals[JTETAPTMAX]) != std::string::npos){
      if(!ProcessUFloat(valStr, jtEtaPtMax, JTETAPTMAX)) continue;
    }



    if(tempStr.substr(0, validConfigVals[DOJTETAPTLOGBINS].size()).find(validConfigVals[DOJTETAPTLOGBINS]) != std::string::npos){
      if(!ProcessBool(valStr, doJtEtaPtLogBins, DOJTETAPTLOGBINS)){
	doJtEtaPtLogBins = false;
	continue;
      }
    }  

    if(tempStr.substr(0, validConfigVals[DOJTETAPTCUSTOMBINS].size()).find(validConfigVals[DOJTETAPTCUSTOMBINS]) != std::string::npos){
      if(!ProcessBool(valStr, doJtEtaPtCustomBins, DOJTETAPTCUSTOMBINS)){
	doJtEtaPtCustomBins = false;
	continue;
      }
    }  

    if(tempStr.substr(0, validConfigVals[JTETAPTCUSTOMBINS].size()).find(validConfigVals[JTETAPTCUSTOMBINS]) != std::string::npos){
      valStr = valStr + ",";
      while(valStr.find(",,") != std::string::npos){
	valStr.replace(valStr.find(",,"), 2, ",");
      }

      std::string tempValStr = valStr;
      bool isGoodString = true;
      while(tempValStr.find(",") != std::string::npos){
        if(!StringIsGoodFloat(tempValStr.substr(0, tempValStr.find(",")))){
	  std::cout << tempStr << " value \'" << valStr << "\', specifically \'" << tempValStr.substr(0, tempValStr.find(",")) << "\' is invalid. Return empty." << std::endl;
          isGoodString = false;
          break;
        }
        else tempValStr.replace(0, tempValStr.find(",")+1, "");
      }
      if(!isGoodString) continue;

      // setting
      jtEtaPtCustomBins.clear();
      configInputs[JTETAPTCUSTOMBINS] = valStr;
      while(valStr.find(",") != std::string::npos){
	jtEtaPtCustomBins.push_back(std::stof(valStr.substr(0, valStr.find(","))));
	valStr.replace(0, valStr.find(",")+1, "");
      }
      if(valStr.size() != 0) jtEtaPtCustomBins.push_back(std::stof(valStr));

      nConfigInputs[JTETAPTCUSTOMBINS]++;      
    }



    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(tempStr.substr(0, validConfigVals[NJTPTETABINS].size()).find(validConfigVals[NJTPTETABINS]) != std::string::npos){
      if(!ProcessUInt(valStr, nJtPtEtaBins, NJTPTETABINS)) continue;
    }

    if(tempStr.substr(0, validConfigVals[JTPTETAMIN].size()).find(validConfigVals[JTPTETAMIN]) != std::string::npos){
      if(!ProcessFloat(valStr, jtPtEtaMin, JTPTETAMIN)) continue;
    }

    if(tempStr.substr(0, validConfigVals[JTPTETAMAX].size()).find(validConfigVals[JTPTETAMAX]) != std::string::npos){
      if(!ProcessFloat(valStr, jtPtEtaMax, JTPTETAMAX)) continue;
    }

    if(tempStr.substr(0, validConfigVals[DOJTPTETAABS].size()).find(validConfigVals[DOJTPTETAABS]) != std::string::npos){
      if(!ProcessBool(valStr, doJtPtEtaAbs, DOJTPTETAABS)){
	doJtPtEtaAbs = false;
	continue;
      }
    }

    if(tempStr.substr(0, validConfigVals[DOJTPTETACUSTOMBINS].size()).find(validConfigVals[DOJTPTETACUSTOMBINS]) != std::string::npos){
	if(!ProcessBool(valStr, doJtPtEtaCustomBins, DOJTPTETACUSTOMBINS)){
	  doJtPtEtaCustomBins = false;
	  continue;
	}
      }

    if(tempStr.substr(0, validConfigVals[JTPTETACUSTOMBINS].size()).find(validConfigVals[JTPTETACUSTOMBINS]) != std::string::npos){
      valStr = valStr + ",";
      while(valStr.find(",,") != std::string::npos){
	valStr.replace(valStr.find(",,"), 2, ",");
      }

      std::string tempValStr = valStr;
      bool isGoodString = true;
      while(tempValStr.find(",") != std::string::npos){
        if(!StringIsGoodFloat(tempValStr.substr(0, tempValStr.find(",")))){
	  std::cout << tempStr << " value \'" << valStr << "\', specifically \'" << tempValStr.substr(0, tempValStr.find(",")) << "\' is invalid. Return empty." << std::endl;
          isGoodString = false;
          break;
        }
        else tempValStr.replace(0, tempValStr.find(",")+1, "");
      }
      if(!isGoodString) continue;

      // setting
      jtPtEtaCustomBins.clear();
      configInputs[JTPTETACUSTOMBINS] = valStr;
      while(valStr.find(",") != std::string::npos){
	jtPtEtaCustomBins.push_back(std::stof(valStr.substr(0, valStr.find(","))));
	valStr.replace(0, valStr.find(",")+1, "");
      }
      if(valStr.size() != 0) jtPtEtaCustomBins.push_back(std::stof(valStr));

      nConfigInputs[JTPTETACUSTOMBINS]++;      
    }

    if(tempStr.substr(0, validConfigVals[FITACCEPTPROBABILITY].size()).find(validConfigVals[FITACCEPTPROBABILITY]) != std::string::npos){
      if(!ProcessUFloat(valStr, fitAcceptProbability, FITACCEPTPROBABILITY)) continue;
    }

    if(tempStr.substr(0, validConfigVals[DOITERATIVEFIT].size()).find(validConfigVals[DOITERATIVEFIT]) != std::string::npos){
      if(!ProcessBool(valStr, doIterativeFit, DOITERATIVEFIT)){
	doIterativeFit = false;
	continue;
      }
    }

    if(tempStr.substr(0, validConfigVals[FITITERATIONS].size()).find(validConfigVals[FITITERATIONS]) != std::string::npos){
      if(!ProcessUInt(valStr, fitIterations, FITITERATIONS)) continue;
    }
    
    if(tempStr.substr(0, validConfigVals[FITITERATIONINTERVAL].size()).find(validConfigVals[FITITERATIONINTERVAL]) != std::string::npos){
      if(!ProcessUFloat(valStr, fitIterationInterval, FITITERATIONINTERVAL)) continue;
    }
    

    if(tempStr.substr(0, validConfigVals[DOWEIGHTS].size()).find(validConfigVals[DOWEIGHTS]) != std::string::npos){
      if(!ProcessBool(valStr, doWeights, DOWEIGHTS)){
        doWeights = false;
        continue;
      }
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(tempStr.substr(0, validConfigVals[DOWEIGHTTRUNC].size()).find(validConfigVals[DOWEIGHTTRUNC]) != std::string::npos){
      if(!ProcessBool(valStr, doWeightTrunc, DOWEIGHTTRUNC)){
        doWeightTrunc = false;
        continue;
      }
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(tempStr.substr(0, validConfigVals[PTHATWEIGHTS].size()).find(validConfigVals[PTHATWEIGHTS]) != std::string::npos){
      valStr = valStr + ",";
      while(valStr.find(",,") != std::string::npos){
	valStr.replace(valStr.find(",,"), 2, ",");
      }

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      std::string tempValStr = valStr;
      bool isGoodString = true;
      while(tempValStr.find(",") != std::string::npos){
	std::cout << tempValStr << std::endl;
        if(!StringIsGoodUFloat(tempValStr.substr(0, tempValStr.find(",")))){
	  std::cout << tempStr << " value \'" << valStr << "\', specifically \'" << tempValStr.substr(0, tempValStr.find(",")) << "\' is invalid. Return empty." << std::endl;
          isGoodString = false;
          break;
        }
        else tempValStr.replace(0, tempValStr.find(",")+1, "");
      }
      if(!isGoodString) continue;

      if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

      // setting
      configInputs[PTHATWEIGHTS] = valStr;
      nConfigInputs[PTHATWEIGHTS]++;

      while(valStr.find(",") != std::string::npos){
	pthatWeights.push_back(std::stof(valStr.substr(0, valStr.find(","))));
	valStr.replace(0, valStr.find(",")+1, "");
      }
      if(valStr.size() != 0) pthatWeights.push_back(std::stof(valStr));
    }

    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

    if(tempStr.substr(0, validConfigVals[DOPTHATSTAGGER].size()).find(validConfigVals[DOPTHATSTAGGER]) != std::string::npos){
      if(!ProcessBool(valStr, doPthatStagger, DOPTHATSTAGGER)){
        doPthatStagger = false;
        continue;
      }
    }

    if(tempStr.substr(0, validConfigVals[STAGGEROFFSET].size()).find(validConfigVals[STAGGEROFFSET]) != std::string::npos){
      if(!ProcessUFloat(valStr, staggerOffset, STAGGEROFFSET)) continue;
    }

    if(tempStr.substr(0, validConfigVals[KEEPFLAVORLESSJETS].size()).find(validConfigVals[KEEPFLAVORLESSJETS]) != std::string::npos){
      if(!ProcessBool(valStr, keepFlavorlessJets, KEEPFLAVORLESSJETS)){
        keepFlavorlessJets = true;
        continue;
      }
    }

    if(tempStr.substr(0, validConfigVals[SETFLAVORLESSJETSASGLUON].size()).find(validConfigVals[SETFLAVORLESSJETSASGLUON]) != std::string::npos){
      if(!ProcessBool(valStr, setFlavorlessJetsAsGluon, SETFLAVORLESSJETSASGLUON)){
        setFlavorlessJetsAsGluon = false;
        continue;
      }
    }

    if(tempStr.substr(0, validConfigVals[DOQGREWEIGHT].size()).find(validConfigVals[DOQGREWEIGHT]) != std::string::npos){
      if(!ProcessBool(valStr, doQGReweight, DOQGREWEIGHT)){
        doQGReweight = false;
        continue;
      }
    }

    if(tempStr.substr(0, validConfigVals[QGREWEIGHTFILENAME].size()).find(validConfigVals[QGREWEIGHTFILENAME]) != std::string::npos){
      qgReweightFileName = valStr;
      configInputs[QGREWEIGHTFILENAME] = qgReweightFileName;
      nConfigInputs[QGREWEIGHTFILENAME]++;
    }

    if(tempStr.substr(0, validConfigVals[DOCORRECTIONS].size()).find(validConfigVals[DOCORRECTIONS]) != std::string::npos){
      if(!ProcessBool(valStr, doCorrections, DOCORRECTIONS)){
        doCorrections = false;
        continue;
      }
    }

    if(tempStr.substr(0, validConfigVals[NCENTBINS].size()).find(validConfigVals[NCENTBINS]) != std::string::npos){
      if(!ProcessUInt(valStr, nCentBins, NCENTBINS)) continue;
    }

    if(tempStr.substr(0, validConfigVals[CENTBINS].size()).find(validConfigVals[CENTBINS]) != std::string::npos){
      valStr = valStr + ",";
      while(valStr.find(",,") != std::string::npos){
	valStr.replace(valStr.find(",,"), 2, ",");
      }
      std::string tempValStr = valStr;
      bool isGoodString = true;
      while(tempValStr.find(",") != std::string::npos){
        if(!StringIsGoodUFloat(tempValStr.substr(0, tempValStr.find(",")))){
	  std::cout << tempStr << " value \'" << valStr << "\', specifically \'" << tempValStr.substr(0, tempValStr.find(",")) << "\' is invalid. Return empty." << std::endl;
          isGoodString = false;
          break;
        }
        else tempValStr.replace(0, tempValStr.find(",")+1, "");
      }
      if(!isGoodString) continue;

      // setting
      centBins.clear();
      configInputs[CENTBINS] = valStr;
      while(valStr.find(",") != std::string::npos){
	centBins.push_back((unsigned int)std::stoi(valStr.substr(0, valStr.find(","))));
	valStr.replace(0, valStr.find(",")+1, "");
      }
      if(valStr.size() != 0) centBins.push_back((unsigned int)std::stoi(valStr));

      nConfigInputs[CENTBINS]++;      
    }
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
    
    if(tempStr.substr(0, validConfigVals[GAMMAPTHATSTAGGER].size()).find(validConfigVals[GAMMAPTHATSTAGGER]) != std::string::npos){
      if(!ProcessUFloat(valStr, gammaPtHatStagger, GAMMAPTHATSTAGGER)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MINGAMMAPT].size()).find(validConfigVals[MINGAMMAPT]) != std::string::npos){
      if(!ProcessUFloat(valStr, minGammaPt, MINGAMMAPT)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MAXGAMMAPT].size()).find(validConfigVals[MAXGAMMAPT]) != std::string::npos){
      if(!ProcessUFloat(valStr, maxGammaPt, MAXGAMMAPT)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MINGAMMAETA].size()).find(validConfigVals[MINGAMMAETA]) != std::string::npos){
      if(!ProcessFloat(valStr, minGammaEta, MINGAMMAETA)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MAXGAMMAETA].size()).find(validConfigVals[MAXGAMMAETA]) != std::string::npos){
      if(!ProcessFloat(valStr, maxGammaEta, MAXGAMMAETA)) continue;
    }

    if(tempStr.substr(0, validConfigVals[DOGAMMAJTDPHICUT].size()).find(validConfigVals[DOGAMMAJTDPHICUT]) != std::string::npos){
      if(!ProcessBool(valStr, doGammaJtDPhiCut, DOGAMMAJTDPHICUT)){
        doGammaJtDPhiCut = false;
        continue;
      }
    }

    if(tempStr.substr(0, validConfigVals[GAMMAJTDPHICUT].size()).find(validConfigVals[GAMMAJTDPHICUT]) != std::string::npos){
      if(!ProcessUFloat(valStr, gammaJtDPhiCut, GAMMAJTDPHICUT)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MINLEPTONPT].size()).find(validConfigVals[MINLEPTONPT]) != std::string::npos){
      if(!ProcessUFloat(valStr, minLeptonPt, MINLEPTONPT)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MINLEPTONETA].size()).find(validConfigVals[MINLEPTONETA]) != std::string::npos){
      if(!ProcessFloat(valStr, minLeptonEta, MINLEPTONETA)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MAXLEPTONETA].size()).find(validConfigVals[MAXLEPTONETA]) != std::string::npos){
      if(!ProcessFloat(valStr, maxLeptonEta, MAXLEPTONETA)) continue;
    }

    if(tempStr.substr(0, validConfigVals[DOLEPTONFLAVORCUTS].size()).find(validConfigVals[DOLEPTONFLAVORCUTS]) != std::string::npos){
      if(!ProcessBool(valStr, doLeptonFlavorCuts, DOLEPTONFLAVORCUTS)){
        doLeptonFlavorCuts = false;
        continue;
      }
    }

    if(tempStr.substr(0, validConfigVals[MINELECTRONPT].size()).find(validConfigVals[MINELECTRONPT]) != std::string::npos){
      if(!ProcessUFloat(valStr, minElectronPt, MINELECTRONPT)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MINELECTRONETA].size()).find(validConfigVals[MINELECTRONETA]) != std::string::npos){
      if(!ProcessFloat(valStr, minElectronEta, MINELECTRONETA)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MAXELECTRONETA].size()).find(validConfigVals[MAXELECTRONETA]) != std::string::npos){
      if(!ProcessFloat(valStr, maxElectronEta, MAXELECTRONETA)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MINMUONPT].size()).find(validConfigVals[MINMUONPT]) != std::string::npos){
      if(!ProcessUFloat(valStr, minMuonPt, MINMUONPT)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MINMUONETA].size()).find(validConfigVals[MINMUONETA]) != std::string::npos){
      if(!ProcessFloat(valStr, minMuonEta, MINMUONETA)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MAXMUONETA].size()).find(validConfigVals[MAXMUONETA]) != std::string::npos){
      if(!ProcessFloat(valStr, maxMuonEta, MAXMUONETA)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MINZPT].size()).find(validConfigVals[MINZPT]) != std::string::npos){
      if(!ProcessUFloat(valStr, minZPt, MINZPT)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MAXZPT].size()).find(validConfigVals[MAXZPT]) != std::string::npos){
      if(!ProcessUFloat(valStr, maxZPt, MAXZPT)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MINZM].size()).find(validConfigVals[MINZM]) != std::string::npos){
      if(!ProcessUFloat(valStr, minZM, MINZM)) continue;
    }

    if(tempStr.substr(0, validConfigVals[MAXZM].size()).find(validConfigVals[MAXZM]) != std::string::npos){
      if(!ProcessUFloat(valStr, maxZM, MAXZM)) continue;
    }

    if(tempStr.substr(0, validConfigVals[DOZJTDPHICUT].size()).find(validConfigVals[DOZJTDPHICUT]) != std::string::npos){
      if(!ProcessBool(valStr, doZJtDPhiCut, DOZJTDPHICUT)){
        doZJtDPhiCut = false;
        continue;
      }
    }

    if(tempStr.substr(0, validConfigVals[ZJTDPHICUT].size()).find(validConfigVals[ZJTDPHICUT]) != std::string::npos){
      if(!ProcessUFloat(valStr, zJtDPhiCut, ZJTDPHICUT)) continue;
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  bool isGoodEventType = false;
  for(unsigned int iter = 0; iter < nEventTypes; iter++){
    if(eventTypeStr.size() == validEventTypes[iter].size() && eventTypeStr.find(validEventTypes[iter]) != std::string::npos){
      isGoodEventType = true;
      break;
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(!isGoodEventType){
    std::cout << "EVENTTYPE, \'" << eventTypeStr << "\', is not a valid EVENTTYPE. Please pick from: ";
    for(unsigned int iter = 0; iter < nEventTypes; iter++){
      if(iter < nEventTypes-1) std::cout << validEventTypes[iter] << ", ";
      else std::cout << validEventTypes[iter] << ".";
    }
    std::cout << " Return false" << std::endl;

    ResetConfigParser();
    return false;
  }
  else{
    if(eventTypeStr.find("DIJET") != std::string::npos) isDijet = true;
    else if(eventTypeStr.find("GAMMAJET") != std::string::npos) isGammaJet = true;
    else if(eventTypeStr.find("ZJET") != std::string::npos) isZJet = true;
  }

  if(outNameStr.size() == 0){
    std::cout << "OUTNAME \'" << outNameStr << "\' is empty. return false" << std::endl;
    ResetConfigParser();
    return false;
  }
  else if(outNameStr.find(".") != std::string::npos && outNameStr.substr(outNameStr.size()-5, 5).find(".root") == std::string::npos){
    std::cout << "OUTNAME \'" << outNameStr << "\' is not valid string, must end in \'.root\'. return false" << std::endl;
    ResetConfigParser();
    return false;
  }
  else if(outNameStr.substr(outNameStr.size()-5, 5).find(".root") == std::string::npos){
    outNameStr = outNameStr + ".root";
  }
  configInputs[OUTNAME] = outNameStr;

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(nPthats == 0){
    ResetConfigParser();
    return false;
  }
  if(pthats.size() != nPthats){
    std::cout << "NPTHAT, \'" << nPthats << "\', does not match PTHAT size, \'" << pthats.size() << "\'. Return false." << std::endl;
    
    ResetConfigParser();

    return false;
  }
  
  if(inputStrings.size() != nPthats){
    std::cout << "NPTHAT, \'" << nPthats << "\', does not match INPUT size, \'" << inputStrings.size() << "\'. Return false." << std::endl;

    ResetConfigParser();

    return false;
  }

  for(unsigned int iter = 0; iter < nValidConfigVals; iter++){
    if(iter == INPUT) continue;
    else if(nConfigInputs[iter] > 1){
      std::cout << "Input \'" << validConfigVals[iter] << "\' has multiple input values in given config file. Return false." << std::endl;
      
      ResetConfigParser();

      return false;
    }
  }
  

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  unsigned int tempIter = 0;
  while(tempIter < pthats.size()){
    bool doSwap = false;
    
    for(unsigned int iter2 = tempIter+1; iter2 < pthats.size(); iter2++){
      if(pthats.at(tempIter) > pthats.at(iter2)){
	unsigned int tempPthat = pthats.at(tempIter);
	pthats.at(tempIter) = pthats.at(iter2);
	pthats.at(iter2) = tempPthat;

	doSwap = true;
      }
    }

    if(!doSwap) tempIter++;
  }
  
  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(unsigned int pthatIter = 0; pthatIter < nPthats; pthatIter++){
    bool pthatIsFound = false;
    std::string findStr = std::to_string(pthats.at(pthatIter)) + ",";

    for(unsigned int pthatIter2 = 0; pthatIter2 < nPthats; pthatIter2++){
      if(inputStrings.at(pthatIter2).find(findStr) != std::string::npos){
	std::string tempStr = inputStrings.at(pthatIter2);
	inputStrings.at(pthatIter2) = inputStrings.at(pthatIter); 
	inputStrings.at(pthatIter) = tempStr;
	pthatIsFound = true;
      }
    }

    if(!pthatIsFound){
      std::cout << "No INPUT found for PTHAT \'" << pthats.at(pthatIter) << "\'. The following are given:" << std::endl;
      PrintInputs();
      std::cout << "Please give input in form \'PTHAT,FILENAME\'. Return false" << std::endl;
      
      ResetConfigParser();

      return false;
    }
  }

  for(unsigned int pthatIter = 0; pthatIter < nPthats; pthatIter++){
    std::cout << "Checking inputs... " << pthatIter+1 << "/" << nPthats << std::endl;

    inputStrings.at(pthatIter).replace(0, inputStrings.at(pthatIter).find(",")+1, "");
    
    std::string tempInputStr = inputStrings.at(pthatIter);

    while(tempInputStr.size() != 0){

      std::string tempInputStr2 = "";
      if(tempInputStr.find(",") != std::string::npos){
	tempInputStr2 = tempInputStr.substr(0, tempInputStr.find(","));
	tempInputStr.replace(0, tempInputStr.find(",")+1, "");
      }
      else{
	tempInputStr2 = tempInputStr;
	tempInputStr = "";
      }

      inputFileList.push_back(tempInputStr2);
      inputFilePtHats.push_back(pthats.at(pthatIter));

      if(!checkFile(tempInputStr2)){
	std::cout << "INPUT \'" << tempInputStr2 << "\' given for PTHAT==" << pthats.at(pthatIter) << " is not a valid file. Return false" << std::endl;
	
	ResetConfigParser();
	
	return false;
      }
      else{
	if(pthatIter == 0){
	  TFile* tempJetInFile_p = new TFile(tempInputStr2.c_str(), "READ");
	  jetTypesFinal = returnRootFileContentsList(tempJetInFile_p, "TTree", "JetAnalyzer");
	  tempJetInFile_p->Close();
	  delete tempJetInFile_p;
	  
	  unsigned int jetIterFinal = 0;
	  while(jetTypesFinal.size() > jetIterFinal){
	    bool keepBool = true;
	    
	    for(unsigned int jetIterKeep = 0; jetIterKeep < jetTypesKeep.size(); jetIterKeep++){
	      if(jetTypesFinal.at(jetIterFinal).find(jetTypesKeep.at(jetIterKeep)) == std::string::npos){
		jetTypesFinal.erase(jetTypesFinal.begin()+jetIterFinal);
		keepBool = false;
		break;
	      }
	    }

	    for(unsigned int jetIterRemove = 0; jetIterRemove < jetTypesRemove.size(); jetIterRemove++){
	      if(jetTypesFinal.at(jetIterFinal).find(jetTypesRemove.at(jetIterRemove)) != std::string::npos){
		jetTypesFinal.erase(jetTypesFinal.begin()+jetIterFinal);
		keepBool = false;
		break;
	      }
	    }
	    
	    if(keepBool) jetIterFinal++;
	  }


	  if(jetTypesFinal.size() == 0){
	    std::cout << "Given JETTYPES, \'" << jetTypes << "\', return no collections in file \'" << tempInputStr2 << "\'. return false." << std::endl;
	    
	    ResetConfigParser();
	  
	    return false;
	  }
	}
	else{
	  std::vector<std::string> tempJetTypes;
	  
	  TFile* tempJetInFile_p = new TFile(tempInputStr2.c_str(), "READ");
	  tempJetTypes = returnRootFileContentsList(tempJetInFile_p, "TTree", "JetAnalyzer");
	  tempJetInFile_p->Close();
	  delete tempJetInFile_p;
	  
	  unsigned int jetIterFinal = 0;
	  while(jetTypesFinal.size() > jetIterFinal){
	    bool keepBool = false;
	    
	    for(unsigned int jetIterTemp = 0; jetIterTemp < tempJetTypes.size(); jetIterTemp++){
	      if(jetTypesFinal.at(jetIterFinal).find(tempJetTypes.at(jetIterTemp)) != std::string::npos && jetTypesFinal.at(jetIterFinal).size() == tempJetTypes.at(jetIterTemp).size()){
		keepBool = true;
	      break;
	      }
	    }
	    
	    if(keepBool) jetIterFinal++;
	    else{
	      jetTypesFinal.erase(jetTypesFinal.begin()+jetIterFinal);
	    }
	  }
	
	  if(jetTypesFinal.size() == 0){
	    std::cout << "Given JETTYPES, \'" << jetTypes << "\', return no collections in file \'" << tempInputStr2 << "\', after filtering on first file. return false." << std::endl;

	    ResetConfigParser();
	    
	    return false;
	  }
	}
      }
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(0 == nJtPtBins){
    std::cout << "NJTPTBINS, \'" << nJtPtBins << "\', is not a valid value. Return false." << std::endl;

    ResetConfigParser();

    return false;
  }

  if(jtPtMin >= jtPtMax){
    std::cout << "JTPTMIN, \'" << jtPtMin << "\', is greater than or equal to JTPTMAX, \'" << jtPtMax << "\'. Return false" << std::endl;

    ResetConfigParser();
    return false;
  }

  if(doJtPtCustomBins && doJtPtLogBins){
    std::cout << "JTPTCUSTOMBINS AND DOJTPTLOGBINGS both TRUE. Return false" << std::endl;
    
    ResetConfigParser();
    return false;
  }

  if(doJtPtCustomBins){
    if(jtPtCustomBins.size()-1 != nJtPtBins){
      std::cout << "DOJTPTCUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTPTCUSTOMBINS size-1, \'" << jtPtCustomBins.size() << " - 1\' not equal to NJTPTBINS, \'" << nJtPtBins << "\'. Return false" << std::endl;
      
      ResetConfigParser();
      return false;
    }

    bool areCustomBinsOrdered = true;
    
    for(unsigned int jtPtBinIter = 0; jtPtBinIter < nJtPtBins; jtPtBinIter++){
      if(jtPtCustomBins.at(jtPtBinIter) > jtPtCustomBins.at(jtPtBinIter+1)){
	areCustomBinsOrdered = false;
	break;
      }
    }

    if(!areCustomBinsOrdered){
      std::cout << "DOJTPTCUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTPTCUSTOMBINS, \'";
      
      for(unsigned int jtPtBinIter = 0; jtPtBinIter < nJtPtBins; jtPtBinIter++){
	std::cout << jtPtCustomBins.at(jtPtBinIter) << ",";
      }

      std::cout << jtPtCustomBins.at(nJtPtBins) << "\' not ordered. Return false" << std::endl;
      
      ResetConfigParser();
      return false;
    }

    if(jtPtCustomBins.at(0) != jtPtMin){
      std::cout << "DOJTPTCUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTPTCUSTOMBINS min, \'" << jtPtCustomBins.at(0) << "\' not equal to JTPTMIN, \'" << jtPtMin << "\'. Return false" << std::endl;

      ResetConfigParser();
      return false;
    }

    if(jtPtCustomBins.at(nJtPtBins) != jtPtMax){
      std::cout << "DOJTPTCUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTPTCUSTOMBINS max, \'" << jtPtCustomBins.at(nJtPtBins) << "\' not equal to JTPTMAX, \'" << jtPtMax << "\'. Return false" << std::endl;

      ResetConfigParser();
      return false;
    }
  }


  if(jtEtaMin >= jtEtaMax){
    std::cout << "JTETAMIN, \'" << jtEtaMin << "\', is greater than or equal to JTETAMAX, \'" << jtEtaMax << "\'. Return false" << std::endl;

    ResetConfigParser();
    return false;
  }

  if(doJtEtaCustomBins){
    if(jtEtaCustomBins.size()-1 != nJtEtaBins){
      std::cout << "DOJTETACUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTETACUSTOMBINS size-1, \'" << jtEtaCustomBins.size() << " - 1\' not equal to nJtEtaBins, \'" << nJtEtaBins << "\'. Return false" << std::endl;
      
      ResetConfigParser();
      return false;
    }

    bool areCustomBinsOrdered = true;
    
    for(unsigned int jtEtaBinIter = 0; jtEtaBinIter < nJtEtaBins; jtEtaBinIter++){
      if(jtEtaCustomBins.at(jtEtaBinIter) > jtEtaCustomBins.at(jtEtaBinIter+1)){
	areCustomBinsOrdered = false;
	break;
      }
    }

    if(!areCustomBinsOrdered){
      std::cout << "DOJTETACUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTETACUSTOMBINS, \'";
      
      for(unsigned int jtEtaBinIter = 0; jtEtaBinIter < nJtEtaBins; jtEtaBinIter++){
	std::cout << jtEtaCustomBins.at(jtEtaBinIter) << ",";
      }

      std::cout << jtEtaCustomBins.at(nJtEtaBins) << "\' not ordered. Return false" << std::endl;
      
      ResetConfigParser();
      return false;
    }

    if(jtEtaCustomBins.at(0) != jtEtaMin){
      std::cout << "DOJTETACUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTETACUSTOMBINS min, \'" << jtEtaCustomBins.at(0) << "\' not equal to jtEtaMin, \'" << jtEtaMin << "\'. Return false" << std::endl;

      ResetConfigParser();
      return false;
    }

    if(jtEtaCustomBins.at(nJtEtaBins) != jtEtaMax){
      std::cout << "DOJTETACUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTETACUSTOMBINS max, \'" << jtEtaCustomBins.at(nJtEtaBins) << "\' not equal to jtEtaMax, \'" << jtEtaMax << "\'. Return false" << std::endl;

      ResetConfigParser();
      return false;
    }
  }

  if(jtEtaPtMin >= jtEtaPtMax){
    std::cout << "JTETAPTMIN, \'" << jtEtaPtMin << "\', is greater than or equal to JTETAPTMAX, \'" << jtEtaPtMax << "\'. Return false" << std::endl;

    ResetConfigParser();
    return false;
  }

  
  if(doJtEtaPtCustomBins && doJtEtaPtLogBins){
    std::cout << "JTETAPTCUSTOMBINS AND DOJTETAPTLOGBINGS both TRUE. Return false" << std::endl;

    ResetConfigParser();
    return false;
  }

  if(doJtEtaPtCustomBins){
    if(jtEtaPtCustomBins.size()-1 != nJtEtaPtBins){
      std::cout << "DOJTETAPTCUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTETAPTCUSTOMBINS size-1, \'" << jtEtaPtCustomBins.size() << " - 1\' not equal to nJtEtaPtBins, \'" << nJtEtaPtBins << "\'. Return false" << std::endl;
      
      ResetConfigParser();
      return false;
    }

    bool areCustomBinsOrdered = true;
    
    for(unsigned int jtEtaPtBinIter = 0; jtEtaPtBinIter < nJtEtaPtBins; jtEtaPtBinIter++){
      if(jtEtaPtCustomBins.at(jtEtaPtBinIter) > jtEtaPtCustomBins.at(jtEtaPtBinIter+1)){
	areCustomBinsOrdered = false;
	break;
      }
    }

    if(!areCustomBinsOrdered){
      std::cout << "DOJTETAPTCUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTETAPTCUSTOMBINS, \'";
      
      for(unsigned int jtEtaPtBinIter = 0; jtEtaPtBinIter < nJtEtaPtBins; jtEtaPtBinIter++){
	std::cout << jtEtaPtCustomBins.at(jtEtaPtBinIter) << ",";
      }

      std::cout << jtEtaPtCustomBins.at(nJtEtaPtBins) << "\' not ordered. Return false" << std::endl;
      
      ResetConfigParser();
      return false;
    }

    if(jtEtaPtCustomBins.at(0) != jtEtaPtMin){
      std::cout << "DOJTETAPTCUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTETAPTCUSTOMBINS min, \'" << jtEtaPtCustomBins.at(0) << "\' not equal to jtEtaPtMin, \'" << jtEtaPtMin << "\'. Return false" << std::endl;

      ResetConfigParser();
      return false;
    }

    if(jtEtaPtCustomBins.at(nJtEtaPtBins) != jtEtaPtMax){
      std::cout << "DOJTETAPTCUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTETAPTCUSTOMBINS max, \'" << jtEtaPtCustomBins.at(nJtEtaPtBins) << "\' not equal to jtEtaPtMax, \'" << jtEtaPtMax << "\'. Return false" << std::endl;

      ResetConfigParser();
      return false;
    }
  }

  if(jtPtEtaMin >= jtPtEtaMax){
    std::cout << "JTPTETAMIN, \'" << jtPtEtaMin << "\' is greater than or equal to JTPTETAMAX, \'" << jtPtEtaMax << "\'. Return false." << std::endl;

    ResetConfigParser();
    return false;
  }

  if(doJtPtEtaAbs){
    if(jtPtEtaMin < 0){
      std::cout << "JTPTETAMIN, \'" << jtPtEtaMin << "\' is less than zero and DOJTPTETAABS is TRUE. Return false." << std::endl;
      
      ResetConfigParser();
      return false;
    }

    if(jtPtEtaMax < 0){
      std::cout << "JTPTETAMAX, \'" << jtPtEtaMax << "\' is less than zero and DOJTPTETAABS is TRUE. Return false." << std::endl;
      
      ResetConfigParser();
      return false;
    }
  }

  if(doJtPtEtaCustomBins){
    if(jtPtEtaCustomBins.size()-1 != nJtPtEtaBins){
      std::cout << "DOJTPTETACUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTPTETACUSTOMBINS size-1, \'" << jtPtEtaCustomBins.size() << " - 1\' not equal to nJtPtEtaBins, \'" << nJtPtEtaBins << "\'. Return false" << std::endl;
      
      ResetConfigParser();
      return false;
    }

    bool areCustomBinsOrdered = true;
    
    for(unsigned int jtPtEtaBinIter = 0; jtPtEtaBinIter < nJtPtEtaBins; jtPtEtaBinIter++){
      if(jtPtEtaCustomBins.at(jtPtEtaBinIter) > jtPtEtaCustomBins.at(jtPtEtaBinIter+1)){
	areCustomBinsOrdered = false;
	break;
      }
    }

    if(!areCustomBinsOrdered){
      std::cout << "DOJTPTETACUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTPTETACUSTOMBINS, \'";
      
      for(unsigned int jtPtEtaBinIter = 0; jtPtEtaBinIter < nJtPtEtaBins; jtPtEtaBinIter++){
	std::cout << jtPtEtaCustomBins.at(jtPtEtaBinIter) << ",";
      }

      std::cout << jtPtEtaCustomBins.at(nJtPtEtaBins) << "\' not ordered. Return false" << std::endl;
      
      ResetConfigParser();
      return false;
    }

    if(jtPtEtaCustomBins.at(0) != jtPtEtaMin){
      std::cout << "DOJTPTETACUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTPTETACUSTOMBINS min, \'" << jtPtEtaCustomBins.at(0) << "\' not equal to jtPtEtaMin, \'" << jtPtEtaMin << "\'. Return false" << std::endl;

      ResetConfigParser();
      return false;
    }

    if(jtPtEtaCustomBins.at(nJtPtEtaBins) != jtPtEtaMax){
      std::cout << "DOJTPTETACUSTOMBINS == TRUE" << std::endl;
      std::cout << "JTPTETACUSTOMBINS max, \'" << jtPtEtaCustomBins.at(nJtPtEtaBins) << "\' not equal to jtPtEtaMax, \'" << jtPtEtaMax << "\'. Return false" << std::endl;

      ResetConfigParser();
      return false;
    }
  }
  
  if(fitAcceptProbability > 1){
    std::cout << "FITACCEPTPROBABILITY value, \'" << fitAcceptProbability << "\' is greater than 1. Return false." << std::endl;
    ResetConfigParser();
    return false;
  }

  if(doIterativeFit && fitIterations*fitIterationInterval > 1){
    std::cout << "DOITERATIVEFIT == TRUE" << std::endl;
    std::cout << "FITITERATIONS, \'" << fitIterations << "\', times FITITERATIONINTERVAL, \'" << fitIterationInterval << "\', is greater than 1. Return false." << std::endl;
    ResetConfigParser();
    return false;
  }

  if((doWeights && doPthatStagger) || (!doWeights && !doPthatStagger)){
    std::cout << "DOWEIGHTS, \'" << doWeights << "\', and DOPTHATSTAGGER, \'" << doPthatStagger << "\', both have same value. Please choose one. Return false." << std::endl;

    ResetConfigParser();

    return false;
  }
  
  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(doWeights){
    if(pthatWeights.size() != pthats.size()){
      std::cout << "If DOWEIGHTS option selected, PTHATS size, \'" << pthats.size() << "\', must match PTHATWEIGHTS size, \'" << pthatWeights.size() << "\'. return false." << std::endl;
      
      ResetConfigParser();
      
      return false;
    }
    else{
      unsigned int weightPos = 0;
      while(weightPos < pthatWeights.size()){
	bool doSwapWeight = false;

	for(unsigned int weightIter = weightPos+1; weightIter < pthatWeights.size(); weightIter++){
	  if(pthatWeights.at(weightIter) > pthatWeights.at(weightPos)){
	    double tempWeight = pthatWeights.at(weightPos);
	    pthatWeights.at(weightPos) = pthatWeights.at(weightIter);
	    pthatWeights.at(weightIter) = tempWeight;
	    
	    doSwapWeight = true;
	  }
	}
	
	if(!doSwapWeight) weightPos++;
      }
    }
  }
  
  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(0 == nCentBins && isPbPb){
    std::cout << "NCENTBINS, \'" << nCentBins << "\', is not a valid value. Return false." << std::endl;

    ResetConfigParser();

    return false;
  }

  if(centBins.size() != nCentBins+1 && isPbPb){
    std::cout << "NCENTBINS, \'" << nCentBins << "\', is not equal to CENTBINS size \'" << centBins.size() << "\'. Return false." << std::endl;

    ResetConfigParser();

    return false;
  }

  if(isPbPb){
    for(unsigned int centIter = 0; centIter < centBins.size(); centIter++){
      if(centBins.at(centIter) > 100){
	std::cout << "Input centBin, \'" << centBins.at(centIter) << "\', is not in 0-100 range. Return false" << std::endl;

	ResetConfigParser();

	return false;
      }
    }
    unsigned int centPos = 0; 
    while(centPos < centBins.size()){
      Bool_t doSwap = false;
      for(unsigned int centIter = centPos+1; centIter < centBins.size(); centIter++){
	if(centBins.at(centPos) > centBins.at(centIter)){
	  int tempCentBin = centBins.at(centPos);
	  centBins.at(centPos) = centBins.at(centIter);
	  centBins.at(centIter) = tempCentBin;
	}
      }

      if(!doSwap) centPos++;
    }
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(minGammaPt >= maxGammaPt){
    std::cout << "MINGAMMAPT, \'" << minGammaPt << "\', is greater than or equal to MAXGAMMAPT, \'" << maxGammaPt << "\'. Return false" << std::endl;

    ResetConfigParser();
    return false;
  }

  if(minGammaEta >= maxGammaEta){
    std::cout << "MINGAMMAETA, \'" << minGammaEta << "\', is greater than or equal to MAXGAMMAETA, \'" << maxGammaEta << "\'. Return false" << std::endl;

    ResetConfigParser();
    return false;
  }


  if(minLeptonEta >= maxLeptonEta){
    std::cout << "MINLEPTONETA, \'" << minLeptonEta << "\', is greater than or equal to MAXLEPTONETA, \'" << maxLeptonEta << "\'. Return false" << std::endl;

    ResetConfigParser();
    return false;
  }

  if(minElectronEta >= maxElectronEta){
    std::cout << "MINELECTRONETA, \'" << minElectronEta << "\', is greater than or equal to MAXELECTRONETA, \'" << maxElectronEta << "\'. Return false" << std::endl;

    ResetConfigParser();
    return false;
  }

  if(minMuonEta >= maxMuonEta){
    std::cout << "MINMUONETA, \'" << minMuonEta << "\', is greater than or equal to MAXMUONETA, \'" << maxMuonEta << "\'. Return false" << std::endl;

    ResetConfigParser();
    return false;
  }

  if(minZPt >= maxZPt){
    std::cout << "MINZPT, \'" << minZPt << "\', is greater than or equal to MAXZPT, \'" << maxZPt << "\'. Return false" << std::endl;

    ResetConfigParser();
    return false;
  }

  if(minZM >= maxZM){
    std::cout << "MINZM, \'" << minZM << "\', is greater than or equal to MAXZM, \'" << maxZM << "\'. Return false" << std::endl;

    ResetConfigParser();
    return false;
  }

  return true;
}

void jecConfigParser::ResetConfigParser()
{
  configFileName = "";
  eventTypeStr = ""; //0
  outNameStr = ""; //1
  isDijet = false; 
  isGammaJet = false;
  isZJet = false;
  nPthats = 0; //2
  pthats.clear(); //3
  inputStrings.clear(); //4
  inputFileList.clear();
  inputFilePtHats.clear();
  isPbPb = false; //5
  jetTypes = ""; //6
  jetTypesKeep.clear();
  jetTypesRemove.clear();
  jetTypesFinal.clear();
  nJtPtBins = 0; //7
  jtPtMin = 30.; //8
  jtPtMax = 100.; //9
  doJtPtLogBins = false; //10
  doJtPtCustomBins = false; //11
  jtPtCustomBins.clear(); //12
  nJtEtaBins = 16; //13
  jtEtaMin = -1.6; //14
  jtEtaMax = 1.6; //15
  doJtEtaCustomBins = false; //16
  jtEtaCustomBins.clear(); //17
  nJtEtaPtBins = 3; //18
  jtEtaPtMin = 35; //19
  jtEtaPtMax = 10000; //20
  doJtEtaPtLogBins = false; //21
  doJtEtaPtCustomBins = false; //22
  jtEtaPtCustomBins.clear(); //23
  nJtPtEtaBins = 3; //24
  jtPtEtaMin = 0.0; //25
  jtPtEtaMax = 1.6; //26
  doJtPtEtaAbs = false; //27
  doJtPtEtaCustomBins = false; //28
  jtPtEtaCustomBins.clear(); //29
  fitAcceptProbability = 0.05; //30
  doIterativeFit = true; //31
  fitIterations = 3; //32
  fitIterationInterval = 0.05; //33
  doWeights = false; //34
  doWeightTrunc = false; //35
  pthatWeights.clear(); //36
  doPthatStagger = true; //37
  staggerOffset = 20.; //38
  keepFlavorlessJets = true; //39
  setFlavorlessJetsAsGluon = false; //40
  doQGReweight = false; //41
  qgReweightFileName = ""; //42
  doCorrections = false; //43
  nCentBins = 2; //44
  centBins.clear(); //45
  centBins.push_back(100);
  centBins.push_back(30);
  centBins.push_back(0);
  minGammaPt = 40.; //46
  maxGammaPt = 10000.; //47
  minGammaEta = -1.44; //48
  maxGammaEta = 1.44; //49
  gammaPtHatStagger = 0.; //50
  doGammaJtDPhiCut = false; //51
  gammaJtDPhiCut = 7.*TMath::Pi()/8.; //52
  minLeptonPt = 15.; //53
  minLeptonEta = -1.44; //54
  maxLeptonEta = 1.44; //55
  doLeptonFlavorCuts = false; //56
  minElectronPt = 15.; //57
  minElectronEta = -1.44; //58
  maxElectronEta = 1.44; //59
  minMuonPt = 10.; //60
  minMuonEta = -2.5; //61
  maxMuonEta = 2.5; //62
  minZPt = 20.; //63
  maxZPt = 10000.; //64
  minZM = 50.; //65
  maxZM = 130.; //66
  doZJtDPhiCut = false; //67
  zJtDPhiCut = 7.*TMath::Pi()/8.; //68

  for(unsigned int iter = 0; iter < nValidConfigVals; iter++){
    nConfigInputs[iter] = 0;
    configInputs[iter] = defaultConfigInputs[iter];
  }
  
  return;
}

std::string jecConfigParser::GetConfigFileName(){return configFileName;}
std::string jecConfigParser::GetConfigFileNameNoExt()
{
  std::string tempConfigFileName = configFileName;
  while(tempConfigFileName.find(".") != std::string::npos){
    tempConfigFileName.replace(tempConfigFileName.find("."), tempConfigFileName.size() - tempConfigFileName.find("."), "");
  }
  return tempConfigFileName;
}
std::string jecConfigParser::GetEventType(){return eventTypeStr;}
std::string jecConfigParser::GetOutName(){return outNameStr;}
bool jecConfigParser::GetIsDijet(){return isDijet;}
bool jecConfigParser::GetIsGammaJet(){return isGammaJet;}
bool jecConfigParser::GetIsZJet(){return isZJet;}
unsigned int jecConfigParser::GetNPthats(){return nPthats;}
unsigned int jecConfigParser::GetNInputs(){return inputFileList.size();}

unsigned int jecConfigParser::GetPthat(unsigned int pos)
{
  if(pos >= nPthats){
    std::cout << "Requested input position, \'" << pos << "\', is >= nPthats, \'" << nPthats << "\'. Return 0" << std::endl;
    return 0;
  }
  
  return pthats.at(pos);
}


void jecConfigParser::PrintPthats()
{
  std::cout << "Pthats: ";

  for(unsigned int iter = 0; iter < nPthats; iter++){
    if(iter < nPthats-1) std::cout << pthats.at(iter) << ", ";
    else std::cout << pthats.at(nPthats-1) << ".";
  }

  std::cout << std::endl;

  return;
}

void jecConfigParser::PrintInputs()
{
  std::cout << "Inputs:" << std::endl;

  for(unsigned int iter = 0; iter < nPthats; iter++){
    std::cout << " Input " << iter << "/" << nPthats << ": " << inputStrings.at(iter) << std::endl;
  }

  return;
}

std::string jecConfigParser::GetInput(unsigned int pos)
{
  if(pos >= inputFileList.size()){
    std::cout << "Requested input position, \'" << pos << "\', is >= fileList.size(), \'" << inputFileList.size() << "\'. Return empty string" << std::endl;
    return "";
  }
  
  return inputFileList.at(pos);
}

unsigned int jecConfigParser::GetInputPtHat(unsigned int pos)
{
  if(pos >= inputFilePtHats.size()){
    std::cout << "Requested input pthat, \'" << pos << "\', is >= inputFilePtHats.size(), \'" << inputFilePtHats.size() << "\'. Return 10000" << std::endl;
    return 10000;
  }
  
  return inputFilePtHats.at(pos);
}

unsigned int jecConfigParser::GetInputPtHatPos(unsigned int pos)
{
  if(pos >= inputFilePtHats.size()){
    std::cout << "Requested input pthat, \'" << pos << "\', is >= inputFilePtHats.size(), \'" << inputFilePtHats.size() << "\'. Return 10000" << std::endl;
    return 10000;
  }

  unsigned int retPos = 100000;
  
  for(unsigned int iter = 0; iter < pthats.size(); iter++){
    if(pthats.at(iter) == inputFilePtHats.at(pos)) retPos = iter;
  }

  return retPos;
}

bool jecConfigParser::GetIsPbPb(){return isPbPb;}
std::string jecConfigParser::GetJetTypes(){return jetTypes;}
void jecConfigParser::PrintJetTypesKeep()
{
  std::cout << "Jet types to keep filter: ";
  for(unsigned int iter = 0; iter < jetTypesKeep.size()-1; iter++){
    std::cout << "\'" << jetTypesKeep.at(iter) << "\', ";
  }
  if(jetTypesKeep.size() != 0) std::cout << "\'" << jetTypesKeep.at(jetTypesKeep.size()-1) << "\'." << std::endl;
  
  return;
}

void jecConfigParser::PrintJetTypesRemove()
{
  std::cout << "Jet types to keep filter: ";
  for(unsigned int iter = 0; iter < jetTypesRemove.size()-1; iter++){
    std::cout << "\'" << jetTypesRemove.at(iter) << "\', ";
  }
  if(jetTypesRemove.size() != 0) std::cout << "\'" << jetTypesRemove.at(jetTypesRemove.size()-1) << "\'." << std::endl;
  
  return;
}

void jecConfigParser::PrintJetTypesFinal()
{
  std::cout << "Jet types to keep filter: ";
  for(unsigned int iter = 0; iter < jetTypesFinal.size()-1; iter++){
    std::cout << "\'" << jetTypesFinal.at(iter) << "\', ";
  }
  if(jetTypesFinal.size() != 0) std::cout << "\'" << jetTypesFinal.at(jetTypesFinal.size()-1) << "\'." << std::endl;
  
  return;
}

std::vector<std::string> jecConfigParser::GetJetTypesFinal(){return jetTypesFinal;}

unsigned int jecConfigParser::GetNJtPtBins(){return nJtPtBins;}
float jecConfigParser::GetJtPtMin(){return jtPtMin;}
float jecConfigParser::GetJtPtMax(){return jtPtMax;}
bool jecConfigParser::GetDoJtPtLogBins(){return doJtPtLogBins;}
bool jecConfigParser::GetDoJtPtCustomBins(){return doJtPtCustomBins;}

void jecConfigParser::FillJtPtCustomBins(Float_t jtPtBinArr[])
{
  for(unsigned int jtPtBinIter = 0; jtPtBinIter < nJtPtBins+1; jtPtBinIter++){
    jtPtBinArr[jtPtBinIter] = jtPtCustomBins.at(jtPtBinIter);
  }

  return;
}

void jecConfigParser::FillJtPtCustomBins(Double_t jtPtBinArr[])
{
  for(unsigned int jtPtBinIter = 0; jtPtBinIter < nJtPtBins+1; jtPtBinIter++){
    jtPtBinArr[jtPtBinIter] = jtPtCustomBins.at(jtPtBinIter);
  }

  return;
}

unsigned int jecConfigParser::GetNJtEtaBins(){return nJtEtaBins;}
float jecConfigParser::GetJtEtaMin(){return jtEtaMin;}
float jecConfigParser::GetJtEtaMax(){return jtEtaMax;}
bool jecConfigParser::GetDoJtEtaCustomBins(){return doJtEtaCustomBins;}

void jecConfigParser::FillJtEtaCustomBins(Float_t jtEtaBinArr[])
{
  for(unsigned int jtEtaBinIter = 0; jtEtaBinIter < nJtEtaBins+1; jtEtaBinIter++){
    jtEtaBinArr[jtEtaBinIter] = jtEtaCustomBins.at(jtEtaBinIter);
  }

  return;
}

void jecConfigParser::FillJtEtaCustomBins(Double_t jtEtaBinArr[])
{
  for(unsigned int jtEtaBinIter = 0; jtEtaBinIter < nJtEtaBins+1; jtEtaBinIter++){
    jtEtaBinArr[jtEtaBinIter] = jtEtaCustomBins.at(jtEtaBinIter);
  }

  return;
}

unsigned int jecConfigParser::GetNJtEtaPtBins(){return nJtEtaPtBins;}
float jecConfigParser::GetJtEtaPtMin(){return jtEtaPtMin;}
float jecConfigParser::GetJtEtaPtMax(){return jtEtaPtMax;}
bool jecConfigParser::GetDoJtEtaPtLogBins(){return doJtEtaPtLogBins;}
bool jecConfigParser::GetDoJtEtaPtCustomBins(){return doJtEtaPtCustomBins;}

void jecConfigParser::FillJtEtaPtCustomBins(Float_t jtEtaPtBinArr[])
{
  for(unsigned int jtEtaPtBinIter = 0; jtEtaPtBinIter < nJtEtaPtBins+1; jtEtaPtBinIter++){
    jtEtaPtBinArr[jtEtaPtBinIter] = jtEtaPtCustomBins.at(jtEtaPtBinIter);
  }

  return;
}

void jecConfigParser::FillJtEtaPtCustomBins(Double_t jtEtaPtBinArr[])
{
  for(unsigned int jtEtaPtBinIter = 0; jtEtaPtBinIter < nJtEtaPtBins+1; jtEtaPtBinIter++){
    jtEtaPtBinArr[jtEtaPtBinIter] = jtEtaPtCustomBins.at(jtEtaPtBinIter);
  }

  return;
}

int jecConfigParser::GetJtEtaPtBinPos(const float jtPt)
{
  int binPos = -1;

  double tempJtEtaPtBins[nJtEtaPtBins+1];
  if(doJtEtaPtLogBins) getLogBins(jtEtaPtMin, jtEtaPtMax, nJtEtaPtBins, tempJtEtaPtBins);
  else if(doJtEtaPtCustomBins) FillJtEtaPtCustomBins(tempJtEtaPtBins);
  else getLinBins(jtEtaPtMin, jtEtaPtMax, nJtEtaPtBins, tempJtEtaPtBins);

  for(unsigned int iter = 0; iter < nJtEtaPtBins; iter++){
    if(tempJtEtaPtBins[iter] <= jtPt && tempJtEtaPtBins[iter+1] > jtPt){
      binPos = iter;
      break;
    }
  }

  return binPos;
}

unsigned int jecConfigParser::GetNJtPtEtaBins(){return nJtPtEtaBins;}
float jecConfigParser::GetJtPtEtaMin(){return jtPtEtaMin;}
float jecConfigParser::GetJtPtEtaMax(){return jtPtEtaMax;}
bool jecConfigParser::GetDoJtPtEtaAbs(){return doJtPtEtaAbs;}
bool jecConfigParser::GetDoJtPtEtaCustomBins(){return doJtPtEtaCustomBins;}

void jecConfigParser::FillJtPtEtaCustomBins(Float_t jtPtEtaBinArr[])
{
  for(unsigned int jtPtEtaBinIter = 0; jtPtEtaBinIter < nJtPtEtaBins+1; jtPtEtaBinIter++){
    jtPtEtaBinArr[jtPtEtaBinIter] = jtPtEtaCustomBins.at(jtPtEtaBinIter);
  }

  return;
}

void jecConfigParser::FillJtPtEtaCustomBins(Double_t jtPtEtaBinArr[])
{
  for(unsigned int jtPtEtaBinIter = 0; jtPtEtaBinIter < nJtPtEtaBins+1; jtPtEtaBinIter++){
    jtPtEtaBinArr[jtPtEtaBinIter] = jtPtEtaCustomBins.at(jtPtEtaBinIter);
  }

  return;
}

float jecConfigParser::GetFitAcceptProbability(){return fitAcceptProbability;}
bool jecConfigParser::GetDoIterativeFit(){return doIterativeFit;}
unsigned int jecConfigParser::GetFitIterations(){return fitIterations;}
float jecConfigParser::GetFitIterationInterval(){return fitIterationInterval;}

int jecConfigParser::GetJtPtEtaBinPos(const float jtEta)
{
  int binPos = -1;

  double tempJtPtEtaBins[nJtPtEtaBins+1];
  if(doJtPtEtaCustomBins) FillJtPtEtaCustomBins(tempJtPtEtaBins);
  else getLinBins(jtPtEtaMin, jtPtEtaMax, nJtPtEtaBins, tempJtPtEtaBins);

  float tempJtEta = jtEta;
  if(doJtPtEtaAbs) tempJtEta = TMath::Abs(tempJtEta);

  for(unsigned int iter = 0; iter < nJtPtEtaBins; iter++){
    if(tempJtPtEtaBins[iter] <= tempJtEta && tempJtPtEtaBins[iter+1] > tempJtEta){
      binPos = iter;
      break;
    }
  }

  return binPos;
}


bool jecConfigParser::GetDoWeights(){return doWeights;}
bool jecConfigParser::GetDoWeightTrunc(){return doWeightTrunc;}
bool jecConfigParser::GetDoPthatStagger(){return doPthatStagger;}

bool jecConfigParser::GetKeepFlavorlessJets(){return keepFlavorlessJets;}
bool jecConfigParser::GetSetFlavorlessJetsAsGluon(){return setFlavorlessJetsAsGluon;}

bool jecConfigParser::GetDoQGReweight(){return doQGReweight;}
std::string jecConfigParser::GetQGReweightFileName(){return qgReweightFileName;}
bool jecConfigParser::GetDoCorrections(){return doCorrections;}

bool jecConfigParser::KeepEventGamma(const float gammaPt, const float gammaEta)
{
  if(!isGammaJet){
    std::cout << "Warning: jecConfigParser::KeepEventGamma called for ISGAMMAJET == FALSE. Return false." << std::endl;
    return false;
  }

  if(gammaPt < minGammaPt) return false;
  else if(gammaPt > maxGammaPt) return false; 
  else if(gammaEta > maxGammaEta) return false;
  else if(gammaEta < minGammaEta) return false;
  else return true;
}


bool jecConfigParser::KeepEventLeptons(const float leptonPt1, const float leptonEta1, const float leptonPt2, const float leptonEta2, const int flavorVal)
{
  if(doLeptonFlavorCuts && TMath::Abs(flavorVal) != 11 && TMath::Abs(flavorVal) != 13){
    std::cout << "Warning: jecConfigParser::KeepEventLeptons called w/ DOLEPTONFLAVORCUTS == TRUE but invalid flavorVal, \'" << flavorVal << "\'. Please use 11,-11,13,-13. Return false." << std::endl;
    return false;
  }

  if(!isZJet){
    std::cout << "Warning: jecConfigParser::KeepEventLeptons called for ISZJET == FALSE. Return false." << std::endl;
    return false;
  }

  bool keepLeptons = true;

  if(!doLeptonFlavorCuts){
    if(leptonPt1 < minLeptonPt) keepLeptons = false;
    else if(leptonEta1 < minLeptonEta || leptonEta1 > maxLeptonEta) keepLeptons = false;
    else if(leptonPt2 < minLeptonPt) keepLeptons = false;
    else if(leptonEta2 < minLeptonEta || leptonEta2 > maxLeptonEta) keepLeptons = false;
  }
  else if(TMath::Abs(flavorVal) == 11){
    if(leptonPt1 < minElectronPt) keepLeptons = false;
    else if(leptonEta1 < minElectronEta || leptonEta1 > maxElectronEta) keepLeptons = false;
    else if(leptonPt2 < minElectronPt) keepLeptons = false;
    else if(leptonEta2 < minElectronEta || leptonEta2 > maxElectronEta) keepLeptons = false;
  }
  else if(TMath::Abs(flavorVal) == 13){
    if(leptonPt1 < minMuonPt) keepLeptons = false;
    else if(leptonEta1 < minMuonEta || leptonEta1 > maxMuonEta) keepLeptons = false;
    else if(leptonPt2 < minMuonPt) keepLeptons = false;
    else if(leptonEta2 < minMuonEta || leptonEta2 > maxMuonEta) keepLeptons = false;
  }

  return keepLeptons;
}

bool jecConfigParser::KeepEventZ(const float zPt, const float zM)
{
  if(!isZJet){
    std::cout << "Warning: jecConfigParser::KeepEventZ called for ISZJET == FALSE. Return false." << std::endl;
    return false;
  }

  if(zPt < minZPt || zPt > maxZPt) return false;
  else if(zM < minZM || zM > maxZM) return false;
  else return true;
}

bool jecConfigParser::PassesZJetDPhiCut(const float zPhi, const float jtPhi)
{
  if(!isZJet){
    std::cout << "Warning: jecConfigParser::PassesZJetDPhiCut called for ISZJET == FALSE. Return false." << std::endl;
    return false;
  }

  if(!doZJtDPhiCut){
    std::cout << "Warning: jecConfigParser::PassesZJetDPhiCut callsed for DOZJTDPHICUT == FALSE. Return false." << std::endl;
    return false;
  }

  float dPhi = TMath::Abs(getDPHI(zPhi, jtPhi));
  bool retVal = false;
  if(dPhi >= zJtDPhiCut) retVal = true;

  return retVal;
}


//EDIT Requires modification for multiple binnings
float jecConfigParser::GetJtWeight(const unsigned int pthatPos, const float jtPt, const float jtEta)
{
  if(jtPt < jtPtMin || jtPt > jtPtMax) return 0.;

  if(jtEta > jtEtaMax) return 0.;
  if(jtEta < jtEtaMin) return 0.;

  int truePtHatPos = -1;
  for(unsigned int iter = 0; iter < nPthats; iter++){
    if(pthats.at(iter) == inputFilePtHats.at(pthatPos)) truePtHatPos = iter;
  }

  if(truePtHatPos == -1){
    std::cout << "Pthat not found at jecConfigParser line " << __LINE__ << ". Return 0 weight" << std::endl;
    return 0;
  }

  float weight = 1.;

  if(isZJet) return weight;

  /*
  if(truePtHatPos == 0) return 0.999565;
  else if(truePtHatPos == 1) return 0.187379;
  else if(truePtHatPos == 2) return 0.0434138;
  else if(truePtHatPos == 3) return 0.0094069;
  else if(truePtHatPos == 4) return 0.00211448;			      
  */

  int jtPtBinPthatPos[nJtPtBins+1];
  double jtPtBins[nJtPtBins+1];
  if(doJtPtLogBins) getLogBins(jtPtMin, jtPtMax, nJtPtBins, jtPtBins);
  else if(doJtPtCustomBins){
    for(unsigned int jtPtBinIter = 0; jtPtBinIter < nJtPtBins+1; jtPtBinIter++){
      jtPtBins[jtPtBinIter] = jtPtCustomBins.at(jtPtBinIter);
    }
  }
  else getLinBins(jtPtMin, jtPtMax, nJtPtBins, jtPtBins);

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(unsigned int iter = 0; iter < nJtPtBins; iter++){
    jtPtBinPthatPos[iter] = -1;
    
    for(unsigned iter2 = 0; iter2 < nPthats-1; iter2++){
      if(pthats.at(iter2)+staggerOffset <= jtPtBins[iter] && pthats.at(iter2+1)+staggerOffset > jtPtBins[iter]){
	jtPtBinPthatPos[iter] = (int)iter2;
	break;
      }
    }
    if(jtPtBinPthatPos[iter] == -1 && pthats.at(nPthats-1)+staggerOffset <= jtPtBins[iter]) jtPtBinPthatPos[iter] = nPthats-1;
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
  
  int jtPtPos = -1;

  if(doWeights) weight = 1.;
  else if(doPthatStagger){
    for(unsigned int iter = 0; iter < nJtPtBins; iter++){
      if(jtPt > jtPtBins[iter] && jtPt < jtPtBins[iter+1]){
	jtPtPos = iter;
	break;
      }
    }

    if(jtPtBinPthatPos[jtPtPos] < 0) weight = 0.;
    else if(jtPtBinPthatPos[jtPtPos] != truePtHatPos) weight = 0.;
    else weight = 1.;    
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  return weight;
}


double jecConfigParser::GetPtHatWeight(const float inPtHat)
{
  double weight = -1.;

  if(!doWeights) return 1;

  for(unsigned int iter = 0; iter < nPthats-1; iter++){
    if(inPtHat < pthats.at(iter+1)){
      weight = pthatWeights.at(iter);
      break;
    } 
  }
  if(inPtHat >= pthats.at(nPthats-1)) weight = pthatWeights.at(nPthats-1);

  if(weight < 0){
    std::cout << "No pthat found for pthat \'" << inPtHat << "\', weight returned is 0" << std::endl;
    std::cout << "Available pthat bins: ";

    for(unsigned int iter = 0; iter < pthats.size()-1; iter++){
      std::cout << pthats.at(iter) << "-" << pthats.at(iter+1) << ", ";
    }
    std::cout << pthats.at(nPthats-2) << "-" << pthats.at(nPthats-1) << "." << std::endl;

    weight = 0;
  }

  return weight;
}


double jecConfigParser::GetTruncPtHatWeight(const float inPtHat, const float inJtPt)
{
  double weight = -1.;

  if(!doWeights) return 1;

  if(!doWeightTrunc) return GetPtHatWeight(inPtHat);

  int pthatPos = -1;

  for(unsigned int iter = 0; iter < nPthats-1; iter++){
    if(inPtHat < pthats.at(iter+1)){
      pthatPos = iter;
      break;
    } 
  }
  if(inPtHat >= pthats.at(nPthats-1)) pthatPos = nPthats-1;

  int jtPos = -1;

  for(unsigned int iter = 0; iter < nPthats-1; iter++){
    if(inJtPt < pthats.at(iter+1)){
      jtPos = iter;
      break;
    } 
  }
  if(inJtPt >= pthats.at(nPthats-1)) jtPos = nPthats-1;

  if(jtPos - pthatPos >= 1) weight = 0;
  else weight = GetPtHatWeight(inPtHat);

  return weight;
}


unsigned int jecConfigParser::GetNCentBins()
{
  if(isPbPb) return nCentBins;
  else return 1;
}

std::vector<unsigned int> jecConfigParser::GetCentBins()
{
  std::vector<unsigned int> dummyVect;
  if(isPbPb) return centBins;
  else return dummyVect;
}

unsigned int jecConfigParser::GetCentBinFromPos(const unsigned int centPos)
{
  if(!isPbPb) return 0;
  else{
    if(centPos > centBins.size()){
      std::cout << "Requested pos, \'" << centPos << "\', greater than number of bins. return 0." << std::endl;
      return 0;
    }
    else return centBins.at(centPos);
  }
}

unsigned int jecConfigParser::GetCentBinFromCent(const unsigned int cent)
{
  if(cent > 100){
    std::cout << "Cent given, \'" << cent << "\', is not in allowed range 0-100%. Return 0"  << std::endl;
    return 0;
  }
  
  unsigned int centPos = 0;

  if(isPbPb){
    for(unsigned int centIter = 0; centIter < nCentBins; centIter++){
      if(cent >= centBins[centIter] && cent < centBins[centIter+1]){
	centPos = centIter;
	break;
      }
    }
  }

  return centPos;
}


int jecConfigParser::GetCentBinFromHiBin(const unsigned int cent)
{
  if(!isPbPb) return 0;

  if(cent/2. < centBins.at(0) || cent/2. >= centBins.at(nCentBins)) return -1;
  
  int centPos = 0;

  for(unsigned int centIter = 0; centIter < nCentBins; centIter++){
    if(cent/2 >= centBins[centIter] && cent/2 < centBins[centIter+1]){
      centPos = centIter;
      break;
    }
  }

  return centPos;
}


void jecConfigParser::PrintCentBins()
{
  std::cout << "Centbins: ";

  for(unsigned int iter = 0; iter < centBins.size(); iter++){
    if(iter < centBins.size()-1) std::cout << centBins.at(iter) << ", ";
    else std::cout << centBins.at(iter) << ".";
  }

  std::cout << std::endl;

  return;
}


float jecConfigParser::GetMinGammaPt(){return minGammaPt;}
float jecConfigParser::GetMaxGammaPt(){return maxGammaPt;}
float jecConfigParser::GetMinGammaEta(){return minGammaEta;}
float jecConfigParser::GetMaxGammaEta(){return maxGammaEta;}
float jecConfigParser::GetGammaPtHatStagger(){return gammaPtHatStagger;}
bool jecConfigParser::GetDoGammaJtDPhiCut(){return doGammaJtDPhiCut;}
float jecConfigParser::GetGammaJtDPhiCut(){return gammaJtDPhiCut;}

float jecConfigParser::GetMinLeptonPt(){return minLeptonPt;}
float jecConfigParser::GetMinLeptonEta(){return minLeptonEta;}
float jecConfigParser::GetMaxLeptonEta(){return maxLeptonEta;}

bool jecConfigParser::GetDoLeptonFlavorCuts(){return doLeptonFlavorCuts;}

float jecConfigParser::GetMinElectronPt(){return minElectronPt;}
float jecConfigParser::GetMinElectronEta(){return minElectronEta;}
float jecConfigParser::GetMaxElectronEta(){return maxElectronEta;}
float jecConfigParser::GetMinMuonPt(){return minMuonPt;}
float jecConfigParser::GetMinMuonEta(){return minMuonEta;}
float jecConfigParser::GetMaxMuonEta(){return maxMuonEta;}

float jecConfigParser::GetMinZPt(){return minZPt;}
float jecConfigParser::GetMaxZPt(){return maxZPt;}
float jecConfigParser::GetMinZM(){return minZM;}
float jecConfigParser::GetMaxZM(){return maxZM;}
bool jecConfigParser::GetDoZJtDPhiCut(){return doZJtDPhiCut;}
float jecConfigParser::GetZJtDPhiCut(){return zJtDPhiCut;}

void jecConfigParser::WriteConfigParamsToRootFile(TFile* writeFile_p)
{
  writeFile_p->cd();
  
  TDirectory* dir_p = writeFile_p->GetDirectory("configParamsDir");
  if(dir_p){
    dir_p->cd();
  }
  else{
    dir_p = writeFile_p->mkdir("configParamsDir");
    dir_p->cd();
  }
  
  TNamed nameStr("ConfigFileName", configFileName.c_str());
  nameStr.Write("", TObject::kOverwrite);

  for(unsigned int iter = 0; iter < nValidConfigVals; iter++){
    if(validConfigVals[iter].find("PTHAT") != std::string::npos && validConfigVals[iter].size() == 5){
      std::string firstStr = "PTHAT";
      std::string secondStr = "";
      for(unsigned int iter2 = 0; iter2 < nPthats-1; iter2++){
	secondStr = secondStr + std::to_string(pthats.at(iter2)) + ", ";
      }
      if(nPthats != 0) secondStr = secondStr + std::to_string(pthats.at(nPthats-1));

      TNamed configStr(firstStr.c_str(), secondStr.c_str());
      configStr.Write("", TObject::kOverwrite);
    }
    else if(validConfigVals[iter].find("INPUT") != std::string::npos && validConfigVals[iter].size() == 5){
      for(unsigned int iter2 = 0; iter2 < nPthats; iter2++){
	std::string firstStr = "INPUT_PTHAT" + std::to_string(pthats.at(iter2));
	std::string secondStr = inputStrings.at(iter2);

	TNamed configStr(firstStr.c_str(), secondStr.c_str());
	configStr.Write("", TObject::kOverwrite);
      }
    }
    else{
      std::string configInputToWrite = configInputs[iter];
      if(configInputToWrite.size() == 0) configInputToWrite = defaultConfigInputs[iter];
      TNamed configStr(validConfigVals[iter].c_str(), configInputToWrite.c_str());
      configStr.Write("", TObject::kOverwrite);
    }
  }

  writeFile_p->cd();

  return;
}


//begin private functions
std::string jecConfigParser::returnLowerStr(std::string inStr)
{
  for(unsigned int iter = 0; iter < inStr.size(); iter++){
    for(unsigned int iter2 = 0; iter2 < alphaUpperStr.size(); iter2++){
      if(inStr.substr(iter, 1).find(alphaUpperStr.at(iter2)) != std::string::npos){
        inStr.replace(iter, 1, alphaLowerStr.substr(iter2, 1));
        break;
      }
    }
  }

  return inStr;
}

bool jecConfigParser::isTrueFalseStr(std::string trueFalseStr)
{
  trueFalseStr = returnLowerStr(trueFalseStr);

  for(unsigned int iter = 0; iter < nValidTrueFalse; iter++){
    if(validTrue[iter].size() == trueFalseStr.size() && validTrue[iter].find(trueFalseStr) != std::string::npos) return true;
    if(validFalse[iter].size() == trueFalseStr.size() && validFalse[iter].find(trueFalseStr) != std::string::npos) return true;
  }

  return false;
}

bool jecConfigParser::parseTrueFalseStr(std::string trueFalseStr)
{
  trueFalseStr = returnLowerStr(trueFalseStr);
  for(unsigned int iter = 0; iter < nValidTrueFalse; iter++){
    if(validTrue[iter].size() == trueFalseStr.size() && validTrue[iter].find(trueFalseStr) != std::string::npos) return true;
    if(validFalse[iter].size() == trueFalseStr.size() && validFalse[iter].find(trueFalseStr) != std::string::npos)return false;
  }

  std::cout << "Input \'" << trueFalseStr << "\' is invalid. Auto-return false. Please debug." << std::endl;
  return false;
}

#endif
