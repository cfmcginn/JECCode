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

//include headers
#include "include/doGlobalDebug.h"
#include "include/checkMakeDir.h"
#include "include/getLogBins.h"
#include "include/getLinBins.h"
#include "include/returnRootFileContentsList.h"


class jecConfigParser{
 private:
  const static unsigned int nEventTypes = 3;
  const std::string validEventTypes[nEventTypes] = {"DIJET", "GAMMAJET", "ZJET"};
  const static unsigned int nMaxPtHat = 10;
  const std::string txtStr = ".txt";
  const std::string numStr = "0123456789";
  const std::string commaStr = ",";
  const std::string dotStr = ".";
  const std::string minusStr = "-";
  const std::string alphaUpperStr = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  const std::string alphaLowerStr = "abcdefghijklmnopqrstuvwxyz";
  const static unsigned int nValidTrueFalse = 2;
  const std::string validTrue[nValidTrueFalse] = {"true", "1"};
  const std::string validFalse[nValidTrueFalse] = {"false", "0"};


  const static unsigned int nValidConfigVals = 24;
  enum configIter {EVENTTYPE, //0
		   OUTNAME, //1
		   NPTHAT, //2
		   PTHAT, //3
		   INPUT, //4
		   ISPBPB, //5
		   JETTYPES, //6
		   NJTPTBINS, //7
		   JTPTLOW, //8
		   JTPTHI, //9
		   JTETAMAX, //10
		   JTETAPTTHRESH, //11
		   JTETABINS, //12
		   JTETAPTBINS, //13
		   DOJTPTLOGBINS, //14
		   DOWEIGHTS, //15
		   DOWEIGHTTRUNC, //16
		   PTHATWEIGHTS, //17
		   DOPTHATSTAGGER, //18
		   STAGGEROFFSET, //19
		   NCENTBINS, //20
		   CENTBINS, //21
		   MINGAMMAPT, //22
		   GAMMAPTHATSTAGGER}; //23
 
  const std::string validConfigVals[nValidConfigVals] = {"EVENTTYPE", //0
							 "OUTNAME", //1
							 "NPTHAT", //2
							 "PTHAT", //3
							 "INPUT", //4
							 "ISPBPB", //5
							 "JETTYPES", //6
							 "NJTPTBINS", //7
							 "JTPTLOW", //8
							 "JTPTHI", //9
							 "JTETAMAX", //10
							 "JTETAPTTHRESH", //11
							 "JTETABINS", //12
							 "JTETAPTBINS", //13
							 "DOJTPTLOGBINS", //14
							 "DOWEIGHTS", //15
							 "DOWEIGHTTRUNC", //16
							 "PTHATWEIGHTS", //17
							 "DOPTHATSTAGGER", //18
							 "STAGGEROFFSET", //19
							 "NCENTBINS", //20
							 "CENTBINS", //21
							 "MINGAMMAPT", //22
							 "GAMMAPTHATSTAGGER"}; //23


  const std::string configTypes[nValidConfigVals] = {"std::string", //0
						     "std::string", //1
						     "unsigned int", //2
						     "std::vector<unsigned int>", //3
						     "std::vector<std::string>", //4
						     "bool", //5
						     "std::string", //6
						     "unsigned int", //7
						     "float", //8
						     "float", //9
						     "float", //10
						     "float", //11
						     "unsigned int", //12
						     "unsigned int", //13
						     "bool", //14
						     "bool", //15
						     "bool", //16
						     "std::vector<float>", //17
						     "bool", //18
						     "float", //19
						     "unsgined int", //20
						     "std::vector<unsigned int>", //21
						     "float", //22
						     "float"}; //23
  

  const std::string defaultConfigInputs[nValidConfigVals] = {"", //0
							     "", //1
							     "", //2
							     "", //3
							     "", //4
							     "FALSE", //5
							     "", //6
							     "", //7
							     "30", //8
							     "100", //9
							     "1.6", //10
							     "35", //11
							     "16", //12
							     "3", //13
							     "FALSE", //14
							     "FALSE", //15
							     "FALSE", //16
							     "", //17
							     "TRUE", //18
							     "20", //19
							     "2", //20
							     "100,30,0", //21
							     "40.", //22
							     "0."}; //23

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
						  0}; //23
  

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
						""}; //23

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
  float jtPtLow = 30.;
  float jtPtHi = 100.;
  bool doJtPtLogBins = false;

  float jtEtaMax = 1.6;
  float jtEtaPtThresh = 35;
  unsigned int jtEtaBins = 16;
  unsigned int jtPtEtaBins = 3;

  bool doWeights = false;
  bool doWeightTrunc = false;
  std::vector<double> pthatWeights;
  bool doPthatStagger = true;
  float staggerOffset = 20.;

  int nCentBins = 2;
  std::vector<unsigned int> centBins = {100, 30, 0};

  float minGammaPt = 40.;
  float gammaPtHatStagger = 0.;

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
  bool SetConfigParser(const std::string);
  void ResetConfigParser();
  std::string GetConfigFileName();
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
  float GetJtPtLow();
  float GetJtPtHi();
  bool GetDoJtPtLogBins();
  float GetJtEtaMax();
  float GetJtEtaPtThresh();
  unsigned int GetJtEtaBins();
  unsigned int GetJtPtEtaBins();
  bool GetDoWeights();
  bool GetDoWeightTrunc();
  bool GetDoPthatStagger();
  float GetJtWeight(const unsigned int, const float, const float, const float);
  double GetPtHatWeight(const float);
  double GetTruncPtHatWeight(const float, const float);
  unsigned int GetNCentBins();
  std::vector<unsigned int> GetCentBins();
  unsigned int GetCentBinFromPos(const unsigned int);
  unsigned int GetCentBinFromCent(const unsigned int);
  int GetCentBinFromHiBin(const unsigned int);
  void PrintCentBins();
  float GetMinGammaPt();
  float GetGammaPtHatStagger();

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
    std::cout << floatString << "\' is invalid, must be number from \'" << numStr << "\' or \'" << dotStr << "\' or \'" << minusStr << "\'. If \'" << minusStr << "\', must be first char. return false" << std::endl;
  }

  return isNum;
}

bool jecConfigParser::StringIsGoodUFloat(const std::string floatString)
{
  bool isUFloat = StringIsGoodFloat(floatString);
  if(!isUFloat) return isUFloat;

  if(floatString.find(minusStr) != std::string::npos){
    isUFloat = false;
    std::cout << floatString << "\' is invalid, must be number from \'" << numStr << "\' or \'" << dotStr << "\', non-negative. return false" << std::endl;
  }

  return isUFloat;
}

bool jecConfigParser::StringIsGoodInt(const std::string intString)
{
  bool isInt = StringIsGoodFloat(intString);
  if(!isInt) return isInt;

  if(intString.find(dotStr) != std::string::npos){
    std::cout << intString << "\' is invalid, must be number from \'" << numStr << "\' or \'" << minusStr << "\'. If \'" << minusStr << "\', must be first char. return false" << std::endl;
    isInt = false;
  }
  return isInt;
}

bool jecConfigParser::StringIsGoodUInt(std::string intString)
{
  bool isUInt = StringIsGoodUFloat(intString);
  if(!isUInt) return isUInt;

  if(intString.find(dotStr) != std::string::npos && intString.find(dotStr) != intString.size()-1){
    std::cout << intString << "\' is invalid, must be number from \'" << numStr << "\', non-negative and non-decimal. return false" << std::endl;
    isUInt = false;
  }
  return isUInt;
}

bool jecConfigParser::SetConfigParser(const std::string inConfigFile)
{
  std::vector<std::string> lines;
  if(!checkFile(inConfigFile)){
    std::cout << "Input jecConfig txt file, \'" << inConfigFile << "\', is not valid file. return false" << std::endl;
    return false;
  }  
  else if(inConfigFile.size() < txtStr.size()){
    std::cout << "Input jecConfig txt file, \'" << inConfigFile << "\', doesn't end in \'.txt\'. return false" << std::endl;
    return false;
  }
  else if(inConfigFile.substr(inConfigFile.size()-txtStr.size(), txtStr.size()).find(txtStr) == std::string::npos){
    std::cout << "Input jecConfig txt file, \'" << inConfigFile << "\', doesn't end in \'.txt\'. return false" << std::endl;
    return false;
  }

  configFileName = inConfigFile;

  std::ifstream inFile(inConfigFile.c_str());
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

  //begin setting class params
  const unsigned int nLines = lines.size();  
  for(unsigned int iter = 0; iter < nLines; iter++){
    tempStr = lines.at(iter);
    std::string valStr = tempStr.substr(tempStr.find("=")+1, tempStr.size() - tempStr.find("=")+1);
    tempStr.replace(tempStr.find("="), tempStr.size() - tempStr.find("="), "");

    for(unsigned int iter = 0; iter < tempStr.size(); iter++){
      if(alphaLowerStr.find(tempStr.at(iter)) != std::string::npos) tempStr.replace(iter, 1, alphaUpperStr.substr(alphaLowerStr.find(tempStr.at(iter)), 1));
    }

    if(tempStr.substr(0, validConfigVals[EVENTTYPE].size()).find(validConfigVals[EVENTTYPE]) != std::string::npos){
      eventTypeStr = valStr;
    }

    if(tempStr.substr(0, validConfigVals[OUTNAME].size()).find(validConfigVals[OUTNAME]) != std::string::npos){
      outNameStr = valStr;
    }

    if(tempStr.substr(0, validConfigVals[NPTHAT].size()).find(validConfigVals[NPTHAT]) != std::string::npos){
      if(!StringIsGoodUInt(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Setting to 0." << std::endl;
        nPthats = 0;
        continue;
      }

      // setting
      nPthats = (unsigned int)std::stoi(valStr);

      // finally check its below class cap
      if(nPthats > nMaxPtHat){
	std::cout << tempStr << " value \'" << valStr << "\' is greater than class cap \'" << nMaxPtHat << "\'. Consider raising cap or using fewer pthats. Setting to 0." << std::endl;
	nPthats = 0;
	valStr = "";
      }
      else{
	configInputs[NPTHAT] = valStr;
	nConfigInputs[NPTHAT]++;
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

    if(tempStr.substr(0, validConfigVals[INPUT].size()).find(validConfigVals[INPUT]) != std::string::npos){
      inputStrings.push_back(valStr);
      nConfigInputs[INPUT]++;      
    }

    if(tempStr.substr(0, validConfigVals[ISPBPB].size()).find(validConfigVals[ISPBPB]) != std::string::npos){
      if(!isTrueFalseStr(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Please give \'TRUE\' or \'FALSE\'. Return false" << std::endl;
	
	ResetConfigParser();
	
	return false;
      }
      
      configInputs[ISPBPB] = valStr;
      nConfigInputs[ISPBPB]++;

      if(parseTrueFalseStr(valStr)) isPbPb = true;
      else isPbPb = false;
    }

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

    if(tempStr.substr(0, validConfigVals[NJTPTBINS].size()).find(validConfigVals[NJTPTBINS]) != std::string::npos){
      if(!StringIsGoodUInt(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Setting to 0." << std::endl;
        nJtPtBins = 0;
        continue;
      }
      // setting
      configInputs[NJTPTBINS] = valStr;
      nJtPtBins = (unsigned int)std::stoi(valStr);
      nConfigInputs[NJTPTBINS]++;
    }

    if(tempStr.substr(0, validConfigVals[JTPTLOW].size()).find(validConfigVals[JTPTLOW]) != std::string::npos){
      if(!StringIsGoodUFloat(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Setting to 30." << std::endl;
        jtPtLow = 30.;
        continue;
      }

      // setting
      configInputs[JTPTLOW] = valStr;
      nConfigInputs[JTPTLOW]++;
      jtPtLow = std::stof(valStr);
    }

    if(tempStr.substr(0, validConfigVals[JTPTHI].size()).find(validConfigVals[JTPTHI]) != std::string::npos){
      if(!StringIsGoodUFloat(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Setting to 100." << std::endl;
        jtPtHi = 100.;
        continue;
      }
      // setting
      configInputs[JTPTHI] = valStr;
      nConfigInputs[JTPTHI]++;
      jtPtHi = std::stof(valStr);
    }

    if(tempStr.substr(0, validConfigVals[DOJTPTLOGBINS].size()).find(validConfigVals[DOJTPTLOGBINS]) != std::string::npos){
      if(!isTrueFalseStr(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Defaulting to false. Return false" << std::endl;
	doJtPtLogBins = false;
	continue;
      }
      
      configInputs[DOJTPTLOGBINS] = valStr;
      nConfigInputs[DOJTPTLOGBINS]++;
      
      if(parseTrueFalseStr(valStr)) doJtPtLogBins = true;
      else doJtPtLogBins = false;
    }  

    if(tempStr.substr(0, validConfigVals[JTETAMAX].size()).find(validConfigVals[JTETAMAX]) != std::string::npos){
      if(!StringIsGoodUFloat(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Setting to 1.6." << std::endl;
        jtEtaMax = 1.6;
        continue;
      }
      // setting
      configInputs[JTETAMAX] = valStr;
      nConfigInputs[JTETAMAX]++;
      jtEtaMax = std::stof(valStr);
    }



    if(tempStr.substr(0, validConfigVals[JTETAPTTHRESH].size()).find(validConfigVals[JTETAPTTHRESH]) != std::string::npos){
      if(!StringIsGoodUFloat(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Setting to 35." << std::endl;
        jtEtaPtThresh = 35.;
        continue;
      }
      // setting
      configInputs[JTETAPTTHRESH] = valStr;
      nConfigInputs[JTETAPTTHRESH]++;
      jtEtaPtThresh = std::stof(valStr);
    }


    if(tempStr.substr(0, validConfigVals[JTETABINS].size()).find(validConfigVals[JTETABINS]) != std::string::npos){
      if(!StringIsGoodUInt(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Setting to 16." << std::endl;
        jtEtaBins = 16;
        continue;
      }
      // setting
      configInputs[JTETABINS] = valStr;
      nConfigInputs[JTETABINS]++;
      jtEtaBins = std::stof(valStr);
    }


    if(tempStr.substr(0, validConfigVals[JTETAPTBINS].size()).find(validConfigVals[JTETAPTBINS]) != std::string::npos){
      if(!StringIsGoodUInt(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Setting to 3." << std::endl;
        jtPtEtaBins = 3;
        continue;
      }
      // setting
      configInputs[JTETAPTBINS] = valStr;
      nConfigInputs[JTETAPTBINS]++;
      jtPtEtaBins = std::stof(valStr);
    }

    if(tempStr.substr(0, validConfigVals[DOWEIGHTS].size()).find(validConfigVals[DOWEIGHTS]) != std::string::npos){
      if(!isTrueFalseStr(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Defaulting to false. Return false" << std::endl;
	doWeights = false;
	continue;
      }

      configInputs[DOWEIGHTS] = valStr;
      nConfigInputs[DOWEIGHTS]++;

      if(parseTrueFalseStr(valStr)) doWeights = true;
      else doWeights = false;
    }

    if(tempStr.substr(0, validConfigVals[DOWEIGHTTRUNC].size()).find(validConfigVals[DOWEIGHTTRUNC]) != std::string::npos){
      if(!isTrueFalseStr(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Defaulting to false. Return false" << std::endl;	
	doWeightTrunc = false;
	continue;
      }

      configInputs[DOWEIGHTTRUNC] = valStr;
      nConfigInputs[DOWEIGHTTRUNC]++;

      if(parseTrueFalseStr(valStr)) doWeightTrunc = true;
      else doWeightTrunc = false;
    }

    if(tempStr.substr(0, validConfigVals[PTHATWEIGHTS].size()).find(validConfigVals[PTHATWEIGHTS]) != std::string::npos){
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
      configInputs[PTHATWEIGHTS] = valStr;
      nConfigInputs[PTHATWEIGHTS]++;

      while(valStr.find(",") != std::string::npos){
	pthatWeights.push_back(std::stof(valStr.substr(0, valStr.find(","))));
	valStr.replace(0, valStr.find(",")+1, "");
      }
      if(valStr.size() != 0) pthatWeights.push_back(std::stof(valStr));
    }

    if(tempStr.substr(0, validConfigVals[DOPTHATSTAGGER].size()).find(validConfigVals[DOPTHATSTAGGER]) != std::string::npos){
      if(!isTrueFalseStr(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Defaulting to true. Return true" << std::endl;
	doPthatStagger = true;
	continue;
      }

      configInputs[DOPTHATSTAGGER] = valStr;	
      nConfigInputs[DOPTHATSTAGGER]++;

      if(parseTrueFalseStr(valStr)) doPthatStagger = true;
      else doPthatStagger = false;
    }

    if(tempStr.substr(0, validConfigVals[STAGGEROFFSET].size()).find(validConfigVals[STAGGEROFFSET]) != std::string::npos){
      if(!StringIsGoodUInt(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Setting to 20." << std::endl;
	staggerOffset = 20.;
        continue;
      }
      // setting
      configInputs[STAGGEROFFSET] = valStr;
      nConfigInputs[STAGGEROFFSET]++;
      staggerOffset = std::stof(valStr);
    }

    if(tempStr.substr(0, validConfigVals[NCENTBINS].size()).find(validConfigVals[NCENTBINS]) != std::string::npos){
      if(!StringIsGoodUInt(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Setting to 0." << std::endl;
        nCentBins = 0;
        continue;
      }
      // setting                                                                                                   
      configInputs[NCENTBINS] = valStr;
      nCentBins = (unsigned int)std::stoi(valStr);
      nConfigInputs[NCENTBINS]++;
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
    
    if(tempStr.substr(0, validConfigVals[GAMMAPTHATSTAGGER].size()).find(validConfigVals[GAMMAPTHATSTAGGER]) != std::string::npos){
      if(!StringIsGoodUInt(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Setting to 0." << std::endl;
        gammaPtHatStagger = 0.;
        continue;
      }

      // setting                                                                                                   
      configInputs[GAMMAPTHATSTAGGER] = valStr;
      gammaPtHatStagger = std::stof(valStr);
      nConfigInputs[GAMMAPTHATSTAGGER]++;
    }

    if(tempStr.substr(0, validConfigVals[MINGAMMAPT].size()).find(validConfigVals[MINGAMMAPT]) != std::string::npos){
      if(!StringIsGoodUFloat(valStr)){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid. Setting to 40." << std::endl;
	minGammaPt = 40.;
        continue;
      }
      // setting                                                                                                   
      configInputs[MINGAMMAPT] = valStr;
      minGammaPt = std::stof(valStr);
      nConfigInputs[MINGAMMAPT]++;
    }
  }

  bool isGoodEventType = false;
  for(unsigned int iter = 0; iter < nEventTypes; iter++){
    if(eventTypeStr.size() == validEventTypes[iter].size() && eventTypeStr.find(validEventTypes[iter]) != std::string::npos){
      isGoodEventType = true;
      break;
    }
  }

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

  if(0 == nJtPtBins){
    std::cout << "NJTPTBINS, \'" << nJtPtBins << "\', is not a valid value. Return false." << std::endl;

    ResetConfigParser();

    return false;
  }

  if(jtPtLow >= jtPtHi){
    std::cout << "JTPTLOW, \'" << jtPtLow << "\', is greater than JTPTHI, \'" << jtPtHi << "\'. Return false" << std::endl;

    ResetConfigParser();
    return false;
  }

  if((doWeights && doPthatStagger) || (!doWeights && !doPthatStagger)){
    std::cout << "DOWEIGHTS, \'" << doWeights << "\', and DOPTHATSTAGGER, \'" << doPthatStagger << "\', both have same value. Please choose one. Return false." << std::endl;

    ResetConfigParser();

    return false;
  }


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

  if(0 == nCentBins && isPbPb){
    std::cout << "NCENTBINS, \'" << nCentBins << "\', is not a valid value. Return false." << std::endl;

    ResetConfigParser();

    return false;
  }

  if((int)centBins.size() != nCentBins+1 && isPbPb){
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

  return true;
}

void jecConfigParser::ResetConfigParser()
{
  configFileName = "";
  eventTypeStr = "";
  outNameStr = "";
  isDijet = false;
  isGammaJet = false;
  isZJet = false;
  nPthats = 0;
  pthats.clear();
  inputStrings.clear();
  inputFileList.clear();
  inputFilePtHats.clear();
  isPbPb = false;
  jetTypes = "";
  jetTypesKeep.clear();
  jetTypesRemove.clear();
  jetTypesFinal.clear();
  nJtPtBins = 0;
  jtPtLow = 30.;
  jtPtHi = 100.;
  jtEtaMax = 1.6;
  jtEtaPtThresh = 35;
  jtEtaBins = 16;
  jtPtEtaBins = 3;
  doJtPtLogBins = false;
  doWeights = false;
  doWeightTrunc = false;
  pthatWeights.clear();
  doPthatStagger = true;
  staggerOffset = 20.;
  nCentBins = 2;
  centBins.clear();
  centBins.push_back(100);
  centBins.push_back(30);
  centBins.push_back(0);
  minGammaPt = 40.;
  gammaPtHatStagger = 0.;

  for(unsigned int iter = 0; iter < nValidConfigVals; iter++){
    nConfigInputs[iter] = 0;
    configInputs[iter] = defaultConfigInputs[iter];
  }
  
  return;
}

std::string jecConfigParser::GetConfigFileName(){return configFileName;}
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
float jecConfigParser::GetJtPtLow(){return jtPtLow;}
float jecConfigParser::GetJtPtHi(){return jtPtHi;}
bool jecConfigParser::GetDoJtPtLogBins(){return doJtPtLogBins;}
float jecConfigParser::GetJtEtaMax(){return jtEtaMax;}
float jecConfigParser::GetJtEtaPtThresh(){return jtEtaPtThresh;}
unsigned int jecConfigParser::GetJtEtaBins(){return jtEtaBins;}
unsigned int jecConfigParser::GetJtPtEtaBins(){return jtPtEtaBins;}
bool jecConfigParser::GetDoWeights(){return doWeights;}
bool jecConfigParser::GetDoWeightTrunc(){return doWeightTrunc;}
bool jecConfigParser::GetDoPthatStagger(){return doPthatStagger;}

float jecConfigParser::GetJtWeight(const unsigned int pthatPos, const float jtPt, const float jtEta, const float bosonPt = 0.)
{
  if(jtPt < jtPtLow || jtPt > jtPtHi) return 0.;

  if(TMath::Abs(jtEta) > jtEtaMax) return 0.;

  if(isGammaJet && bosonPt < minGammaPt) return 0.;

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
  if(doJtPtLogBins) getLogBins(jtPtLow, jtPtHi, nJtPtBins, jtPtBins);
  else getLinBins(jtPtLow, jtPtHi, nJtPtBins, jtPtBins);

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
    for(int centIter = 0; centIter < nCentBins; centIter++){
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

  for(int centIter = 0; centIter < nCentBins; centIter++){
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
float jecConfigParser::GetGammaPtHatStagger(){return gammaPtHatStagger;}


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
