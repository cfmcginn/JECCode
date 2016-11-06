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
  const static unsigned int nMaxPtHat = 10;
  const std::string txtStr = ".txt";
  const std::string numStr = "0123456789";
  const std::string commaStr = ",";
  const std::string dotStr = ".";
  const std::string alphaUpperStr = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  const std::string alphaLowerStr = "abcdefghijklmnopqrstuvwxyz";
  const static unsigned int nValidTrueFalse = 2;
  const std::string validTrue[nValidTrueFalse] = {"true", "1"};
  const std::string validFalse[nValidTrueFalse] = {"false", "0"};


  const static unsigned int nValidConfigVals = 12;
  enum configIter {NPTHAT, PTHAT, INPUT, ISPBPB, JETTYPES, NJTPTBINS, JTPTLOW, JTPTHI, DOJTPTLOGBINS, DOWEIGHTS, DOPTHATSTAGGER, STAGGEROFFSET};
  const std::string validConfigVals[nValidConfigVals] = {"NPTHAT", "PTHAT", "INPUT", "ISPBPB", "JETTYPES", "NJTPTBINS", "JTPTLOW", "JTPTHI", "DOJTPTLOGBINS", "DOWEIGHTS", "DOPTHATSTAGGER", "STAGGEROFFSET"};

  const std::string defaultConfigInputs[nValidConfigVals] = {"", "", "", "FALSE", "", "", "30", "100", "FALSE", "FALSE", "TRUE", "20"};
  unsigned int nConfigInputs[nValidConfigVals] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::string configInputs[nValidConfigVals] = {"", "", "", "", "", "", "", "", "", "", "", ""};
  std::string configFileName = "";
  

  unsigned int nPthats = 0;
  std::vector<unsigned int> pthats;
  std::vector<std::string> inputStrings;
  bool isPbPb = false;

  std::string jetTypes = "";
  std::vector<std::string> jetTypesKeep;
  std::vector<std::string> jetTypesRemove;
  std::vector<std::string> jetTypesFinal;

  unsigned int nJtPtBins = 0;
  float jtPtLow = 30.;
  float jtPtHi = 100.;
  bool doJtPtLogBins = false;

  bool doWeights = false;
  bool doPthatStagger = true;
  float staggerOffset = 20.;

  std::string returnLowerStr(std::string);
  bool isTrueFalseStr(std::string);
  bool parseTrueFalseStr(std::string);

 public:
  jecConfigParser();
  jecConfigParser(const std::string);
  bool SetConfigParser(const std::string);
  void ResetConfigParser();
  std::string GetConfigFileName();
  unsigned int GetNPthats();
  unsigned int GetPthat(const unsigned int);
  void PrintPthats();
  void PrintInputs();
  std::string GetInput(const unsigned int);
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
  bool GetDoWeights();
  bool GetDoPthatStagger();
  float GetJtWeight(const unsigned int, const float);
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


    if(tempStr.substr(0, validConfigVals[NPTHAT].size()).find(validConfigVals[NPTHAT]) != std::string::npos){
      if(valStr.find("-") != std::string::npos){ // check if negative
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be non-negative. Setting to 0" << std::endl;
	nPthats = 0;
	continue;
      }
      else if(valStr.find(".") != std::string::npos){ // check if integer
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be integer (no \'.\'). Setting to 0" << std::endl;
        nPthats = 0;
	continue;
      }
      
      bool isNum = true; // check if number
      for(unsigned int iter2 = 0; iter2 < valStr.size(); iter2++){
	if(numStr.find(valStr.at(iter2)) == std::string::npos){
	  isNum = false;
	  break;
	}
      }
      if(!isNum){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be pure number from \'" << numStr << "\'. Setting to 0" << std::endl;
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
    if(tempStr.substr(0, validConfigVals[PTHAT].size()).find(validConfigVals[PTHAT]) != std::string::npos){
      while(valStr.find(",,") != std::string::npos){
	valStr.replace(valStr.find(",,"), 2, ",");
      }

      if(valStr.find("-") != std::string::npos){ // check if negative
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be non-negative. Return empty." << std::endl;
	continue;
      }
      else if(valStr.find(".") != std::string::npos){ // check if integer
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be integer (no \'.\'). Return empty" << std::endl;
	continue;
      }
      
      bool isNum = true; // check if number
      for(unsigned int iter2 = 0; iter2 < valStr.size(); iter2++){
	if(numStr.find(valStr.at(iter2)) == std::string::npos && commaStr.find(valStr.at(iter2)) == std::string::npos){
	  isNum = false;
	  break;
	}
      }
      if(!isNum){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be pure number from \'" << numStr << "\' or \',\'. Return empty" << std::endl;
	continue;
      }
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
      if(valStr.find("-") != std::string::npos){ // check if negative
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be non-negative. Setting to 0" << std::endl;
	nJtPtBins = 0;
	continue;
      }
      else if(valStr.find(".") != std::string::npos){ // check if integer
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be integer (no \'.\'). Setting to 0" << std::endl;
        nJtPtBins = 0;
	continue;
      }
      
      bool isNum = true; // check if number
      for(unsigned int iter2 = 0; iter2 < valStr.size(); iter2++){
	if(numStr.find(valStr.at(iter2)) == std::string::npos){
	  isNum = false;
	  break;
	}
      }
      if(!isNum){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be pure number from \'" << numStr << "\'. Setting to 0" << std::endl;
	nJtPtBins = 0;
	continue;
      }
      // setting
      configInputs[NJTPTBINS] = valStr;
      nJtPtBins = (unsigned int)std::stoi(valStr);
      nConfigInputs[NJTPTBINS]++;
    }


    if(tempStr.substr(0, validConfigVals[JTPTLOW].size()).find(validConfigVals[JTPTLOW]) != std::string::npos){
      if(valStr.find("-") != std::string::npos){ // check if negative
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be non-negative. Setting to default" << std::endl;
	jtPtLow = 30.;
	continue;
      }
      
      bool isNum = true; // check if number
      for(unsigned int iter2 = 0; iter2 < valStr.size(); iter2++){
	if(numStr.find(valStr.at(iter2)) == std::string::npos && dotStr.find(valStr.at(iter2)) == std::string::npos){
	  isNum = false;
	  break;
	}
      }
      if(!isNum){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be pure number from \'" << numStr << "\' or \'.\'. Setting to default" << std::endl;
	jtPtLow = 30.;
	continue;
      }
      // setting
      configInputs[JTPTLOW] = valStr;
      nConfigInputs[JTPTLOW]++;
      jtPtLow = std::stof(valStr);
    }

    if(tempStr.substr(0, validConfigVals[JTPTHI].size()).find(validConfigVals[JTPTHI]) != std::string::npos){
      if(valStr.find("-") != std::string::npos){ // check if negative
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be non-negative. Setting to default" << std::endl;
	jtPtHi = 100.;
	continue;
      }
      
      bool isNum = true; // check if number
      for(unsigned int iter2 = 0; iter2 < valStr.size(); iter2++){
	if(numStr.find(valStr.at(iter2)) == std::string::npos && dotStr.find(valStr.at(iter2)) == std::string::npos){
	  isNum = false;
	  break;
	}
      }
      if(!isNum){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be pure number from \'" << numStr << "\' or \'.\'. Setting to default" << std::endl;
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
      if(valStr.find("-") != std::string::npos){ // check if negative
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be non-negative. Setting to default" << std::endl;
	staggerOffset = 20.;
	continue;
      }
      
      bool isNum = true; // check if number
      for(unsigned int iter2 = 0; iter2 < valStr.size(); iter2++){
	if(numStr.find(valStr.at(iter2)) == std::string::npos && dotStr.find(valStr.at(iter2)) == std::string::npos){
	  isNum = false;
	  break;
	}
      }
      if(!isNum){
	std::cout << tempStr << " value \'" << valStr << "\' is invalid, must be pure number from \'" << numStr << "\' or \'.\'. Setting to default" << std::endl;
	staggerOffset = 20.;
	continue;
      }
      // setting
      configInputs[STAGGEROFFSET] = valStr;
      nConfigInputs[STAGGEROFFSET]++;
      staggerOffset = std::stof(valStr);
    }
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
    inputStrings.at(pthatIter).replace(0, inputStrings.at(pthatIter).find(",")+1, "");
    if(!checkFile(inputStrings.at(pthatIter))){
      std::cout << "INPUT \'" << inputStrings.at(pthatIter) << "\' given for PTHAT==" << pthats.at(pthatIter) << " is not a valid file. Return false" << std::endl;

      ResetConfigParser();

      return false;
    }
    else{
      //EDIT HERE
      if(pthatIter == 0){
	TFile* tempJetInFile_p = new TFile(inputStrings.at(pthatIter).c_str(), "READ");
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
	  std::cout << "Given JETTYPES, \'" << jetTypes << "\', return no collections in file \'" << inputStrings.at(pthatIter) << "\'. return false." << std::endl;

	  ResetConfigParser();

	  return false;
	}
      }
      else{
	std::vector<std::string> tempJetTypes;

        TFile* tempJetInFile_p = new TFile(inputStrings.at(pthatIter).c_str(), "READ");
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
	  std::cout << "Given JETTYPES, \'" << jetTypes << "\', return no collections in file \'" << inputStrings.at(pthatIter) << "\', after filtering on first file. return false." << std::endl;

          ResetConfigParser();

          return false;
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

  return true;
}

void jecConfigParser::ResetConfigParser()
{
  configFileName = "";
  nPthats = 0;
  pthats.clear();
  inputStrings.clear();
  isPbPb = false;
  jetTypes = "";
  jetTypesKeep.clear();
  jetTypesRemove.clear();
  jetTypesFinal.clear();
  nJtPtBins = 0;
  jtPtLow = 30.;
  jtPtHi = 100.;
  doJtPtLogBins = false;
  doWeights = false;
  doPthatStagger = true;
  staggerOffset = 20.;

  for(unsigned int iter = 0; iter < nValidConfigVals; iter++){
    nConfigInputs[iter] = 0;
    configInputs[iter] = defaultConfigInputs[iter];
  }
  
  return;
}

std::string jecConfigParser::GetConfigFileName(){return configFileName;}
unsigned int jecConfigParser::GetNPthats(){return nPthats;}

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
  if(pos >= nPthats){
    std::cout << "Requested input position, \'" << pos << "\', is >= nPthats, \'" << nPthats << "\'. Return empty string" << std::endl;
    return "";
  }
  
  return inputStrings.at(pos);
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
bool jecConfigParser::GetDoWeights(){return doWeights;}
bool jecConfigParser::GetDoPthatStagger(){return doPthatStagger;}

float jecConfigParser::GetJtWeight(const unsigned int pthatPos, const float jtPt)
{
  if(jtPt < jtPtLow || jtPt > jtPtHi) return 0.;

  float weight = 1.;

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
    else if((unsigned int)jtPtBinPthatPos[jtPtPos] != pthatPos) weight = 0.;
    else weight = 1.;    
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  return weight;
}

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
