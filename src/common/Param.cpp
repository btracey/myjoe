
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Param.h"



std::ostream& operator<< (std::ostream &out, Param &param)
{  
  for (std::vector<string>::iterator it=param.token.begin(); it<param.token.end(); it++)
    out << *it << " ";  
  return out;  
}


void ParamMap::Remove(const string &str) {

  map<string, vector<Param> >::iterator it = paraMap.find(str);
  if (it == paraMap.end()) {
    cerr << "Error: remove Param from ParamMap failed, could not find: " << str << endl;
  }
  else
    paraMap.erase (it);
}

int ParamMap::getIntParam(const string &str) {

  int iParamVec = 0;

  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end()) {
    cerr << "Error: could not find integer param: " << str << endl;
    throw(-1);
  }

  return it->second.begin()->getInt();    // return always first interger of parameter in vector
                                          // there is nothing like a multiple defined variable,
                                          // only commands can be called several times with different arguments
}

int ParamMap::getIntParam(const string &str, const int defaultValue) {

  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end())
    return defaultValue;

  return it->second.begin()->getInt();
}

int ParamMap::getIntParam(const string &str, const string &numberStr) {

  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end())
  {
    add(str+"="+numberStr);
    map<string, vector<Param> >::iterator it = paraMap.find(str);
    return it->second.begin()->getInt();
  }

  return it->second.begin()->getInt();
}

bool ParamMap::getBoolParam(const string &str) {

  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end()) {
    cerr << "Error: could not find bool param: " << str << endl;
    throw(-1);
  }

  return it->second.begin()->getBool();
}

bool ParamMap::getBoolParam(const string &str,const bool defaultValue) {

  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end())
    return defaultValue;

  return it->second.begin()->getBool();
}

double ParamMap::getDoubleParam(const string &str) {

  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end()) {
    cerr << "Error: could not find double param: " << str << endl;
    throw(-1);
  }

  return it->second.begin()->getDouble();
}

double ParamMap::getDoubleParam(const string &str,const double defaultValue) {

  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end())
    return defaultValue;

  return it->second.begin()->getDouble();
}

double ParamMap::getDoubleParam(const string &str, const string &numberStr) {

  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end())
  {
    add(str+"="+numberStr);
    map<string, vector<Param> >::iterator it = paraMap.find(str);
    return it->second.begin()->getDouble();
  }

  return it->second.begin()->getDouble();
}

string ParamMap::getStringParam(const string &str) {

  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end()) {
    cerr << "Error: could not find string param: "<< str << endl;
    throw(-1);
  }

  return it->second.begin()->getString();
}

string ParamMap::getStringParam(const string &str, const string &strValue) {

  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end())
  {
    add(str+"="+strValue);
    map<string, vector<Param> >::iterator it = paraMap.find(str);
    return it->second.begin()->getString();
  }

  return it->second.begin()->getString();
}

bool ParamMap::checkParam(const string &str) {

  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end())    return false;
  else                        return true;
}

/*
bool ParamMap::getParam(Param &p, const string &str, size_t iVecPara) {

  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end())
    return false;
  else {
    p = it->second.at(iVecPara);
    return true;
  }
}*/

bool ParamMap::getParam(vector<Param> *&v_ptr, const string &str) {

  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end())
    return false;
  else {
    v_ptr = &it->second;
    return true;
  }
}

bool ParamMap::getParam(Param *&v_ptr, const string &str, size_t iVecPara) {

  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end())
    return false;
  else {
    if (it->second.size()-1 < iVecPara)
    {
      cerr << "parameter vector access out of range in map: " << str << endl;
      throw(-1);
    }

    v_ptr = &it->second.at(iVecPara);
    return true;
  }
}


Param * ParamMap::getParam(const string &str, size_t iVecPara) {

  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end())
  {
    return NULL;
    /*
    cerr << "ERROR: parameter: " << str << " in ParamMap not found!" << endl;
    throw(-1);
    */
  }
  else
    return &it->second.at(iVecPara);
}
/*
Param * ParamMap::operator[](const string &str) {

  map<string, Param>::iterator it = paraMap.find(str);

  if (it == paraMap.end()) {
    if (mpi_rank == 0)
    cout << "parameter \""<< str << "\" not found"<< endl;
    throw(-1);
  }
  else
  return &it->second;
}*/

vector<Param> * ParamMap::getParamVector(const string &str)
{
  map<string, vector<Param> >::iterator it = paraMap.find(str);

  if (it == paraMap.end())
  {
    cerr << " could not get ParamVector from map: " << str << " not present!"<< endl;
    throw(-1);
  }
  else
    return &it->second;
}

void ParamMap::add(const string &str) {

  int hash = checkHash(str);

  if (hash != 1)
  {
    Param p(str.substr(0, hash));             // cut the string until a hash appears, otherwise use the whole string (hash = str.length())

    if (!checkParam(p.getParamName()))        // if name not present -> create list with first Param entry
    {
      vector<Param> paramVec;
      paramVec.push_back(p);
      paraMap[p.getParamName()] = paramVec;
    }
    else                                      // if name already exits push_back Param to STL vector in map
    {
      vector<Param> *paramVec = getParamVector(p.getParamName());
      paramVec->push_back(p);
    }
  }
}

void ParamMap::addParamsFromFile(const char *name,const bool verbose) {

  // note that verbose has default "false" in Param.h prototype

  ifstream infile;

  infile.open(name, ifstream::in);
  if (!infile.is_open()) {
    if (verbose)
      cout << " > Warning: param file \"" << name << "\" is not present." << endl;
    return;
  }
  
  char str[MAX_LINE_LENGTH];

  if (verbose)
    cout << "************** param file ****************" << endl;

  while (infile.good()) {
    infile.getline(str, MAX_LINE_LENGTH);
    if (verbose)
      cout << str << endl;
    add(str);
  }

  if (verbose)
    cout << "*********** end of param file ************" << endl;

  infile.close();
}
 
void ParamMap::addParamsFromArgs(int argc,char * argv[],const bool verbose) {

  // note that verbose has default "false" in Param.h prototype
    
  string str;
  assert( str.length() == 0 );

  if (verbose)
    cout << "************** param args ****************" << endl;
      
  int iargc = 1;
  while (iargc < argc) {
    
    // look for args that start with --, and skip anything before this...
    if ( (strlen(argv[iargc]) >= 3) && (argv[iargc][0] == '-') && (argv[iargc][1] == '-') ) {
      
      if (str.length() > 0) {
	if (verbose)
	  cout << str << endl;
	add(str);
      }

      // start the next...
      str = &(argv[iargc][2]);
      
    }
    else {
      
      str.append(" ");
      str.append(argv[iargc]);
      
    }
    
    iargc++;
      
  }
    
  // last one could be out of the loop...
  if (str.length() > 0) {
    if (verbose)
      cout << str << endl;
    add(str);
  }

  if (verbose)
    cout << "*********** end of param args ************" << endl;
  
}
  
void ParamMap::writeParamToFileASCIIstyleC(char fname[]) {

  FILE *fp;

  if ((fp = fopen(fname, "wt")) == NULL) {
    printf("could not open file %s to write PARAMETERS\n", fname);
    return;
  }

  for (map<string, vector<Param> >::iterator it=paraMap.begin(); it!=paraMap.end(); it++)
    for (vector<Param>::iterator itVec=it->second.begin(); itVec!=it->second.end(); itVec++)
      {
	itVec->print();
	itVec->writeToFile(fp);
      }
  
  fclose(fp);
}

void ParamMap::dump() {
  
  cout << "************** start of params ****************" << endl;
  
  for (map<string, vector<Param> >::iterator it=paraMap.begin(); it!=paraMap.end(); it++)
    for (vector<Param>::iterator itVec=it->second.begin(); itVec!=it->second.end(); itVec++)
      itVec->print();

  cout << "**************  end of params  ****************" << endl;

}




