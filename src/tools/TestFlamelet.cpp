/**
 * \brief Tool that check if the flamelets do intersect (to be then included in the createchemtable tool)
 *
 * \author Vincent Terrapon 
 * \date August 2009
 * \version 1.0
 */

#include "Flamelet.h"

/* Global variables */
int             nFiles;
Flamelet        *myFlamelets;
vector <string> VarProg;
vector <string> FlameletList;
vector <double> WeightProg;
string          VarComp;
string          CombustionModel;
int             CombustionRegime;


void TestFlameletInitialize(ParamMap *myInputPtr)
/* Read input file and flamelets */
{
  cout << endl << "***** Flamelet initialization *****" << endl << endl;
  
  Param    *p = NULL;  
  string   FlameletListFilename, file, buffer_str;
  ifstream fin;
  
  CombustionModel = myInputPtr->getStringParam("COMBUSTION_MODEL", "FPVA");
  
  if (CombustionModel == "FPVA")
  {
    CombustionRegime = 2;
    p = myInputPtr->getParam("VAR_PROG");
    if (p == NULL)
    {
      cerr << "### Missing input VAR_PROG in input file! ###" << endl;
      throw(-1);
    }
    int nVarProg = p->getSize() - 1;
    if (nVarProg <= 0)
    {
      cerr << "### Missing variables for progress variable in input file! ###" << endl;
      throw(-1);
    }
    for (int i = 0; i < nVarProg; i++)
      VarProg.push_back(p->getString(i + 1));
    
    p = myInputPtr->getParam("WEIGHT_PROG");
    if (p == NULL)
    {
      cout << "No weights found in input file for progress variable, set to default value 1.0" << endl;
      for (int i = 0; i < nVarProg; i++)
        WeightProg.push_back(1.0);
    }
    else
    {
      if (nVarProg != p->getSize()-1)
      {
        cerr << "### Number of variables for progress variable does not match the number of weights in input file! ###" << endl;
        throw(-1);
      }
      for (int i = 0; i < nVarProg; i++)
        WeightProg.push_back(p->getDouble(i + 1));
    }
  }
  else if (CombustionModel == "STEADY")
  {
    CombustionRegime = 1;
    VarProg.push_back("chi");
    WeightProg.push_back(1.0);
  }
  else
  {
    cerr << "### Combustion model not recognized! ###" << endl;
    throw(-1);
  }

  // Read variable to use for comparison
  VarComp = myInputPtr->getStringParam("VAR_COMPARISON", "PROG");

  // Read list of flamelet files
  FlameletListFilename = myInputPtr->getStringParam("FLAMELET_LIST_FILENAME", "FlameletList.txt");
  char filename[kStrLong];
  size_t dum = FlameletListFilename.copy(&filename[0], kStrLong);
  dum = FlameletListFilename.size();
  if (dum < kStrLong)
    strcpy(&filename[dum], "\0");
  fin.open(filename);
  if (fin.fail())
  {
    cerr << "### Cannot open flamelet list input file " << filename << " ! ###" << endl;
    throw(-1);
  }
  while(getline(fin, buffer_str))
  {
    istringstream buf(buffer_str);
    buf >> file;
    if (buffer_str.empty())
      break;
    FlameletList.push_back(file);
  }
  fin.close();

  // Load the flamelets
  nFiles = FlameletList.size();
  myFlamelets = new Flamelet[nFiles];
  for (int fl=0; fl<nFiles; fl++)
    myFlamelets[fl].LoadFlamelet(CombustionRegime, myInputPtr, FlameletList[fl], VarProg, WeightProg);
  
  cout << endl << "***** Flamelet initialized *****" << endl << endl;
  return;
}

/**********************************************************************************************/

void TestFlameletFinalize()
/* Clean memory */
{
  if (myFlamelets != NULL)  delete [] myFlamelets;
  
  cout << endl << "***** Flamelet finalized *****" << endl << endl;
  
  return;
}

/**********************************************************************************************/

void CheckFlameletIntersection()
/* Loop over all flamelets and check if they intersect (order n^2/2 curve interpolations!) 
   Algorithm is not optimized and could be done better, but easier that way...
   Assumes that flamelets variable to compare starts at 0 and ends at 0. i.e., Z=0 and Z=1 are the
   only intersection points.
   The algorithm uses the fact that the linear interpolation between two points is a convex function
 */
{
  int    nPtot1, nPtot2, np1, np2;
  int    isLarger;
  double dVarComp, Y1, Y2, Xintersec;
  int    J1, J2;
    
  // Take one flamelet ...
  for (int nfl1=0; nfl1<nFiles-1; nfl1++)
  {
    cout << "Flamelet chi_st=" << myFlamelets[nfl1].GetFlameletChiSt() << ": " << endl;
    nPtot1 = myFlamelets[nfl1].GetFlameletSize();
    J1 = myFlamelets[nfl1].GetFlameletVariableIndex(VarComp);
    // ... and compare it with all other flamelets
    for (int nfl2=nfl1+1; nfl2<nFiles; nfl2++)
    {
      nPtot2 = myFlamelets[nfl2].GetFlameletSize();
      J2 = myFlamelets[nfl2].GetFlameletVariableIndex(VarComp);
      
      // First find out which curve is larger by comparing maximum value
      double max1 = myFlamelets[nfl1].GetFlameletVariableMax(J1);
      double max2 = myFlamelets[nfl2].GetFlameletVariableMax(J2);
      if (max1 > max2)
        isLarger = 1;
      else if (max2 > max1)
        isLarger = 2;
      else
        isLarger = 0;
      
//      cout << "Flamelet 1 max: " << max1 << "    Flamelet 2 max: " << max2 << endl;
      
      // Go over all interior points of the flamelet 1 and ...
      np2 = 0;
      for (np1=1; np1<nPtot1-1; np1++)
      {
        // ... find the 2 points of flamelet 2 around given point of flamelet 1 and check if isLarger changes
        // np2 is the largest point still smaller that np1 so that np2+1 is greater or equal than np1 (in Z coordinates)
        while ((myFlamelets[nfl2].GetFlameletZ(np2+1) < myFlamelets[nfl1].GetFlameletZ(np1)) && (np2+1 < nPtot2-1))
        {
          // Increment np2 and check each of them where they lie relative to flamelet 1 
          // Compare point value of flamelet 2 to linear interpolation of flamelet 1 between np1-1 and np1
          np2++;
          Y1 = myFlamelets[nfl1].GetFlameletVariable(np1-1, J1) + 
               (myFlamelets[nfl1].GetFlameletVariable(np1, J1) - myFlamelets[nfl1].GetFlameletVariable(np1-1, J1)) / 
               (myFlamelets[nfl1].GetFlameletZ(np1) - myFlamelets[nfl1].GetFlameletZ(np1-1)) * 
               (myFlamelets[nfl2].GetFlameletZ(np2) - myFlamelets[nfl1].GetFlameletZ(np1-1));
          Y2 = myFlamelets[nfl2].GetFlameletVariable(np2, J2);
          dVarComp = Y1 - Y2;
          
          if (((isLarger == 1) || (isLarger == -1)) && (dVarComp < 0.0))  // intersection
          {
            isLarger = 0;
            Xintersec = myFlamelets[nfl2].GetFlameletZ(np2);
//            cout << "A " << myFlamelets[nfl1].GetFlameletZ(np1-1) << " " << myFlamelets[nfl2].GetFlameletZ(np2) << 
//                     " " << myFlamelets[nfl1].GetFlameletZ(np1) << endl;
//            cout << "A " << myFlamelets[nfl1].GetFlameletVariable(np1-1, J1) << " " << myFlamelets[nfl2].GetFlameletVariable(np2, J2) << 
//                     " " << myFlamelets[nfl1].GetFlameletVariable(np1, J1) << endl;
//            cout << "A " << Y1 << " " << Y2 << endl;
//            cout << "A isLarger " << isLarger << endl;
          }
          else if (((isLarger == 2) || (isLarger == -2)) && (dVarComp > 0.0))  // intersection
          {
            isLarger = 0;
            Xintersec = myFlamelets[nfl2].GetFlameletZ(np2);
//            cout << "A " << myFlamelets[nfl1].GetFlameletZ(np1-1) << " " << myFlamelets[nfl2].GetFlameletZ(np2) << 
//                     " " << myFlamelets[nfl1].GetFlameletZ(np1) << endl;
//            cout << "A " << myFlamelets[nfl1].GetFlameletVariable(np1-1, J1) << " " << myFlamelets[nfl2].GetFlameletVariable(np2, J2) << 
//                    " " << myFlamelets[nfl1].GetFlameletVariable(np1, J1) << endl;
//            cout << "A " << Y1 << " " << Y2 << endl;
//            cout << "A isLarger " << isLarger << endl;
          }
          else if ((isLarger == 2) && (dVarComp == 0.0))  // overlap
          {
            isLarger = -2;
            Xintersec = myFlamelets[nfl2].GetFlameletZ(np2);
          }
          else if ((isLarger == 1) && (dVarComp == 0.0))  // overlap
          {
            isLarger = -1;
            Xintersec = myFlamelets[nfl2].GetFlameletZ(np2);
          }
        }
        // now np2 is just before np1, lets compute the value of flamelet 2 variable at point np1
        Y1 = myFlamelets[nfl1].GetFlameletVariable(np1, J1);
        Y2 = myFlamelets[nfl2].GetFlameletVariable(np2, J2) + 
             (myFlamelets[nfl2].GetFlameletVariable(np2+1, J2) - myFlamelets[nfl2].GetFlameletVariable(np2, J2)) / 
             (myFlamelets[nfl2].GetFlameletZ(np2+1) - myFlamelets[nfl2].GetFlameletZ(np2)) * 
             (myFlamelets[nfl1].GetFlameletZ(np1) - myFlamelets[nfl2].GetFlameletZ(np2));
        dVarComp = Y1 - Y2;
        if (((isLarger == 1) || (isLarger == -1)) && (dVarComp < 0.0))  // intersection
        {
          isLarger = 0;
          Xintersec = myFlamelets[nfl1].GetFlameletZ(np1);
//          cout << myFlamelets[nfl2].GetFlameletZ(np2) << " " << myFlamelets[nfl1].GetFlameletZ(np1) << 
//           " " << myFlamelets[nfl2].GetFlameletZ(np2+1) << endl;
//          cout << myFlamelets[nfl2].GetFlameletVariable(np2, J2) << " " << myFlamelets[nfl1].GetFlameletVariable(np1, J1) << 
//           " " << myFlamelets[nfl2].GetFlameletVariable(np2+1, J2) << endl;
//          cout << Y1 << " " << Y2 << endl;
//          cout <<"isLarger " << isLarger << endl;
        }
        else if (((isLarger == 2) || (isLarger == -2)) && (dVarComp > 0.0))  // intersection
        {
          isLarger = 0;
          Xintersec = myFlamelets[nfl1].GetFlameletZ(np1);
//          cout << myFlamelets[nfl2].GetFlameletZ(np2) << " " << myFlamelets[nfl1].GetFlameletZ(np1) << 
//           " " << myFlamelets[nfl2].GetFlameletZ(np2+1) << endl;
//          cout << myFlamelets[nfl2].GetFlameletVariable(np2, J2) << " " << myFlamelets[nfl1].GetFlameletVariable(np1, J1) << 
//           " " << myFlamelets[nfl2].GetFlameletVariable(np2+1, J2) << endl;
//          cout << Y1 << " " << Y2 << endl;
//          cout <<"isLarger " << isLarger << endl;
        }
        else if ((isLarger == 2) && (dVarComp == 0.0))  // overlap
        {
          isLarger = -2;
          Xintersec = myFlamelets[nfl1].GetFlameletZ(np1);
        }
        else if ((isLarger == 1) && (dVarComp == 0.0))  // overlap
        {
          isLarger = -1;
          Xintersec = myFlamelets[nfl1].GetFlameletZ(np1);
        }
      }

      // Flamelet 2 checked agains flamelet 1
      if (isLarger == 0)
        cout << "Error: flamelet chi_st=" << myFlamelets[nfl1].GetFlameletChiSt() << " intersects with flamelet chi_st=" << myFlamelets[nfl2].GetFlameletChiSt() 
             << "   at Z=" << Xintersec << endl;
      if (isLarger < 0)
        cout << "Error: flamelet chi_st=" << myFlamelets[nfl1].GetFlameletChiSt() << " overlaps with flamelet chi_st=" << myFlamelets[nfl2].GetFlameletChiSt() 
             << "   at Z=" << Xintersec << endl;
    }
    cout << endl;
  }
  return;
}


/**********************************************************************************************/

int main()
{
  char InputFileName[kStrLong] = "./TestFlamelet.in";
  ParamMap myInput(InputFileName);

  TestFlameletInitialize(&myInput);
  
  CheckFlameletIntersection();
  
  TestFlameletFinalize();
  
  return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////
