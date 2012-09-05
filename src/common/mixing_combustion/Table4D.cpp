#include "Table4D.h"

Table4D::Table4D()
{
  Init();
  interp = new InterpolationIndex(4,2);
  Allocated = false;
}

Table4D::Table4D(int i1, int i2, int i3, int i4, int i5)
{
  Init();
  Allocate(i1, i2, i3, i4, i5);
  interp = new InterpolationIndex(4,2);
}
//------------------------------//
//         Destructor           //
//------------------------------// 

Table4D::~Table4D() { DestroyTable(); }
  
//------------------------------//
//    De-/Allocation Methods    //
//------------------------------//
void Table4D::Allocate(int i1, int i2, int i3, int i4, int i5)
{
  SetDimensions(i1, i2, i3, i4, i5);
  Allocate();
}

void Table4D::Allocate()
{
  getMem1D(&x1, 0, n1-1, "Chemtable::Load x1", true);
  getMem1D(&x2, 0, n2-1, "Chemtable::Load x2", true);
  getMem1D(&x3, 0, n3-1, "Chemtable::Load x3", true);
  getMem1D(&x4, 0, n4-1, "Chemtable::Load x4", true);
  getMem5D(&Data, 0, n1-1, 0, n2-1, 0, n3-1, 0, n4-1, 0, nvar-1, "Chemtable::Load Data", true);
  SetAllocated();
}

void Table4D::DestroyTable()
{
  if (x1 != NULL)       {freeMem1D(x1, 0, n1-1); x1 = NULL;}
  if (x2 != NULL)       {freeMem1D(x2, 0, n2-1); x2 = NULL;}
  if (x3 != NULL)       {freeMem1D(x3, 0, n3-1); x3 = NULL;}
  if (x4 != NULL)       {freeMem1D(x4, 0, n3-1); x4 = NULL;}
  if (Data != NULL)     {freeMem5D(Data, 0, n1-1, 0, n2-1, 0, n3-1, 0, n4-1, 0, nvar-1);   Data = NULL;}
}
//------------------------------//
//     Check Bound Methods      //
//------------------------------//
void Table4D::SetAllocated() { if ( (x1 != NULL)&&
                                    (x2 != NULL)&&
                                    (x3 != NULL)&&
                                    (x4 != NULL)&&
                                    (Data != NULL) ) Allocated = true;}

bool Table4D::IsAllocated()
{
  SetAllocated();
  return Allocated;
}

void Table4D::Init(void) {
  n1 = 0;
  n2 = 0;
  n3 = 0;
  n4 = 0;
  nvar =0;

  x1   = NULL;
  x2   = NULL;
  x3   = NULL;
  x4   = NULL;
  Data = NULL;

  point[0] = -1.0;
  point[1] = -1.0;
  point[2] = -1.0;
  point[3] = -1.0;

  interp_mode[0] = 0; // Binary Serach
  interp_mode[1] = 0; // Linear mesh followed by linear stretching for Zm
  interp_mode[2] = 0; // Power law for Zvar which is Sz in fact
  interp_mode[3] = 0; // linear law for C
}


void Table4D::CopyCoordinate(int n, double* x, double* copy) 
{
  for (int i=0; i<n; ++i) { copy[i] = x[i]; }
}

void Table4D::CopyCoordinate1(double* x)  { CopyCoordinate(n1,x1,x); }
void Table4D::CopyCoordinate2(double* x)  { CopyCoordinate(n2,x2,x); }
void Table4D::CopyCoordinate3(double* x)  { CopyCoordinate(n3,x3,x); }
void Table4D::CopyCoordinate4(double* x)  { CopyCoordinate(n4,x4,x); }


// GetArrayInfo
void Table4D::GetArrayInfo(int ivar, float &minval, float &maxval,
                           int &mi1, int &mi2, int &mi3, int &mi4,
                           int &ma1, int &ma2, int &ma3, int &ma4)
{
  GetArrayInfo4D(Data, ivar, n1, n2, n3, n4,
                 minval, maxval,
                 mi1, mi2, mi3, mi4,
                 ma1, ma2, ma3, ma4);

}

void Table4D::InterpolatePoint_old(void) {
  // Dimension 1
  BinarySearch(interp->index[0], interp->weight[0][0], x1, 0, n1-1, point[0]);
  interp->weight[0][1] = 1.0 -interp->weight[0][0];
  // Dimension 2
  BinarySearch(interp->index[1], interp->weight[1][0], x2, 0, n2-1, point[1]);
  interp->weight[1][1] = 1.0 -interp->weight[1][0];
  // Dimension 3
  BinarySearch(interp->index[2], interp->weight[2][0], x3, 0, n3-1, point[2]);
  interp->weight[2][1] = 1.0 -interp->weight[2][0];
  // Dimension 4
  BinarySearch(interp->index[3], interp->weight[3][0], x4, 0, n4-1, point[3]);
  interp->weight[3][1] = 1.0 -interp->weight[3][0];
}

void Table4D::InterpolatePoint(void) {
  // Dimension 1
  ComputeIndexAndWeight(interp->index[0], interp->weight[0][0], x1, n1, point[0], interp_mode[0]);
  interp->weight[0][1] = 1.0 -interp->weight[0][0];
  // Dimension 2
  ComputeIndexAndWeight(interp->index[1], interp->weight[1][0], x2, n2, point[1], interp_mode[1]);
  interp->weight[1][1] = 1.0 -interp->weight[1][0];
  // Dimension 3
  ComputeIndexAndWeight(interp->index[2], interp->weight[2][0], x3, n3, point[2], interp_mode[2]);
  interp->weight[2][1] = 1.0 -interp->weight[2][0];
  // Dimension 4
  ComputeIndexAndWeight(interp->index[3], interp->weight[3][0], x4, n4, point[3], interp_mode[3]);
  interp->weight[3][1] = 1.0 -interp->weight[3][0];
}

void Table4D::InterpolatePoint(double coor1, double coor2, double coor3, double coor4) {
  if ((coor1 != point[0])||
      (coor2 != point[1])||
      (coor3 != point[2])||
      (coor4 != point[3])) {

    SetInterpPoint(coor1, coor2, coor3, coor4);
    InterpolatePoint();

  }
}

void Table4D::InterpolatePoint_old(double coor1, double coor2, double coor3, double coor4) {
  if ((coor1 != point[0])||
      (coor2 != point[1])||
      (coor3 != point[2])||
      (coor4 != point[3])) {

    SetInterpPoint(coor1, coor2, coor3, coor4);
    InterpolatePoint_old();

  }
}

double Table4D::Interpolate(double coor1, double coor2, double coor3, double coor4, int ivar) {
  InterpolatePoint(coor1, coor2, coor3, coor4);
  if ((ivar<0)||(ivar>=nvar)){
    cerr << "Failed to access variable number " << ivar<< " , nvar= "<<nvar<<endl;
    throw(-1);
  }
  return Interpolate(ivar);
}

double Table4D::Interpolate(int ivar) {
  double g[2], f[2], h[2];
  double val = 0.0;
  for (int l=0; l<2; ++l) {
    h[l] = 0.0;
    for (int k=0; k < 2; ++k) {
      g[k] = 0.0;
      for (int j= 0; j<2; ++j) {
        f[j] = 0.0;
        for (int i = 0; i < 2; ++i)
          f[j] += interp->weight[0][i] * (double) Data[interp->index[0]+i]
                                                       [interp->index[1]+j]
                                                        [interp->index[2]+k]
                                                         [interp->index[3]+l][ivar];
        g[k] += interp->weight[1][j] * f[j];
      }
      h[l] += interp->weight[2][k] * g[k];
    }
    val += interp->weight[3][l] * h[l];
  }

  return val;
}

void Table4D::ComputeIndexAndWeight(int &index, double &weight, double *vec, int size, double coor, int mode) {
  switch (mode) {
    case 0 : // Look for it
      BinarySearch(index, weight, vec, 0, size-1, coor);
      break;
    case 1 : // The mesh law is linear
      ComputeIndexAndWeightLin(index, weight, vec, size, coor);
      break;
    case 2 : // The mesh law is a power law
      ComputeIndexAndWeightPow(index, weight, vec, size, coor, 2.7); // default power law value is 2.7
      break;
    case 3:
      ComputeIndexAndWeightTwoLin(index, weight, vec, size, coor, 5.0e4, 5.0e5, 1.5e5 ,11); // default P_middle = 1.5e5 ; iz_fin = 21
      break;
    case 4:
      ComputeIndexAndWeightLinAndStretched(index, weight, vec, size, coor, 0.03 , size/3); // default Zst = 0.03 ; icut = nz/3
      break;
    default:
      cout << "Wrong interp mode" << endl;
      throw(-1);
  }
}

void Table4D::ComputeIndexAndWeightLin(int &index, double &weight, double *vec, int size, double coor) {
  if (coor <= vec[0]) {
    index = 0;
    weight = 1.0;
    return;
  }
  if (coor >= vec[size-1]) {
    index = size - 2;
    weight = 0.0;
    return;
  }

  // Now coor is strictly inside the domain
  index = int ( coor*( double (size)-1.0) );
  weight = (vec[index+1]-coor)/(vec[index+1]-vec[index]);
}

void Table4D::ComputeIndexAndWeightTwoLin(int &index, double &weight, double *vec, int size, double coor,
                                          double z_min, double z_max, double z_cut, int i_cut) {
  if (coor <= vec[0]) {
    index = 0;
    weight = 1.0;
    return;
  }
  if (coor >= vec[size-1]) {
    index = size - 2;
    weight = 0.0;
    return;
  }

  // Now coor is strictly inside the domain
  if (coor <= z_cut) {
    // Linear mesh for [0,Zcut]
    //index = int ( coor/z_cut*(double (i_cut)-1.0));
    index = int ( (coor-z_cut)/(z_min-z_cut)*(double (i_cut)-1.0) );
  } else {
    // Linear mesh for ]Zcut, 1]
    //index = int ((coor-z_cut)/(1.0-z_cut)*(double (size-i_cut)) ) + i_cut -1;
    index = int ( (coor-z_cut)/(z_max-z_cut)*(double (size - i_cut)-1.0) );
  }
  weight = (vec[index+1]-coor)/(vec[index+1]-vec[index]);
}

void Table4D::ComputeIndexAndWeightLinAndStretched(int &index, double &weight, double *vec, int size, double coor,
                                                   double z_cut, int i_cut) {
  if (coor <= vec[0]) {
    index = 0;
    weight = 1.0;
    return;
  }
  if (coor >= vec[size-1]) {
    index = size - 2;
    weight = 0.0;
    return;
  }

  // Now coor is strictly inside the domain
  if (coor <= z_cut) { // Linear mesh for [0,Zcut]
    index = int ( coor/z_cut*(double (i_cut)) );
  } else {
    // Mesh with linear growth to reach Z=1
    // Z(i) = (i-i_cut)*dz + 0.5*(i-i_cut)(i-i_cut-1)*alpha + Z_cut
    // alpha from: Zm(size-1) = 1.0
    double dz = z_cut/ (double (i_cut));
    double m = (double) (size - 1 - i_cut);
    double r = 1.0 - z_cut;
    double alpha = 2.0 * (r - m * dz) / (m * (m - 1.0));
    double delta = (dz-0.5*alpha)*(dz-0.5*alpha)-2.0*alpha*(z_cut-coor); // b^2 - 4ac
    double a = (0.5*alpha-dz+sqrt(delta))/alpha;// (-b + sqrt(delta))/(2a)
    index = int (a) + i_cut;
  }
  weight = (vec[index+1]-coor)/(vec[index+1]-vec[index]);
}

void Table4D::ComputeIndexAndWeightPow(int &index, double &weight, double *vec, int size, double coor, double power) {
  if (coor <= vec[0]) {
    index = 0;
    weight = 1.0;
    return;
  }
  if (coor >= vec[size-1]) {
    index = size - 2;
    weight = 0.0;
    return;
  }

  // Now coor is strictly inside the domain
  index = int ( (pow(coor,1.0/power)*( double (size)-1.0)) );
  weight = (vec[index+1]-coor)/(vec[index+1]-vec[index]);
}

void Table4D::BCAST(int mpi_root, MPI_Comm mpi_comm) {
  MPI_Bcast(&n1,1,MPI_INT,mpi_root,mpi_comm);
  MPI_Bcast(&n2,1,MPI_INT,mpi_root,mpi_comm);
  MPI_Bcast(&n3,1,MPI_INT,mpi_root,mpi_comm);
  MPI_Bcast(&n4,1,MPI_INT,mpi_root,mpi_comm);
  MPI_Bcast(&nvar,1,MPI_INT,mpi_root,mpi_comm);
  if (mpi_rank != 0) Allocate(n1, n2, n3, n4, nvar);
  // Wait for everyone to finish allocation
  MPI_Barrier(mpi_comm);
  // Send coordinate arrays
  MPI_Bcast(x1,n1,MPI_DOUBLE,mpi_root,mpi_comm);
  MPI_Bcast(x2,n2,MPI_DOUBLE,mpi_root,mpi_comm);
  MPI_Bcast(x3,n3,MPI_DOUBLE,mpi_root,mpi_comm);
  MPI_Bcast(x4,n4,MPI_DOUBLE,mpi_root,mpi_comm);
  // Send Multi-dimensional array
  for (int i1=0; i1<n1; ++i1)
    for (int i2=0; i2<n2; ++i2)
      for (int i3=0; i3<n3; ++i3)
        for (int i4=0; i4<n4; ++i4)
          MPI_Bcast(Data[i1][i2][i3][i4],nvar,MPI_FLOAT,mpi_root,mpi_comm);
}
