#include "Table3D.h"

Table3D::Table3D()
{
  Init();
  interp = new InterpolationIndex(3,2);
  Allocated = false;
}

Table3D::Table3D(int i1, int i2, int i3, int i4)
{
  Init();
  Allocate(i1, i2, i3, i4);
  interp = new InterpolationIndex(3,2);
}
//------------------------------//
//         Destructor           //
//------------------------------// 

Table3D::~Table3D() { DestroyTable(); }
  
//------------------------------//
//    De-/Allocation Methods    //
//------------------------------//
void Table3D::Allocate(int i1, int i2, int i3, int i4)
{
  SetDimensions(i1, i2, i3, i4);
  Allocate();
}

void Table3D::Allocate()
{
  getMem1D(&x1, 0, n1-1, "Chemtable::Load x1", true);
  getMem1D(&x2, 0, n2-1, "Chemtable::Load x2", true);
  getMem1D(&x3, 0, n3-1, "Chemtable::Load x3", true);
  getMem4D(&Data, 0, n1-1, 0, n2-1, 0, n3-1, 0, nvar-1, "Chemtable::Load Data", true);
  SetAllocated();
}

void Table3D::DestroyTable()
{
  if (x1 != NULL)       {freeMem1D(x1, 0, n1-1); x1 = NULL;}
  if (x2 != NULL)       {freeMem1D(x2, 0, n2-1); x2 = NULL;}
  if (x3 != NULL)       {freeMem1D(x3, 0, n3-1); x3 = NULL;}
  if (Data != NULL)     {freeMem4D(Data, 0, n1-1, 0, n2-1, 0, n3-1, 0, nvar-1);   Data = NULL;}
}
//------------------------------//
//     Check Bound Methods      //
//------------------------------//
void Table3D::SetAllocated() { if ( (x1 != NULL)&&
                                    (x2 != NULL)&&
                                    (x3 != NULL)&&
                                    (Data != NULL) ) Allocated = true;}

bool Table3D::IsAllocated()
{
  SetAllocated();
  return Allocated;
}

void Table3D::Init(void) {
  n1 = 0;
  n2 = 0;
  n3 = 0;
  nvar =0;

  x1   = NULL;
  x2   = NULL;
  x3   = NULL;
  Data = NULL;

  point[0] = -1.0;
  point[1] = -1.0;
  point[2] = -1.0;

  interp_mode[0] = 0; // Linear mesh followed by linear stretching for Zm
  interp_mode[1] = 0; // Power law for Zvar which is Sz in fact
  interp_mode[2] = 0; // linear law for C
}


void Table3D::CopyCoordinate(int n, double* x, double* copy) 
{
  for (int i=0; i<n; ++i) { copy[i] = x[i]; }
}

void Table3D::CopyCoordinate1(double* x)  { CopyCoordinate(n1,x1,x); }
void Table3D::CopyCoordinate2(double* x)  { CopyCoordinate(n2,x2,x); }
void Table3D::CopyCoordinate3(double* x)  { CopyCoordinate(n3,x3,x); }


// GetArrayInfo
void Table3D::GetArrayInfo(int ivar, float &minval, float &maxval,
                           int &mi1, int &mi2, int &mi3,
                           int &ma1, int &ma2, int &ma3)
{
  GetArrayInfo3D(Data, ivar, n1, n2, n3,
                 minval, maxval,
                 mi1, mi2, mi3,
                 ma1, ma2, ma3);

}

void Table3D::InterpolatePoint_old(void) {
  // Dimension 1
  BinarySearch(interp->index[0], interp->weight[0][0], x1, 0, n1-1, point[0]);
  interp->weight[0][1] = 1.0 -interp->weight[0][0];
  // Dimension 2
  BinarySearch(interp->index[1], interp->weight[1][0], x2, 0, n2-1, point[1]);
  interp->weight[1][1] = 1.0 -interp->weight[1][0];
  // Dimension 3
  BinarySearch(interp->index[2], interp->weight[2][0], x3, 0, n3-1, point[2]);
  interp->weight[2][1] = 1.0 -interp->weight[2][0];
}

void Table3D::InterpolatePoint(void) {
  // Dimension 1
  ComputeIndexAndWeight(interp->index[0], interp->weight[0][0], x1, n1, point[0], interp_mode[0]);
  interp->weight[0][1] = 1.0 -interp->weight[0][0];
  // Dimension 2
  ComputeIndexAndWeight(interp->index[1], interp->weight[1][0], x2, n2, point[1], interp_mode[1]);
  interp->weight[1][1] = 1.0 -interp->weight[1][0];
  // Dimension 3
  ComputeIndexAndWeight(interp->index[2], interp->weight[2][0], x3, n3, point[2], interp_mode[2]);
  interp->weight[2][1] = 1.0 -interp->weight[2][0];
}

void Table3D::InterpolatePoint(double coor1, double coor2, double coor3) {
  if ((coor1 != point[0])||
      (coor2 != point[1])||
      (coor3 != point[2])) {

    SetInterpPoint(coor1, coor2, coor3);
    InterpolatePoint();

  }
}

void Table3D::InterpolatePoint_old(double coor1, double coor2, double coor3) {
  if ((coor1 != point[0])||
      (coor2 != point[1])||
      (coor3 != point[2])) {

    SetInterpPoint(coor1, coor2, coor3);
    InterpolatePoint_old();

  }
}

double Table3D::Interpolate(double coor1, double coor2, double coor3, int ivar) {
  InterpolatePoint(coor1, coor2, coor3);
  if ((ivar<0)||(ivar>=nvar)){
    cerr << "Failed to access variable number " << ivar<< " , nvar= "<<nvar<<endl;
    throw(-1);
  }
  return Interpolate(ivar);
}

double Table3D::Interpolate(int ivar) {
  double g[2], f[2];
  double val = 0.0;
  for (int k=0; k < 2; ++k) {
    g[k] = 0.0;
    for (int j= 0; j<2; ++j) {
      f[j] = 0.0;
      for (int i = 0; i < 2; ++i)
        f[j] += interp->weight[0][i] * (double) Data[interp->index[0]+i]
                                                     [interp->index[1]+j]
                                                      [interp->index[2]+k][ivar];
      g[k] += interp->weight[1][j] * f[j];
    }
    val += interp->weight[2][k] * g[k];
  }

  return val;
}

void Table3D::ComputeIndexAndWeight(int &index, double &weight, double *vec, int size, double coor, int mode) {
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

void Table3D::ComputeIndexAndWeightLin(int &index, double &weight, double *vec, int size, double coor) {
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

void Table3D::ComputeIndexAndWeightTwoLin(int &index, double &weight, double *vec, int size, double coor,
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

void Table3D::ComputeIndexAndWeightLinAndStretched(int &index, double &weight, double *vec, int size, double coor,
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

void Table3D::ComputeIndexAndWeightPow(int &index, double &weight, double *vec, int size, double coor, double power) {
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

void Table3D::BCAST(int mpi_root, MPI_Comm mpi_comm) {
  MPI_Bcast(&n1,1,MPI_INT,mpi_root,mpi_comm);
  MPI_Bcast(&n2,1,MPI_INT,mpi_root,mpi_comm);
  MPI_Bcast(&n3,1,MPI_INT,mpi_root,mpi_comm);
  MPI_Bcast(&nvar,1,MPI_INT,mpi_root,mpi_comm);
  if (mpi_rank != 0) Allocate(n1, n2, n3, nvar);
  // Wait for everyone to finish allocation
  MPI_Barrier(mpi_comm);
  // Send coordinate arrays
  MPI_Bcast(x1,n1,MPI_DOUBLE,mpi_root,mpi_comm);
  MPI_Bcast(x2,n2,MPI_DOUBLE,mpi_root,mpi_comm);
  MPI_Bcast(x3,n3,MPI_DOUBLE,mpi_root,mpi_comm);
  // Send Multi-dimensional array
  for (int i1=0; i1<n1; ++i1)
    for (int i2=0; i2<n2; ++i2)
      for (int i3=0; i3<n3; ++i3)
        MPI_Bcast(Data[i1][i2][i3],nvar,MPI_FLOAT,mpi_root,mpi_comm);
}

void Table3D::CheckNan(void) {
	for (int i=0; i< n1; ++i) {if (x1[i] != x1[i]) { cout << "Nan in x1 at i=" <<i<<endl; throw(-1); } }
	for (int i=0; i< n2; ++i) {if (x2[i] != x2[i]) { cout << "Nan in x2 at i=" <<i<<endl; throw(-1); } }
  for (int i=0; i< n3; ++i) {if (x3[i] != x3[i]) { cout << "Nan in x3 at i=" <<i<<endl; throw(-1); } }
  for (int i1=0; i1<n1; ++i1)
    for (int i2=0; i2<n2; ++i2)
      for (int i3=0; i3<n3; ++i3)
        for (int ivar=0; ivar<nvar; ++ivar) {
      	  if (Data[i1][i2][i3][ivar] != Data[i1][i2][i3][ivar]) {
      	  	cout << "Nan in Data at (" <<i1 << "," << i2 << "," << i3<<"," << ivar<<")"<<endl;
      	  	throw(-1);
      	  }
        }
}
