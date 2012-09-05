#include "UgpWithCvFakeOp.h"
#include <math.h>

// The following routine is required to build the face-based operators.
// It is part of the lapack distribution -- the location of this library
// must be set in the LAPACK_LIB definition in Makefile.in

/*
  SUBROUTINE DGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK, INFO )

  *  DGGLSE solves the linear equality-constrained least squares (LSE)
  *  problem:
  *
  *          minimize || c - A*x ||_2   subject to   B*x = d
  *
  *  where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
  *  M-vector, and d is a given P-vector. It is assumed that
  *  P <= N <= M+P, and
  *
  *           rank(B) = P and  rank( (A) ) = N.
  *                                ( (B) )
  *
  *  These conditions ensure that the LSE problem has a unique solution,
  *  which is obtained using a generalized RQ factorization of the
  *  matrices (B, A) given by
  *
  *     B = (0 R)*Q,   A = Z*T*Q.

*/

extern "C" void dgglse_(int * M,int * N,int * P,double * A,int * LDA,
			double * B,int * LDB,double * C,double * D,
			double * X,double * WORK,int * LWORK,int * INFO);

void UgpWithCvFakeOp::buildOperators() {
  
  if (mpi_rank == 0)
    cout << "buildOperators()" << endl;

  /*
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // locate ifa_debug...
  ifa_debug = -1;
  // prism 64 grid...
  //double xdebug[3] = { -1.6238, 1.80, 0.0 };
  // cart 64 grid...
  double xdebug[3] = { 0.0, 0.078125, 0.0 };
  FOR_IFA {
  double d2 = 0.0;
  FOR_I3 {
  double dx = xdebug[i] - x_fa[ifa][i];
  d2 += dx*dx;
  }
  if (sqrt(d2) < 1.0E-2) {
  assert( ifa_debug == -1 );
  ifa_debug = ifa;
  }
  }
  if (ifa_debug != -1)
  cout << "found ifa_debug on rank: " << mpi_rank << " location: " << 
  x_fa[ifa_debug][0] << " " << x_fa[ifa_debug][1] << " " << x_fa[ifa_debug][2] << endl;
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  */

  cvofa0_i = new int[nfa+1];
  cvofa1_i = new int[nfa+1];

  for (int iter = 0; iter < 2; ++iter) {
    cvofa0_i[0] = 0;
    cvofa1_i[0] = 0;
    FOR_IFA {

      const int icv0 = cvofa[ifa][0];
      assert( (icv0 >= 0)&&(icv0 < ncv) );
      const int icv1 = cvofa[ifa][1];
      assert( (icv1 >= 0)&&(icv1 < ncv_gf) );
      
      {

	// use icv0 and its immediate nbrs to build the operator, PLUS
	// any fake cells associated with boundary faces of icv0...
	const int noc_f = nbocv_i[icv0];
	const int noc_l = nbocv_i[icv0+1]-1;
	cvofa0_i[ifa+1] = cvofa0_i[ifa] + noc_l - noc_f + 1;
	int cof = cvofa0_i[ifa];
	if (iter == 1) {
	  // icv0 is the diagonal - add it first
	  assert( nbocv_v[noc_f] == icv0 ); 
	  cvofa0_v[cof++] = icv0;
	  // icv1 is second...
	  cvofa0_v[cof++] = icv1;
	  // then the rest...
	  for (int noc = noc_f+1; noc <= noc_l; ++noc) { // skip diag
	    const int icv = nbocv_v[noc];
	    if (icv != icv1)
	      cvofa0_v[cof++] = icv;
	  }
	}
	
	// also loop on faces of icv0 and, for any boundary faces,
	// include the fake cv in the operator...
	const int foc_f = faocv_i[icv0];
	const int foc_l = faocv_i[icv0+1]-1;
	for (int foc = foc_f; foc <= foc_l; ++foc) {
	  const int ifa_nbr = faocv_v[foc];
	  if (ifa_nbr < nfa_b) {
	    assert( cvofa[ifa_nbr][0] == icv0 );
	    assert( (cvofa[ifa_nbr][1] >= ncv_g)&&(cvofa[ifa_nbr][1] < ncv_gf) );
	    cvofa0_i[ifa+1] += 1;
	    if (iter == 1) {
	      const int icv = cvofa[ifa_nbr][1];
	      // might already be added as icv1 for face ifa...
	      if (icv != icv1)
		cvofa0_v[cof++] = cvofa[ifa_nbr][1];
	    }
	  }
	}
      }
	
      if ( (icv1 >= 0)&&(icv1 < ncv) ) {
	assert( ifa >= nfa_bpi );
	const int noc_f = nbocv_i[icv1];
	const int noc_l = nbocv_i[icv1+1]-1;
	cvofa1_i[ifa+1] = cvofa1_i[ifa] + noc_l - noc_f + 1;
	int cof = cvofa1_i[ifa];
	if (iter == 1) {
	  // icv1 is the diagonal - add it first
	  assert( nbocv_v[noc_f] == icv1 ); 
	  cvofa1_v[cof++] = icv1;
	  // icv0 is second...
	  cvofa1_v[cof++] = icv0;
	  // then the rest - skip the diag... 
	  for (int noc = noc_f+1; noc <= noc_l; ++noc) {
	    const int icv = nbocv_v[noc];
	    if (icv != icv0)
	      cvofa1_v[cof++] = icv;
	  }
	}
	// also loop on faces of icv1 and, for any boundary faces,
	// include the fake cv in the operator...
	const int foc_f = faocv_i[icv1];
	const int foc_l = faocv_i[icv1+1]-1;
	for (int foc = foc_f; foc <= foc_l; ++foc) {
	  const int ifa_nbr = faocv_v[foc];
	  if (ifa_nbr < nfa_b) {
	    assert( cvofa[ifa_nbr][0] == icv1 );
	    assert( (cvofa[ifa_nbr][1] >= ncv_g)&&(cvofa[ifa_nbr][1] < ncv_gf) );
	    cvofa1_i[ifa+1] += 1;
	    if (iter == 1) {
	      const int icv = cvofa[ifa_nbr][1];
	      assert( icv != icv0 );
	      cvofa1_v[cof++] = cvofa[ifa_nbr][1];
	    }
	  }
	}
      }
      else {
	// there is no 1_ operator for faces that do not have their
	// icv1 local...
	assert( ifa < nfa_bpi );
	cvofa1_i[ifa+1] = cvofa1_i[ifa];
	assert( cvofa1_i[ifa+1] == 0 ); // should all be at the start of the list
      }
    }
    
    if (iter == 0) {
      cvofa0_s = cvofa0_i[nfa];
      cvofa0_v = new int[cvofa0_s];
      cvofa1_s = cvofa1_i[nfa];
      cvofa1_v = new int[cvofa1_s];
    }
    
  }
  
  // now the value reconstruction operators...
  
  // to build these operators so that they are constrained by the finite
  // volume assumtions that the polynomial integrated over the cell is
  // the same as the mean value, we need a way to perform these integrations.
  // For cells we own, this is straight-forward, because we have all the data
  // associated with cell nodes and can use sub-tet integration. For ghost
  // (and eventually fake) cells, we do not know, and we need to push this data
  // over in temporary arrays...

  // --------------------------------
  // faces of ghost cells...
  // --------------------------------

  // put -1 in non-boundary cv's...
  for (int icv = ncv_bp; icv < ncv; ++icv) {      
    cv_flag[icv] = -1;
  }
  for (int icv = ncv; icv < ncv_g; icv++) {
    cv_flag[icv] = -1;
  }
    
  // cv's 0..ncv_bp-1 are potentially ghosts (interprocessor/periodic)... 
  for (int icv = 0; icv < ncv_bp; ++icv) {
    cv_flag[icv] = faocv_i[icv+1] - faocv_i[icv];
    assert( cv_flag[icv] >= 4 ); // should be atleast 4 face nbrs (i.e. tet)...
  }
    
  updateCvData(cv_flag,REPLACE_DATA);
    
  int * faocv_g_i = new int[ncv_g-ncv+1];
  faocv_g_i[0] = 0;
    
  // the number of faces per ghost cell...
  int my_nfoc_max = 0;

  for (int icv = ncv; icv < ncv_g; icv++) {
    assert( cv_flag[icv] >= 4 ); // should be atleast 4 face nbrs (i.e. tet)...
    faocv_g_i[icv-ncv+1] =  faocv_g_i[icv-ncv] + cv_flag[icv];
    my_nfoc_max = max( my_nfoc_max, cv_flag[icv] );
  }
      
  int nfoc_max;
  MPI_Allreduce(&my_nfoc_max,&nfoc_max,1,MPI_INT,MPI_MAX,mpi_comm);
    
  //if (mpi_rank == 0)
  //  cout << " > nfoc_max: " << nfoc_max << endl;
    
  int nfa_g = faocv_g_i[ncv_g-ncv];

  // --------------------------------
  // all ghost face coordinates are exchanged
  // thru ghost cells...
  // --------------------------------
   
  double (*x_fa_g)[3] = new double[nfa_g][3];
  int * noofa_g_i = new int[nfa_g+1];
    
  double (*r2_tmp)[3] = new double[ncv_g][3];
    
  for (int foc = 0; foc < nfoc_max; ++foc) {
      
    // step 1: populate r2_tmp with the relevant face centroid... 
    for (int icv = 0; icv < ncv_bp; ++icv) {
      const int nfoc = faocv_i[icv+1]-faocv_i[icv];
      if (foc < nfoc) {
	const int ifa = faocv_v[faocv_i[icv]+foc];
	FOR_I3 r2_tmp[icv][i] = x_fa[ifa][i];
	// put the nof count into cv_flag...
	cv_flag[icv] = noofa_i[ifa+1]-noofa_i[ifa];
	assert( cv_flag[icv] >= 3 ); // must be atleast 3 nodes
      }
    }

    // step 2: exchange...
    updateCvData(r2_tmp,REPLACE_TRANSLATE_DATA); // apply any translation to this coordinate data
    updateCvData(cv_flag,REPLACE_DATA);
      
    // step 3: populate x_fa_g...
    for (int icv = ncv; icv < ncv_g; icv++) {
      const int nfoc = faocv_g_i[icv-ncv+1]-faocv_g_i[icv-ncv];
      if (foc < nfoc) {
	const int ifa = faocv_g_i[icv-ncv]+foc; // "ghost" faces are just contiguous
	FOR_I3 x_fa_g[ifa][i] = r2_tmp[icv][i];
	noofa_g_i[ifa+1] = cv_flag[icv]; // put into ifa+1
	assert( cv_flag[icv] >= 3 ); // must be atleast 3 nodes
      }
    }
      
  }

  // complete noofa_g_i...
  int my_nnof_max = 0;
  noofa_g_i[0] = 0;
  for (int ifa = 0; ifa < nfa_g; ifa++) {
    my_nnof_max = max( my_nnof_max, noofa_g_i[ifa+1] );
    noofa_g_i[ifa+1] += noofa_g_i[ifa];
  }
  int noofa_g_s = noofa_g_i[nfa_g];
  int nnof_max;
  MPI_Allreduce(&my_nnof_max,&nnof_max,1,MPI_INT,MPI_MAX,mpi_comm);

  //if (mpi_rank == 0)
  //  cout << " > nnof_max: " << nnof_max << endl;
    
  // --------------------------------
  // nodes of ghost cells...
  // --------------------------------
    
  for (int icv = 0; icv < ncv_bp; ++icv) {
    cv_flag[icv] = noocv_i[icv+1] - noocv_i[icv];
    assert( cv_flag[icv] >= 4 ); // should be atleast 4 nodes (i.e. tet)...
  }
  updateCvData(cv_flag,REPLACE_DATA);
    
  int * noocv_g_i = new int[ncv_g-ncv+1];
  noocv_g_i[0] = 0;
    
  // the number of nodes per ghost cell...
  int my_nnoc_max = 0;
    
  for (int icv = ncv; icv < ncv_g; icv++) {
    assert( cv_flag[icv] >= 4 ); // should be atleast 4 nodes (i.e. tet)...
    noocv_g_i[icv-ncv+1] =  noocv_g_i[icv-ncv] + cv_flag[icv];
    my_nnoc_max = max( my_nnoc_max, cv_flag[icv] );
  }
      
  int nnoc_max;
  MPI_Allreduce(&my_nnoc_max,&nnoc_max,1,MPI_INT,MPI_MAX,mpi_comm);

  //if (mpi_rank == 0)
  //  cout << " > nnoc_max: " << nnoc_max << endl;
    
  int nno_g = noocv_g_i[ncv_g-ncv];

  // --------------------------------
  // all ghost node coordinates are exchanged
  // thru ghost cells...
  // --------------------------------
   
  double (*x_no_g)[3] = new double[nno_g][3];
    
  for (int noc = 0; noc < nnoc_max; ++noc) {
      
    // step 1: populate r2_tmp with the relevant node coord...
    for (int icv = 0; icv < ncv_bp; ++icv) {
      const int nnoc = noocv_i[icv+1]-noocv_i[icv];
      if (noc < nnoc) {
	const int ino = noocv_v[noocv_i[icv]+noc];
	FOR_I3 r2_tmp[icv][i] = x_no[ino][i];
      }
    }
      
    // step 2: exchange...
    updateCvData(r2_tmp,REPLACE_TRANSLATE_DATA); // apply any translation to this coordinate data
      
    // step 3: populate x_no_g...
    for (int icv = ncv; icv < ncv_g; icv++) {
      const int nnoc = noocv_g_i[icv-ncv+1]-noocv_g_i[icv-ncv];
      if (noc < nnoc) {
	const int ino = noocv_g_i[icv-ncv]+noc; // "ghost" nodes are just contiguous
	FOR_I3 x_no_g[ino][i] = r2_tmp[icv][i];
      }
    }
      
  }

  delete[] r2_tmp;

  // --------------------------------
  // finally we need the node-of-face 
  // connectivity for this temporary ghost 
  // data...
  // --------------------------------

  int * noofa_g_v = new int[noofa_g_s];

  for (int foc = 0; foc < nfoc_max; ++foc) {
    for (int nof = 0; nof < nnof_max; ++nof) {
	
      for (int icv = 0; icv < ncv_bp; ++icv) {
	const int nfoc = faocv_i[icv+1]-faocv_i[icv];
	if (foc < nfoc) {
	    
	  const int ifa = faocv_v[faocv_i[icv]+foc];
	  const int nnof = noofa_i[ifa+1]-noofa_i[ifa];
	  if (nof < nnof) {

	    // populate the no_flag with the local node number...
	    const int noc_f = noocv_i[icv];
	    const int noc_l = noocv_i[icv+1]-1;
	    for (int noc = noc_f; noc <= noc_l; ++noc)
	      no_flag[noocv_v[noc]] = noc-noc_f; // 0,1,2...
	      
	    if (cvofa[ifa][0] == icv) {
	      // face is already outward pointing...
	      const int ino = noofa_v[noofa_i[ifa]+nof];
	      cv_flag[icv] = no_flag[ino];
	    }
	    else {
	      assert( cvofa[ifa][1] == icv );
	      // flip face...
	      const int ino = noofa_v[noofa_i[ifa+1]-1-nof];
	      cv_flag[icv] = no_flag[ino];
	    }
		
	  }
	}
      }

      updateCvData(cv_flag,REPLACE_DATA);
	
      for (int icv = ncv; icv < ncv_g; icv++) {
	const int nfoc = faocv_g_i[icv-ncv+1]-faocv_g_i[icv-ncv];
	if (foc < nfoc) {
	  const int ifa = faocv_g_i[icv-ncv]+foc; // "ghost" faces are just contiguous
	  const int nnof = noofa_g_i[ifa+1]-noofa_g_i[ifa];
	  if (nof < nnof) {
	    noofa_g_v[noofa_g_i[ifa]+nof] = noocv_g_i[icv-ncv]+cv_flag[icv];
	  }
	}
      }

    }
  }
    
  // done - now build least squares problem. This requires integration
  // of the polynomials in each arbitrary volume - here we use 4-point
  // sub-tet quadrature...

  //cout << "1" << endl;

  // 4-point symmetric tet integration...
  // each point has weight 1/4...
    
  double tet_int_wgt_4[4][4];
    
  tet_int_wgt_4[0][0] = 0.5854101966249685;
  tet_int_wgt_4[0][1] = 0.1381966011250105;
  tet_int_wgt_4[0][2] = 0.1381966011250105;
  tet_int_wgt_4[0][3] = 0.1381966011250105;
    
  tet_int_wgt_4[1][1] = 0.5854101966249685;
  tet_int_wgt_4[1][2] = 0.1381966011250105;
  tet_int_wgt_4[1][3] = 0.1381966011250105;
  tet_int_wgt_4[1][0] = 0.1381966011250105;
    
  tet_int_wgt_4[2][2] = 0.5854101966249685;
  tet_int_wgt_4[2][3] = 0.1381966011250105;
  tet_int_wgt_4[2][0] = 0.1381966011250105;
  tet_int_wgt_4[2][1] = 0.1381966011250105;
    
  tet_int_wgt_4[3][3] = 0.5854101966249685;
  tet_int_wgt_4[3][0] = 0.1381966011250105;
  tet_int_wgt_4[3][1] = 0.1381966011250105;
  tet_int_wgt_4[3][2] = 0.1381966011250105;

  const int NCOF_MAX = 25; // 1 + 4*6

  // [A]{x} = {C} stores the least squares problem
  double A_q[(NCOF_MAX-2)*5];
  double A_l1[(NCOF_MAX-1)*4];
  double A_l2[(NCOF_MAX-2)*4];
  
  double C_q[NCOF_MAX-2];
  double C_l1[NCOF_MAX-1];
  double C_l2[NCOF_MAX-2];
  
  // to back out the cof for any row of A...
  int cof_of_A_q[NCOF_MAX-2];
  int cof_of_A_l1[NCOF_MAX-1];
  int cof_of_A_l2[NCOF_MAX-2];

  // [B]{x} = {D} stores the constraints...
  double B_q[2*5];
  double B_l1[1*4];
  double B_l2[2*4];

  double D_q[2];
  double D_l1[1];
  double D_l2[2];
  
  // to back out the cof for any row of B...
  int cof_of_B_q[2];
  int cof_of_B_l1[1];
  int cof_of_B_l2[2];

  // tmp arrays corresponding to the largest...
  double Atmp[NCOF_MAX*5];
  double Ctmp[NCOF_MAX];
  double Btmp[2*5];
  double Dtmp[2];

  // X array for answer...
  double X[5];

  // and Work array req'd by the DGGLSE routine... 
  const int LWORK_MAX = NCOF_MAX*(1+NCOF_MAX);
  double WORK[LWORK_MAX];
  int LWORK = LWORK_MAX;

  /*
   *  DGGLSE solves the linear equality-constrained least squares (LSE)
   *  problem:
   *
   *          minimize || c - A*x ||_2   subject to   B*x = d
   *
   *  where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
   *  M-vector, and d is a given P-vector. It is assumed that
   *  P <= N <= M+P, and
   *
   *           rank(B) = P and  rank( (A) ) = N.
   *                                ( (B) )
   *
   *  These conditions ensure that the LSE problem has a unique solution,
   *  which is obtained using a generalized RQ factorization of the
   *  matrices (B, A) given by
   *
   *     B = (0 R)*Q,   A = Z*T*Q.
   */
  
  // these are used to store the final operators...
  
  cvofa0_value = new double[cvofa0_s];
  cvofa1_value = new double[cvofa1_s];
  
  cvofa0_grad = new double[cvofa0_s][3];
  cvofa1_grad = new double[cvofa1_s][3];

  // but build 3 different operators for now...
  
  // the _q (quadratic) operator has 2 constraints - one for each cell
  // on either side of the face, and weighted LS elsewhere... 
  double * cvofa0_value_q = cvofa0_value; // point into the final operator for convenience.
  double * cvofa1_value_q = cvofa1_value;
  
  // the _l1 (linear 1) operator have a single constraint only on the 
  // relevant biased cell...
  double * cvofa0_value_l1 = new double[cvofa0_s];
  double * cvofa1_value_l1 = new double[cvofa1_s];

  // the _l2 (linear 2) operator has 2 constraints - one for each cell
  // on either side of the face, and weighted LS elsewhere. On a uniform
  // grid, this does not introduce any bias, so is used for the gradients...
  double * cvofa0_value_l2 = new double[cvofa0_s];
  double * cvofa1_value_l2 = new double[cvofa1_s];

  double (*cvofa0_grad_l2)[3] = cvofa0_grad;
  double (*cvofa1_grad_l2)[3] = cvofa1_grad;

  // this is A LOT of memory, but it allows us to optimize the
  // operators...

  //int my_op_err = 0; // uncomment to get really bad parts of the grid
  for (int ifa = 0; ifa < nfa; ++ifa) {
    
    // for recording operator errors...
    fa_flag[ifa] = 0;
    
    /*
      cout << "ifa: " << ifa << " cvofa[ifa]: " << cvofa[ifa][0] << " " << cvofa[ifa][1] << endl;
      cout << "XXXXX: " << x_fa[ifa][0] << " " << x_fa[ifa][1] << " " << x_fa[ifa][2] << endl; 
      getchar();
    */

    // build a local orthogonal basis, with the first
    // vector aligned with "n", XXXXX or perhaps "s"...
      
    // n-version...
    // USE THIS ONE?
    double e0[3] = { fa_normal[ifa][0], fa_normal[ifa][1], fa_normal[ifa][2] };
    {
      double mag = sqrt( e0[0]*e0[0] + e0[1]*e0[1] + e0[2]*e0[2] );
      FOR_I3 e0[i] /= mag;
    }

    /*
    // s-version...
    // DO NOT USE THIS VERSION - it was for testing only, and somewhat
    // counter-intuitively, the n-version does better, even on skewed 
    // hex meshes.
    double e0[3];
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];
      FOR_I3 e0[i] = x_cv[icv1][i]-x_cv[icv0][i];
      double mag = sqrt( e0[0]*e0[0] + e0[1]*e0[1] + e0[2]*e0[2] );
      FOR_I3 e0[i] /= mag;
    }
    */

    // the second and third vectors don't matter too much - they just have
    // to form a right-handed orthogonal set with e0... 
    // for the second vector, cross the first face edge with e0...
      
    double dx_ed[3]; 
    FOR_I3 dx_ed[i] = 
      x_no[noofa_v[noofa_i[ifa]+1]][i] - 
      x_no[noofa_v[noofa_i[ifa]]][i];
    double e1[3] = {
      dx_ed[1]*e0[2] - dx_ed[2]*e0[1],
      dx_ed[2]*e0[0] - dx_ed[0]*e0[2],
      dx_ed[0]*e0[1] - dx_ed[1]*e0[0] };
    {
      double mag = sqrt( e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2] );
      FOR_I3 e1[i] /= mag;
    }
      
    // and the third is just e0 x e1...

    double e2[3] = {
      e0[1]*e1[2] - e0[2]*e1[1],
      e0[2]*e1[0] - e0[0]*e1[2],
      e0[0]*e1[1] - e0[1]*e1[0] };
    {
      double mag = sqrt( e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2] );
      // should not require normalization...
      assert( fabs(mag - 1.0) < 1.0E-12 );
      FOR_I3 e2[i] /= mag;
    }

    // ===============================================
    // left side operator cvofa0_value...
    // ===============================================
      
    {
	
      const int cof_f = cvofa0_i[ifa];
      const int cof_l = cvofa0_i[ifa+1]-1;

      // figure out the scaling associated with each vector to keep the system 
      // reasonably conditioned...
      double delta0 = 0.0;
      double delta1 = 0.0;
      double delta2 = 0.0;
      for (int cof = cof_f; cof <= cof_l; ++cof) {
	const int icv = cvofa0_v[cof];
	double dx = x_cv[icv][0] - x_fa[ifa][0];
	double dy = x_cv[icv][1] - x_fa[ifa][1];
	double dz = x_cv[icv][2] - x_fa[ifa][2];
	double dxp = dx*e0[0] + dy*e0[1] + dz*e0[2];
	double dyp = dx*e1[0] + dy*e1[1] + dz*e1[2];
	double dzp = dx*e2[0] + dy*e2[1] + dz*e2[2];
	delta0 = max(delta0,fabs(dxp));
	delta1 = max(delta1,fabs(dyp));
	delta2 = max(delta2,fabs(dzp));
      }
      assert( delta0 > 0.0 );
      assert( delta1 > 0.0 );
      assert( delta2 > 0.0 );

      int N4 = 4;
      int N5 = 5;    
      
      int P_q = 2; // 2 constraints...
      int M_q = cof_l - cof_f + 1 - P_q; // remainder goes into M	
      int ia_q = 0; // row location in the least squares matrix A or vector C
      int ib_q = 0; // row location in the constraint matrix B or vector D
      
      int P_l1 = 1; // 1 constraints...
      int M_l1 = cof_l - cof_f + 1 - P_l1; // remainder goes into M	
      int ia_l1 = 0; // row location in the least squares matrix A or vector C
      int ib_l1 = 0; // row location in the constraint matrix B or vector D

      int P_l2 = 2; // 2 constraints...
      int M_l2 = cof_l - cof_f + 1 - P_l2; // remainder goes into M	
      int ia_l2 = 0; // row location in the least squares matrix A or vector C
      int ib_l2 = 0; // row location in the constraint matrix B or vector D

      const int ncof = cof_l - cof_f + 1; 
      assert( ncof >= 5 ); // we always need atleast 5 cells to close the quadratic polynomial 
      assert( ncof <= NCOF_MAX );
	
      for (int cof = cof_f; cof <= cof_l; ++cof) {
	const int icv = cvofa0_v[cof];
	  
	/*
	  cout << "ifa has icv: " << icv << endl;
	  cout << "XXXXX: " << x_cv[icv][0] << " " << x_cv[icv][1] << " " << x_cv[icv][2] << endl; 
	*/
	  
	// the primed polynomial is phi = c0 + c1*xp + c2*yp + c3*zp + c4*xp*xp
	double coeff[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
	double weight = 1.0;
	  
	if (icv < ncv) {
	    
	  // local cv...
	  const int foc_f = faocv_i[icv];
	  const int foc_l = faocv_i[icv+1]-1;
	  for (int foc = foc_f; foc <= foc_l; ++foc) {
	    const int ifa_nbr = faocv_v[foc];
	    double dx_fa[3] = {
	      x_fa[ifa_nbr][0] - x_cv[icv][0],
	      x_fa[ifa_nbr][1] - x_cv[icv][1],
	      x_fa[ifa_nbr][2] - x_cv[icv][2] };
	    int nof_begin,nof_end,nof_inc;
	    if (cvofa[ifa_nbr][0] == icv) {
	      nof_begin = noofa_i[ifa_nbr];
	      nof_end = noofa_i[ifa_nbr+1];
	      nof_inc = 1;
	    }
	    else {
	      assert(cvofa[ifa_nbr][1] == icv);
	      nof_begin = noofa_i[ifa_nbr+1]-1;
	      nof_end = noofa_i[ifa_nbr]-1;
	      nof_inc = -1;
	    }
	    int ino1 = noofa_v[nof_end-nof_inc];
	    for (int nof = nof_begin; nof != nof_end; nof += nof_inc) {
	      int ino0 = ino1;
	      ino1 = noofa_v[nof];
	      double dx_no0[3] = {
		x_no[ino0][0] - x_cv[icv][0],
		x_no[ino0][1] - x_cv[icv][1],
		x_no[ino0][2] - x_cv[icv][2] };
	      double dx_no1[3] = {
		x_no[ino1][0] - x_cv[icv][0],
		x_no[ino1][1] - x_cv[icv][1],
		x_no[ino1][2] - x_cv[icv][2] };
	      double this_vol = 
		dx_fa[0]*(dx_no0[1]*dx_no1[2]-dx_no0[2]*dx_no1[1]) +
		dx_fa[1]*(dx_no0[2]*dx_no1[0]-dx_no0[0]*dx_no1[2]) +
		dx_fa[2]*(dx_no0[0]*dx_no1[1]-dx_no0[1]*dx_no1[0]);
	      assert( this_vol > 0.0 );
	      // cycle through the 4 integration points asssociated with this cell...
	      for (int ip = 0; ip < 4; ip++) {
		// the integration point vector relative to the flux face centroid is...
		double dx[3];
		FOR_I3 dx[i] = 
		  tet_int_wgt_4[ip][0]*x_no[ino0][i] + 
		  tet_int_wgt_4[ip][1]*x_no[ino1][i] + 
		  tet_int_wgt_4[ip][2]*x_fa[ifa_nbr][i] + 
		  tet_int_wgt_4[ip][3]*x_cv[icv][i] - 
		  x_fa[ifa][i];
		// convert this vector into the "primed" coordinate system...
		double xp = (dx[0]*e0[0] + dx[1]*e0[1] + dx[2]*e0[2])/delta0;
		double yp = (dx[0]*e1[0] + dx[1]*e1[1] + dx[2]*e1[2])/delta1;
		double zp = (dx[0]*e2[0] + dx[1]*e2[1] + dx[2]*e2[2])/delta2;
		// the polynomial coefficients...
		coeff[0] += this_vol;
		coeff[1] += this_vol*xp;
		coeff[2] += this_vol*yp;
		coeff[3] += this_vol*zp;
		coeff[4] += this_vol*xp*xp;
	      }
	    }
	  }
	    
	}
	else if (icv < ncv_g) {
	    
	  //cout << "ghost" << endl; 

	  // a ghost cv...
	  const int ifa_nbr_f = faocv_g_i[icv-ncv];
	  const int ifa_nbr_l = faocv_g_i[icv-ncv+1]-1;
	  for (int ifa_nbr = ifa_nbr_f; ifa_nbr <= ifa_nbr_l; ++ifa_nbr) {
	    double dx_fa[3] = {
	      x_fa_g[ifa_nbr][0] - x_cv[icv][0],
	      x_fa_g[ifa_nbr][1] - x_cv[icv][1],
	      x_fa_g[ifa_nbr][2] - x_cv[icv][2] };
	    const int nof_f = noofa_g_i[ifa_nbr];
	    const int nof_l = noofa_g_i[ifa_nbr+1]-1;
	    // loop on edges ino0-ino1...
	    int ino1 = noofa_g_v[nof_l];
	    for (int nof = nof_f; nof <= nof_l; ++nof) {
	      int ino0 = ino1;
	      ino1 = noofa_g_v[nof];
	      double dx_no0[3] = {
		x_no_g[ino0][0] - x_cv[icv][0],
		x_no_g[ino0][1] - x_cv[icv][1],
		x_no_g[ino0][2] - x_cv[icv][2] };
	      double dx_no1[3] = {
		x_no_g[ino1][0] - x_cv[icv][0],
		x_no_g[ino1][1] - x_cv[icv][1],
		x_no_g[ino1][2] - x_cv[icv][2] };
	      double this_vol = 
		dx_fa[0]*(dx_no0[1]*dx_no1[2]-dx_no0[2]*dx_no1[1]) +
		dx_fa[1]*(dx_no0[2]*dx_no1[0]-dx_no0[0]*dx_no1[2]) +
		dx_fa[2]*(dx_no0[0]*dx_no1[1]-dx_no0[1]*dx_no1[0]);
	      assert( this_vol > 0.0 );
	      // cycle through the 4 integration points asssociated with this cell...
	      for (int ip = 0; ip < 4; ip++) {
		// the integration point vector relative to the flux face centroid is...
		double dx[3];
		FOR_I3 dx[i] = 
		  tet_int_wgt_4[ip][0]*x_no_g[ino0][i] + 
		  tet_int_wgt_4[ip][1]*x_no_g[ino1][i] + 
		  tet_int_wgt_4[ip][2]*x_fa_g[ifa_nbr][i] + 
		  tet_int_wgt_4[ip][3]*x_cv[icv][i] - 
		  x_fa[ifa][i];
		// convert this vector into the "primed" coordinate system...
		double xp = (dx[0]*e0[0] + dx[1]*e0[1] + dx[2]*e0[2])/delta0;
		double yp = (dx[0]*e1[0] + dx[1]*e1[1] + dx[2]*e1[2])/delta1;
		double zp = (dx[0]*e2[0] + dx[1]*e2[1] + dx[2]*e2[2])/delta2;
		// the polynomial coefficients...
		// no need to include the factor of 0.25 here - we normalize
		// the coefficients below...
		coeff[0] += this_vol;
		coeff[1] += this_vol*xp;
		coeff[2] += this_vol*yp;
		coeff[3] += this_vol*zp;
		coeff[4] += this_vol*xp*xp;
	      }
	    }
	  }
	}
	else {

	  //cout << "fake" << endl;

	  // this is a fake cell that must be reflected through its boundary
	  // face normal. boundary faces and fake cells are ordered in the same way, 
	  // so we can get the relevant face as follows...
	    
	  assert( (icv >= ncv_g)&&(icv < ncv_gf) );
	  const int ifa_b = icv-ncv_g;
	  assert( cvofa[ifa_b][1] == icv );
	    
	  // get the unit normal - it is used to relect the integration points...
	  double area = sqrt( fa_normal[ifa_b][0]*fa_normal[ifa_b][0] + 
			      fa_normal[ifa_b][1]*fa_normal[ifa_b][1] + 
			      fa_normal[ifa_b][2]*fa_normal[ifa_b][2] );
	  double unitNormal[3] = { fa_normal[ifa_b][0]/area,
				   fa_normal[ifa_b][1]/area,
				   fa_normal[ifa_b][2]/area };
	    
	  // the near-boundary cv is...
	  int icv0 = cvofa[ifa_b][0];

	  // as a minimum geometric requirement on mesh quality, make sure the icv0-ifa_b vector
	  // dotted with the normal is positive...
	  assert( unitNormal[0]*(x_fa[ifa_b][0]-x_cv[icv0][0]) +
		  unitNormal[1]*(x_fa[ifa_b][1]-x_cv[icv0][1]) +
		  unitNormal[2]*(x_fa[ifa_b][2]-x_cv[icv0][2]) > 0.0 );
		    	    
	  /*
	    cout << "icv0 is: " << icv0 << endl;
	    cout << "n is: " << unitNormal[0] << " " << unitNormal[1] << " " << unitNormal[2] << endl;
	  */

	  const int foc_f = faocv_i[icv0];
	  const int foc_l = faocv_i[icv0+1]-1;
	  for (int foc = foc_f; foc <= foc_l; ++foc) {
	    const int ifa_nbr = faocv_v[foc];
	    double dx_fa[3] = {
	      x_fa[ifa_nbr][0] - x_cv[icv0][0],
	      x_fa[ifa_nbr][1] - x_cv[icv0][1],
	      x_fa[ifa_nbr][2] - x_cv[icv0][2] };
	    int nof_begin,nof_end,nof_inc;
	    if (cvofa[ifa_nbr][0] == icv0) {
	      nof_begin = noofa_i[ifa_nbr];
	      nof_end = noofa_i[ifa_nbr+1];
	      nof_inc = 1;
	    }
	    else {
	      assert(cvofa[ifa_nbr][1] == icv0);
	      nof_begin = noofa_i[ifa_nbr+1]-1;
	      nof_end = noofa_i[ifa_nbr]-1;
	      nof_inc = -1;
	    }
	    int ino1 = noofa_v[nof_end-nof_inc];
	    for (int nof = nof_begin; nof != nof_end; nof += nof_inc) {
	      int ino0 = ino1;
	      ino1 = noofa_v[nof];
	      double dx_no0[3] = {
		x_no[ino0][0] - x_cv[icv0][0],
		x_no[ino0][1] - x_cv[icv0][1],
		x_no[ino0][2] - x_cv[icv0][2] };
	      double dx_no1[3] = {
		x_no[ino1][0] - x_cv[icv0][0],
		x_no[ino1][1] - x_cv[icv0][1],
		x_no[ino1][2] - x_cv[icv0][2] };
	      double this_vol = 
		dx_fa[0]*(dx_no0[1]*dx_no1[2]-dx_no0[2]*dx_no1[1]) +
		dx_fa[1]*(dx_no0[2]*dx_no1[0]-dx_no0[0]*dx_no1[2]) +
		dx_fa[2]*(dx_no0[0]*dx_no1[1]-dx_no0[1]*dx_no1[0]);
	      assert( this_vol > 0.0 );
	      // cycle through the 4 integration points asssociated with this cell...
	      for (int ip = 0; ip < 4; ip++) {
		// the integration point vector relative to the reflected face centroid is...
		double dx[3];
		FOR_I3 dx[i] = 
		  tet_int_wgt_4[ip][0]*x_no[ino0][i] + 
		  tet_int_wgt_4[ip][1]*x_no[ino1][i] + 
		  tet_int_wgt_4[ip][2]*x_fa[ifa_nbr][i] + 
		  tet_int_wgt_4[ip][3]*x_cv[icv0][i] - 
		  x_fa[ifa_b][i]; // relative to boundary face initially.
		// relect each quadrature point through the boundary face normal...
		double delta = dx[0]*unitNormal[0] + dx[1]*unitNormal[1] + dx[2]*unitNormal[2];
		// this should be negative, unless the face is really twisted...
		/*
		  if ( delta >= 0.0 ) {
		  fa_flag[ifa] = 1;
		  my_op_err += 1;
		  }
		*/
		FOR_I3 dx[i] += x_fa[ifa_b][i] - x_fa[ifa][i] - 2.0*delta*unitNormal[i];
		// convert this vector into the "primed" coordinate system...
		double xp = (dx[0]*e0[0] + dx[1]*e0[1] + dx[2]*e0[2])/delta0;
		double yp = (dx[0]*e1[0] + dx[1]*e1[1] + dx[2]*e1[2])/delta1;
		double zp = (dx[0]*e2[0] + dx[1]*e2[1] + dx[2]*e2[2])/delta2;
		// the polynomial coefficients...
		coeff[0] += this_vol;
		coeff[1] += this_vol*xp;
		coeff[2] += this_vol*yp;
		coeff[3] += this_vol*zp;
		coeff[4] += this_vol*xp*xp;
	      }
	    }
	  }
	    
	}

	// normalize out the volume and multiply by the weight...
	// this gives an area-weighting to the imposed constraints,
	// which recovers the 1D behaviour for certain cases of 3D grids
	// with local refinement...
	  
	coeff[1] *= weight/coeff[0];
	coeff[2] *= weight/coeff[0];
	coeff[3] *= weight/coeff[0];
	coeff[4] *= weight/coeff[0];
	coeff[0] = weight;

	// _q and _l2 have the same cv's constrained...
	if ( (icv == cvofa[ifa][0])||(icv == cvofa[ifa][1]) ) {
	  // this is one of the constrained cv's next to the face, so populate B...
	  // _q...
	  for (int j = 0; j < 5; j++) B_q[ib_q+j*P_q] = coeff[j];
	  D_q[ib_q] = weight;
	  cof_of_B_q[ib_q++] = cof;
	  // _l2...
	  for (int j = 0; j < 4; j++) B_l2[ib_l2+j*P_l2] = coeff[j];
	  D_l2[ib_l2] = weight;
	  cof_of_B_l2[ib_l2++] = cof;
	}
	else {
	  // this is one of the least squares cv's, so populate A...
	  // _q
	  for (int j = 0; j < 5; j++) A_q[ia_q+j*M_q] = coeff[j];
	  C_q[ia_q] = weight;
	  cof_of_A_q[ia_q++] = cof;
	  // _l2...
	  for (int j = 0; j < 4; j++) A_l2[ia_l2+j*M_l2] = coeff[j];
	  C_l2[ia_l2] = weight;
	  cof_of_A_l2[ia_l2++] = cof;
	}
	
	// _l1 has only icv0 constrained for this left operator...
	if ( icv == cvofa[ifa][0] ) {
	  // this is one of the constrained cv's next to the face, so populate B...
	  // _l1...
	  for (int j = 0; j < 4; j++) B_l1[ib_l1+j*P_l1] = coeff[j];
	  D_l1[ib_l1] = weight;
	  cof_of_B_l1[ib_l1++] = cof;
	}
	else {
	  // this is one of the least squares cv's, so populate A...
	  // _l1...
	  for (int j = 0; j < 4; j++) A_l1[ia_l1+j*M_l1] = coeff[j];
	  C_l1[ia_l1] = weight;
	  cof_of_A_l1[ia_l1++] = cof;
	}


      }
	
      assert( ia_q == M_q );
      assert( ia_l1 == M_l1 );
      assert( ia_l2 == M_l2 );

      assert( ib_q == P_q );
      assert( ib_l1 == P_l1 );
      assert( ib_l2 == P_l2 );
	
      // _q quadratic reconstruction in normal direction, linear in other...

      for (int ii = 0; ii < M_q; ++ii) {
	for (int i = 0; i < M_q; ++i) Ctmp[i] = 0.0; 
	for (int i = 0; i < P_q; ++i) Dtmp[i] = 0.0;
	Ctmp[ii] = C_q[ii]; 
	for (int ij = 0; ij < M_q*5; ++ij) Atmp[ij] = A_q[ij];
	for (int ij = 0; ij < P_q*5; ++ij) Btmp[ij] = B_q[ij];
	int info;
	dgglse_(&M_q,&N5,&P_q,Atmp,&M_q,Btmp,&P_q,Ctmp,Dtmp,X,WORK,&LWORK,&info);
	assert( info == 0 );
	cvofa0_value_q[cof_of_A_q[ii]] = X[0];
      }
      for (int ii = 0; ii < P_q; ++ii) {
	for (int i = 0; i < M_q; ++i) Ctmp[i] = 0.0; 
	for (int i = 0; i < P_q; ++i) Dtmp[i] = 0.0;
	Dtmp[ii] = D_q[ii]; 
	for (int ij = 0; ij < M_q*5; ++ij) Atmp[ij] = A_q[ij];
	for (int ij = 0; ij < P_q*5; ++ij) Btmp[ij] = B_q[ij];
	int info;
	dgglse_(&M_q,&N5,&P_q,Atmp,&M_q,Btmp,&P_q,Ctmp,Dtmp,X,WORK,&LWORK,&info);
	assert( info == 0 );
	cvofa0_value_q[cof_of_B_q[ii]] = X[0];
      }
      
      // _l1 single constraint linear reconstruction in all directions...

      for (int ii = 0; ii < M_l1; ++ii) {
	for (int i = 0; i < M_l1; ++i) Ctmp[i] = 0.0; 
	for (int i = 0; i < P_l1; ++i) Dtmp[i] = 0.0;
	Ctmp[ii] = C_l1[ii]; 
	for (int ij = 0; ij < M_l1*4; ++ij) Atmp[ij] = A_l1[ij];
	for (int ij = 0; ij < P_l1*4; ++ij) Btmp[ij] = B_l1[ij];
	int info;
	dgglse_(&M_l1,&N4,&P_l1,Atmp,&M_l1,Btmp,&P_l1,Ctmp,Dtmp,X,WORK,&LWORK,&info);
	assert( info == 0 );
	cvofa0_value_l1[cof_of_A_l1[ii]] = X[0];
      }
      for (int ii = 0; ii < P_l1; ++ii) {
	for (int i = 0; i < M_l1; ++i) Ctmp[i] = 0.0; 
	for (int i = 0; i < P_l1; ++i) Dtmp[i] = 0.0;
	Dtmp[ii] = D_l1[ii]; 
	for (int ij = 0; ij < M_l1*4; ++ij) Atmp[ij] = A_l1[ij];
	for (int ij = 0; ij < P_l1*4; ++ij) Btmp[ij] = B_l1[ij];
	int info;
	dgglse_(&M_l1,&N4,&P_l1,Atmp,&M_l1,Btmp,&P_l1,Ctmp,Dtmp,X,WORK,&LWORK,&info);
	assert( info == 0 );
	cvofa0_value_l1[cof_of_B_l1[ii]] = X[0];
      }
      
      // _l2 double constraint linear reconstruction in all directions...
      // this introduces minimal bias, and is used for the gradient... 
      
      for (int ii = 0; ii < M_l2; ++ii) {
	for (int i = 0; i < M_l2; ++i) Ctmp[i] = 0.0; 
	for (int i = 0; i < P_l2; ++i) Dtmp[i] = 0.0;
	Ctmp[ii] = C_l2[ii]; 
	for (int ij = 0; ij < M_l2*4; ++ij) Atmp[ij] = A_l2[ij];
	for (int ij = 0; ij < P_l2*4; ++ij) Btmp[ij] = B_l2[ij];
	int info;
	dgglse_(&M_l2,&N4,&P_l2,Atmp,&M_l2,Btmp,&P_l2,Ctmp,Dtmp,X,WORK,&LWORK,&info);
	assert( info == 0 );
	cvofa0_value_l2[cof_of_A_l2[ii]] = X[0];
	FOR_I3 cvofa0_grad_l2[cof_of_A_l2[ii]][i] = X[1]*e0[i]/delta0 + X[2]*e1[i]/delta1 + X[3]*e2[i]/delta2;
      }
      for (int ii = 0; ii < P_l2; ++ii) {
	for (int i = 0; i < M_l2; ++i) Ctmp[i] = 0.0; 
	for (int i = 0; i < P_l2; ++i) Dtmp[i] = 0.0;
	Dtmp[ii] = D_l2[ii]; 
	for (int ij = 0; ij < M_l2*4; ++ij) Atmp[ij] = A_l2[ij];
	for (int ij = 0; ij < P_l2*4; ++ij) Btmp[ij] = B_l2[ij];
	int info;
	dgglse_(&M_l2,&N4,&P_l2,Atmp,&M_l2,Btmp,&P_l2,Ctmp,Dtmp,X,WORK,&LWORK,&info);
	assert( info == 0 );
	cvofa0_value_l2[cof_of_B_l2[ii]] = X[0];
	FOR_I3 cvofa0_grad_l2[cof_of_B_l2[ii]][i] = X[1]*e0[i]/delta0 + X[2]*e1[i]/delta1 + X[3]*e2[i]/delta2;
      }

    }

    // ==================================================================
    // for the right-side operator, we only build for internal faces...
    // ==================================================================

    if ( ifa >= nfa_bpi ) {
	
      const int cof_f = cvofa1_i[ifa];
      const int cof_l = cvofa1_i[ifa+1]-1;

      // figure out the scaling associated with each vector to keep the system 
      // reasonably conditioned...
      double delta0 = 0.0;
      double delta1 = 0.0;
      double delta2 = 0.0;
      for (int cof = cof_f; cof <= cof_l; ++cof) {
	const int icv = cvofa1_v[cof];
	double dx = x_cv[icv][0] - x_fa[ifa][0];
	double dy = x_cv[icv][1] - x_fa[ifa][1];
	double dz = x_cv[icv][2] - x_fa[ifa][2];
	double dxp = dx*e0[0] + dy*e0[1] + dz*e0[2];
	double dyp = dx*e1[0] + dy*e1[1] + dz*e1[2];
	double dzp = dx*e2[0] + dy*e2[1] + dz*e2[2];
	delta0 = max(delta0,fabs(dxp));
	delta1 = max(delta1,fabs(dyp));
	delta2 = max(delta2,fabs(dzp));
      }
      assert( delta0 > 0.0 );
      assert( delta1 > 0.0 );
      assert( delta2 > 0.0 );

      int N4 = 4;
      int N5 = 5;    
      
      int P_q = 2; // 2 constraints...
      int M_q = cof_l - cof_f + 1 - P_q; // remainder goes into M	
      int ia_q = 0; // row location in the least squares matrix A or vector C
      int ib_q = 0; // row location in the constraint matrix B or vector D
      
      int P_l1 = 1; // 1 constraints...
      int M_l1 = cof_l - cof_f + 1 - P_l1; // remainder goes into M	
      int ia_l1 = 0; // row location in the least squares matrix A or vector C
      int ib_l1 = 0; // row location in the constraint matrix B or vector D

      int P_l2 = 2; // 2 constraints...
      int M_l2 = cof_l - cof_f + 1 - P_l2; // remainder goes into M	
      int ia_l2 = 0; // row location in the least squares matrix A or vector C
      int ib_l2 = 0; // row location in the constraint matrix B or vector D
      
      const int ncof = cof_l - cof_f + 1; 
      assert( ncof >= 5 ); // we always need atleast 5 cells to close the polynomial
      assert( ncof <= NCOF_MAX );

      for (int cof = cof_f; cof <= cof_l; ++cof) {
	const int icv = cvofa1_v[cof];

	// the primed polynomial is phi = c0 + c1*xp + c2*yp + c3*zp + c4*xp*xp
	double coeff[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
	double weight = 1.0;
	  
	if (icv < ncv) {
	    
	  // local cv...
	  const int foc_f = faocv_i[icv];
	  const int foc_l = faocv_i[icv+1]-1;
	  for (int foc = foc_f; foc <= foc_l; ++foc) {
	    const int ifa_nbr = faocv_v[foc];
	    double dx_fa[3] = {
	      x_fa[ifa_nbr][0] - x_cv[icv][0],
	      x_fa[ifa_nbr][1] - x_cv[icv][1],
	      x_fa[ifa_nbr][2] - x_cv[icv][2] };
	    int nof_begin,nof_end,nof_inc;
	    if (cvofa[ifa_nbr][0] == icv) {
	      nof_begin = noofa_i[ifa_nbr];
	      nof_end = noofa_i[ifa_nbr+1];
	      nof_inc = 1;
	    }
	    else {
	      assert(cvofa[ifa_nbr][1] == icv);
	      nof_begin = noofa_i[ifa_nbr+1]-1;
	      nof_end = noofa_i[ifa_nbr]-1;
	      nof_inc = -1;
	    }
	    int ino1 = noofa_v[nof_end-nof_inc];
	    for (int nof = nof_begin; nof != nof_end; nof += nof_inc) {
	      int ino0 = ino1;
	      ino1 = noofa_v[nof];
	      double dx_no0[3] = {
		x_no[ino0][0] - x_cv[icv][0],
		x_no[ino0][1] - x_cv[icv][1],
		x_no[ino0][2] - x_cv[icv][2] };
	      double dx_no1[3] = {
		x_no[ino1][0] - x_cv[icv][0],
		x_no[ino1][1] - x_cv[icv][1],
		x_no[ino1][2] - x_cv[icv][2] };
	      double this_vol = 
		dx_fa[0]*(dx_no0[1]*dx_no1[2]-dx_no0[2]*dx_no1[1]) +
		dx_fa[1]*(dx_no0[2]*dx_no1[0]-dx_no0[0]*dx_no1[2]) +
		dx_fa[2]*(dx_no0[0]*dx_no1[1]-dx_no0[1]*dx_no1[0]);
	      assert( this_vol > 0.0 );
	      // cycle through the 4 integration points asssociated with this cell...
	      for (int ip = 0; ip < 4; ip++) {
		// the integration point vector relative to the flux face centroid is...
		double dx[3];
		FOR_I3 dx[i] = 
		  tet_int_wgt_4[ip][0]*x_no[ino0][i] + 
		  tet_int_wgt_4[ip][1]*x_no[ino1][i] + 
		  tet_int_wgt_4[ip][2]*x_fa[ifa_nbr][i] + 
		  tet_int_wgt_4[ip][3]*x_cv[icv][i] - 
		  x_fa[ifa][i];
		// convert this vector into the "primed" coordinate system...
		double xp = (dx[0]*e0[0] + dx[1]*e0[1] + dx[2]*e0[2])/delta0;
		double yp = (dx[0]*e1[0] + dx[1]*e1[1] + dx[2]*e1[2])/delta1;
		double zp = (dx[0]*e2[0] + dx[1]*e2[1] + dx[2]*e2[2])/delta2;
		// the polynomial coefficients...
		coeff[0] += this_vol;
		coeff[1] += this_vol*xp;
		coeff[2] += this_vol*yp;
		coeff[3] += this_vol*zp;
		coeff[4] += this_vol*xp*xp;
	      }
	    }
	  }
	    
	}
	else if (icv < ncv_g) {
	    
	  // a ghost cv...
	  const int ifa_nbr_f = faocv_g_i[icv-ncv];
	  const int ifa_nbr_l = faocv_g_i[icv-ncv+1]-1;
	  for (int ifa_nbr = ifa_nbr_f; ifa_nbr <= ifa_nbr_l; ++ifa_nbr) {
	    double dx_fa[3] = {
	      x_fa_g[ifa_nbr][0] - x_cv[icv][0],
	      x_fa_g[ifa_nbr][1] - x_cv[icv][1],
	      x_fa_g[ifa_nbr][2] - x_cv[icv][2] };
	    const int nof_f = noofa_g_i[ifa_nbr];
	    const int nof_l = noofa_g_i[ifa_nbr+1]-1;
	    // loop on edges ino0-ino1...
	    int ino1 = noofa_g_v[nof_l];
	    for (int nof = nof_f; nof <= nof_l; ++nof) {
	      int ino0 = ino1;
	      ino1 = noofa_g_v[nof];
	      double dx_no0[3] = {
		x_no_g[ino0][0] - x_cv[icv][0],
		x_no_g[ino0][1] - x_cv[icv][1],
		x_no_g[ino0][2] - x_cv[icv][2] };
	      double dx_no1[3] = {
		x_no_g[ino1][0] - x_cv[icv][0],
		x_no_g[ino1][1] - x_cv[icv][1],
		x_no_g[ino1][2] - x_cv[icv][2] };
	      double this_vol = 
		dx_fa[0]*(dx_no0[1]*dx_no1[2]-dx_no0[2]*dx_no1[1]) +
		dx_fa[1]*(dx_no0[2]*dx_no1[0]-dx_no0[0]*dx_no1[2]) +
		dx_fa[2]*(dx_no0[0]*dx_no1[1]-dx_no0[1]*dx_no1[0]);
	      assert( this_vol > 0.0 );
	      // cycle through the 4 integration points asssociated with this cell...
	      for (int ip = 0; ip < 4; ip++) {
		// the integration point vector relative to the flux face centroid is...
		double dx[3];
		FOR_I3 dx[i] = 
		  tet_int_wgt_4[ip][0]*x_no_g[ino0][i] + 
		  tet_int_wgt_4[ip][1]*x_no_g[ino1][i] + 
		  tet_int_wgt_4[ip][2]*x_fa_g[ifa_nbr][i] + 
		  tet_int_wgt_4[ip][3]*x_cv[icv][i] - 
		  x_fa[ifa][i];
		// convert this vector into the "primed" coordinate system...
		double xp = (dx[0]*e0[0] + dx[1]*e0[1] + dx[2]*e0[2])/delta0;
		double yp = (dx[0]*e1[0] + dx[1]*e1[1] + dx[2]*e1[2])/delta1;
		double zp = (dx[0]*e2[0] + dx[1]*e2[1] + dx[2]*e2[2])/delta2;
		// the polynomial coefficients...
		// no need to include the factor of 0.25 here - we normalize
		// the coefficients below...
		coeff[0] += this_vol;
		coeff[1] += this_vol*xp;
		coeff[2] += this_vol*yp;
		coeff[3] += this_vol*zp;
		coeff[4] += this_vol*xp*xp;
	      }
	    }
	  }
	}
	else {
	    
	  // this is a fake cell that must be reflected through its boundary
	  // face normal. boundary faces and fake cells are ordered in the same way, 
	  // so we can get the relevant face as follows...
	    
	  assert( (icv >= ncv_g)&&(icv < ncv_gf) );
	  const int ifa_b = icv-ncv_g;
	  assert( cvofa[ifa_b][1] == icv );
	    
	  // get the unit normal - it is used to relect the integration points...
	  double area = sqrt( fa_normal[ifa_b][0]*fa_normal[ifa_b][0] + 
			      fa_normal[ifa_b][1]*fa_normal[ifa_b][1] + 
			      fa_normal[ifa_b][2]*fa_normal[ifa_b][2] );
	  double unitNormal[3] = { fa_normal[ifa_b][0]/area,
				   fa_normal[ifa_b][1]/area,
				   fa_normal[ifa_b][2]/area };
	    
	  // the near-boundary cv is...
	  int icv0 = cvofa[ifa_b][0];

	  // as a minimum geometric requirement...
	  assert( unitNormal[0]*(x_fa[ifa_b][0]-x_cv[icv0][0]) +
		  unitNormal[1]*(x_fa[ifa_b][1]-x_cv[icv0][1]) +
		  unitNormal[2]*(x_fa[ifa_b][2]-x_cv[icv0][2]) > 0.0 );

	  const int foc_f = faocv_i[icv0];
	  const int foc_l = faocv_i[icv0+1]-1;
	  for (int foc = foc_f; foc <= foc_l; ++foc) {
	    const int ifa_nbr = faocv_v[foc];
	    double dx_fa[3] = {
	      x_fa[ifa_nbr][0] - x_cv[icv0][0],
	      x_fa[ifa_nbr][1] - x_cv[icv0][1],
	      x_fa[ifa_nbr][2] - x_cv[icv0][2] };
	    int nof_begin,nof_end,nof_inc;
	    if (cvofa[ifa_nbr][0] == icv0) {
	      nof_begin = noofa_i[ifa_nbr];
	      nof_end = noofa_i[ifa_nbr+1];
	      nof_inc = 1;
	    }
	    else {
	      assert(cvofa[ifa_nbr][1] == icv0);
	      nof_begin = noofa_i[ifa_nbr+1]-1;
	      nof_end = noofa_i[ifa_nbr]-1;
	      nof_inc = -1;
	    }
	    int ino1 = noofa_v[nof_end-nof_inc];
	    for (int nof = nof_begin; nof != nof_end; nof += nof_inc) {
	      int ino0 = ino1;
	      ino1 = noofa_v[nof];
	      double dx_no0[3] = {
		x_no[ino0][0] - x_cv[icv0][0],
		x_no[ino0][1] - x_cv[icv0][1],
		x_no[ino0][2] - x_cv[icv0][2] };
	      double dx_no1[3] = {
		x_no[ino1][0] - x_cv[icv0][0],
		x_no[ino1][1] - x_cv[icv0][1],
		x_no[ino1][2] - x_cv[icv0][2] };
	      double this_vol = 
		dx_fa[0]*(dx_no0[1]*dx_no1[2]-dx_no0[2]*dx_no1[1]) +
		dx_fa[1]*(dx_no0[2]*dx_no1[0]-dx_no0[0]*dx_no1[2]) +
		dx_fa[2]*(dx_no0[0]*dx_no1[1]-dx_no0[1]*dx_no1[0]);
	      assert( this_vol > 0.0 );
	      // cycle through the 3 integration points asssociated with this cell...
	      for (int ip = 0; ip < 4; ip++) {
		// the integration point vector relative to the flux face centroid is...
		double dx[3];
		FOR_I3 dx[i] = 
		  tet_int_wgt_4[ip][0]*x_no[ino0][i] + 
		  tet_int_wgt_4[ip][1]*x_no[ino1][i] + 
		  tet_int_wgt_4[ip][2]*x_fa[ifa_nbr][i] + 
		  tet_int_wgt_4[ip][3]*x_cv[icv0][i] - 
		  x_fa[ifa_b][i]; // relative to boundary face initially.
		// relect each quadrature point through the boundary face normal...
		double delta = dx[0]*unitNormal[0] + dx[1]*unitNormal[1] + dx[2]*unitNormal[2];
		// this should be negative, unless the face is really twisted...
		/*
		  if ( delta >= 0.0 ) {
		  fa_flag[ifa] = 1;
		  my_op_err += 1;
		  }
		*/
		FOR_I3 dx[i] += x_fa[ifa_b][i] - x_fa[ifa][i] - 2.0*delta*unitNormal[i];
		// convert this vector into the "primed" coordinate system...
		double xp = (dx[0]*e0[0] + dx[1]*e0[1] + dx[2]*e0[2])/delta0;
		double yp = (dx[0]*e1[0] + dx[1]*e1[1] + dx[2]*e1[2])/delta1;
		double zp = (dx[0]*e2[0] + dx[1]*e2[1] + dx[2]*e2[2])/delta2;
		// the polynomial coefficients...
		coeff[0] += this_vol;
		coeff[1] += this_vol*xp;
		coeff[2] += this_vol*yp;
		coeff[3] += this_vol*zp;
		coeff[4] += this_vol*xp*xp;
	      }
	    }
	  }

	}
	  
	// normalize out the volume and multiply by the weight...

	coeff[1] *= weight/coeff[0];
	coeff[2] *= weight/coeff[0];
	coeff[3] *= weight/coeff[0];
	coeff[4] *= weight/coeff[0];
	coeff[0] = weight;

	// _q and _l2 have the same cv's constrained...
	if ( (icv == cvofa[ifa][0])||(icv == cvofa[ifa][1]) ) {
	  // this is one of the constrained cv's next to the face, so populate B...
	  // _q...
	  for (int j = 0; j < 5; j++) B_q[ib_q+j*P_q] = coeff[j];
	  D_q[ib_q] = weight;
	  cof_of_B_q[ib_q++] = cof;
	  // _l2...
	  for (int j = 0; j < 4; j++) B_l2[ib_l2+j*P_l2] = coeff[j];
	  D_l2[ib_l2] = weight;
	  cof_of_B_l2[ib_l2++] = cof;
	}
	else {
	  // this is one of the least squares cv's, so populate A...
	  // _q
	  for (int j = 0; j < 5; j++) A_q[ia_q+j*M_q] = coeff[j];
	  C_q[ia_q] = weight;
	  cof_of_A_q[ia_q++] = cof;
	  // _l2...
	  for (int j = 0; j < 4; j++) A_l2[ia_l2+j*M_l2] = coeff[j];
	  C_l2[ia_l2] = weight;
	  cof_of_A_l2[ia_l2++] = cof;
	}
	
	// _l1 has only icv1 constrained for this right operator...
	if ( icv == cvofa[ifa][1] ) {
	  // this is one of the constrained cv's next to the face, so populate B...
	  // _l1...
	  for (int j = 0; j < 4; j++) B_l1[ib_l1+j*P_l1] = coeff[j];
	  D_l1[ib_l1] = weight;
	  cof_of_B_l1[ib_l1++] = cof;
	}
	else {
	  // this is one of the least squares cv's, so populate A...
	  // _l1...
	  for (int j = 0; j < 4; j++) A_l1[ia_l1+j*M_l1] = coeff[j];
	  C_l1[ia_l1] = weight;
	  cof_of_A_l1[ia_l1++] = cof;
	}

      }

      assert( ia_q == M_q );
      assert( ia_l1 == M_l1 );
      assert( ia_l2 == M_l2 );

      assert( ib_q == P_q );
      assert( ib_l1 == P_l1 );
      assert( ib_l2 == P_l2 );
	
      // _q quadratic reconstruction in normal direction, linear in other...

      for (int ii = 0; ii < M_q; ++ii) {
	for (int i = 0; i < M_q; ++i) Ctmp[i] = 0.0; 
	for (int i = 0; i < P_q; ++i) Dtmp[i] = 0.0;
	Ctmp[ii] = C_q[ii]; 
	for (int ij = 0; ij < M_q*5; ++ij) Atmp[ij] = A_q[ij];
	for (int ij = 0; ij < P_q*5; ++ij) Btmp[ij] = B_q[ij];
	int info;
	dgglse_(&M_q,&N5,&P_q,Atmp,&M_q,Btmp,&P_q,Ctmp,Dtmp,X,WORK,&LWORK,&info);
	assert( info == 0 );
	cvofa1_value_q[cof_of_A_q[ii]] = X[0];
      }
      for (int ii = 0; ii < P_q; ++ii) {
	for (int i = 0; i < M_q; ++i) Ctmp[i] = 0.0; 
	for (int i = 0; i < P_q; ++i) Dtmp[i] = 0.0;
	Dtmp[ii] = D_q[ii]; 
	for (int ij = 0; ij < M_q*5; ++ij) Atmp[ij] = A_q[ij];
	for (int ij = 0; ij < P_q*5; ++ij) Btmp[ij] = B_q[ij];
	int info;
	dgglse_(&M_q,&N5,&P_q,Atmp,&M_q,Btmp,&P_q,Ctmp,Dtmp,X,WORK,&LWORK,&info);
	assert( info == 0 );
	cvofa1_value_q[cof_of_B_q[ii]] = X[0];
      }
      
      // _l1 single constraint linear reconstruction in all directions...

      for (int ii = 0; ii < M_l1; ++ii) {
	for (int i = 0; i < M_l1; ++i) Ctmp[i] = 0.0; 
	for (int i = 0; i < P_l1; ++i) Dtmp[i] = 0.0;
	Ctmp[ii] = C_l1[ii]; 
	for (int ij = 0; ij < M_l1*4; ++ij) Atmp[ij] = A_l1[ij];
	for (int ij = 0; ij < P_l1*4; ++ij) Btmp[ij] = B_l1[ij];
	int info;
	dgglse_(&M_l1,&N4,&P_l1,Atmp,&M_l1,Btmp,&P_l1,Ctmp,Dtmp,X,WORK,&LWORK,&info);
	assert( info == 0 );
	cvofa1_value_l1[cof_of_A_l1[ii]] = X[0];
      }
      for (int ii = 0; ii < P_l1; ++ii) {
	for (int i = 0; i < M_l1; ++i) Ctmp[i] = 0.0; 
	for (int i = 0; i < P_l1; ++i) Dtmp[i] = 0.0;
	Dtmp[ii] = D_l1[ii]; 
	for (int ij = 0; ij < M_l1*4; ++ij) Atmp[ij] = A_l1[ij];
	for (int ij = 0; ij < P_l1*4; ++ij) Btmp[ij] = B_l1[ij];
	int info;
	dgglse_(&M_l1,&N4,&P_l1,Atmp,&M_l1,Btmp,&P_l1,Ctmp,Dtmp,X,WORK,&LWORK,&info);
	assert( info == 0 );
	cvofa1_value_l1[cof_of_B_l1[ii]] = X[0];
      }
      
      // _l2 double constraint linear reconstruction in all directions...
      // this introduces minimal bias, and is used for the gradient... 
      
      for (int ii = 0; ii < M_l2; ++ii) {
	for (int i = 0; i < M_l2; ++i) Ctmp[i] = 0.0; 
	for (int i = 0; i < P_l2; ++i) Dtmp[i] = 0.0;
	Ctmp[ii] = C_l2[ii]; 
	for (int ij = 0; ij < M_l2*4; ++ij) Atmp[ij] = A_l2[ij];
	for (int ij = 0; ij < P_l2*4; ++ij) Btmp[ij] = B_l2[ij];
	int info;
	dgglse_(&M_l2,&N4,&P_l2,Atmp,&M_l2,Btmp,&P_l2,Ctmp,Dtmp,X,WORK,&LWORK,&info);
	assert( info == 0 );
	cvofa1_value_l2[cof_of_A_l2[ii]] = X[0];
	FOR_I3 cvofa1_grad_l2[cof_of_A_l2[ii]][i] = X[1]*e0[i]/delta0 + X[2]*e1[i]/delta1 + X[3]*e2[i]/delta2;
      }
      for (int ii = 0; ii < P_l2; ++ii) {
	for (int i = 0; i < M_l2; ++i) Ctmp[i] = 0.0; 
	for (int i = 0; i < P_l2; ++i) Dtmp[i] = 0.0;
	Dtmp[ii] = D_l2[ii]; 
	for (int ij = 0; ij < M_l2*4; ++ij) Atmp[ij] = A_l2[ij];
	for (int ij = 0; ij < P_l2*4; ++ij) Btmp[ij] = B_l2[ij];
	int info;
	dgglse_(&M_l2,&N4,&P_l2,Atmp,&M_l2,Btmp,&P_l2,Ctmp,Dtmp,X,WORK,&LWORK,&info);
	assert( info == 0 );
	cvofa1_value_l2[cof_of_B_l2[ii]] = X[0];
	FOR_I3 cvofa1_grad_l2[cof_of_B_l2[ii]][i] = X[1]*e0[i]/delta0 + X[2]*e1[i]/delta1 + X[3]*e2[i]/delta2;
      }

    }
      
  }

  // cleanup...
  delete[] faocv_g_i;
  delete[] x_fa_g;
  delete[] noofa_g_i;
  delete[] noocv_g_i;
  delete[] x_no_g;
  delete[] noofa_g_v;

  // error checking...
  /*
    int op_err;
    MPI_Allreduce(&my_op_err,&op_err,1,MPI_INT,MPI_SUM,mpi_comm);
    if (op_err > 0) {
    if (mpi_rank == 0)
    cerr << "Error: ops failed at " << op_err << " faces. Writing bad faces for diagnostics." << endl;
    FOR_ICV cv_flag[icv] = 0;
    updateFaI1(fa_flag,MAX_DATA);
    for (int ifa = 0; ifa < nfa_bpi; ++ifa)
    if (fa_flag[ifa] != 0)
    cv_flag[cvofa[ifa][0]] = 1;
    for (int ifa = nfa_bpi; ifa < nfa; ++ifa)
    if (fa_flag[ifa] != 0)
    cv_flag[cvofa[ifa][0]] = cv_flag[cvofa[ifa][1]] = 1;
    writeFlaggedFacesTecplot("failed_faces.dat");
    writeFlaggedCvsTecplot("failed_cvs.dat");
    MPI_Barrier(mpi_comm);
    throw(-1);
    }
  */

  // reconstruction operator stats...
  {
    double my_diag_min[3] = { 1.0E+20, 1.0E+20, 1.0E+20 };    
    double my_coeff_max[3] = { 0.0, 0.0, 0.0 };
    FOR_IFA {
      {
	const int cof_f = cvofa0_i[ifa];
	const int cof_l = cvofa0_i[ifa+1]-1;
	// diag is the first...
	my_diag_min[0] = min( my_diag_min[0], cvofa0_value_q[cof_f] );
	my_diag_min[1] = min( my_diag_min[1], cvofa0_value_l1[cof_f] );
	my_diag_min[2] = min( my_diag_min[2], cvofa0_value_l2[cof_f] );
	for (int cof = cof_f+1; cof <= cof_l; ++cof) { // skip diag (first)
	  my_coeff_max[0] = max( my_coeff_max[0], fabs(cvofa0_value_q[cof]) );
	  my_coeff_max[1] = max( my_coeff_max[1], fabs(cvofa0_value_l1[cof]) );
	  my_coeff_max[2] = max( my_coeff_max[2], fabs(cvofa0_value_l2[cof]) );
	}
      }
      {
	const int cof_f = cvofa1_i[ifa];
	const int cof_l = cvofa1_i[ifa+1]-1;
	my_diag_min[0] = min( my_diag_min[0], cvofa1_value_q[cof_f] );
	my_diag_min[1] = min( my_diag_min[1], cvofa1_value_l1[cof_f] );
	my_diag_min[2] = min( my_diag_min[2], cvofa1_value_l2[cof_f] );
	for (int cof = cof_f+1; cof <= cof_l; ++cof) { // skip diag (first)
	  my_coeff_max[0] = max( my_coeff_max[0], fabs(cvofa1_value_q[cof]) );
	  my_coeff_max[1] = max( my_coeff_max[1], fabs(cvofa1_value_l1[cof]) );
	  my_coeff_max[2] = max( my_coeff_max[2], fabs(cvofa1_value_l2[cof]) );
	}
      }
    }
    double diag_min[3];
    MPI_Reduce(my_diag_min,diag_min,3,MPI_DOUBLE,MPI_MIN,0,mpi_comm);
    double coeff_max[3];
    MPI_Reduce(my_coeff_max,coeff_max,3,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
    if (mpi_rank == 0) {
      cout << " > min diag max coeff (_q operator):  " << diag_min[0] << " " << coeff_max[0] << endl;
      cout << " > min diag max coeff (_l1 operator): " << diag_min[1] << " " << coeff_max[1] << endl;
      cout << " > min diag max coeff (_l2 operator): " << diag_min[2] << " " << coeff_max[2] << endl;
    }
  }

  /*
  // take a look...
  int ifa_debug = nfa_bpi;
  if (ifa_debug >= 0) {
    cout << "ifa_debug: " << ifa_debug << endl;
    {
      cout << "left operator: " << endl;
      const int cof_f = cvofa0_i[ifa_debug];
      const int cof_l = cvofa0_i[ifa_debug+1]-1;
      cout << " > _q operator:   ";
      for (int cof = cof_f; cof <= cof_l; ++cof)
	cout << cvofa0_value_q[cof] << " ";
      cout << endl;
      cout << " > _l1 operator:  ";
      for (int cof = cof_f; cof <= cof_l; ++cof)
	cout << cvofa0_value_l1[cof] << " ";
      cout << endl;
      cout << " > _l2 operator:  ";
      for (int cof = cof_f; cof <= cof_l; ++cof)
	cout << cvofa0_value_l2[cof] << " ";
      cout << endl;
    }
    {
      cout << "right operator: " << endl;
      const int cof_f = cvofa1_i[ifa_debug];
      const int cof_l = cvofa1_i[ifa_debug+1]-1;
      cout << " > _q operator:   ";
      for (int cof = cof_f; cof <= cof_l; ++cof)
	cout << cvofa1_value_q[cof] << " ";
      cout << endl;
      cout << " > _l1 operator:  ";
      for (int cof = cof_f; cof <= cof_l; ++cof)
	cout << cvofa1_value_l1[cof] << " ";
      cout << endl;
      cout << " > _l2 operator:  ";
      for (int cof = cof_f; cof <= cof_l; ++cof)
	cout << cvofa1_value_l2[cof] << " ";
      cout << endl;
    }
    cout << endl;
  }
  MPI_Pause("ifa_debug");
  */
  
  // ===============================================
  // set the final operators...
  // ===============================================

  // this is an attempt to ensure the one-sided operators have a 
  // weighting that is dominated by the upwind cell...

  /*
  combineOperatorsGR(cvofa0_value,cvofa0_grad,cvofa1_value,cvofa1_grad,
		     cvofa0_value_q,cvofa1_value_q,
		     cvofa0_value_l1,cvofa1_value_l1,
		     cvofa0_grad_l2,cvofa1_grad_l2);
  */

  /*
  combineOperatorsLinear(cvofa0_value,cvofa0_grad,cvofa1_value,cvofa1_grad,
			 cvofa0_value_q,cvofa1_value_q,
			 cvofa0_value_l1,cvofa1_value_l1,
			 cvofa0_grad_l2,cvofa1_grad_l2);
  */
  
  // my experiments...
  combineOperatorsHam(cvofa0_value,cvofa0_grad,cvofa1_value,cvofa1_grad,
		      cvofa0_value_q,cvofa1_value_q,
		      cvofa0_value_l1,cvofa1_value_l1,
		      cvofa0_grad_l2,cvofa1_grad_l2);

  // cleanup...
  
  delete[] cvofa0_value_l1;
  delete[] cvofa1_value_l1;
  
  delete[] cvofa0_value_l2;
  delete[] cvofa1_value_l2;
  
}

void UgpWithCvFakeOp::combineOperatorsGR(double * cvofa0_value,double (*cvofa0_grad)[3],double * cvofa1_value,double (*cvofa1_grad)[3],
					 double * cvofa0_value_q,double * cvofa1_value_q,
					 double * cvofa0_value_l1,double * cvofa1_value_l1,
					 double (*cvofa0_grad_l2)[3],double (*cvofa1_grad_l2)[3]) {

  if (mpi_rank == 0)
    cout << " > combineOperatorsGR" << endl;
    
    // we need to take some combination of the operators here...
    // 1. the result has to be sufficiently biased to allow the dissipation mechanism to work...
    // 2. accuracy...
    
    const double gr_max = 0.75;
    int my_op_count[3] = { 0, 0, 0 };
    
    // --------------------------
    // left, including boundary faces...
    // --------------------------
    
  for (int ifa = 0; ifa < nfa; ++ifa) {

    const int cof_f = cvofa0_i[ifa];
    const int cof_l = cvofa0_i[ifa+1]-1;
   
    // compute sum_{j,i!=j}(|a_ij|)/a_ii, the "Greshgorin ratio" for the 
    // quadratic scheme...

    double gr_q = 0.0;
    for (int cof = cof_f+1; cof <= cof_l; ++cof)
      gr_q += fabs(cvofa0_value_q[cof]);
    // normalize by the diagonal...
    assert( cvofa0_value_q[cof_f] > 0.0 );
    gr_q /= cvofa0_value_q[cof_f];
    
    // if this is less than the target, just use the quadratic operator. On
    // uniform hex grids, this number is 0.6 ((1/6+1/3)/(5/6))...

    if ( gr_q <= gr_max ) {
      // fine - use quadratic operator...
      my_op_count[0] += 1;
      for (int cof = cof_f; cof <= cof_l; ++cof)
	cvofa0_value[cof] = cvofa0_value_q[cof];
    }
    else {
      
      // compute the "Greshgorin ratio" for the linear scheme...
      double gr_l1 = 0.0;
      for (int cof = cof_f+1; cof <= cof_l; ++cof)
	gr_l1 += fabs(cvofa0_value_l1[cof]);
      assert( cvofa0_value_l1[cof_f] > 0.0 );
      gr_l1 /= cvofa0_value_l1[cof_f];
      
      if ( gr_l1 <= gr_max ) {
	// the linear scheme has "Greshgorin ratio" below the target,
	// so we can use a linear combination of quadratic and linear
	// to achieve the target...
	my_op_count[1] += 1;
	for (int cof = cof_f; cof <= cof_l; ++cof)
	  cvofa0_value[cof] = ((gr_max-gr_l1)*cvofa0_value_q[cof] + (gr_q-gr_max)*cvofa0_value_q[cof])/(gr_q-gr_l1);
      }
      else {
	// both linear and quadratic are above the max - linearly combine linear
	// with 1st-order upwind (gr=0) to produce final operator...
	my_op_count[2] += 1;
	/*
	for (int cof = cof_f; cof <= cof_l; ++cof)
	  cvofa0_value[cof] = (gr_max-0.0)*cvofa0_value_l1[cof]/(gr_l1-0.0);
	// upwind has 1 on the diag...
	cvofa0_value[cof_f] += (gr_l1-gr_max)*1.0/(gr_l1-0.0);
	*/
	// or, just linear...
 	for (int cof = cof_f; cof <= cof_l; ++cof)
	  cvofa0_value[cof] = cvofa0_value_l1[cof];

      }

    }
    
    // for grad, just use the non-biased _l2... 
    for (int cof = cof_f; cof <= cof_l; ++cof)
      FOR_I3 cvofa0_grad[cof][i] = cvofa0_grad_l2[cof][i];
    
  }

  // --------------------------
  // right...
  // --------------------------
  
  for (int ifa = nfa_bpi; ifa < nfa; ++ifa) {

    const int cof_f = cvofa1_i[ifa];
    const int cof_l = cvofa1_i[ifa+1]-1;

    // compute sum_{j,i!=j}(|a_ij|)/a_ii, the "Greshgorin ratio" for the 
    // quadratic scheme...
    
    double gr_q = 0.0;
    for (int cof = cof_f+1; cof <= cof_l; ++cof)
      gr_q += fabs(cvofa1_value_q[cof]);
    // normalize by the diagonal...
    assert( cvofa1_value_q[cof_f] > 0.0 );
    gr_q /= cvofa1_value_q[cof_f];
    
    // if this is less than the target, just use the quadratic operator. On
    // uniform hex grids, this number is 0.6 ((1/6+1/3)/(5/6))...
    
    if ( gr_q <= gr_max ) {
      // fine - use quadratic operator...
      my_op_count[0] += 1;
      for (int cof = cof_f; cof <= cof_l; ++cof)
	cvofa1_value[cof] = cvofa1_value_q[cof];
    }
    else {
      
      // compute the "Greshgorin ratio" for the linear scheme...
      double gr_l1 = 0.0;
      for (int cof = cof_f+1; cof <= cof_l; ++cof)
	gr_l1 += fabs(cvofa1_value_l1[cof]);
      assert( cvofa1_value_l1[cof_f] > 0.0 );
      gr_l1 /= cvofa1_value_l1[cof_f];
      
      if ( gr_l1 <= gr_max ) {
	// the linear scheme has "Greshgorin ratio" below the target,
	// so we can use a linear combination of quadratic and linear
	// to achieve the target...
	my_op_count[1] += 1;
	for (int cof = cof_f; cof <= cof_l; ++cof)
	  cvofa1_value[cof] = ((gr_max-gr_l1)*cvofa1_value_q[cof] + (gr_q-gr_max)*cvofa1_value_q[cof])/(gr_q-gr_l1);
      }
      else {
	// both linear and quadratic are above the max - linearly combine linear
	// with 1st-order upwind (gr=0) to produce final operator...
	my_op_count[2] += 1;
	/*
	for (int cof = cof_f; cof <= cof_l; ++cof)
	  cvofa1_value[cof] = (gr_max-0.0)*cvofa1_value_l1[cof]/(gr_l1-0.0);
	// upwind has 1 on the diag...
	cvofa1_value[cof_f] += (gr_l1-gr_max)*1.0/(gr_l1-0.0);
	*/
	// or just linear...
	for (int cof = cof_f; cof <= cof_l; ++cof)
	  cvofa1_value[cof] = cvofa1_value_l1[cof];
	
      }
      
    }
    
    // for grad, just use the non-biased _l2... 
    for (int cof = cof_f; cof <= cof_l; ++cof)
      FOR_I3 cvofa1_grad[cof][i] = cvofa1_grad_l2[cof][i];
    
  }

  int op_count[3];
  MPI_Reduce(my_op_count,op_count,3,MPI_INT,MPI_SUM,0,mpi_comm);
  if (mpi_rank == 0) {
    //int op_count_sum = op_count[0] + op_count[1] + op_count[2];
    cout << " > op_count: fully _q:    " << op_count[0] << endl;
    cout << " > op_count: q _l1 blend: " << op_count[1] << endl;
    cout << " > op_count: fully _l1:   " << op_count[2] << endl;
  }

}

void UgpWithCvFakeOp::combineOperatorsLinear(double * cvofa0_value,double (*cvofa0_grad)[3],double * cvofa1_value,double (*cvofa1_grad)[3],
					     double * cvofa0_value_q,double * cvofa1_value_q,
					     double * cvofa0_value_l1,double * cvofa1_value_l1,
					     double (*cvofa0_grad_l2)[3],double (*cvofa1_grad_l2)[3]) {

  if (mpi_rank == 0)
    cout << " > combineOperatorsLinear" << endl;
    
  // --------------------------
  // left, including boundary faces...
  // --------------------------
  
  for (int ifa = 0; ifa < nfa; ++ifa) {    
    const int cof_f = cvofa0_i[ifa];
    const int cof_l = cvofa0_i[ifa+1]-1;
  
    for (int cof = cof_f; cof <= cof_l; ++cof)
      cvofa0_value[cof] = cvofa0_value_l1[cof];
   
    // for grad, just use the non-biased _l2... 
    for (int cof = cof_f; cof <= cof_l; ++cof)
      FOR_I3 cvofa0_grad[cof][i] = cvofa0_grad_l2[cof][i];
    
  }

  // --------------------------
  // right...
  // --------------------------
  
  for (int ifa = nfa_bpi; ifa < nfa; ++ifa) {

    const int cof_f = cvofa1_i[ifa];
    const int cof_l = cvofa1_i[ifa+1]-1;

    for (int cof = cof_f; cof <= cof_l; ++cof)
      cvofa1_value[cof] = cvofa1_value_l1[cof];
    
    // for grad, just use the non-biased _l2... 
    for (int cof = cof_f; cof <= cof_l; ++cof)
      FOR_I3 cvofa1_grad[cof][i] = cvofa1_grad_l2[cof][i];

  }

}

void UgpWithCvFakeOp::combineOperatorsHam(double * cvofa0_value,double (*cvofa0_grad)[3],double * cvofa1_value,double (*cvofa1_grad)[3],
					  double * cvofa0_value_q,double * cvofa1_value_q,
					  double * cvofa0_value_l1,double * cvofa1_value_l1,
					  double (*cvofa0_grad_l2)[3],double (*cvofa1_grad_l2)[3]) {

  if (mpi_rank == 0)
    cout << " > combineOperatorsHam" << endl;

  int my_op_count[2] = { 0, 0 }; // store the count of quadratic and linear operators
  
  // --------------------------
  // left, including boundary faces...
  // --------------------------
  
  for (int ifa = 0; ifa < nfa; ++ifa) {    
    const int cof_f = cvofa0_i[ifa];
    const int cof_l = cvofa0_i[ifa+1]-1;
    
    // select the operator with the smaller coeff rms (including the diag).
    // on uniform hex meshes, this will preference the quadratic operator, 
    // but on bad grids, it may preference the linear...

    double rms_q = 0.0;
    double rms_l1 = 0.0;
    for (int cof = cof_f; cof <= cof_l; ++cof) {
      rms_q += cvofa0_value_q[cof]*cvofa0_value_q[cof];
      rms_l1 += cvofa0_value_l1[cof]*cvofa0_value_l1[cof];
    }
    
    if ((rms_q <= rms_l1)&&(cvofa0_value_q[cof_f] > cvofa0_value_l1[cof_f]*0.75)) {
      // take quadratic...
      my_op_count[0] += 1;
      for (int cof = cof_f; cof <= cof_l; ++cof)
	cvofa0_value[cof] = cvofa0_value_q[cof];
    }
    else {
      // take linear...
      my_op_count[1] += 1;
      for (int cof = cof_f; cof <= cof_l; ++cof)
	cvofa0_value[cof] = cvofa0_value_l1[cof];
    }
      
    // for grad, just use the non-biased _l2... 
    for (int cof = cof_f; cof <= cof_l; ++cof)
      FOR_I3 cvofa0_grad[cof][i] = cvofa0_grad_l2[cof][i];
    
  }

  // --------------------------
  // right...
  // --------------------------
  
  for (int ifa = nfa_bpi; ifa < nfa; ++ifa) {

    const int cof_f = cvofa1_i[ifa];
    const int cof_l = cvofa1_i[ifa+1]-1;

    double rms_q = 0.0;
    double rms_l1 = 0.0;
    for (int cof = cof_f; cof <= cof_l; ++cof) {
      rms_q += cvofa1_value_q[cof]*cvofa1_value_q[cof];
      rms_l1 += cvofa1_value_l1[cof]*cvofa1_value_l1[cof];
    }

    if ((rms_q <= rms_l1)&&(cvofa1_value_q[cof_f] > cvofa1_value_l1[cof_f]*0.75)) {
      // take quadratic...
      my_op_count[0] += 1;
      for (int cof = cof_f; cof <= cof_l; ++cof)
	cvofa1_value[cof] = cvofa1_value_q[cof];
    }
    else {
      // take linear...
      my_op_count[1] += 1;
      for (int cof = cof_f; cof <= cof_l; ++cof)
	cvofa1_value[cof] = cvofa1_value_l1[cof];
    }

    // for grad, just use the non-biased _l2... 
    for (int cof = cof_f; cof <= cof_l; ++cof)
      FOR_I3 cvofa1_grad[cof][i] = cvofa1_grad_l2[cof][i];
    
  }

  // report stats...
  int op_count[2];
  MPI_Reduce(my_op_count,op_count,2,MPI_INT,MPI_SUM,0,mpi_comm);
  if (mpi_rank == 0) {
    //int op_count_sum = op_count[0] + op_count[1] + op_count[2];
    cout << " > op_count: _q:  " << op_count[0] << endl;
    cout << " > op_count: _l1: " << op_count[1] << endl;
  }

}

void UgpWithCvFakeOp::checkOperators() {
    
  if (mpi_rank == 0)
    cout << "checkOperators()" << endl;

  // do a test...
  // populate phi field including ghost and fake with a linear function...

  double * phi = new double[ncv_gf];
  for (int icv = 0; icv < ncv_gf; ++icv)
    phi[icv] = 1.1 + 1.21*x_cv[icv][0] + 2.321*x_cv[icv][1] + 3.4321*x_cv[icv][2];
    
  // now check reconstructions...
    
  // left side operators...
    
  double my_buf[2] = { 0.0, 0.0 };

  /*
  std::stringstream ss;
  ss << "example." << mpi_rank;

  ofstream myfile;
  myfile.open (ss.str().c_str());
  */

  for (int ifa = 0; ifa < nfa; ++ifa) {
    double phi_f = 0.0;
    double grad_phi_f[3] = { 0.0, 0.0, 0.0 };
    const int cof_f = cvofa0_i[ifa];
    const int cof_l = cvofa0_i[ifa+1]-1;
    for (int cof = cof_f; cof <= cof_l; ++cof) {
      const int icv = cvofa0_v[cof];
      phi_f += cvofa0_value[cof]*phi[icv];
      FOR_I3 grad_phi_f[i] += cvofa0_grad[cof][i]*phi[icv];
    }
    double phi_f_exact = 1.1 + 1.21*x_fa[ifa][0] + 2.321*x_fa[ifa][1] + 3.4321*x_fa[ifa][2];
    my_buf[0] = max( fabs(phi_f - phi_f_exact), my_buf[0] );

    /*
    // take a look...
    if (fabs(phi_f - phi_f_exact) > 1.0E-8) {
      myfile << "got large error in left operator of face: " << ifa << endl;
      for (int cof = cof_f; cof <= cof_l; ++cof) {
	const int icv = cvofa0_v[cof];
	myfile << " > x,y,z,coeff: " << x_cv[icv][0] << " " << x_cv[icv][1] << " " << x_cv[icv][2] << " " << cvofa0_value[cof] << endl;
      }
    }
    */

    my_buf[1] = max( fabs(grad_phi_f[0] - 1.21), my_buf[1] );
    my_buf[1] = max( fabs(grad_phi_f[1] - 2.321), my_buf[1] );
    my_buf[1] = max( fabs(grad_phi_f[2] - 3.4321), my_buf[1] );
  }
    
  // right side operators...
    
  for (int ifa = nfa_bpi; ifa < nfa; ++ifa) {
    double phi_f = 0.0;      
    double grad_phi_f[3] = { 0.0, 0.0, 0.0 };
    const int cof_f = cvofa1_i[ifa];
    const int cof_l = cvofa1_i[ifa+1]-1;
    for (int cof = cof_f; cof <= cof_l; ++cof) {
      const int icv = cvofa1_v[cof];
      phi_f += cvofa1_value[cof]*phi[icv];
      FOR_I3 grad_phi_f[i] += cvofa1_grad[cof][i]*phi[icv];
    }
    double phi_f_exact = 1.1 + 1.21*x_fa[ifa][0] + 2.321*x_fa[ifa][1] + 3.4321*x_fa[ifa][2];
    my_buf[0] = max( fabs(phi_f - phi_f_exact), my_buf[0] );

    /*
    // take a look...
    if (fabs(phi_f - phi_f_exact) > 1.0E-8) {
      myfile << "got large error in right operator of face: " << ifa << endl;
      for (int cof = cof_f; cof <= cof_l; ++cof) {
	const int icv = cvofa1_v[cof];
	myfile << " > x,y,z,coeff: " << x_cv[icv][0] << " " << x_cv[icv][1] << " " << x_cv[icv][2] << " " << cvofa1_value[cof] << endl;
      }
    }
    */

    my_buf[1] = max( fabs(grad_phi_f[0] - 1.21), my_buf[1] );
    my_buf[1] = max( fabs(grad_phi_f[1] - 2.321), my_buf[1] );
    my_buf[1] = max( fabs(grad_phi_f[2] - 3.4321), my_buf[1] );
  }
  
  //myfile.close();
    
  double buf[2];
  MPI_Reduce(my_buf,buf,2,MPI_DOUBLE,MPI_MAX,0,mpi_comm);
  if (mpi_rank == 0) {
    cout << " > max value error (should be small): " << buf[0] << endl;
    cout << " > max gradient error (should be small): " << buf[1] << endl;
  }
    
  delete[] phi;
    
}
  
// these are some structs and a sorting operator req'd by the skewNorm calc...
  
typedef struct {
  int i_global,j_global;
  double coeff[3];
} MatrixEntry;
  
class MatrixEntryCompare : public std::binary_function<MatrixEntry&, MatrixEntry&, bool> {
public: 
  bool operator()(const MatrixEntry& a,const MatrixEntry& b) {
    // sort based on i_global...
    return( a.i_global < b.i_global );
  }
};
  

void UgpWithCvFakeOp::calcSkewNormStuff(const double fa_alpha_coeff) {
  
  // the differencing matrix should be skew-symmetric. Compute the 
  // row norm of the Dx + Dx^T as an indication of how non-skew-symmetric 
  // the differencing matrix actually is...
  
  if (mpi_rank == 0)
    cout << "calcSkewNormStuff(), fa_alpha_coeff: " << fa_alpha_coeff << endl;
  
  // computing the row-norm of "D + D^T" where
  // D is the differencing matrix integrated over volume, i.e.:
  // D_i phi = sum( 0.5*( phi_l + phi_r)*n_i*A_face )

  // ==============================================================
  // start by storing the number of nb's of cv in each cv_flag...
  // ==============================================================

  // for fakes mainly, but not really neccessary, because the nbocv_i/v does
  // not reference the faces - they should already be eliminated...
  FOR_ICV_GF cv_flag[icv] = -1; 
  FOR_ICV cv_flag[icv] = nbocv_i[icv+1] - nbocv_i[icv];
  updateCvData(cv_flag,REPLACE_DATA);
    
  // now every cv need to build a nbr list that includes ALL its immediate
  // nbrs AND their nbrs. Since we do not know exactly how big this will be, 
  // include space for all nbrs and their nbrs - some of which will be common...
    
  int * nb2ocv_i = new int [ncv+1];
  nb2ocv_i[0] = 0;
  FOR_ICV {
    // diagonal...
    nb2ocv_i[icv+1] = nb2ocv_i[icv] + 1;
    // cycle through our immediate nbrs and add their nbr counts minus 1 (we
    // are already in the list)...
    const int noc_f = nbocv_i[icv];
    const int noc_l = nbocv_i[icv+1]-1;
    // skip diagonal...
    for (int noc = noc_f+1; noc <= noc_l; ++noc) {
      const int icv_nbr = nbocv_v[noc];
      assert( (icv_nbr >= 0)&&(icv_nbr < ncv_g) );
      assert( cv_flag[icv_nbr] >= 2 );
      nb2ocv_i[icv+1] += cv_flag[icv_nbr]-1;
    }
  }
  int nb2ocv_s = nb2ocv_i[ncv];
  int * nb2ocv_v = new int[nb2ocv_s];
  for (int n2oc = 0; n2oc < nb2ocv_s; ++n2oc)
    nb2ocv_v[n2oc] = -1;

  {
    double my_buf[2];
    my_buf[0] = (double)nb2ocv_s; // use double here to avoid overflow
    my_buf[1] = (double)ncv;
    double buf[2];
    MPI_Reduce(my_buf,buf,2,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0)
      cout << " > average nb2ocv count per cv: " << buf[0]/buf[1] << endl;
  }

  // put a global numbering in cv_flag. Here we use cvora, which should 
  // already be initialized...
  assert( cvora != NULL );
  assert( cvora[0] == 0 );
  assert( cvora[mpi_rank+1]-cvora[mpi_rank] == ncv ); 
  FOR_ICV_GF cv_flag[icv] = -1; 
  FOR_ICV cv_flag[icv] = icv + cvora[mpi_rank];
  updateCvData(cv_flag,REPLACE_DATA);
    
  // cycle through the faces once and count the number of messages to go to
  // nbrs...
  // FOR NOW, ignore non-cartesian transforms, but in the future, we should 
  // include them to properly rotate the operators at these interfaces...

  // build the full DpDt "D plus D transpose" differencing operator...
  double (*DpDt)[3] = new double[nb2ocv_s][3];
  for (int n2oc = 0; n2oc < nb2ocv_s; ++n2oc)
    FOR_I3 DpDt[n2oc][i] = 0.0;

  // we also need a vector of non-local matrix entries...
  vector<MatrixEntry> matrixEntry;

  // do local stuff in two passes - on the first, simply count the non-local matrix entries,
  // and on the second pack them, and fill the local matrix entries...
    
  for (int iter = 0; iter < 2; ++iter) {
      
    int meCount = 0; // matrix entry count
      
    FOR_IFA {
	
      // this face based operator, combined with the normal, is used to build the
      // differencing operator associated with two cells, icv0 and icv1...
	
      int i0_global = cv_flag[cvofa[ifa][0]];
	
      // i0_global should ALWAYS be local...
      assert( (i0_global >= cvora[mpi_rank])&&(i0_global < cvora[mpi_rank+1]) );
	
      int i1_global = cv_flag[cvofa[ifa][1]];
      if ( i1_global < 0) {
	// this must be a "fake" cv. Assume Neumann closure in building this
	// operator for now...
	assert( cvofa[ifa][1] >= ncv_g );
	i1_global = i0_global;
      }
	
      // i1 global may not be local...
      int i1_nonlocal = 0;
      if ( (i1_global < cvora[mpi_rank])||(i1_global >= cvora[mpi_rank+1]) )
	i1_nonlocal = 1;
	
      {
	  
	// left-side reconstruction...
	  
	const int cof_f = cvofa0_i[ifa];
	const int cof_l = cvofa0_i[ifa+1]-1;
	for (int cof = cof_f; cof <= cof_l; ++cof) {
	  const int icv = cvofa0_v[cof];
	  int j_global = cv_flag[icv];
	  if ( j_global < 0) {
	    // this must be a "fake" cv. Assume Neumann closure in building this
	    // operator for now...
	    assert( icv >= ncv_g );
	    j_global = i0_global;
	  }
	  int j_nonlocal = 0;
	  if ( (j_global < cvora[mpi_rank])||(j_global >= cvora[mpi_rank+1]) )
	    j_nonlocal = 1;
	  if (iter == 0) {
	    // element [i0,j] is always local...
	    // element [j,i0] may by non-local...
	    if (j_nonlocal) meCount++; 
	    // element [i1,j] may be non-local...
	    if (i1_nonlocal) meCount++;
	    // element [j,i1] may be non-local...	  
	    if (j_nonlocal) meCount++; 
	  }
	  else {
	    // 0.5*coeff*phi[j]*face_normal[ifa] gets added to i0_global's operator
	    // and subtracted from i1_global's...
	    double coeff[3];
	    FOR_I3 coeff[i] = 0.5*cvofa0_value[cof]*fa_normal[ifa][i];
	    // element [i0,j] is always local...	      
	    addMatrixEntry(DpDt,coeff,i0_global,j_global,nb2ocv_i,nb2ocv_v);
	    // element [j,i0] may by non-local...
	    if (j_nonlocal) {
	      matrixEntry[meCount].i_global = j_global;
	      matrixEntry[meCount].j_global = i0_global;
	      FOR_I3 matrixEntry[meCount].coeff[i] = coeff[i];
	      meCount++; 
	    }
	    else {
	      addMatrixEntry(DpDt,coeff,j_global,i0_global,nb2ocv_i,nb2ocv_v);
	    }
	
	    // the other 2 entries get the negative coeff...
	    FOR_I3 coeff[i] = -coeff[i];

	    // element [i1,j] may be non-local...
	    if (i1_nonlocal) {
	      matrixEntry[meCount].i_global = i1_global;
	      matrixEntry[meCount].j_global = j_global;
	      FOR_I3 matrixEntry[meCount].coeff[i] = coeff[i];
	      meCount++; 
	    }
	    else {
	      addMatrixEntry(DpDt,coeff,i1_global,j_global,nb2ocv_i,nb2ocv_v);
	    }
	    // element [j,i1] may be non-local...	  
	    if (j_nonlocal) {
	      matrixEntry[meCount].i_global = j_global;
	      matrixEntry[meCount].j_global = i1_global;
	      FOR_I3 matrixEntry[meCount].coeff[i] = coeff[i];
	      meCount++; 
	    }
	    else {		
	      addMatrixEntry(DpDt,coeff,j_global,i1_global,nb2ocv_i,nb2ocv_v);
	    }
	  }
	}
      }

      if (ifa >= nfa_bpi) {
	  
	// right side reconstruction...
	  
	const int cof_f = cvofa1_i[ifa];
	const int cof_l = cvofa1_i[ifa+1]-1;
	for (int cof = cof_f; cof <= cof_l; ++cof) {
	  const int icv = cvofa1_v[cof];
	  int j_global = cv_flag[icv];
	  if ( j_global < 0) {
	    // this must be a "fake" cv. Assume Neumann closure in building this
	    // operator for now...
	    assert( icv >= ncv_g );
	    j_global = i1_global;
	  }
	  int j_nonlocal = 0;
	  if ( (j_global < cvora[mpi_rank])||(j_global >= cvora[mpi_rank+1]) )
	    j_nonlocal = 1;
	  if (iter == 0) {
	    // element [i0,j] is always local...
	    // element [j,i0] may by non-local...
	    if (j_nonlocal) meCount++; 
	    // element [i1,j] may be non-local...
	    if (i1_nonlocal) meCount++;
	    // element [j,i1] may be non-local...	  
	    if (j_nonlocal) meCount++; 
	  }
	  else {
	    // 0.5*coeff*phi[j]*face_normal[ifa] gets added to i0_global's operator
	    // and subtracted from i1_global's...
	    double coeff[3];
	    FOR_I3 coeff[i] = 0.5*cvofa1_value[cof]*fa_normal[ifa][i];
	    // element [i0,j] is always local...	      
	    addMatrixEntry(DpDt,coeff,i0_global,j_global,nb2ocv_i,nb2ocv_v);
	    // element [j,i0] may by non-local...
	    if (j_nonlocal) {
	      matrixEntry[meCount].i_global = j_global;
	      matrixEntry[meCount].j_global = i0_global;
	      FOR_I3 matrixEntry[meCount].coeff[i] = coeff[i];
	      meCount++; 
	    }
	    else {
	      addMatrixEntry(DpDt,coeff,j_global,i0_global,nb2ocv_i,nb2ocv_v);
	    }

	    // the other 2 entries get the negative coeff...
	    FOR_I3 coeff[i] = -coeff[i];

	    // element [i1,j] may be non-local...
	    if (i1_nonlocal) {
	      matrixEntry[meCount].i_global = i1_global;
	      matrixEntry[meCount].j_global = j_global;
	      FOR_I3 matrixEntry[meCount].coeff[i] = coeff[i];
	      meCount++; 
	    }
	    else {
	      addMatrixEntry(DpDt,coeff,i1_global,j_global,nb2ocv_i,nb2ocv_v);
	    }
	    // element [j,i1] may be non-local...	  
	    if (j_nonlocal) {
	      matrixEntry[meCount].i_global = j_global;
	      matrixEntry[meCount].j_global = i1_global;
	      FOR_I3 matrixEntry[meCount].coeff[i] = coeff[i];
	      meCount++; 
	    }
	    else {		
	      addMatrixEntry(DpDt,coeff,j_global,i1_global,nb2ocv_i,nb2ocv_v);
	    }
	      
	  }
	}
      }
	
    }
      
    //cout << "meCount: " << meCount << endl;
      
    // on the first time, resize the non-local matrixEntry vector...
    if (iter == 0)
      matrixEntry.resize(meCount);
    else
      assert( matrixEntry.size() == meCount );
    
  }

  // now we have the packed vector matrixEntry, which we need to sort based on 
  // global row index so we can send it out to the relevant processors...
    
  std::sort(matrixEntry.begin(),matrixEntry.end(),MatrixEntryCompare());
    
  // the above sort allows us to figure out how many we are going to send 
  // to each processor in linear time... 

  int * send_count = new int[mpi_size];
  for (int i = 0; i < mpi_size; ++i)
    send_count[i] = 0;
  int irank = 0;
  for (int ime = 0; ime < matrixEntry.size(); ++ime) {
    while (matrixEntry[ime].i_global >= cvora[irank+1]) {
      irank++;
      assert( irank < mpi_size );
    }
    send_count[irank] += 1;
  }
    
  // send buffers - could save some space here and reuse the matrixEntry vector, but
  // clear it instead...
    
  int * send_ij = new int[matrixEntry.size()*2];
  double * send_coeff = new double[matrixEntry.size()*3];
    
  // populate the send buffers...
    
  for (int ime = 0; ime < matrixEntry.size(); ++ime) {
    send_ij[2*ime] = matrixEntry[ime].i_global;
    send_ij[2*ime+1] = matrixEntry[ime].j_global;
    FOR_I3 send_coeff[3*ime+i] = matrixEntry[ime].coeff[i];
  }
    
  // not needed any more...
    
  matrixEntry.clear();
    
  // complete send stuff and setup recv stuff...
    
  int * send_disp = new int[mpi_size];
  send_disp[0] = 0;
  for (int i = 1; i < mpi_size; ++i)
    send_disp[i] = send_count[i-1] + send_disp[i-1];
    
  int * recv_count = new int[mpi_size];
  MPI_Alltoall(send_count,1,MPI_INT,recv_count,1,MPI_INT,mpi_comm);
    
  int * recv_disp = new int[mpi_size];
  recv_disp[0] = 0;
  for (int i = 1; i < mpi_size; i++)
    recv_disp[i] = recv_count[i-1] + recv_disp[i-1];
    
  int recv_count_sum = recv_disp[mpi_size-1] + recv_count[mpi_size-1];
    
  int * recv_ij = new int[recv_count_sum*2];
  double * recv_coeff = new double[recv_count_sum*3];
    
  // multiply send_count and send_disp by 2...
  for (int i = 0; i < mpi_size; ++i) {
    send_count[i] *= 2;
    recv_count[i] *= 2;
    send_disp[i] *= 2;
    recv_disp[i] *= 2;
  }
    
  MPI_Alltoallv(send_ij,send_count,send_disp,MPI_INT,
		recv_ij,recv_count,recv_disp,MPI_INT,mpi_comm);
    
  // multiply send_count and send_disp by 3...
  for (int i = 0; i < mpi_size; ++i) {
    send_count[i] /= 2; send_count[i] *= 3;
    recv_count[i] /= 2; recv_count[i] *= 3;
    send_disp[i] /= 2; send_disp[i] *= 3;
    recv_disp[i] /= 2; recv_disp[i] *= 3;
  }
    
  MPI_Alltoallv(send_coeff,send_count,send_disp,MPI_DOUBLE,
		recv_coeff,recv_count,recv_disp,MPI_DOUBLE,mpi_comm);

  delete[] send_count; 
  delete[] send_disp;
  delete[] send_ij;
  delete[] send_coeff;
    
  for (int ime = 0; ime < recv_count_sum; ++ime) {
    int i_global = recv_ij[2*ime];
    assert( (i_global >= cvora[mpi_rank])&&(i_global < cvora[mpi_rank+1]) );
    int j_global = recv_ij[2*ime+1];
    double coeff[3];
    FOR_I3 coeff[i] = recv_coeff[3*ime+i];
    addMatrixEntry(DpDt,coeff,i_global,j_global,nb2ocv_i,nb2ocv_v);
  }

  delete[] recv_count; 
  delete[] recv_disp;
  delete[] recv_ij;
  delete[] recv_coeff;
    
  // we can now compute the norm of the symmetric part of the difference
  // operator in the direction of each face...

  // on boundaries, only the one side is valid...
  assert( fa_alpha == NULL );
  fa_alpha = new double[nfa];
  FOR_IFA fa_alpha[ifa] = 0.0;

  FOR_ICV {
      
    const int n2oc_f = nb2ocv_i[icv];
    const int n2oc_l = nb2ocv_i[icv+1]-1;

    const int foc_f = faocv_i[icv];
    const int foc_l = faocv_i[icv+1]-1;
      
    // cycle through the faces of this cv...
    for (int foc = foc_f; foc <= foc_l; ++foc) {
      const int ifa = faocv_v[foc];
      // area and unit normal...
      double area = sqrt( fa_normal[ifa][0]*fa_normal[ifa][0] +
			  fa_normal[ifa][1]*fa_normal[ifa][1] +
			  fa_normal[ifa][2]*fa_normal[ifa][2] );
      double nx = fa_normal[ifa][0]/area;
      double ny = fa_normal[ifa][1]/area;
      double nz = fa_normal[ifa][2]/area;
      // for each face, cycle through the elements of DpDt...
      for (int n2oc = n2oc_f; n2oc <= n2oc_l; ++n2oc) {
	if (nb2ocv_v[n2oc] == -1)
	  break;
	// compute the norm so there is no component cancellation...
	/*
	  fa_alpha[ifa] += 0.25*( (DpDt[n2oc][0]*nx)*(DpDt[n2oc][0]*nx) +
	  (DpDt[n2oc][1]*ny)*(DpDt[n2oc][1]*ny) +
	  (DpDt[n2oc][2]*nz)*(DpDt[n2oc][2]*nz) );
	*/
	// 
	// actually i think the following is more appropriate - makes norm for a face 
	// independent of its specific Cartiesian orientation...
	// note that 0.25 is 0.5**2, because each face should be visited twice.
	// 
	fa_alpha[ifa] += 0.25*( (DpDt[n2oc][0]*nx + DpDt[n2oc][1]*ny + DpDt[n2oc][2]*nz)*
				(DpDt[n2oc][0]*nx + DpDt[n2oc][1]*ny + DpDt[n2oc][2]*nz) );
	// further note:
	// it might be possible to reduce the dissipation in certain cases by adding the
	// left and right differencing operators at the face before taking the norm, but
	// that would involve alot of extra messaging for faces that don't have both
	// their cv's local, and would not change the fact that the tri meshes still have
	// some inherent symmetric part even in the prosm (perfect equilateral) mesh case.
      }
    }
      
  }

  // now boundary faces have only been visited once, internal faces twice, so...
    
  updateFaR1(fa_alpha,ADD_DATA);

  // double the boundary values...

  for (int ifa = 0; ifa < nfa_b; ++ifa) 
    fa_alpha[ifa] *= 2.0;

  // now normalize using the area...

  FOR_IFA {

    double area2 = 
      fa_normal[ifa][0]*fa_normal[ifa][0] +
      fa_normal[ifa][1]*fa_normal[ifa][1] +
      fa_normal[ifa][2]*fa_normal[ifa][2];
      
    // this has an arbitrary order (1) coeff in front...
      
    fa_alpha[ifa] = fa_alpha_coeff*sqrt( fa_alpha[ifa]/area2 );
      
    // the limiting value for fa_alpha should be 1.0. If it is greater than 1, then
    // we need to reduce the coeff corrections towards first order to yield 1...  

    if (fa_alpha[ifa] > 1.0) {
	
      /*
      // left operator always exists...
      {
	const int cof_f = cvofa0_i[ifa];
	const int cof_l = cvofa0_i[ifa+1]-1;
	cvofa0_value[cof_f] = 1.0;
	for (int cof = cof_f+1; cof <= cof_l; ++cof) { // skip first
	  double coeff = cvofa0_value[cof];
	  cvofa0_value[cof] /= fa_alpha[ifa];
	  cvofa0_value[cof_f] -= cvofa0_value[cof];
	}
      }
	
      // right operators for faces with both cv's local...
      if ( ifa >= nfa_bpi ) {
	const int cof_f = cvofa1_i[ifa];
	const int cof_l = cvofa1_i[ifa+1]-1;
	cvofa1_value[cof_f] = 1.0;
	for (int cof = cof_f+1; cof <= cof_l; ++cof) { // skip first
	  double coeff = cvofa1_value[cof];
	  cvofa1_value[cof] /= fa_alpha[ifa];
	  cvofa1_value[cof_f] -= cvofa1_value[cof];
	}
      }
      */
	
      // and set to 1.0...
      fa_alpha[ifa] = 1.0;
	
    }
      
  }
    
  /*
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if (ifa_debug >= 0)
  cout << "got fa_alpha for fa_debug: " << fa_alpha[ifa_debug] << endl;
  MPI_Pause("ifa_debug");
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  */
    
  // HACK: zero all fa_alpha's...
  //FOR_IFA fa_alpha[ifa] = 0.0;
    
  dumpScalarRange(fa_alpha,nfa,"FA_ALPHA");

  // take a look...
  /*
    FOR_IFA {
    cout << "ifa, normal, fa_alpha: " << ifa << " " << 
    fa_normal[ifa][0] << " " <<
    fa_normal[ifa][1] << " " <<
    fa_normal[ifa][2] << " " <<
    fa_alpha[ifa] << endl;
    }
  */
    
  delete[] DpDt;
  delete[] nb2ocv_i;
  delete[] nb2ocv_v;

}

void UgpWithCvFakeOp::simpleOperators() {

  // makes the reconstruction operator 1/2 + 1/2...
    
  if (mpi_rank == 0)
    cout << "simpleOperators()" << endl;
    
  // left side operators...
    
  for (int ifa = 0; ifa < nfa; ++ifa) {
    const int cof_f = cvofa0_i[ifa];
    const int cof_l = cvofa0_i[ifa+1]-1;
    for (int cof = cof_f; cof <= cof_l; ++cof) {
      const int icv = cvofa0_v[cof];
      if ( (icv == cvofa[ifa][0])||(icv == cvofa[ifa][1]) ) {
	cvofa0_value[cof] = 0.5;
      }
      else {
	cvofa0_value[cof] = 0.0;
      }
    }
  }

  // right side operators...
    
  for (int ifa = nfa_bpi; ifa < nfa; ++ifa) {
    const int cof_f = cvofa1_i[ifa];
    const int cof_l = cvofa1_i[ifa+1]-1;
    for (int cof = cof_f; cof <= cof_l; ++cof) {
      const int icv = cvofa1_v[cof];
      if ( (icv == cvofa[ifa][0])||(icv == cvofa[ifa][1]) ) {
	cvofa1_value[cof] = 0.5;
      }
      else {
	cvofa1_value[cof] = 0.0;
      }
    }
  }
    
}
  
void UgpWithCvFakeOp::firstOrderOperators() {
  
  // makes the reconstruction operator first order...
  
  if (mpi_rank == 0)
    cout << "firstOrderOperators()" << endl;

  FOR_IFA fa_alpha[ifa] = 1.0;
    
  // left side operators...
    
  for (int ifa = 0; ifa < nfa; ++ifa) {
    const int cof_f = cvofa0_i[ifa];
    const int cof_l = cvofa0_i[ifa+1]-1;
    for (int cof = cof_f; cof <= cof_l; ++cof) {
      const int icv = cvofa0_v[cof];
      if ( icv == cvofa[ifa][0] ) {
	assert( cof == cof_f );
	cvofa0_value[cof] = 1.0;
      }
      else {
	cvofa0_value[cof] = 0.0;
      }
    }
  }
    
  // right side operators...
    
  for (int ifa = nfa_bpi; ifa < nfa; ++ifa) {
    const int cof_f = cvofa1_i[ifa];
    const int cof_l = cvofa1_i[ifa+1]-1;
    for (int cof = cof_f; cof <= cof_l; ++cof) {
      const int icv = cvofa1_v[cof];
      if ( icv == cvofa[ifa][1] ) {
	assert( cof == cof_f );
	cvofa1_value[cof] = 1.0;
      }
      else {
	cvofa1_value[cof] = 0.0;
      }
    }
  }

  
    
}

void UgpWithCvFakeOp::firstOrderOperatorsFlaggedCvs() {
  
  // makes the reconstruction operator first order...
  
  if (mpi_rank == 0)
    cout << "firstOrderOperatorsFlaggedCvs()" << endl;

  FOR_IFA fa_flag[ifa] = 0;
  
  FOR_ICV {
    if (cv_flag[icv] != 0) {
      const int foc_f = faocv_i[icv];
      const int foc_l = faocv_i[icv+1]-1;
      for (int foc = foc_f; foc <= foc_l; ++foc) {
	const int ifa = faocv_v[foc];
	fa_flag[ifa] = 2; // use a 2 to indicate face next to a flagged cv
      }
    }
  }
  updateFaI1(fa_flag,MAX_DATA);
  
  // also set neighboring fa_alpha to 1.0...
  
  FOR_ICV {
    if (cv_flag[icv] == 0) {
      const int foc_f = faocv_i[icv];
      const int foc_l = faocv_i[icv+1]-1;
      int foc;
      for (foc = foc_f; foc <= foc_l; ++foc) {
	const int ifa = faocv_v[foc];
	if (fa_flag[ifa] == 2)
	  break;
      }
      if (foc <= foc_l) {
	// this means we broke out, so set any zero fa_flag's to 1
	for (foc = foc_f; foc <= foc_l; ++foc) {
	  const int ifa = faocv_v[foc];
	  if (fa_flag[ifa] == 0)
	    fa_flag[ifa] = 1;
	}
      }
    }
  }
  updateFaI1(fa_flag,MAX_DATA);

  // any non-zero fa_flag gets fa_alpha == 1...
  
  FOR_IFA {
    if (fa_flag[ifa] != 0) 
      fa_alpha[ifa] = 1.0;
  }

  // left side operators...
  
  for (int ifa = 0; ifa < nfa; ++ifa) {
    if (fa_flag[ifa] == 2) {
      const int cof_f = cvofa0_i[ifa];
      const int cof_l = cvofa0_i[ifa+1]-1;
      for (int cof = cof_f; cof <= cof_l; ++cof) {
	const int icv = cvofa0_v[cof];
	if ( icv == cvofa[ifa][0] ) {
	  assert( cof == cof_f );
	  cvofa0_value[cof] = 1.0;
	}
	else {
	  cvofa0_value[cof] = 0.0;
	}
      }
    }
  }
    
  // right side operators...
    
  for (int ifa = nfa_bpi; ifa < nfa; ++ifa) {
    if (fa_flag[ifa] == 2) {
      const int cof_f = cvofa1_i[ifa];
      const int cof_l = cvofa1_i[ifa+1]-1;
      for (int cof = cof_f; cof <= cof_l; ++cof) {
	const int icv = cvofa1_v[cof];
	if ( icv == cvofa[ifa][1] ) {
	  assert( cof == cof_f );
	  cvofa1_value[cof] = 1.0;
	}
	else {
	  cvofa1_value[cof] = 0.0;
	}
      }
    }
  }
    
}


