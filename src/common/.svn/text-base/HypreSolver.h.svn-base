#ifndef HYPRE_SOLVER_H
#define HYPRE_SOLVER_H


#ifdef WITH_HYPRE_VER_2_4
#define WITH_HYPRE
#endif


#ifdef WITH_HYPRE

#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

class HypreSolver {
  
private:

  int init_flag;

  HYPRE_IJMatrix Aij;
  HYPRE_IJVector bij;
  HYPRE_IJVector xij;

  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_ParVector par_b;
  HYPRE_ParVector par_x;

  HYPRE_Solver solver;
  HYPRE_Solver precon;

  double zeroAbs_; 
  double zeroRel_;
  int maxIter_;
  
  int currIter;

public:
  
  HypreSolver() : init_flag(0), zeroAbs_(1.0e-10), zeroRel_(1.0e-4), maxIter_(50), currIter(0) {
    
  }
  
  HypreSolver(const double zeroAbs, const double zeroRel, const int maxIter) 
    : init_flag(0), zeroAbs_(zeroAbs), zeroRel_(zeroRel), maxIter_(maxIter), currIter(0) {
    
  }

  void solveAMG(double * phi,
                double * A,double * rhs,
                int * noora,int * nbono_i,int * nbono_v,
                const double zero) {

    // Algebraic Multigrid interface to hypre
    //
    // although this routine is written with nodal notation (nbono_i/v etc), it
    // can be used for any parallel CSR (compact-storage-row) matrix and connectivity
    // (e.g. CV-based connectivity - i.e. nbocv_i/v structure), although you may need 
    // to eliminate duplications due to periodicity in some cases.
    //
    // noora[mpi_size+1] - "node-of-rank" an integer array describing the data layout (row
    //                     layout) across the processors. mpi_rank owns nodes 
    //
    //                          ino_global_first = noora[mpi_rank] up to and including 
    // 
    //                          ino_global_last = noora[mpi_rank+1]-1. 
    //
    //                     Note that noora[0] == 0, and and the global number of nodes will 
    //                     
    //                          nno_global = noora[mpi_size]
    //
    //                     Our local node size (nno_local, say) will be 
    //
    //                          nno_local = noora[mpi_rank+1]-noora[mpi_rank]. 
    //
    // phi[nno_local]      solution to [A]{phi} = {rhs}, also used as an initial guess 
    // A[nbono_i[nno_local]] CSR matrix 
    // rhs[nno_local]      right-hand side
    // noora               descibed above
    // nbono_i/v           connectivity structure, HOWEVER entries in nbono_v (i.e. the nodal
    //                     values) are stored in terms of global node number
    //                     On any processor, the global node number is related to the local 
    //                     node number as: 
    //
    //                           ino_global = ino_local + noora[mpi_rank]
    //
    // zero is the convergence tol
    //
    
    if (init_flag == 0) {
      
      init_flag = 1;

      // A must exist the first time...
      assert(A != NULL);
      initMatrixAndVectors(phi,A,rhs,noora,nbono_i,nbono_v);
      
      initAMGSolver(zero);
      
    }
    else {
      
      reinitVectors(phi,rhs,noora);
      
    }
    
    // solve...
    
    HYPRE_BoomerAMGSolve(solver,parcsr_A,par_b,par_x);
    
    // copy solution back to phi...

    getSolution(phi,noora);

  }
  
  void solveBCGStab(double * phi,
                double * A,double * rhs,
                int * noora,int * nbono_i,int * nbono_v,
                const double zero = -1.0) {
    
    if (init_flag == 0) {
      
      init_flag = 1;
      
      // A must exist the first time...
      assert(A != NULL);
      initMatrixAndVectors(phi,A,rhs,noora,nbono_i,nbono_v);
      
      initBCGStabSolver(zero);
      
    }
    else {
      
      reinitVectors(phi,rhs,noora);
      
      // HACK Rene: A not constant -> need to reinit A ???
      reinitMatrix(A,noora,nbono_i,nbono_v);
    }
    
    // solve...
    HYPRE_BiCGSTABSolve(solver, (HYPRE_Matrix)parcsr_A, 
      (HYPRE_Vector)par_b, (HYPRE_Vector)par_x); // par_b, par_x?
    
    // copy solution back to phi...
    
    getSolution(phi,noora);
    
  }
  
  void solveGMRESWithPILUTPrecon(double * phi, 
                                 double * A, double * rhs,
                                 int * noora, int * nbono_i, int * nbono_v, 
                                 const double zero = -1.0) {

    if (init_flag == 0) {

      init_flag = 1;

      // A must exist the first time...
      assert(A != NULL);
      initMatrixAndVectors(phi, A, rhs, noora, nbono_i, nbono_v);

      initGMRESWithPILUTPrecon(zero);
    }
    else {

      reinitVectors(phi,rhs,noora);

      // HACK Rene: A not constant -> need to reinit A ??? 
      reinitMatrix(A, noora,nbono_i,nbono_v);
    }

    // solve...
    HYPRE_GMRESSolve(solver, 
        (HYPRE_Matrix) parcsr_A, (HYPRE_Vector) par_b, (HYPRE_Vector) par_x);

    // copy solution back to phi...

    getSolution(phi, noora);

  }
 

// EUCLID precon: had to update to HYPRE version 2.4.0!!!
#ifdef WITH_HYPRE_VER_2_4
  
  void solveGMRESWithEuclidPrecon(double * phi, 
                                  double * A, double * rhs,
                                  int * noora, int * nbono_i, int * nbono_v, 
                                  const double zero = -1.0) {

    if (init_flag == 0) {

      init_flag = 1;

      // A must exist the first time...
      assert(A != NULL);
      initMatrixAndVectors(phi, A, rhs, noora, nbono_i, nbono_v);

      initGMRESWithEuclidPrecon(zero);
    }
    else {

      reinitVectors(phi,rhs,noora);

      // HACK Rene: A not constant -> need to reinit A ???
      reinitMatrix(A, noora, nbono_i, nbono_v);
    }
    
    if (currIter >= 10) 
      HYPRE_GMRESSetup(solver, (HYPRE_Matrix) parcsr_A, (HYPRE_Vector) par_b, (HYPRE_Vector) par_x);
    
    // solve...
    HYPRE_GMRESSolve(solver, 
        (HYPRE_Matrix) parcsr_A, (HYPRE_Vector) par_b, (HYPRE_Vector) par_x);
    
    
    HYPRE_GMRESGetNumIterations(solver, &currIter);

    double norm;
    HYPRE_GMRESGetFinalRelativeResidualNorm(solver, &norm);
    
    if (mpi_rank == 0)
      cout << "  HYPRE iter: " << currIter << ", norm: " << norm << endl;

    // copy solution back to phi...

    getSolution(phi, noora);

  }
#endif
  
  void solvePCG(double * phi,
                double * A,double * rhs,
                int * noora,int * nbono_i,int * nbono_v,
                const double zero) {
    
    if (init_flag == 0) {
      
      init_flag = 1;
      
      // A must exist the first time...
      assert(A != NULL);
      initMatrixAndVectors(phi,A,rhs,noora,nbono_i,nbono_v);
      
      initPCGSolver(zero);
      
    }
    else {
      
      reinitVectors(phi,rhs,noora);
      
    }
    
    // solve...
    
    HYPRE_ParCSRPCGSolve(solver,parcsr_A,par_b,par_x);
    
    // copy solution back to phi...

    getSolution(phi,noora);
    
  }

  void solvePCGWithAmgPrecon(double * phi,
                             double * A,double * rhs,
                             int * noora,int * nbono_i,int * nbono_v,
                             const double zero) {
    
    if (init_flag == 0) {
      
      init_flag = 1;

      // A must exist the first time...
      assert(A != NULL);
      initMatrixAndVectors(phi,A,rhs,noora,nbono_i,nbono_v);
      
      //initPCGWithParaSailsPreconSolver(zero);
      initPCGWithAMGPreconSolver(zero);

    }
    else {
      
      reinitVectors(phi,rhs,noora);
      
    }
    
    // solve...
    
    HYPRE_ParCSRPCGSolve(solver,parcsr_A,par_b,par_x);
    
    // copy solution back to phi...

    getSolution(phi,noora);
    
  }
  
  
private:

  
  void initMatrixAndVectors(double * phi,
                            double * A,double * rhs,
                            int * noora,int * nbono_i,int * nbono_v) {
    
    HYPRE_IJMatrixCreate(mpi_comm,
                         noora[mpi_rank],noora[mpi_rank+1]-1,
                         noora[mpi_rank],noora[mpi_rank+1]-1,&Aij);
    
    HYPRE_IJMatrixSetObjectType(Aij,HYPRE_PARCSR);
    
    
      
    int nno = noora[mpi_rank+1] - noora[mpi_rank];
    int * diag_sizes = new int[nno];
    int * offd_sizes = new int[nno];
    int row_nnz_max = 0;
    for (int ino = 0; ino < nno; ino++) {
      diag_sizes[ino] = offd_sizes[ino] = 0;
      int non_f = nbono_i[ino]; 
      int non_l = nbono_i[ino+1]-1; 
      row_nnz_max = max(row_nnz_max,non_l-non_f+1);
      for (int non = non_f; non <= non_l; non++) {
        int ino_nbr = nbono_v[non];
        if ((ino_nbr >= noora[mpi_rank])&&(ino_nbr < noora[mpi_rank+1]))
          diag_sizes[ino] += 1;
        else
          offd_sizes[ino] += 1;
      }
    }
    
    HYPRE_IJMatrixSetDiagOffdSizes(Aij,
                                   (const int *) diag_sizes,
                                   (const int *) offd_sizes);
    
    delete[] diag_sizes;
    delete[] offd_sizes;

    HYPRE_IJMatrixInitialize(Aij);

    int * cols = new int[row_nnz_max];
    double * values = new double[row_nnz_max];
    for (int ino = 0; ino < nno; ino++) {
      int non_f = nbono_i[ino]; 
      int non_l = nbono_i[ino+1]-1; 
      for (int non = non_f; non <= non_l; non++) {
        int ino_nbr = nbono_v[non];
        cols[non-non_f] = ino_nbr;
        values[non-non_f] = A[non];
      }
      int nnz = non_l-non_f+1;
      int row = ino + noora[mpi_rank];
      // first element of cols should be diagonal...
      // except in certain cases, so forget it
      //assert( cols[0] == row );
      // set the values...
      HYPRE_IJMatrixSetValues( Aij, 1, &nnz, &row,
                               (const int *) cols,
                               (const double *) values );
    }
      
    delete[] cols;
    delete[] values;
      
    HYPRE_IJMatrixAssemble(Aij);
      
    void *object;
    HYPRE_IJMatrixGetObject(Aij,&object);
    parcsr_A = (HYPRE_ParCSRMatrix)object;
      
    // ------------------------
    // and the vectors...
    // ------------------------
      
    // rhs: b...
      
    HYPRE_IJVectorCreate(mpi_comm,noora[mpi_rank],noora[mpi_rank+1]-1,&bij);
    HYPRE_IJVectorSetObjectType(bij,HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(bij);
    HYPRE_IJVectorSetValues(bij,nno,NULL,rhs);
    HYPRE_IJVectorAssemble(bij); // required?
    HYPRE_IJVectorGetObject(bij,(void **) &par_b);
      
    // x...
      
    HYPRE_IJVectorCreate(mpi_comm,noora[mpi_rank],noora[mpi_rank+1]-1,&xij);
    HYPRE_IJVectorSetObjectType(xij,HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(xij);
    HYPRE_IJVectorSetValues(xij,nno,NULL,phi);
    HYPRE_IJVectorAssemble(xij); // required?
    HYPRE_IJVectorGetObject(xij, (void **) &par_x);

  }
  
  

  void reinitVectors(double * phi,double * rhs,int * noora) {
    
    int nno = noora[mpi_rank+1] - noora[mpi_rank];
    
    // on subsequent calls, we do not need to initialize the solver or matrix,
    // just reset the rhs and initial guess...
    HYPRE_IJVectorInitialize(bij);
    HYPRE_IJVectorSetValues(bij,nno,NULL,rhs);
    HYPRE_IJVectorAssemble(bij); // required?
    HYPRE_IJVectorGetObject(bij, (void **) &par_b);


    HYPRE_IJVectorInitialize(xij);
    HYPRE_IJVectorSetValues(xij,nno,NULL,phi);
    HYPRE_IJVectorAssemble(xij); // required?
    HYPRE_IJVectorGetObject(xij, (void **) &par_x);


  }

  void reinitMatrix(double * A, int * noora, int * nbono_i, int * nbono_v) {

    HYPRE_IJMatrixInitialize(Aij);

    static int row_nnz_max = 0;
    int nno = noora[mpi_rank+1] - noora[mpi_rank];
    
    if (row_nnz_max == 0)
    {
      for (int ino = 0; ino < nno; ino++) 
      {
        int non_f = nbono_i[ino]; 
        int non_l = nbono_i[ino+1]-1; 
        row_nnz_max = max(row_nnz_max,non_l-non_f+1);
      }
    }

    static int * cols = new int[row_nnz_max];
    static double * values = new double[row_nnz_max];
    for (int ino = 0; ino < nno; ino++) 
    {
      int non_f = nbono_i[ino]; 
      int non_l = nbono_i[ino+1]-1; 
      for (int non = non_f; non <= non_l; non++) 
      {
        int ino_nbr = nbono_v[non];
        cols[non-non_f] = ino_nbr;
        values[non-non_f] = A[non];
      }
      int nnz = non_l-non_f+1;
      int row = ino + noora[mpi_rank];

      HYPRE_IJMatrixSetValues( Aij, 1, &nnz, &row, (const int *) cols, (const double *) values );
    }

    HYPRE_IJMatrixAssemble(Aij);
    HYPRE_IJMatrixGetObject(Aij, (void **) &parcsr_A);

//    delete[] cols;
//    delete[] values;

  }

  void getSolution(double * phi,int * noora) {
      
    // out local size...
    int nno = noora[mpi_rank+1] - noora[mpi_rank];
    
    // would like to do this, and it appears in atleast one
    // of the hypre example routines (ij_mv.c), but !@#$%
    // arg! passing NULL does not pull out the solution 
    // properly...
    //HYPRE_IJVectorGetValues(xij,nno,NULL,phi);
    
    // So... have to build rows[] containing the global index of
    // our processor's rows...
    int * rows = new int[nno];
    for (int ino = 0; ino < nno; ino++)
      rows[ino] = ino + noora[mpi_rank];
    HYPRE_IJVectorGetValues(xij,nno,rows,phi);
    delete[] rows;
    
  }
  

  void initAMGSolver(const double zero) {
    
    // ------------------------
    // solver...
    // ------------------------
      
    HYPRE_BoomerAMGCreate(&solver);
      
    /* Set some parameters (See Reference Manual for more parameters) */
      
    HYPRE_BoomerAMGSetPrintLevel(solver, 3);  /* 3 - initialization plus convergence every iter, 
                                                 2 - convergence every iter
                                                 1 - just the initialization, quiet during iters
                                              */

    //#ifdef _SKIP_THIS
    HYPRE_BoomerAMGSetCoarsenType(solver, 6); /* Falgout coarsening */
    HYPRE_BoomerAMGSetRelaxType(solver, 3);   /* G-S/Jacobi hybrid relaxation */
    //HYPRE_BoomerAMGSetNumSweeps(solver, 1);   /* Sweeps on each level */
    //#endif

    HYPRE_BoomerAMGSetCycleRelaxType(solver, 3, 3); /* Modify coarsest relaxation so it is NOT GE */
    HYPRE_BoomerAMGSetMaxLevels(solver, 20);  /* maximum number of levels */
    HYPRE_BoomerAMGSetMaxIter(solver, 50);
      
    HYPRE_BoomerAMGSetTol(solver, zero);      /* conv. tolerance */
    
    /* experiment with coarsening parameters... */
      
    HYPRE_BoomerAMGSetStrongThreshold(solver, 0.3); /* 0.5 better for 3D? */
      
    /* Now setup */
      
    HYPRE_BoomerAMGSetup(solver,parcsr_A,par_b,par_x);

  }

  void initGMRESWithPILUTPrecon(const double zero) {

    // ------------------------
    // solver...
    // ------------------------

    HYPRE_ParCSRGMRESCreate(mpi_comm, &solver);
    
    int k_dim = 5;
    HYPRE_GMRESSetKDim(solver, k_dim);
    
    HYPRE_GMRESSetMaxIter(solver, 1000);
    
    if (zero != -1.0)       HYPRE_GMRESSetTol(solver, zero);
    else                    HYPRE_GMRESSetTol(solver, zeroRel_);
    
    HYPRE_GMRESSetLogging(solver, 1);
    
    int ioutdat = 3;
    HYPRE_GMRESSetPrintLevel(solver, ioutdat);

    int ierr = HYPRE_ParCSRPilutCreate(mpi_comm, &precon);
    if (ierr)
    {
      printf("Error in ParPilutCreate\n");
    }

    HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSolve,
                                  (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSetup, precon);

    double drop_tol = -1;

    if (drop_tol >= 0)
      HYPRE_ParCSRPilutSetDropTolerance(precon, drop_tol);

    int nonzeros_to_keep = -1;

    if (nonzeros_to_keep >= 0)
      HYPRE_ParCSRPilutSetFactorRowSize(precon, nonzeros_to_keep);

    HYPRE_Solver precon_gotten;

    HYPRE_GMRESGetPrecond(solver, &precon_gotten);

    if (precon_gotten != precon)
      cout << "HYPRE_GMRESGetPrecond got bad precon" << endl;

    HYPRE_GMRESSetup(solver, 
        (HYPRE_Matrix) parcsr_A, (HYPRE_Vector) par_b, (HYPRE_Vector) par_x);

  }

// EUCLID precon: had to update to HYPRE version 2.4.0!!!
#ifdef WITH_HYPRE_VER_2_4
  
  void initGMRESWithEuclidPrecon(const double zero) {

    
    // ------------------------
    // precon...
    // ------------------------

    /*
    double sai_threshold = 0.1;
    double sai_filter = 0.1;
    int poutdat = 1;
    int max_levels = 1;
    HYPRE_ParaSailsCreate(mpi_comm, &precon);
    
    HYPRE_ParaSailsSetParams(precon, sai_threshold, max_levels);
    HYPRE_ParaSailsSetFilter(precon, sai_filter);
    HYPRE_ParaSailsSetLogging(precon, poutdat);
    HYPRE_ParaSailsSetSym(precon, 0);
    */
    
    
    
    HYPRE_EuclidCreate(mpi_comm, &precon);

    int eu_level = 3;           // default = -1;
    double eu_ilut = 0.0;
    double eu_sparse_A = 0.0;
    int eu_row_scale = 0;
    int eu_bj = 0;

    if (eu_level > -1)    HYPRE_EuclidSetLevel(precon, eu_level);
    if (eu_ilut)          HYPRE_EuclidSetILUT(precon, eu_ilut);
    if (eu_sparse_A)      HYPRE_EuclidSetSparseA(precon, eu_sparse_A);
    if (eu_row_scale)     HYPRE_EuclidSetRowScale(precon, eu_row_scale);
    if (eu_bj)            HYPRE_EuclidSetBJ(precon, eu_bj);

    int eu_stats = 3;
    int eu_mem = 0;

    HYPRE_EuclidSetStats(precon, eu_stats);
    HYPRE_EuclidSetMem(precon, eu_mem);

    
    
    // ------------------------
    // solver...
    // ------------------------

    HYPRE_ParCSRGMRESCreate(mpi_comm, &solver);

    int k_dim = 10;
    HYPRE_GMRESSetKDim(solver, k_dim);

    HYPRE_GMRESSetMaxIter(solver, maxIter_);

    HYPRE_GMRESSetAbsoluteTol(solver, zeroAbs_);

    if (zero != -1.0)       HYPRE_GMRESSetTol(solver, zero);
    else                    HYPRE_GMRESSetTol(solver, zeroRel_);

    HYPRE_GMRESSetLogging(solver, 1);

    int ioutdat = 1;
    HYPRE_GMRESSetPrintLevel(solver, ioutdat);

//    HYPRE_GMRESSetRelChange(solver, 1);

/*    HYPRE_GMRESSetPrecond(solver,
                         (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
                         (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup,
                          precon);*/
    
    HYPRE_GMRESSetPrecond(solver,
                         (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,
                         (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup,
                          precon);

    HYPRE_Solver precon_gotten;

    HYPRE_GMRESGetPrecond(solver, &precon_gotten);

    if (precon_gotten != precon)
      cout << "HYPRE_GMRESGetPrecond got bad precon" << endl;

    HYPRE_GMRESSetup(solver, 
        (HYPRE_Matrix) parcsr_A, (HYPRE_Vector) par_b, (HYPRE_Vector) par_x);
  }
  
#endif

  void initBCGStabSolver(const double zero) {
    
    // ------------------------
    // solver...
    // ------------------------

    HYPRE_ParCSRBiCGSTABCreate(mpi_comm,&solver);

    
    if (zero != -1.0)
    {
      HYPRE_BiCGSTABSetMaxIter(solver,1000);
      HYPRE_BiCGSTABSetTol(solver, zero);
    }
    else
    {
#ifdef WITH_HYPRE_VER_2_4   // only available in hypre version 2_4
      HYPRE_BiCGSTABSetAbsoluteTol(solver,zeroAbs_);
#endif
      HYPRE_BiCGSTABSetTol(solver, zeroRel_);
      HYPRE_BiCGSTABSetMaxIter(solver,maxIter_);
    }
    
    //HYPRE_BiCGSTABSetLogging(solver,ioutdat);
    HYPRE_BiCGSTABSetPrintLevel(solver,3);
    
    
    precon = NULL;
    HYPRE_BiCGSTABSetPrecond(solver,
			     (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScale,
			     (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScaleSetup,
			     precon);
    
    HYPRE_BiCGSTABSetup(solver, (HYPRE_Matrix)parcsr_A, 
			(HYPRE_Vector)par_b, (HYPRE_Vector)par_x); // par_b, par_x?
    
  }

  void initPCGSolver(const double zero) {
    
    // ------------------------
    // solver...
    // ------------------------

    HYPRE_ParCSRPCGCreate(mpi_comm,&solver);
      
    /* Set some parameters (See Reference Manual for more parameters) */
    HYPRE_PCGSetMaxIter(solver,5000); /* max iterations */
    HYPRE_PCGSetTol(solver,zero); /* conv. tolerance */
    HYPRE_PCGSetTwoNorm(solver,1); /* use the two norm as the stopping criteria */
    HYPRE_PCGSetPrintLevel(solver,1); /* print solve info: 1: quiet, 3: every iter */
    /* HYPRE_PCGSetLogging(*solver, 1); */ /* needed to get run info later */
    
    /* diagonal scaling */
    precon = NULL;
    HYPRE_PCGSetPrecond(solver,
                        (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScale,
                        (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScaleSetup,
                        precon);

    /* Now setup (and solve! later) */
    HYPRE_ParCSRPCGSetup(solver,parcsr_A,par_b,par_x);
  
  }


  void initPCGWithAMGPreconSolver(const double zero) {

    HYPRE_ParCSRPCGCreate(mpi_comm,&solver);
    
    /* Set some parameters (See Reference Manual for more parameters) */
    HYPRE_PCGSetMaxIter(solver,1000); /* max iterations */
    HYPRE_PCGSetTol(solver,zero); /* conv. tolerance */
    HYPRE_PCGSetTwoNorm(solver,1); /* use the two norm as the stopping criteria */
    HYPRE_PCGSetPrintLevel(solver,1); /* print solve info: 1 quiet, 3 every iter */
    /* HYPRE_PCGSetLogging(*solver, 1); */ /* needed to get run info later */
  
    /* Now set up the AMG preconditioner and specify any parameters */
    HYPRE_BoomerAMGCreate(&precon);
    /* use 2 here - just the grid hierarchy at the start... */
    HYPRE_BoomerAMGSetPrintLevel(precon, 1); /* print amg solution info */
    HYPRE_BoomerAMGSetCoarsenType(precon, 6);
    HYPRE_BoomerAMGSetRelaxType(precon, 6); /* Sym G.S./Jacobi hybrid */
    HYPRE_BoomerAMGSetNumSweeps(precon, 1);
    HYPRE_BoomerAMGSetTol(precon, 0.0); /* conv. tolerance zero */
    HYPRE_BoomerAMGSetMaxIter(precon, 1); /* do only one iteration! */
  
    /* Set the PCG preconditioner */
    HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
                        (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precon);
    
    /* Now setup (and solve! later) */
    HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
 
  }

  void initPCGWithParaSailsPreconSolver(const double zero) {

    HYPRE_ParCSRPCGCreate(mpi_comm,&solver);
    
    /* Set some parameters (See Reference Manual for more parameters) */
    HYPRE_PCGSetMaxIter(solver,1000); /* max iterations */
    HYPRE_PCGSetTol(solver,zero); /* conv. tolerance */
    HYPRE_PCGSetTwoNorm(solver,1); /* use the two norm as the stopping criteria */
    HYPRE_PCGSetPrintLevel(solver,3); /* print solve info */
    /* HYPRE_PCGSetLogging(*solver, 1); */ /* needed to get run info later */

    int max_levels = 1;
    double   sai_threshold = 0.1;
    double   sai_filter = 0.1;

    HYPRE_ParaSailsCreate(mpi_comm, &precon);
    HYPRE_ParaSailsSetParams(precon, sai_threshold, max_levels);
    HYPRE_ParaSailsSetFilter(precon, sai_filter);
    /* HYPRE_ParaSailsSetLogging(precon, poutdat); */
    
    HYPRE_PCGSetPrecond(solver,
                        (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
                        (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup,
                        precon);

  }



};

#else

// if the user has not built the hypre libraries, provide an empty
// HypreSolver class that errors is instantiated. This means you
// can still explicitly include a hypre solver as a POINTER in your
// solver, but do not instantiate it...

class HypreSolver {
  
public:
  
  HypreSolver() {
    
    std::cerr << "Error: to use hypre you need to compile WITH_HYPRE." << std::endl;
    throw(-1);
 
  }
  
  HypreSolver(const double zeroAbs, const double zeroRel, const int maxIter) {
      
    std::cerr << "Error: to use hypre you need to compile WITH_HYPRE." << std::endl;
    throw(-1);
 
  }

  // only need to provide blank interface to public methods...

  void solveAMG(double * phi,
                double * A,double * rhs,
                int * noora,int * nbono_i,int * nbono_v,
                const double zero) {

    
  }
  
  void solveBCGStab(double * phi,
                    double * A,double * rhs,
                    int * noora,int * nbono_i,int * nbono_v,
                    const double zero = -1.0) {
    
  }

  void solveGMRESWithPILUTPrecon(double * phi, 
                                 double * A, double * rhs,
                                 int * noora, int * nbono_i, int * nbono_v, 
                                 const double zero = -1.0) {
    
  }
  
  void solveGMRESWithEuclidPrecon(double * phi, 
                                  double * A, double * rhs,
                                  int * noora, int * nbono_i, int * nbono_v, 
                                  const double zero = -1.0) {
    
  }
  
  void solvePCG(double * phi,
                double * A,double * rhs,
                int * noora,int * nbono_i,int * nbono_v,
                const double zero) {

  }

  void solvePCGWithAmgPrecon(double * phi,
                             double * A,double * rhs,
                             int * noora,int * nbono_i,int * nbono_v,
                             const double zero) {
  }


};

#endif
  

#endif
