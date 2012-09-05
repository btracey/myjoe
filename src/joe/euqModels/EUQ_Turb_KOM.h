#ifndef EUQ_TURB_KOM_H
#define EUQ_TURB_KOM_H

#include "turbModels/TurbModel_KOM.h"

#include "UgpWithCvCompFlow.h"
#include "JoeWithModels.h"


class EUQ_Turb_KOM :virtual public JoeWithModels, virtual public RansTurbKOm
{
public:   // constructors

	EUQ_Turb_KOM(char *name) : JoeWithModels(name) , UgpWithCvCompFlow(name)
	{
		if (mpi_rank == 0)
		{
			cout << "EUQ_Turb_KOM()" << endl;

			// output what src terms are active
			if ( checkParam("EUQ_SRC_MOMENTUM") ) cout << "EUQ_SRC_MOMENTUM : active" << endl;
			if ( checkParam("EUQ_SRC_ENERGY") )   cout << "EUQ_SRC_ENERGY   : active" << endl;

			// output how eigenvalue perturbations will occur
			if (checkParam("EUQ_EVAL"))
			{
				// optional spit out input file params here
				string eval_marker     = getParam("EUQ_EVAL")->getString("MARKER");
				string eval_interp     = getParam("EUQ_EVAL")->getString("INTERP");
				string eval_mode       = getParam("EUQ_EVAL")->getString("MODE");
				double eval_cut        = getParam("EUQ_EVAL")->getDouble("CUTOFF");
				double eval_delta_max  = getParam("EUQ_EVAL")->getDouble("DELTA_MAX");
				int eval_corner        = getParam("EUQ_EVAL")->getInt("CORNER");

				cout << "EUQ_EVAL " << endl;
				cout << "    mode         : " << eval_mode << endl;
				cout << "    marker       : " << eval_marker << endl;
				cout << "    interpolation: " << eval_interp << endl;
				cout << "    cutoff       : " << eval_cut << endl;
				cout << "    delta_max    : " << eval_delta_max << endl;
				cout << "    corner       : " << eval_corner << endl;


			}
			if (checkParam("EUQ_EVEC"))
			{
				cout << "EUQ_EVEC: active" << endl;
			}
			if (checkParam("EUQ_K"))
			{
				cout << "EUQ_K: active" << endl;
			}

			// enforcment/flagging of realizability via the barycentric map
			if (checkParam("FORCE_REALIZABILITY"))
			{
				string realize_method = getParam("FORCE_REALIZABILITY")->getString("METHOD");
				cout << "FORCE_REALIZABILITY method: " << realize_method << endl;
			}
		}

		// register turbulence stress tensor terms
		rij_diag = NULL;     registerVector(rij_diag,"rij_diag",CV_DATA);
		rij_offdiag = NULL;  registerVector(rij_offdiag,"rij_offdiag",CV_DATA);
		sij_diag = NULL;     registerVector(sij_diag,"sij_diag",CV_DATA);
		sij_offdiag = NULL;  registerVector(sij_offdiag,"sij_offdiag",CV_DATA);

		// reigster turbulence source contribution terms
		srcMomentum = NULL;  registerVector(srcMomentum,"srcMomentum",CV_DATA);
		srcEnergy = NULL;    registerScalar(srcEnergy,"srcEnergy",CV_DATA);

		// register EUQ required variables
		baryCo = NULL;	     registerVector(baryCo,"baryCo",CV_DATA);
		colorC = NULL;	     registerVector(colorC,"colorC",CV_DATA);
		d_bary = NULL;       registerScalar(d_bary,"d_bary",CV_DATA);
		eval_sensor = NULL;  registerScalar(eval_sensor,"eval_sensor",CV_DATA);
		realizable = NULL;   registerScalar(realizable,"realizable",CV_DATA);

		// extras
		PK = NULL;           registerScalar(PK,"PK",CV_DATA);
		P_eps = NULL;        registerScalar(P_eps,"P_eps",CV_DATA);

		/* Define barycentric traingle corner points
		 * 1 comp = corner[0][:], 2 comp = corner[1][:], 3 comp = corner[2][:]
		 * for second index: 0=x, 1=y location
		 */
		corner[0][0] = 1.0;
		corner[0][1] = 0.0;
		corner[1][0] = 0.0;
		corner[1][1] = 0.0;
		corner[2][0] = 0.5;
		corner[2][1] = 0.866025;
	}

	virtual ~EUQ_Turb_KOM() {}

public:   // member vars

	// instanciate turbulence stress tensor terms
	double (*rij_diag)[3];        ///< reynolds stress holders
	double (*rij_offdiag)[3];     ///< reynolds stress holders
	double (*sij_diag)[3];        ///< mean velocity grad holders
	double (*sij_offdiag)[3];     ///< mean velocity grad holders

	// instanciate turbulence source contribution terms
	double (*srcMomentum)[3];     ///< momentum source
	double *srcEnergy;            ///< energy source

	// instanciate turbulent EUQ required variables
	double (*baryCo)[3];		  /// Barycentric map coordinates
	double (*colorC)[3];		  /// Barycentric map colors
	double *d_bary;               ///< barycentric perturbation variable
	double corner[3][2];          ///< barycentric triangle coordinates
	double *eval_sensor;          ///< holder for eigenvalue euq sensor
	double *realizable;			///< flag for whether state is realizable or not

	// extras
	double *PK;                   ///< turbulence production
	double *P_eps;				  ///< production - dissipation

public:   // member functions

	/**
	 * source term called when backward_euler_semincoupled timestepping is used
	 */
	virtual void sourceHookCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
	{
		// generate delta_R_ij, and thus momentum source (srcMomentum) and (eventually) energy source (srcEnergy) terms
		perturbStress();

		if ( checkParam("EUQ_SRC_MOMENTUM") || checkParam("EUQ_SRC_ENERGY") )
		{

			for (int icv = 0; icv < ncv; icv++)
			{
				if ( checkParam("EUQ_SRC_MOMENTUM") )
				{
					for (int i=0; i<3; i++)
					{
						// add momentum source (will be zero if no perturbation)
						rhs[icv][i+1] -= srcMomentum[icv][i] * cv_volume[icv];
					}
				}

				if ( checkParam("EUQ_SRC_ENERGY") )
				{
					// currently not implemented for the energy equation
					// rhs[icv][4] -= srcEnergy[icv] * cv_volume[icv];
				}
			}
		}
	};

	/**
	 * source term called when backward_euler timestepping is used
	 */
	virtual void sourceHook(double *rhs_rho, double (*rhs_rhou)[3], double *rhs_rhoE, double (*A)[5][5])
	{
		// generate delta_R_ij, and thus momentum source (srcMomentum) and (eventually) energy source (srcEnergy) terms
		perturbStress();

		if ( checkParam("EUQ_SRC_MOMENTUM") || checkParam("EUQ_SRC_ENERGY") )
		{

			for (int icv = 0; icv < ncv; icv++)
			{
				if ( checkParam("EUQ_SRC_MOMENTUM") )
				{
					for (int i=0; i<3; i++)
					{
						// add momentum source (will be zero if no perturbation)
						rhs_rhou[icv][i] -= srcMomentum[icv][i] * cv_volume[icv];
					}
				}

				if ( checkParam("EUQ_SRC_ENERGY") )
				{
					// currently not implemented for the energy equation
					// rhs_rhoE[icv] -= srcEnergy[icv] * cv_volume[icv];
				}
			}
		}
	};


	/**
	 *
	 * perturb stresses for epistemic UQ
	 * the perturbation of the Reynolds stresses is only applied in the momentum equations, energy equation is left unaltered
	 *
	 * the UQ momentum equation reads:
	 *
	 *    drui            d                                |                  d                              |
	 *    ---- + conv =  --- (2*mut*Sij - 2/3 rho k dij) + | Tij_perturbed - --- (2*mut*Sij - 2/3 rho k dij) |
	 *     dt            dxj                               |                 dxj                             |
	 *
	 *                   \_____________________________/    \________________________________________________/
	 *                                  |                                           |
	 *                       keep term in viscous fluxes              add only difference to source term
	 *
	 *    Tij_perturbed is then (-rho ui uj)
	 *
	 */
	void perturbStress()
	{

		calcGradVel();

		for (int icv=0; icv<ncv; icv++)
		{

			//####################################################################################
			// build Rate-of-Strain tensor S_ij
			// S_ij = 0.5 * (du_i/dx_j + du_j/dx_i) - 1/3 * du_k/dx_k * d_ij (compressible)
			// S_ij = 0.5 * (du_i/dx_j + du_j/dx_i)                          (incompressible)

			double div = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

			// Diagonal Terms
			// compressible version
			sij_diag[icv][0] = (grad_u[icv][0][0] - 1.0/3.0*div);
			sij_diag[icv][1] = (grad_u[icv][1][1] - 1.0/3.0*div);
			sij_diag[icv][2] = (grad_u[icv][2][2] - 1.0/3.0*div);

			/* incompressible version
     	    sij_diag[icv][0] = grad_u[icv][0][0];
     	    sij_diag[icv][1] = grad_u[icv][1][1];
    	 	sij_diag[icv][2] = grad_u[icv][2][2];
			 */

			// Off-Diagonal Terms
			sij_offdiag[icv][0] = 0.5 * (grad_u[icv][0][1] + grad_u[icv][1][0]);
			sij_offdiag[icv][1] = 0.5 * (grad_u[icv][0][2] + grad_u[icv][2][0]);
			sij_offdiag[icv][2] = 0.5 * (grad_u[icv][1][2] + grad_u[icv][2][1]);

			//####################################################################################
			// build compressible Reynolds stresses from eddy viscosity model (given grad_u, mu_t, kine)
			// <u_i*u_j> = 2/3*k*delta_ij - nu_t (du_i/dx_j + du_j/dx_i - 2/3 du_k/dx_k d_ij)
			// or similarly using S_ij
			// <u_i*u_j> = 2/3*k*delta_ij - 2 * nu_t * S_ij


			double term1 = (2.0/3.0) * kine[icv];

			// Diagonal Terms
			rij_diag[icv][0] = term1 - (muT[icv]/rho[icv]) * 2.0 * sij_diag[icv][0];
			rij_diag[icv][1] = term1 - (muT[icv]/rho[icv]) * 2.0 * sij_diag[icv][1];
			rij_diag[icv][2] = term1 - (muT[icv]/rho[icv]) * 2.0 * sij_diag[icv][2];

			// Off-Diagonal Terms
			rij_offdiag[icv][0] =     -(muT[icv]/rho[icv]) * 2.0 * sij_offdiag[icv][0];
			rij_offdiag[icv][1] =     -(muT[icv]/rho[icv]) * 2.0 * sij_offdiag[icv][1];
			rij_offdiag[icv][2] =     -(muT[icv]/rho[icv]) * 2.0 * sij_offdiag[icv][2];

			//####################################################################################
			// Compute anisotropy tensor

			double aij[3][3],eigv[3][3],eigvt[3][3],eigs[3],eigsnew[3][3];

			for (int i=0; i < 3; i++)
			{
				for (int j=0; j < 3; j++)
				{
					aij[i][j] = 0.0;
					eigv[i][j] = 0.0;
					eigvt[i][j] = 0.0;
					eigsnew[i][j] = 0.0;
				}
				eigs[i] = 0.0;
			}

			// store current value of kine based on r_ij
			double tempk = 0.5 * (rij_diag[icv][0] + rij_diag[icv][1] + rij_diag[icv][2]);

			// build normalized turbulence anisotropy tensor
			// a_ij = r_ij / (2 * k) - d_ij / 3


			// Diagonal Terms
			aij[0][0] = 0.5 * rij_diag[icv][0] / tempk - 1.0/3.0;
			aij[1][1] = 0.5 * rij_diag[icv][1] / tempk - 1.0/3.0;
			aij[2][2] = 0.5 * rij_diag[icv][2] / tempk - 1.0/3.0;

			// Off-Diagonal Terms
			aij[0][1] = 0.5 * rij_offdiag[icv][0] / tempk;
			aij[0][2] = 0.5 * rij_offdiag[icv][1] / tempk;
			aij[1][2] = 0.5 * rij_offdiag[icv][2] / tempk;
			aij[1][0] = aij[0][1];
			aij[2][0] = aij[0][2];
			aij[2][1] = aij[1][2];

			//####################################################################################
			//Compute barycentric map location

			//compute eigenvalues of aij
			eigen_decomposition(aij,eigv,eigs);

			// generate transposeof eigenvalues
			for (int i=0; i < 3; i++)
				for (int j=0; j < 3; j++)
					eigvt[j][i] = eigv[i][j];

			/* eigenvalues are automatically sorted by QR algorithm
			 * such that eigs[2] > eigs[1] > eigs[0]
			 */

			// compute convex combination coefficients
			double c1c = eigs[2] - eigs[1];
			double c2c = 2.0 * (eigs[1] - eigs[0]);
			double c3c = 3.0 * eigs[0] + 1.0;

			// save barycentric coordinates
			// baryCo[icv][0] = 1.0 * c1c + 0.0 * c2c + 0.5 * c3c;
			// baryCo[icv][1] = 0.0 * c1c + 0.0 * c2c + 0.866025 * c3c;
			baryCo[icv][0] = corner[0][0] * c1c + corner[1][0] * c2c + corner[2][0] * c3c;
			baryCo[icv][1] = corner[0][1] * c1c + corner[1][1] * c2c + corner[2][1] * c3c;
			baryCo[icv][2] = 0.0;


			// calculate turbulence production for each cell
			// PK = -<u_i u_j> * d<U_i>/dx_j = -r_ij * grad_u[i][j]
			PK[icv] = 0.0;
			for (int i=0; i < 3; i++)
			{
				// diagonal terms
				PK[icv] -= rij_diag[icv][i] * grad_u[icv][i][i];
			}

			// off-diagonal terms (counted twice due to symmetry)
			PK[icv] -= 2 * rij_offdiag[icv][0] * grad_u[icv][0][1];
			PK[icv] -= 2 * rij_offdiag[icv][1] * grad_u[icv][0][2];
			PK[icv] -= 2 * rij_offdiag[icv][2] * grad_u[icv][1][2];

			//      PK[icv] = muT[icv] * (grad_u[icv][0][1] - grad_u[icv][1][0]) * (grad_u[icv][0][1] - grad_u[icv][1][0]);

			P_eps[icv] = PK[icv] - ( omega[icv] * kine[icv] );

			//####################################################################################
			// Introduce perturbations as specified by Joe.in File

			// build holder for new eigenvalues, new eigenvectors, new k, temporary anisotropy
			double eigs_new[3][3],eigv_new[3][3],eigvt_new[3][3],k_new, aijt[3][3];

			for (int i=0; i < 3; i++)
			{
				for (int j=0; j < 3; j++)
				{
					eigs_new[i][j] = 0.0;
					aijt[i][j] = 0.0;
					// initialize eigenvector to current values, in case of no perturbation
					eigv_new[i][j] = eigv[i][j];
				}
				// initialize eigenvalues to current values, in case of no perturbation
				eigs_new[i][i] = eigs[i];
			}
			// initialize k to current value, in case of no perturbation
			k_new = tempk;

			if (checkParam("FORCE_REALIZABILITY"))
			{
				forceRealizability(icv,eigs_new);
			}

			if (checkParam("EUQ_EVAL"))
			{
				perturbEvals(icv,c1c,c2c,c3c,eigs_new);
			}

			if (checkParam("EUQ_EVEC"))
			{
				perturbEvectors(icv);
			}

			if (checkParam("EUQ_K"))
			{
				perturbK(icv);
			}
			// generate transpose of new eigenvalues
			for (int i=0; i < 3; i++)
				for (int j=0; j < 3; j++)
					eigvt_new[j][i] = eigv_new[i][j];

			// build aij from new eigenvalues, new eigenvectors, new k
			// aij = eigv * eigs * eigv'
			matMult(aijt,eigv_new,eigs_new);
			matMult(aij,aijt,eigvt_new);

			// build reynolds stresses from reconstructed aij, k_new
			rij_diag[icv][0]    = 2.0 * k_new * (aij[0][0] + 1.0/3.0);
			rij_diag[icv][1]    = 2.0 * k_new * (aij[1][1] + 1.0/3.0);
			rij_diag[icv][2]    = 2.0 * k_new * (aij[2][2] + 1.0/3.0);
			rij_offdiag[icv][0] = 2.0 * k_new * aij[0][1];
			rij_offdiag[icv][1] = 2.0 * k_new * aij[0][2];
			rij_offdiag[icv][2] = 2.0 * k_new * aij[1][2];

			// store convex combination coefficients for each cell value (post perturbation)
			colorC[icv][0] = c1c;
			colorC[icv][1] = c2c;
			colorC[icv][2] = c3c;

		} // end icv loop

		// open memory for stresses as scalars, open with static so the memory will be declared only once, don't delete memory at the end!
		static double *rhouu = new double [ncv_g];   static double *rhouu_fa = new double [nfa_b];
		static double *rhouv = new double [ncv_g];   static double *rhouv_fa = new double [nfa_b];
		static double *rhouw = new double [ncv_g];   static double *rhouw_fa = new double [nfa_b];
		static double *rhovv = new double [ncv_g];   static double *rhovv_fa = new double [nfa_b];
		static double *rhovw = new double [ncv_g];   static double *rhovw_fa = new double [nfa_b];
		static double *rhoww = new double [ncv_g];   static double *rhoww_fa = new double [nfa_b];

		// ------------------------------------------------------------------
		// compute (Tau_ij_epistemic - Tau_ij_eddyVisc)
		// ------------------------------------------------------------------
		for (int icv = 0; icv < ncv; icv++)
		{
			double div = grad_u[icv][0][0] + grad_u[icv][1][1] + grad_u[icv][2][2];

			double term1 = (2.0/3.0) * kine[icv];

			/* subtract turbulence model momentum contribution
			 * which has already been applied as a source from the newly
			 * computed reynolds stress. This gives a delta Reynolds stress
			 * value which we will use as an additional source.
			 */
			rij_diag[icv][0]    -= term1 - (muT[icv]/rho[icv]) * 2.0 * sij_diag[icv][0];
			rij_diag[icv][1]    -= term1 - (muT[icv]/rho[icv]) * 2.0 * sij_diag[icv][1];
			rij_diag[icv][2]    -= term1 - (muT[icv]/rho[icv]) * 2.0 * sij_diag[icv][2];
			rij_offdiag[icv][0] -=       - (muT[icv]/rho[icv]) * 2.0 * sij_offdiag[icv][0];
			rij_offdiag[icv][1] -=       - (muT[icv]/rho[icv]) * 2.0 * sij_offdiag[icv][1];
			rij_offdiag[icv][2] -=       - (muT[icv]/rho[icv]) * 2.0 * sij_offdiag[icv][2];

			// build rho*<u'u'>
			rhouu[icv] = rij_diag[icv][0]    * rho[icv];
			rhovv[icv] = rij_diag[icv][1]    * rho[icv];
			rhoww[icv] = rij_diag[icv][2]    * rho[icv];
			rhouv[icv] = rij_offdiag[icv][0] * rho[icv];
			rhouw[icv] = rij_offdiag[icv][1] * rho[icv];
			rhovw[icv] = rij_offdiag[icv][2] * rho[icv];

		} // end icv loop

		// update terms for source computation
		updateCvData(rhouu, REPLACE_DATA);
		updateCvData(rhouv, REPLACE_DATA);
		updateCvData(rhouw, REPLACE_DATA);
		updateCvData(rhovv, REPLACE_DATA);
		updateCvData(rhovw, REPLACE_DATA);
		updateCvData(rhoww, REPLACE_DATA);

		// update computed perturbation values
		updateCvData(rij_diag, REPLACE_ROTATE_DATA);
		updateCvData(rij_offdiag, REPLACE_ROTATE_DATA);
		updateCvData(sij_diag, REPLACE_ROTATE_DATA);
		updateCvData(sij_offdiag, REPLACE_ROTATE_DATA);

		updateCvData(baryCo, REPLACE_ROTATE_DATA);
		updateCvData(colorC, REPLACE_ROTATE_DATA);
		updateCvData(d_bary, REPLACE_DATA);
		updateCvData(eval_sensor, REPLACE_DATA);

		updateCvData(PK, REPLACE_DATA);
		updateCvData(P_eps, REPLACE_DATA);

		// ------------------------------------------------------------------
		// new turbulent stresses are saved as cell values, but have not influenced solution (at faces)
		// Do we need to put loop at faces for FV consistency?!?
		// ------------------------------------------------------------------
		// set stresses at boundaries, set at wall=zero, otherwise "Von Neumann BC"
		for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
			if (zone->getKind() == FA_ZONE_BOUNDARY)
			{
				if (zoneIsWall(zone->getName()))
					for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
					{
						rhouu_fa[ifa] = 0.0;
						rhouv_fa[ifa] = 0.0;
						rhouw_fa[ifa] = 0.0;
						rhovv_fa[ifa] = 0.0;
						rhovw_fa[ifa] = 0.0;
						rhoww_fa[ifa] = 0.0;
					}
				else
					for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
					{
						int icv0 = cvofa[ifa][0];

						rhouu_fa[ifa] = rhouu[icv0];
						rhouv_fa[ifa] = rhouv[icv0];
						rhouw_fa[ifa] = rhouw[icv0];
						rhovv_fa[ifa] = rhovv[icv0];
						rhovw_fa[ifa] = rhovw[icv0];
						rhoww_fa[ifa] = rhoww[icv0];
					}
			}
		// ------------------------------------------
		// calculate momentum source terms
		// ------------------------------------------
		static double *drhouu_dx = new double [ncv_g];  // d(ruu)/dx
		static double *drhouv_dy = new double [ncv_g];  // d(ruv)/dy
		static double *drhouw_dz = new double [ncv_g];  // d(ruw)/dz
		static double *drhovu_dx = new double [ncv_g];  // d(rvu)/dx
		static double *drhovv_dy = new double [ncv_g];  // d(rvv)/dy
		static double *drhovw_dz = new double [ncv_g];  // d(rvw)/dz
		static double *drhowu_dx = new double [ncv_g];  // d(rwu)/dx
		static double *drhowv_dy = new double [ncv_g];  // d(rwv)/dy
		static double *drhoww_dz = new double [ncv_g];  // d(rww)/dz


		static double (*tmp)[3] = new double [ncv_g][3];

		calcCvScalarGrad(tmp, rhouu, rhouu_fa, gradreconstruction, limiterNavierS, kine, epsilonSDWLS);  // calc grad
		for (int icv = 0; icv < ncv; icv++)   drhouu_dx[icv] = tmp[icv][0];   // save only d(ruu)/dx
		updateCvData(drhouu_dx, REPLACE_DATA);

		calcCvScalarGrad(tmp, rhovv, rhovv_fa, gradreconstruction, limiterNavierS, kine, epsilonSDWLS);  // calc grad
		for (int icv = 0; icv < ncv; icv++)   drhovv_dy[icv] = tmp[icv][1];   // save only d(rvv)/dy
		updateCvData(drhovv_dy, REPLACE_DATA);

		calcCvScalarGrad(tmp, rhoww, rhoww_fa, gradreconstruction, limiterNavierS, kine, epsilonSDWLS);  // calc grad
		for (int icv = 0; icv < ncv; icv++)   drhoww_dz[icv] = tmp[icv][2];   // save only d(rww)/dz
		updateCvData(drhoww_dz, REPLACE_DATA);

		calcCvScalarGrad(tmp, rhouv, rhouv_fa, gradreconstruction, limiterNavierS, kine, epsilonSDWLS);  // calc grad
		for (int icv = 0; icv < ncv; icv++)
		{
			drhouv_dy[icv] = tmp[icv][1];   // save only d(ruv)/dy
			drhovu_dx[icv] = tmp[icv][0];   // save only d(rvu)/dx
		}
		updateCvData(drhouv_dy, REPLACE_DATA);
		updateCvData(drhovu_dx, REPLACE_DATA);

		calcCvScalarGrad(tmp, rhouw, rhouw_fa, gradreconstruction, limiterNavierS, kine, epsilonSDWLS);  // calc grad
		for (int icv = 0; icv < ncv; icv++)
		{
			drhouw_dz[icv] = tmp[icv][2];   // save only d(ruw)/dz
			drhowu_dx[icv] = tmp[icv][0];   // save only d(rwu)/dx
		}
		updateCvData(drhouw_dz, REPLACE_DATA);
		updateCvData(drhowu_dx, REPLACE_DATA);

		calcCvScalarGrad(tmp, rhovw, rhovw_fa, gradreconstruction, limiterNavierS, kine, epsilonSDWLS);  // calc grad
		for (int icv = 0; icv < ncv; icv++)
		{
			drhovw_dz[icv] = tmp[icv][2];   // save only d(rvw)/dz
			drhowv_dy[icv] = tmp[icv][1];   // save only d(rwv)/dy
		}
		updateCvData(drhovw_dz, REPLACE_DATA);
		updateCvData(drhowv_dy, REPLACE_DATA);

		// -------------------------------------------
		// enhance stability by filtering source terms - eventually phase out
		// -------------------------------------------
		int filter_iters = getIntParam("EUQ_FILTER_ITERS", "10");

		for (int f=0; f<filter_iters; f++) // filter N=filter_iters times
		{
			filterFieldVar(drhouu_dx);
			filterFieldVar(drhouv_dy);
			filterFieldVar(drhouw_dz);

			filterFieldVar(drhovu_dx);
			filterFieldVar(drhovv_dy);
			filterFieldVar(drhovw_dz);

			filterFieldVar(drhowu_dx);
			filterFieldVar(drhowv_dy);
			filterFieldVar(drhoww_dz);
		}


		for (int icv = 0; icv < ncv; icv++)
		{
			// add to momentum source terms
			srcMomentum[icv][0] = drhouu_dx[icv] + drhouv_dy[icv] + drhouw_dz[icv];
			srcMomentum[icv][1] = drhovu_dx[icv] + drhovv_dy[icv] + drhovw_dz[icv];
			srcMomentum[icv][2] = drhowu_dx[icv] + drhowv_dy[icv] + drhoww_dz[icv];

			// add to energy source terms

		} // end icv loop
	} // end perturb stress routine

	/**
	 *
	 * Perturbations of the turbulence anisotropy eigenvalues for EUQ
	 *
	 */
	void perturbEvals(int icv, double c1c, double c2c,double c3c,double eigs_new[3][3])
	{
		// read in eigenvalue perturbation parameters
		string eval_marker     = getParam("EUQ_EVAL")->getString("MARKER");
		string eval_interp     = getParam("EUQ_EVAL")->getString("INTERP");
		double eval_cut        = getParam("EUQ_EVAL")->getDouble("CUTOFF");
		double eval_delta_max  = getParam("EUQ_EVAL")->getDouble("DELTA_MAX");

		// compute marker for cell
		if (eval_marker == "M")
		{
			double s_vec[3], g_vec[3], vel_mag, g_mag, f_mark;

			// compute (u_k u_k)^0.5
			vel_mag = sqrt(vel[icv][0]*vel[icv][0] + vel[icv][1]*vel[icv][1] + vel[icv][2]*vel[icv][2]);

			// compute u_i / (u_k u_k)^0.5
			s_vec[0] = vel[icv][0] / vel_mag;
			s_vec[1] = vel[icv][1] / vel_mag;
			s_vec[2] = vel[icv][2] / vel_mag;

			// compute g_j vector
			for (int j=0; j < 3; j++)
			{
				g_vec[j] = 0.0;
				for (int i=0; i < 3; i++)
				{
					g_vec[j] += s_vec[i]*grad_u[icv][i][j];
				}

			}

			g_mag = sqrt(g_vec[0]*g_vec[0] + g_vec[1]*g_vec[1] + g_vec[2]*g_vec[2]);


			// compute f = |g_j s_j| / (g_k g_k)^0.5
			f_mark = 0.0;
			for (int i=0; i < 3; i++)
				f_mark += fabs(g_vec[i]*s_vec[i]);

			f_mark /= g_mag;

			// marker definition
			eval_sensor[icv] = f_mark * kine[icv] / (vel_mag * vel_mag);
		}
		else if (eval_marker == "ZERO")
		{
			eval_sensor[icv] = 0.0;
		}
		else {
			if (mpi_rank == 0)
				cerr << "ERROR: EUQ_EVAL MARKER (" << eval_marker <<
				") not recognized, choose MARKER=<M, ZERO>" << endl;
			throw(-1);
		}
		//####################################################################################
		// Eigenvalue perturbation
		//####################################################################################

		if ( eval_interp == "LINEAR") {

			if ( eval_sensor[icv] > eval_cut )
			{
				d_bary[icv] = eval_delta_max;
			}
			else
			{
				d_bary[icv] = (eval_sensor[icv] / eval_cut) * eval_delta_max;
			}

		}
		else if ( eval_interp == "COSINE") {

			if ( eval_sensor[icv] > eval_cut )
			{
				d_bary[icv] = eval_delta_max;
			}
			else
			{
				d_bary[icv] = eval_delta_max * (0.5 - 0.5 * cos(3.14159265 * eval_sensor[icv] / eval_cut)) ;
			}

		}
		// implement hypertangent here too
		else {
			if (mpi_rank == 0)
				cerr << "ERROR: EUQ_EVAL PERTURB (" << eval_interp <<
				") not recognized, choose PERTURB=<LINEAR, COSINE>" << endl;
			throw(-1);
		}

		int eval_corner = getParam("EUQ_EVAL")->getInt("CORNER") - 1;

		// Update barycentric coordinates to perturbed state only if mode is active, otherwise do no changes
		if (getParam("EUQ_EVAL")->getString("MODE") == "ACTIVE")
		{
			baryCo[icv][1] = baryCo[icv][1] + d_bary[icv] * (corner[eval_corner][1] - baryCo[icv][1]);
			baryCo[icv][0] = baryCo[icv][0] + d_bary[icv] * (corner[eval_corner][0] - baryCo[icv][0]);
		}
		// rebuild c1c,c2c,c3c based on new barycentric coordinates
		c3c = baryCo[icv][1] / corner[2][1];
		c1c = baryCo[icv][0] - corner[2][0] * c3c;
		c2c = 1 - c1c - c3c;


		// Build new eigenvalues from new coordinates
		eigs_new[0][0] = (c3c - 1) / 3.0;
		eigs_new[1][1] = 0.5 *c2c + eigs_new[0][0];
		eigs_new[2][2] = c1c + eigs_new[1][1];
	}


	/**
	 *
	 * Perturbations of the turbulence anisotropy eigenvectors for EUQ
	 *
	 */
	void perturbEvectors(int icv)
	{}

	/**
	 *
	 * Perturbations of the turbulence kinetic energy for EUQ
	 *
	 */
	void perturbK(int icv)
	{}

	/**
	 *
	 * Force barycentric map location to realizable states of turbulence (not completed yet)
	 *
	 */
	void forceRealizability(int icv,double eigs_new[3][3])
	{
		// first check if cell location is realizable
		double x0 = baryCo[icv][0];
		double y0 = baryCo[icv][1];

		// if NOT realizable
		if ( ( x0 <= 0.5 && ( y0 <= 0.0 || y0 >= 2.0 * corner[2][1] * x0 )) || ( x0 > 0.5 && ( y0 <= 0.0 || y0 >= -2.0 * corner[2][1] * x0 + 2.0 * corner[2][1] )) )
		{
			/*			// find nearest neighbor in map (location on map boundary)

			// Line C2c-C1c (bottom boundary)
			// perpendicular distance is just y0 coordinate
			double d_21 = fabs(y0);

			double x1 = corner[1][0];
			double y1 = corner[1][1];
			double x2 = corner[0][0];
			double y2 = corner[0][1];

			double slope_p = -(x2-x1)/(y2-y1);
			 */
			realizable[icv] = -1;
		}
		else
		{
			realizable[icv] = 0;
		}

	}

	/**
	 *
	 * Filter reynolds stress terms for stibility purposes
	 *
	 */
	void filterFieldVar(double *var)
	{
		/* Rene's filtering routine
		 * is this necessary given that we use smoothing on the sensor?
		 */

		static double *tmpVar = new double [ncv_g];
		for (int icv = 0; icv < ncv_g; icv++)
			tmpVar[icv] = var[icv];

		for (int icv = 0; icv < ncv; icv++)
		{
			int noc_f = nbocv_i[icv];
			int noc_l = nbocv_i[icv + 1] - 1;

			int count = 1;
			for (int noc = noc_f + 1; noc <= noc_l; noc++)
			{
				int icv_nbr = nbocv_v[noc];
				var[icv] += tmpVar[icv_nbr];
				count++;
			}
			var[icv] /= count;
		}
		updateCvData(var, REPLACE_DATA);
	}


	/**
	 *
	 * Multiply two 3x3 matrices
	 *
	 */
	void matMult(double (*result)[3], double (*m1)[3], double (*m2)[3])  // m1*m2 matrices
	{
		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
			{
				result[i][j] = 0.0;
				for (int k=0; k<3; k++)
					result[i][j] += m1[i][k] * m2[k][j];
			}
	}


	//these routiens are for computing the eigenvalue decomposition
	static double hypot2(double x, double y) {
		return sqrt(x*x+y*y);
	}

	/**
	 *
	 * Symmetric Tridiagonal QL algorithm: tred2
	 *
	 */
	static void tred2(double V[3][3], double d[3], double e[3]) {

		//  This is derived from the Algol procedures tred2 by
		//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
		//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
		//  Fortran subroutine in EISPACK.

		for (int j = 0; j < 3; j++) {
			d[j] = V[3-1][j];
		}

		// Householder reduction to tridiagonal form.

		for (int i = 3-1; i > 0; i--) {

			// Scale to avoid under/overflow.

			double scale = 0.0;
			double h = 0.0;
			for (int k = 0; k < i; k++) {
				scale = scale + fabs(d[k]);
			}
			if (scale == 0.0) {
				e[i] = d[i-1];
				for (int j = 0; j < i; j++) {
					d[j] = V[i-1][j];
					V[i][j] = 0.0;
					V[j][i] = 0.0;
				}
			} else {

				// Generate Householder vector.

				for (int k = 0; k < i; k++) {
					d[k] /= scale;
					h += d[k] * d[k];
				}
				double f = d[i-1];
				double g = sqrt(h);
				if (f > 0) {
					g = -g;
				}
				e[i] = scale * g;
				h = h - f * g;
				d[i-1] = f - g;
				for (int j = 0; j < i; j++) {
					e[j] = 0.0;
				}

				// Apply similarity transformation to remaining columns.

				for (int j = 0; j < i; j++) {
					f = d[j];
					V[j][i] = f;
					g = e[j] + V[j][j] * f;
					for (int k = j+1; k <= i-1; k++) {
						g += V[k][j] * d[k];
						e[k] += V[k][j] * f;
					}
					e[j] = g;
				}
				f = 0.0;
				for (int j = 0; j < i; j++) {
					e[j] /= h;
					f += e[j] * d[j];
				}
				double hh = f / (h + h);
				for (int j = 0; j < i; j++) {
					e[j] -= hh * d[j];
				}
				for (int j = 0; j < i; j++) {
					f = d[j];
					g = e[j];
					for (int k = j; k <= i-1; k++) {
						V[k][j] -= (f * e[k] + g * d[k]);
					}
					d[j] = V[i-1][j];
					V[i][j] = 0.0;
				}
			}
			d[i] = h;
		}

		// Accumulate transformations.

		for (int i = 0; i < 3-1; i++) {
			V[3-1][i] = V[i][i];
			V[i][i] = 1.0;
			double h = d[i+1];
			if (h != 0.0) {
				for (int k = 0; k <= i; k++) {
					d[k] = V[k][i+1] / h;
				}
				for (int j = 0; j <= i; j++) {
					double g = 0.0;
					for (int k = 0; k <= i; k++) {
						g += V[k][i+1] * V[k][j];
					}
					for (int k = 0; k <= i; k++) {
						V[k][j] -= g * d[k];
					}
				}
			}
			for (int k = 0; k <= i; k++) {
				V[k][i+1] = 0.0;
			}
		}
		for (int j = 0; j < 3; j++) {
			d[j] = V[3-1][j];
			V[3-1][j] = 0.0;
		}
		V[3-1][3-1] = 1.0;
		e[0] = 0.0;
	}
	// Symmetric tridiagonal QL algorithm.

	/**
	 *
	 * Symmetric Tridiagonal QL algorithm: tql2
	 *
	 */
	static void tql2(double V[3][3], double d[3], double e[3]) {

		//  This is derived from the Algol procedures tql2, by
		//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
		//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
		//  Fortran subroutine in EISPACK.

		for (int i = 1; i < 3; i++) {
			e[i-1] = e[i];
		}
		e[3-1] = 0.0;

		double f = 0.0;
		double tst1 = 0.0;
		double eps = pow(2.0,-52.0);
		for (int l = 0; l < 3; l++) {

			// Find small subdiagonal element

			tst1 = max(tst1,(fabs(d[l]) + fabs(e[l])));
			int m = l;
			while (m < 3) {
				if (fabs(e[m]) <= eps*tst1) {
					break;
				}
				m++;
			}

			// If m == l, d[l] is an eigenvalue,
			// otherwise, iterate.

			if (m > l) {
				int iter = 0;
				do {
					iter = iter + 1;  // (Could check iteration count here.)

					// Compute implicit shift

					double g = d[l];
					double p = (d[l+1] - g) / (2.0 * e[l]);
					double r = hypot2(p,1.0);
					if (p < 0) {
						r = -r;
					}
					d[l] = e[l] / (p + r);
					d[l+1] = e[l] * (p + r);
					double dl1 = d[l+1];
					double h = g - d[l];
					for (int i = l+2; i < 3; i++) {
						d[i] -= h;
					}
					f = f + h;

					// Implicit QL transformation.

					p = d[m];
					double c = 1.0;
					double c2 = c;
					double c3 = c;
					double el1 = e[l+1];
					double s = 0.0;
					double s2 = 0.0;
					for (int i = m-1; i >= l; i--) {
						c3 = c2;
						c2 = c;
						s2 = s;
						g = c * e[i];
						h = c * p;
						r = hypot2(p,e[i]);
						e[i+1] = s * r;
						s = e[i] / r;
						c = p / r;
						p = c * d[i] - s * g;
						d[i+1] = h + s * (c * g + s * d[i]);

						// Accumulate transformation.

						for (int k = 0; k < 3; k++) {
							h = V[k][i+1];
							V[k][i+1] = s * V[k][i] + c * h;
							V[k][i] = c * V[k][i] - s * h;
						}
					}
					p = -s * s2 * c3 * el1 * e[l] / dl1;
					e[l] = s * p;
					d[l] = c * p;

					// Check for convergence.

				} while (fabs(e[l]) > eps*tst1);
			}
			d[l] = d[l] + f;
			e[l] = 0.0;
		}

		// Sort eigenvalues and corresponding vectors.

		for (int i = 0; i < 3-1; i++) {
			int k = i;
			double p = d[i];
			for (int j = i+1; j < 3; j++) {
				if (d[j] < p) {
					k = j;
					p = d[j];
				}
			}
			if (k != i) {
				d[k] = d[i];
				d[i] = p;
				for (int j = 0; j < 3; j++) {
					p = V[j][i];
					V[j][i] = V[j][k];
					V[j][k] = p;
				}
			}
		}
	}

	/**
	 *
	 * Compute eigenvalue decomposition of a 3x3 tensor into
	 * a 3x3 eigenvector matrix
	 * and a 1x3 vector of the eigenvalues
	 *
	 */
	void eigen_decomposition(double A[3][3], double V[3][3], double d[3]) {
		double e[3];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				V[i][j] = A[i][j];
			}
		}
		tred2(V, d, e);
		tql2(V, d, e);
	}


};



#endif


