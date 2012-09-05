/**
 * \brief Combustion classes for thermodynamics properties calculation of species and mixture.
 *
 * Contains functions to load species parameters, compute thermodynamics quantities 
 * (e.g., mu, lambda, cp, gamma, enthalpy, temperature) as a function of temperature and mixture composition.
 * Based on FlameMaster and adapted.
 * 
 * \author Vincent Terrapon 
 * \date April 2010
 * \version 2.0
 */

#ifndef COMBUSTION_H
#define COMBUSTION_H

//#define mixing_tabulated

#include "CombustionGeneralDefinitions.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <string>
using std::string;
#include <map>
using std::pair;
using std::map;
#include "Param.h"



///////////////////////////////////////////////////////////////////////////////////////////////////////////

/** 
 * Overall constants
 */
#define kR            8314.4            ///< Universal gas constant in [J/(kmol*K)].
#define kTreference   298.15             ///< Reference temperature for which the enthalpy is zero: H(Tref) = 0.
#define kCoefLen      7                  ///< Array size of NASA coefficients for calculation of cp and h.
#define AWT_H         1.008              ///< Atomic weights of H  atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
#define AWT_HE        4.00               ///< Atomic weights of HE atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
#define AWT_C         12.01              ///< Atomic weights of C  atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
#define AWT_N         14.01             ///< Atomic weights of N  atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
#define AWT_O         16.00              ///< Atomic weights of O  atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
#define AWT_F         18.99840           ///< Atomic weights of F  atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
#define AWT_CL        35.453             ///< Atomic weights of CL atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
#define AWT_AR        39.948             ///< Atomic weights of AR atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
#define AWT_BR        79.909             ///< Atomic weights of BR atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
//#define AWT_H         1.00794            ///< Atomic weights of H  atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
//#define AWT_HE        4.00260            ///< Atomic weights of HE atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
//#define AWT_C         12.01070           ///< Atomic weights of C  atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
//#define AWT_N         14.00674           ///< Atomic weights of N  atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
//#define AWT_O         15.99940           ///< Atomic weights of O  atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
//#define AWT_F         18.99840           ///< Atomic weights of F  atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
//#define AWT_CL        35.45270           ///< Atomic weights of CL atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
//#define AWT_AR        39.94800           ///< Atomic weights of AR atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
//#define AWT_BR        79.90400           ///< Atomic weights of BR atom in [kg/(kmol)] used in Species::ComputeSpeciesComposition().
#define om_mu_m1      3.3530622607       ///< Coefficients to compute Species::omega_mu().
#define om_mu_m2      2.53272006         ///< Coefficients to compute Species::omega_mu().
#define om_mu_m3      2.9024238575       ///< Coefficients to compute Species::omega_mu().
#define om_mu_m4      0.11186138893      ///< Coefficients to compute Species::omega_mu().
#define om_mu_m5      0.8662326188       ///< Coefficients to compute Species::omega_mu().
#define om_mu_m6      1.3913958626       ///< Coefficients to compute Species::omega_mu().
#define om_mu_m7      3.158490576        ///< Coefficients to compute Species::omega_mu().
#define om_mu_m8      0.18973411754      ///< Coefficients to compute Species::omega_mu().
#define om_mu_m9      0.00018682962894   ///< Coefficients to compute Species::omega_mu().
#define TMaxClip      5000.0             ///< Maximum temperature at which it is clipped.
#define TMinClip      10.0               ///< Minimum temperature at which it is clipped.
#define TPolyLinMin   300.0              ///< Temperature under which Cp is assumed constant and NASA polynomials not used anymore.
#define TPolyLinMax   3000.0             ///< Temperature above which Cp is assumed constant and NASA polynomials not used anymore.
#define NewtonIterMax 200                ///< Maximum number of iterations in Newton solver for temperature.
#define NewtonEpsilon 1.0e-8             ///< Error tolerance in Newton solver for temperature.
///////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Classes definition */
/**********************/

/*! \brief Contains the composition of different species.
 *
 *  Each variable contains the number of corresponding atoms in the given species.
 *  \warning Only 9 different atoms are defined: C, H, N, O, BR, F, AR, CL, HE.
 */
typedef struct
{
  unsigned char C, H, N, O, BR, F, AR, CL, HE;
} compo;

/***********************************************************************************************************/
/***********************************************************************************************************/

/*! \brief Class containing the quantities associated with a given species and functions to compute the
 *  species properties.
 *
 *  Compute species constants and species properties (e.g., Yk, mu, lambda, cp, h, atomic composition, M).
 */
class Species
{
private:
  /* variables */
  char   name[kNameLen];  ///< Name of species (read from input file).
  double M;               ///< Molecular weight in [kg/kmol]  (computed from atomic composition).
  compo  compo_elems;     ///< Atomic composition (computed from species name).
  double Tmed;            ///< Temperature switch to determine which NASA coefficient to use: cold/hot (read from input file).
  double hot[kCoefLen];   ///< NASA coefficients for T>Tmed [K] (read from input file).
  double cold[kCoefLen];  ///< NASA coefficients for T<Tmed [K] (read from input file).
  int    shape;           ///< Shape of molecule: 0 (atom), 1 (linear), 2 (non linear) (read from input file). NOT USED HERE!
  double k_over_eps;      ///< Temperature factor for evaluating omega_mu. k/eps is in [K], where eps is the characteristic energy and k is the Boltzmann constant (read from input file).
  double sigma;           ///< Collision diameter for evaluation of viscosity mu (read from input file).
  double dipole;          ///< Dipole moment (read from input file). NOT USED HERE!
  double alpha;           ///< Polarizability (read from input file). NOT USED HERE!
  double Zrot;            ///< Rotational relaxation number (read from input file). NOT USED HERE!
  double mu_coeff;        ///< Constant factor for evaluation of mu in [kg/(m*s)] (computed).
  double Yk;              ///< Species mass fraction in [-] (computed from chemistry table or from linear mixing rules).
  double Yk_bc_left;      ///< Species mass fraction at left boundary (Z=0) (used only for mixing)
  double Yk_bc_right;     ///< Species mass fraction at right boundary (Z=1) (used only for mixing)
  double mu;              ///< Species viscosity in [kg/(m*s)] (computed).
  double cp;              ///< Species heat capacity at constant pressure in [J/(kg*K)] (computed with NASA coefficients).
  double lambda;          ///< Species heat diffusion coefficient in [W/(m*K)] (computed).
  double delta;           ///< Working array to compute mu and lambda of mixture
  double Lewis;           ///< Lewis number = lambda/(rho*cp*Dk) = Dth/Dk in [-]. \warning Set to unity!
  double rhoDk;           ///< Diffusion coefficient of species multiplied by density rho * Dk in [kg/m/s].


  /***********************************************************************************************************/
  /* private functions */

  /*! \brief Compute atomic composition of species based on species name.
   *
   * This functions parses the name string of the species and counts the number of C, H, N, O, BR, F, AR, CL and HE atoms. 
   * From FlameMaster but modified.
   * The algorithm currently works only for names that can be described by the following productions: \n
   * name  -> list \n
   * list  -> list desc | desc \n
   * desc  -> id num | id \n
   * id   -> C | H | N | O | F | BR | AR | CL | HE \n
   * num    -> num digit | digit \n
   * digit -> 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 \n
   * \warning Only C, H, N, O, BR, F, AR, CL and HE atoms are considered! Should also be checked for more complex
   * species name (e.g., isomeres).
   */
  void ComputeSpeciesComposition()
  {
    compo *cPtr = &compo_elems;
    char theChar;
    char *ptr = name;
    char *isoPtr = 0;
    int num_of_atoms;
    bool err = false;

    /* initialize everything to zero */
    cPtr->C = (unsigned char) 0, cPtr->H = (unsigned char) 0, 
    cPtr->N = (unsigned char) 0, cPtr->O = (unsigned char) 0, 
    cPtr->F = (unsigned char) 0, cPtr->BR = (unsigned char) 0, 
    cPtr->AR = (unsigned char) 0, cPtr->HE = (unsigned char) 0, 
    cPtr->CL = (unsigned char) 0;

    if (isoPtr = strrchr(ptr, '-'))
    {
      ptr = &isoPtr[1];
    }
    while (*ptr)
    {

      theChar = *ptr++; /* save this character and move to the next one */

      /* if theChar is a 'A', it must be a AR atom or a benzene group (A1=C6H5)
         (at least in this implementation), but let's take a look at the next letter */
      if (theChar == 'A')
      {
        /* if it is 'R', goto next character */
        if (*ptr == 'R')
        {
          ptr++;            
        }
        /* if it is a digit, it is a benzene group, save it as 'X' */
        else if (isdigit(*ptr))
        {
          theChar = 'X'; 
          ptr++;
        }
        /* otherwise this element is not known */
        else
        {
          cPtr->C = cPtr->H = cPtr->N = cPtr->O = cPtr->F = cPtr->BR
          = cPtr->AR = cPtr->HE = cPtr->CL = 0;
          cerr << "### 1 Unknown atom in composition " << theChar << *ptr
          << " ###" << endl;
          throw(-1);
        }
      }
      /* if theChar is a 'B', it must be a BR atom
         (at least in this implementation), but let's take a look at the next letter */
      if (theChar == 'B')
      {
        /* if next letter is not 'R', this element is not known */
        if (*ptr != 'R')
        {
          cPtr->C = cPtr->H = cPtr->N = cPtr->O = cPtr->F = cPtr->BR
          = cPtr->AR = cPtr->HE = cPtr->CL = 0;
          cerr << "### 2 Unknown atom in composition " << theChar << *ptr
          << " ###" << endl;
          throw(-1);
        }
        else
          /* if it is 'R', goto next character */
          ptr++;
      }
      /* if it is 'H', species can be H or HE: check next character */
      if (theChar == 'H')
      {
        /* if next letter is 'E', goto next character */
        if (*ptr == 'E')
        {
          theChar = 'E'; /* Because H is already used in switch below */
          ptr++;
        }
        /* if next letter is not 'E' or one of known atoms, this element is not known */
        else if ((*ptr != 'C') && (*ptr != 'N') && (*ptr != 'F') && (*ptr
            != 'A') && (*ptr != 'H') && (*ptr != 'O') && (*ptr != 'B')
            && (*ptr != '2') && (*ptr != '3') && (*ptr != '4') && (*ptr
                != '5') && (*ptr != '6') && (*ptr != '7') && (*ptr != '8')
                && (*ptr != '9') && (*ptr != ' ') && (*ptr != '\0'))
        {
          cPtr->C = cPtr->H = cPtr->N = cPtr->O = cPtr->F = cPtr->BR
          = cPtr->AR = cPtr->HE = cPtr->CL = 0;
          cerr << "### 3 Unknown atom in composition " << theChar << *ptr
          << " ###" << endl;
          throw(-1);
        }
      }
      /* if it is 'C', species can be C or CL: check next character */
      if (theChar == 'C')
      {
        /* if next letter is 'L', goto next character */
        if (*ptr == 'L')
        {
          theChar = 'L'; /* Because C is already used in switch below */
          ptr++;
        }
        /* if next letter is not 'L' or one of known atoms, this element is not known */
        else if ((*ptr != 'C') && (*ptr != 'N') && (*ptr != 'F') && (*ptr
            != 'A') && (*ptr != 'H') && (*ptr != 'O') && (*ptr != 'B')
            && (*ptr != '2') && (*ptr != '3') && (*ptr != '4') && (*ptr
                != '5') && (*ptr != '6') && (*ptr != '7') && (*ptr != '8')
                && (*ptr != '9') && (*ptr != ' ') && (*ptr != '\0'))
        {
          cPtr->C = cPtr->H = cPtr->N = cPtr->O = cPtr->F = cPtr->BR
          = cPtr->AR = cPtr->HE = cPtr->CL = 0;
          cerr << "### 4 Unknown atom in composition " << theChar << *ptr
          << " ###" << endl;
          throw(-1);
        }
      }

      /* Get the number of atoms */
      num_of_atoms = 0;
      while (*ptr && isdigit(*ptr))
      {
        num_of_atoms *= 10;
        num_of_atoms += *ptr++ - '0';
      }
      if (!num_of_atoms)
        num_of_atoms = 1; /* there is always at least one */

      switch (toupper(theChar))
      {
      case 'C':
        cPtr->C += num_of_atoms;
        break;
      case 'H':
        cPtr->H += num_of_atoms;
        break;
      case 'N':
        cPtr->N += num_of_atoms;
        break;
      case 'O':
        cPtr->O += num_of_atoms;
        break;
      case 'F':
        cPtr->F += num_of_atoms;
        break;
      case 'B':
        cPtr->BR += num_of_atoms;
        break;
      case 'A':
        cPtr->AR += num_of_atoms;
        break;
      case 'E':
        cPtr->HE += num_of_atoms;
        break;
      case 'L':
        cPtr->CL += num_of_atoms;
        break;
      case 'X': /* benzene group */
        cPtr->C += num_of_atoms*6;
        cPtr->H += num_of_atoms*5;
        break;
      default:
        cPtr->C = cPtr->H = cPtr->N = cPtr->O = cPtr->F = cPtr->BR = cPtr->AR
        = cPtr->HE = cPtr->CL = 0;
        err = true;
        break; /* don't do anything right now */
      } /* switch */
      if (err)
      {
        cerr << "### Error during reading of composition, current char is "
        << theChar << " ###" << endl;
        throw(-1);
      }
    } /* while */
    return;
  }

  /***********************************************************************************************************/

  /*! \brief Compute species molar mass M based on composition.
   *
   *  Compute species molar mass M by multiplying the number of atoms present in species with their corresponding atomic mass.
   */
  void ComputeSpeciesMolarMass()
  {
    M = (double) compo_elems.C * AWT_C + (double) compo_elems.H * AWT_H
    + (double) compo_elems.N * AWT_N + (double) compo_elems.O * AWT_O
    + (double) compo_elems.F * AWT_F + (double) compo_elems.BR * AWT_BR
    + (double) compo_elems.CL * AWT_CL + (double) compo_elems.AR * AWT_AR
    + (double) compo_elems.HE * AWT_HE;
    return;
  }

  /***********************************************************************************************************/

  /*! \brief Compute omega_mu as function of dimensionless temperature.
   *
   *  The collision integral for viscosity, omega_mu, is used to compute the viscosity. omega_mu is computed
   *  based on the coefficients #om_mu_m1 to #om_mu_m9.
   *  \param[in] temp Dimensionless temperature T/(eps/k).
   *  \return omega_mu for a given dimensionless temperature in [-].
   */
  double omega_mu(double &temp)
  {
    double num, den;

    num = om_mu_m1 + temp * (om_mu_m2 + temp * (om_mu_m3 + temp * om_mu_m4));
    den = om_mu_m5 + temp * (om_mu_m6 + temp * (om_mu_m7 + temp * (om_mu_m8
        + temp * om_mu_m9)));
    return num / den;
  }

  /***********************************************************************************************************/

public:
  /* Accessors */
  /*! \brief Accessor to species #name. */
  void GetSpeciesName(char * buffer)
  {
    strcpy(buffer, name);
  }
  /*! \brief Accessor to species #name. */
  void SetSpeciesName(char * buffer)
  {
    strcpy(name, buffer);
  }
  /*! \brief Accessor to species molar mass #M. */
  double GetSpeciesMolarMass() const
  {
    return M;
  }
  /*! \brief Accessor to species molar switch temperature #Tmed. */
  double GetSpeciesTmed() const
  {
    return Tmed;
  }
  /*! \brief Accessor to species molar switch temperature #Tmed. */
  void SetSpeciesTmed(double val)
  {
    Tmed = val;
  }
  /*! \brief Accessor to species NASA coefficients #hot. */
  void GetSpeciesHot(double * val)
  {
    for (int i = 0; i < kCoefLen; i++)
      val[i] = hot[i];
  }
  /*! \brief Accessor to species NASA coefficients #hot. */
  void SetSpeciesHot(double * val)
  {
    for (int i = 0; i < kCoefLen; i++)
      hot[i] = val[i];
  }
  /*! \brief Accessor to species NASA coefficients #cold. */
  void GetSpeciesCold(double * val)
  {
    for (int i = 0; i < kCoefLen; i++)
      val[i] = cold[i];
  }
  /*! \brief Accessor to species NASA coefficients #cold. */
  void SetSpeciesCold(double * val)
  {
    for (int i = 0; i < kCoefLen; i++)
      cold[i] = val[i];
  }
  /*! \brief Accessor to species #shape. */
  int GetSpeciesShape() const
  {
    return shape;
  }
  /*! \brief Accessor to species #shape. */
  void SetSpeciesShape(int val)
  {
    shape = val;
  }
  /*! \brief Accessor to species k/eps #k_over_eps. */
  double GetSpeciesKeps() const
  {
    return k_over_eps;
  }
  /*! \brief Accessor to species k/eps #k_over_eps. */
  void SetSpeciesKeps(double val)
  {
    k_over_eps = val;
  }
  /*! \brief Accessor to species collision diameter #sigma. */
  double GetSpeciesSigma() const
  {
    return sigma;
  }
  /*! \brief Accessor to species collision diameter #sigma. */
  void SetSpeciesSigma(double val)
  {
    sigma = val;
  }
  /*! \brief Accessor to species #dipole. */
  double GetSpeciesDipole() const
  {
    return dipole;
  }
  /*! \brief Accessor to species #dipole. */
  void SetSpeciesDipole(double val)
  {
    dipole = val;
  }
  /*! \brief Accessor to species polarizability #alpha. */
  double GetSpeciesAlpha() const
  {
    return alpha;
  }
  /*! \brief Accessor to species polarizability #alpha. */
  void SetSpeciesAlpha(double val)
  {
    alpha = val;
  }
  /*! \brief Accessor to species rotational relaxation number #Zrot. */
  double GetSpeciesZrot() const
  {
    return Zrot;
  }
  /*! \brief Accessor to species rotational relaxation number #Zrot. */
  void SetSpeciesZrot(double val)
  {
    Zrot = val;
  }
  /*! \brief Accessor to species viscosity constant #mu_coeff. */
  double GetSpeciesMucoeff() const
  {
    return mu_coeff;
  }
  /*! \brief Accessor to species mass fraction #Yk. */
  double GetYk() const
  {
    return Yk;
  }
  /*! \brief Accessor to species mass fraction #Yk. */
  void SetYk(double y)
  {
    Yk = y;
  }
  /*! \brief Accessor to species mass fraction at left boundary (Z=0) #Yk_bc_left. */
  double GetYk_bc_left() const
  {
    return Yk_bc_left;
  }
  /*! \brief Accessor to species mass fraction at left boundary (Z=0) #Yk_bc_left. */
  void SetYk_bc_left(double y)
  {
    Yk_bc_left = y;
  }
  /*! \brief Accessor to species mass fraction at right boundary (Z=1) #Yk_bc_right. */
  double GetYk_bc_right() const
  {
    return Yk_bc_right;
  }
  /*! \brief Accessor to species mass fraction at right boundary (Z=1) #Yk_bc_right. */
  void SetYk_bc_right(double y)
  {
    Yk_bc_right = y;
  }
  /*! \brief Accessor to species viscosity #mu. */
  double GetMu() const
  {
    return mu;
  }
  /*! \brief Accessor to species heat capacity #cp. */
  double GetCp() const
  {
    return cp;
  }
  /*! \brief Accessor to species heat diffusion coefficient #lambda. */
  double GetLambda() const
  {
    return lambda;
  }
  /*! \brief Accessor to species working array #delta. */
  double GetDelta() const
  {
    return delta;
  }
  /*! \brief Accessor to species working array #delta. */
  void SetDelta(double x)
  {
    delta = x;
  }
  /*! \brief Accessor to species #Lewis number. */
  double GetLewis() const
  {
    return Lewis;
  }
  /*! \brief Accessor to species #Lewis number. */
  void SetLewis(double y)
  {
    Lewis = y;
  }
  /*! \brief Accessor to species diffusion  coefficient #rhoDk. */
  double GetRhoDk() const
  {
    return rhoDk;
  }
  /*! \brief Accessor to species diffusion coefficicent #rhoDk. */
  void SetRhoDk(double y)
  {
    rhoDk = y;
  }

  /***********************************************************************************************************/
  /* other public functions */

  /*! \brief Initialize species constants.
   *
   *  Computes for the species: the atomic composition #compo_elems, the molar mass #M, the constant part of
   *  the viscosity expression #mu_coeff (µ_coeff = 2.6693e-6 sqrt(M / sigma^2) and set the #Lewis number to 1.
   */
  void ComputeSpeciesConstants()
  {
    ComputeSpeciesComposition();
    ComputeSpeciesMolarMass();
    mu_coeff = 2.6693e-6 * sqrt(M) / (sigma * sigma);
    /* Set Lewis to unity (can be changed later) */
    Lewis = 1.0;
    return;
  }

  /***********************************************************************************************************/

  /*! \brief Compute species enthalpy at given temperature.
   *
   *  Compute species enthalpy based on NASA coefficients. If the temperature is higher than #Tmed, the #hot 
   *  coefficients are used, otherwise the #cold coefficients. R/M is in [J/(kg*K)].
   *  \param[in] temp Temperature in [K].
   *  \return Enthalpy of species in [J/kg] at given temperature.
   */
  double ComputeSpeciesEnthalpy(double temp)
  {
    double *coeff; /* pointer to nasa coefficients */

    if (temp > Tmed)
      coeff = hot;
    else
      coeff = cold;

    return (temp * (coeff[0] + temp * (coeff[1] * 0.5 + temp * (coeff[2] / 3.0
        + temp * (coeff[3] * 0.25 + temp * coeff[4] * 0.2)))) + coeff[5]) * kR
        / M;
  }

  /***********************************************************************************************************/

  /*! \brief Compute species heat capacity at given temperature.
   *
   *  Compute species heat capacity based on NASA coefficients. If the temperature is higher than #Tmed, the #hot 
   *  coefficients are used, otherwise the #cold coefficients. R/M is in [J/(kg*K)].
   *  \param[in] temp Temperature in [K].
   *  \return Heat capacity of species in [J/(kg*K)] at given temperature.
   */
  double ComputeSpeciesCp(double temp)
  {
    double *coeff; /* pointer to nasa coefficients */

    if (temp > Tmed)
      coeff = hot;
    else
      coeff = cold;

    return (coeff[0] + temp * (coeff[1] + temp * (coeff[2] + temp * (coeff[3] + temp * coeff[4])))) * kR / M;
  }

  /***********************************************************************************************************/

  /*! \brief Compute species viscosity at given temperature.
   *
   *  Compute species viscosity (assuming monoatomic gas) using Lennard-Jones parameters: 
   *  #mu = #mu_coeff * sqrt(T) / #omega_mu(T/#k_over_eps). 
   *  See \em Transport \em Phenomena, Bird, Stewart, Lightfoot, 2007, (p.26).
   *  \param[in] temp Temperature in [K].
   *  \return Viscosity of species in [kg/(m*s)] at given temperature.
   */
  double ComputeSpeciesViscosity(double &temp)
  {
    double t = temp * k_over_eps;
    return mu_coeff * sqrt(temp) / omega_mu(t);
  }

  /***********************************************************************************************************/

  /*! \brief Compute species thermal conductivity (heat diffusion coefficient).
   *
   *  Compute species thermal conductivity based on Eucken formula for polyatomic gases: 
   *  #lambda = #mu * (#cp + 1.25 * #kR/#M).
   *  See \em Transport \em Phenomena, Bird, Stewart, Lightfoot, 2007, (p.276).
   *  \return Thermal conductivity of species in [W/(m K)].
   *  \warning Works only if mu and cp have been previously computed!
   */
  double ComputeSpeciesConductivity()
  {
    return mu * (cp + 1.25 * kR / M);
  }

  /***********************************************************************************************************/

  /*! \brief Compute species diffusion coefficient #rhoDk.
   *
   *  Compute species thermal conductivity based on given Lewis number: rhoDk = Lewis * cp / lambda
   *  \return Species diffusion coefficient multiplied by density rho * Dk in [s].
   *  \warning Works only if cp and lambda have been previously computed!
   */
  double ComputeSpeciesDiffusionCoeff()
  {
    return Lewis * cp / lambda;
  }

  /***********************************************************************************************************/

  /*! \brief Compute species properties #mu, #cp, #lambda and #rhoDk at given temperature.
   *
   *  Compute #mu with Species::ComputeSpeciesViscosity(), #cp with Species::ComputeSpeciesCp(),
   *  #lambda with Species::ComputeSpeciesConductivity() and #rhoDk with Species::ComputeSpeciesDiffusionCoeff().
   *  \param[in] temp Temperature in [K].
   */
  void ComputeSpeciesProperties(double &temp)
  {
    mu = ComputeSpeciesViscosity(temp);
    cp = ComputeSpeciesCp(temp);
    lambda = ComputeSpeciesConductivity();
    rhoDk = ComputeSpeciesDiffusionCoeff();
  }

  /***********************************************************************************************************/

  /*! \brief Output species constants and properties on screen.*/
  void OutputSpeciesConstants() const
  {
    int i;

    cout << name << "\t M=" << M << "\t Tmed=" << Tmed << "\t shape=" << shape
    << "\t k_over_eps=" << k_over_eps << "\t sigma=" << sigma
    << "\t mu_coeff=" << mu_coeff << endl;
    for (i = 0; i < kCoefLen; i++)
    {
      cout << "a_h[" << i << "]=" << hot[i] << " \t";
    }
    cout << endl;
    for (i = 0; i < kCoefLen; i++)
    {
      cout << "a_c[" << i << "]=" << cold[i] << " \t";
    }
    cout << endl;
    cout << "C=" << (int) compo_elems.C << "\t H=" << (int) compo_elems.H
    << "\t N=" << (int) compo_elems.N << "\t O=" << (int) compo_elems.O
    << "\t BR=" << (int) compo_elems.BR << "\t F=" << (int) compo_elems.F
    << "\t AR=" << (int) compo_elems.AR << "\t HE=" << (int) compo_elems.HE
    << "\t CL=" << (int) compo_elems.CL << endl;
    cout << "Dipole=" << dipole << "\t alpha=" << alpha << "\t Zrot=" << Zrot
    << endl << endl;
    return;
  }

};

/***********************************************************************************************************/
/***********************************************************************************************************/

/*! \brief Class containing the overall mixture properties and corresponding functions to compute them.
 *
 *  Load species present in mixture, where the mixture composition is read from input file.
 *  Read species constants from thermodynamics input file created with FlameMaster tool \p CreateBinFile. 
 *  2 different versions of thermodynamics input file can be read, controlled by an input flag in Mixture::Load()
 *  read from input file.
 */
class Mixture
{
private:
  /* Variables */
  map<string, Species> mySpecies;                   ///< Map containing the species linking species name (string) to species (Species).
  map<string, Species>::iterator itSp, itSp2;       ///< Iterators to loop over the species.
  pair<map<string, Species>::iterator, bool> retS;  ///< Return value for map insertion to check if species already present.

  int    NkSpecies;                       ///< Number of species in mixture (read from input file).
  char   (*SpeciesName)[kNameLen];        ///< Array containing name of species (read from input file).
  double T;                               ///< Temperature of mixture in [K] (computed).
  double M;                               ///< Molecular mass of mixture in [kg/(kmol)] (computed).
  double R_over_M;                        ///< Mixture gas constant R/M in [J/(kg*K)] (computed).
  double mul;                             ///< Mixture laminar viscosity in [kg/(m*s)] (computed).
  double lambda;                          ///< Mixture thermal conductivity in [W/(m*K)] (computed).
  double gama;                            ///< Heat capacity ratio in [-] (computed).
  double a_sound;                         ///< Speed of sound in mixture in [m/s] (computed).
  double cp;                              ///< Mixture heat capacity in [J/(kg*K)] (computed).
  double H;                               ///< Mixture enthalpy (chemical and sensible) in [J/kg] (computed).
  double Lewis;                           ///< Mixture lewis number (to compute the diffusion coefficient) in [-].
  double rhoDk;                           ///< Mixture diffusion coefficient multiplied by rho in [kg/m/s] (computed assuming Lewis=cst (=lambda/cp/Lewis).
  double pressure;                        ///< Mixture pressure in Pa (computed from equation of state for ideal gas).
  string dominant_species;                ///< Species used to ensure sum of Species::Yk to be 1 by substracting the sum of other species from 1 (read from input file).
  int    CombustionRegime;                ///< Combustion regime: FPVA=2; Steady Flamelet=1;  Mixing=0;  Variable Properties=-1.  
  string OutputToScreen;                  ///< "YES" -> summary of parameters and input data to screen.

public:
  /* accessors */
  /*! \brief Accessor to number of species #NkSpecies in mixture. */
  int GetNkspecies()
  {
    return NkSpecies;
  }
  /*! \brief Accessor to number of species #NkSpecies in mixture. */
  void SetNkspecies(int i)
  {
    NkSpecies = i;
  }
  /*! \brief Accessor to name of species \p i. */
  string GetMixSpecies(int i)
  {
    itSp = mySpecies.begin();
    for (int j = 0; j < i; j++)
    {
      itSp++;
    }
    return itSp->first;
  }
  /*! \brief Accessor to mass fraction of species \p speciesname. */
  double GetMixYk(string speciesname)
  {
    return mySpecies[speciesname].GetYk();
  }
  /*! \brief Accessor to mixture temperature #T. */
  double GetMixT()
  {
    return T;
  }
  /*! \brief Accessor to mixture molecular mass #M. */
  double GetMixM()
  {
    return M;
  }
  /*! \brief Accessor to gas constant #R_over_M (=R/M). */
  double GetMixR_o_M()
  {
    return R_over_M;
  }
  /*! \brief Accessor to mixture laminar viscosity #mul. */
  double GetMixMul()
  {
    return mul;
  }
  /*! \brief Accessor to mixture thermal conductivity #lambda. */
  double GetMixLambda()
  {
    return lambda;
  }
  /*! \brief Accessor to mixture heat capacity ratio #gama. */
  double GetMixGama()
  {
    return gama;
  }
  /*! \brief Accessor to mixture speed of sound #a_sound. */
  double GetMixAsound()
  {
    return a_sound;
  }
  /*! \brief Accessor to mixture heat capacity #cp. */
  double GetMixCp()
  {
    return cp;
  }
  /*! \brief Accessor to mixture enthalpy #H. */
  double GetMixH()
  {
    return H;
  }
  /*! \brief Accessor to mixture Lewis number Le. */
  double GetMixLewis()
  {
    return Lewis;
  }
  /*! \brief Accessor to mixture heat diffusivity coefficient #rhoD_th. */
  double GetMixRhoDk()
  {
    return rhoDk;
  }
  /*! \brief Accessor to mixture #pressure. */
  double GetMixPressure()
  {
    return pressure;
  }
  /*! \brief Accessor to #CombustionRegime (0: Mixing, 1: Combustion). */
  void SetCombustionRegime(int val)
  {
    CombustionRegime = val;
  }
  /*! \brief Accessor to #CombustionRegime (0: Mixing, 1: Combustion). */
  int GetCombustionRegime()
  {
    return CombustionRegime;
  }
  /***********************************************************************************************************/
  /* other public functions */

  /*! \brief Load species contained in the input file from thermodynamics input file, initialize some constants.
   *
   *  The species are stored in a map container with their name (string) as key to retrieve them easily.
   *  Thermodynamics input file (binary) is the same as the input file for ScanMan in FlameMaster, created by
   *  FlameMaster tool \p CreateBinFile. 2 versions of the thermodynamics input file can be read: the official
   *  version from FlameMaster or the version of Guillaume Blanquart (the second version contains additional
   *  species data (e.g., shape, alpha) which are currently not used but still saved. If first version of file
   *  is used, these variables are set to -1. The version of the thermodynamics input file is controled by
   *  a flag in the input file: THERMO_BLANQ_FLAG. 
   *  \param[in] myInputPtr A pointer of type ParamMap pointing on the general input file.
   *  \warning The exact format of the binary thermodynamics input file must match!
   */
  void Load(ParamMap *myInputPtr)
  {
    FILE   *inFp;
    size_t size_read, dum;
    int    i, counter = 0, int_buff, thermo_blanq_flag;
    char   filename[kStrLong], buffer_m[kNameLen];
    Param  *p = NULL;
    Param  *bc_l = NULL;
    Param  *bc_r = NULL;
    char   *q = NULL;
    string key, buffer_str;
    bool   is_dominantsp_here = false;
    double *Yk_left = NULL;
    double *Yk_right = NULL;
    double sum_yk_l = 0.0, sum_yk_r = 0.0;

    /* variables to be read from thermo file. Following struct only used to read thermodynamics input file */
    typedef struct SR_buff_blanquart
    {
      char   name[kNameLen];
      double M;
      compo  composition_elems;
      double Tmed;
      double hot[kCoefLen];
      double cold[kCoefLen];
      int    shape;
      double k_over_eps;
      double sigma;
      double mu;
      double alpha;
      double Zrot;
    } SR_buff_blanquart;
    SR_buff_blanquart sr;

    /* Only used for compatibility with FlameMaster structure (no CL present in this version) */
    typedef struct compo_old
    {
      unsigned char C, H, N, O, BR, F, AR, HE;
    } compo_old;

    /* variables to be read from thermo file. Following struct only used to read thermodynamics input file */
    typedef struct SR_buff_FlameMaster
    {
      char      name[kNameLen];
      double    M;
      compo_old composition_elems;
      double    hot[kCoefLen];
      double    cold[kCoefLen];
      double    k_over_eps;
      double    sigma;
      double    mu_coeff;
      double    D_coeff;
      double    omega_coeff;
      int       id;
    } SR_buff_FlameMaster;
    SR_buff_FlameMaster srfm;

    /* Read input file for thermodynamic input information */
    NkSpecies = myInputPtr->getIntParam("NUMBER_SPECIES");
    if (CombustionRegime > 0) dominant_species = myInputPtr->getStringParam("DOM_SPECIES");
    SpeciesName = new char[NkSpecies][kNameLen];
    Yk_left  = new double[NkSpecies];
    Yk_right = new double[NkSpecies];

    if (CombustionRegime == -1)
      bc_l = myInputPtr->getParam("MIXTURE_COMPO");
    else if (CombustionRegime == 0)
    {
      bc_l = myInputPtr->getParam("MIXING_BC_LEFT");
      bc_r = myInputPtr->getParam("MIXING_BC_RIGHT");
    } 
    p = myInputPtr->getParam("SPECIES");
    if (NkSpecies != p->getSize() - 1)
    {
      cerr << "### Number of species read " << p->getSize() - 1 << " does not match the input number of species " << NkSpecies << " ###" << endl;
      throw(-1);
    }
    for (i = 0; i < NkSpecies; i++)
    {
      buffer_str = p->getString(i + 1);
      if (CombustionRegime <= 0)
        is_dominantsp_here = true;
      else
      {
        if (dominant_species.compare(buffer_str) == 0)
        {
          is_dominantsp_here = true;
        }
      }
      dum = buffer_str.copy(&buffer_m[0], kNameLen);
      dum = buffer_str.size();
      if (dum < kNameLen)
        strcpy(&buffer_m[dum], "\0");
      strcpy(&SpeciesName[i][0], buffer_m);

      if (CombustionRegime == -1)
      {
        Yk_left[i] = bc_l->getDouble(i + 1);
        Yk_right[i] = -1.0;
        sum_yk_l += Yk_left[i];
        sum_yk_r = 1.0;
        if (Yk_left[i] < 0.0)
        {
          cerr << "### Species mass fractions cannot be negative! ###" << endl;
          throw(-1);
        }
      }
      else if (CombustionRegime == 0)
      {
        Yk_left[i]  = bc_l->getDouble(i + 1);
        Yk_right[i] = bc_r->getDouble(i + 1);
        sum_yk_l += Yk_left[i];
        sum_yk_r += Yk_right[i];
        if ((Yk_left[i] < 0.0) || (Yk_right[i] < 0.0))
        {
          cerr << "### Species mass fractions cannot be negative! ###" << endl;
          throw(-1);
        }
      }
      else
      {
        Yk_left[i]  = -1.0;
        Yk_right[i] = -1.0;
        sum_yk_l = 1.0;
        sum_yk_r = 1.0;
      }
    }
    if (!is_dominantsp_here)
    {
      cerr << "### The dominant species " << dominant_species << " is not part of the species list ###" << endl;
      throw(-1);
    }
    if (sum_yk_l != 1.0)
    {
      cerr << "### The sum of species mass fraction at the left boundary is not 1 but " << sum_yk_l << " ###" << endl;
      throw(-1);
    }
    if (sum_yk_r != 1.0)
    {
      cerr << "### The sum of species mass fraction at the right boundary is not 1 but " << sum_yk_r << " ###" << endl;
      throw(-1);
    }

    if ((mpi_rank == 0) && (OutputToScreen == "YES"))
    {
      cout << endl;
      cout << "The combustion regime is: ";
      if (CombustionRegime == -1)
        cout << "Variable Properties" << endl;
      if (CombustionRegime == 0)
        cout << "Mixing" << endl;
      if (CombustionRegime == 1)
        cout << "Steady flamelet" << endl;
      if (CombustionRegime == 2)
        cout << "Flamelet with progress variable FPVA" << endl;
      
      cout << "Combustion part of input file read, the " << NkSpecies << " species are:" << endl;
      for (i = 0; i < NkSpecies; i++)
      {
        strcpy(buffer_m, &SpeciesName[i][0]);
        cout << buffer_m << "\t Yk(Z=0)=" << Yk_left[i] << "\t Yk(Z=1)="  << Yk_right[i] << endl;
      }
      cout << endl;
      if (CombustionRegime > 0) cout << "The dominant species is " << dominant_species << endl;
    }

    buffer_str = myInputPtr->getStringParam("THERMO_INPUT_FILE");
    dum = buffer_str.copy(&filename[0], kStrLong);
    dum = buffer_str.size();
    if (dum < kStrLong)
      strcpy(&filename[dum], "\0");

    /* Read thermodynamic input file */
    if ((mpi_rank == 0) && (OutputToScreen == "YES"))
    {
      cout << endl << "Thermo input file is " << filename << endl;
    }
    if (!(inFp = fopen(filename, "rb")))
    {
      cerr << "### Could not open input file " << filename << " ###" << endl;
      throw(-1);
    }
    thermo_blanq_flag = myInputPtr->getIntParam("THERMO_BLANQ_FLAG", "0");
    /* if input file is from G. Blanquart version of CreateBinFile */
    if (thermo_blanq_flag)
    {

      if ((mpi_rank == 0) && (OutputToScreen == "YES"))
      {
        cout << "Used CreateBinFile version by G. Blanquart for thermo input file" << endl;
      }
      /* Loop over input file, read data and create/load new species if required */
      dum = fread(&int_buff, sizeof(int_buff), 1, inFp);
      if (dum != 1)
      {
        cerr << "### Could not read int_buff in thermo file ###" << endl;
        throw(-1);
      }
      /* if int_buff is not 0, there are diffusion coefficients in file, erase before using CreateBinFile, read not implemented */
      if (int_buff != 0)
      {
        cerr << "### Diffusion coefficients present, but read not yet implemented, change CreateBinFile input file ###" << endl;
        throw(-1);
      }
      while (size_read = fread(&sr, sizeof(SR_buff_blanquart), 1, inFp))
      {
        /* Check if species should be loaded by comparing names */
        for (i = 0; i < NkSpecies; i++)
        {
          strcpy(buffer_m, &SpeciesName[i][0]);
          if (!strcmp(buffer_m, sr.name))
          {
            /* Create new entry in mySpecies map */
            key.assign(buffer_m);
            retS = mySpecies.insert(pair<string, Species> (string(key),Species()));
            if (retS.second == false)
            {
              cout << "### Species " << key << " already in mySpecies (in loop i=" << i << ") ###" << endl;
            }
            else
            {
              /* Populate new species */
              counter++;
              retS.first->second.SetSpeciesName(buffer_m);
              retS.first->second.SetSpeciesTmed(sr.Tmed);
              retS.first->second.SetSpeciesHot(sr.hot);
              retS.first->second.SetSpeciesCold(sr.cold);
              retS.first->second.SetSpeciesShape(sr.shape);
              retS.first->second.SetSpeciesKeps(sr.k_over_eps);
              retS.first->second.SetSpeciesSigma(sr.sigma);
              retS.first->second.SetSpeciesDipole(sr.mu);
              retS.first->second.SetSpeciesAlpha(sr.alpha);
              retS.first->second.SetSpeciesZrot(sr.Zrot);
              retS.first->second.SetYk_bc_left(Yk_left[i]);
              retS.first->second.SetYk_bc_right(Yk_right[i]);
            }

          }
        }
      }
    }

    /* or if it is standard CreateBinFile from FlameMaster */
    else
    {
      if ((mpi_rank == 0) && (OutputToScreen == "YES"))
      {
        cout << "Used standard CreateBinFile version from FlameMaster for thermo input file" << endl;
      }
      /* Loop over input file, read data and create/load new species if required */
      while (size_read = fread(&srfm, sizeof(SR_buff_FlameMaster), 1, inFp))
      {
        /* Check if species should be loaded by comparing names */
        for (i = 0; i < NkSpecies; i++)
        {
          strcpy(buffer_m, &SpeciesName[i][0]);
          if (!strcmp(buffer_m, srfm.name))
          {
            /* Create new entry in mySpecies map */
            key.assign(buffer_m);
            retS = mySpecies.insert(pair<string, Species> (string(key),Species()));
            if (retS.second == false)
            {
              cout << "### Species " << key << " already in mySpecies (in loop i=" << i << ") ###" << endl;
            }
            else
            {
              /* Populate new species */
              counter++;
              retS.first->second.SetSpeciesName(buffer_m);
              retS.first->second.SetSpeciesTmed(1000.0);
              retS.first->second.SetSpeciesHot(srfm.hot);
              retS.first->second.SetSpeciesCold(srfm.cold);
              retS.first->second.SetSpeciesShape(-1);
              retS.first->second.SetSpeciesKeps(srfm.k_over_eps);
              retS.first->second.SetSpeciesSigma(srfm.sigma);
              retS.first->second.SetSpeciesDipole(-1.0);
              retS.first->second.SetSpeciesAlpha(-1.0);
              retS.first->second.SetSpeciesZrot(-1.0);
              retS.first->second.SetYk_bc_left(Yk_left[i]);
              retS.first->second.SetYk_bc_right(Yk_right[i]);
            }

          }
        }
      }
    }

    fclose(inFp);
    delete[] Yk_left;
    delete[] Yk_right;

    /* Check if all species given by input file initialized */
    if (counter != NkSpecies)
    {
      cerr << "### Number of species loaded " << counter << " does not match the input number of species " << NkSpecies << " ###" << endl;
      throw(-1);
    }

    /* Now initialize some species constants */
    for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
    {
      itSp->second.ComputeSpeciesConstants();
    }
    Lewis = 1.0; // HACK: temporary, could be changed later!

    /* Output species data for check */
    if ((mpi_rank == 0) && (OutputToScreen == "YES"))
    {
      cout << endl << "Thermo properties input file read, species properties are:" << endl << endl;
      for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
      {
        itSp->second.OutputSpeciesConstants();
      }
    }
    return;
  }

  /***********************************************************************************************************/

  void Load(int CombReg, ParamMap *myInputPtr, string output = "YES")
  {
    OutputToScreen = output;
    CombustionRegime = CombReg;
    Load(myInputPtr);
  }

  /***********************************************************************************************************/  
  
  /*! \brief Clean mixture and species and deallocate memory */
  void Unload()
  {
    delete [] SpeciesName;
    SpeciesName = 0;
    mySpecies.erase(mySpecies.begin(), mySpecies.end());
    return;
  }

  /***********************************************************************************************************/

  /*! \brief Compute mixture properties based on energy (chemical and sensible) and density.
   *
   *  Compute mixture temperature #T, viscosity #mul, heat conductivity #lambda,
   *  heat capacity #cp, speed of sound #a_sound, #heat capacity ratio #gama, heat diffusivity 
   *  coefficient #rhoD_th, #pressure.\n
   *  The viscosity of the mixture is computed from the mixing rule by C.R. Wilke and the heat conductivity
   *  from the mixing rule by E.A. Mason and S.C. Saxena. \n
   *  See \em Transport \em Phenomena, Bird, Stewart, Lightfoot, 2007, (p.27 and p.276). \n
   *  When temperature #T is computed, the species properties must first be computed before the mixture
   *  properties can be computed (done inside this function).\n
   *  The initial guess used in the Newton iteration to solve for the temperature is set arbitrarily to 1000K,
   *  a better approximation could be used but the Newton solver converges fast so that there is no real need.
   *  \param[in] Energy Chemical and sensible energy of mixture in [J/kg] (from flow solver).
   *  \param[in] rho Mixture density in [kg/m^3] (from flow solver).
   *  \warning The definition of the energy must be consistent with the one used in these routines!
   */
  void ComputeMixProperties(double Energy, double rho)
  {
    double Tinitial = 1000.0;
    double Mi_o_Mj, mui_o_muj, Gij, dummy1, dummy2, deltai;

    T = ComputeMixTemperature(Energy, Tinitial);

    /* Compute species properties mu, cp, lambda */
    for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
    {
      itSp->second.ComputeSpeciesProperties(T);
    }

    /* Compute mu and lamda for mixture using working arrays Gij and Delta_i */
    mul = 0.0;
    lambda = 0.0;
    cp = 0.0;
    for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
    {
      deltai = 0.0;
      Gij = 0.0;
      for (itSp2 = mySpecies.begin(); itSp2 != mySpecies.end(); itSp2++)
      {
        Mi_o_Mj = itSp->second.GetSpeciesMolarMass()
        / itSp2->second.GetSpeciesMolarMass();
        mui_o_muj = itSp->second.GetMu() / itSp2->second.GetMu();
        dummy1 = 1.0 / sqrt(8.0 * (1.0 + Mi_o_Mj));
        dummy2 = 1.0 + sqrt(mui_o_muj * sqrt(1.0 / Mi_o_Mj));
        Gij = dummy1 * dummy2 * dummy2;
        deltai += Gij * Mi_o_Mj * itSp2->second.GetYk();
      }
      itSp->second.SetDelta(deltai);
      dummy1 = itSp->second.GetYk() / deltai;
      mul += dummy1 * itSp->second.GetMu();
      lambda += 1.0 / (1.065 / dummy1 - 0.065) * itSp->second.GetLambda();
      cp += itSp->second.GetYk() * itSp->second.GetCp();
    }

    /* Other properties */
    H = ComputeMixEnthalpy(T);
    gama = cp / (cp - R_over_M);
    a_sound = sqrt(gama * R_over_M * T);
    rhoDk = lambda / cp / Lewis;
    pressure = rho * R_over_M * T;

    return;
  }

  /***********************************************************************************************************/

  /*! \brief Compute mixture viscosity #mul and heat conductivity #lambda.
   *
   *  Compute mixture viscosity #mul and heat conductivity #lambda in function of the temperature.
   *  The viscosity of the mixture is computed from the mixing rule by C.R. Wilke and the heat conductivity
   *  from the mixing rule by E.A. Mason and S.C. Saxena. \n
   *  See \em Transport \em Phenomena, Bird, Stewart, Lightfoot, 2007, (p.27 and p.276). \n
   *  \param[in] temp Temperature in [K].
   *  \warning The species mass fractions must have been first computed!
   */
  void ComputeMixLambdaMu(double temp)
  {
    double Mi_o_Mj, mui_o_muj, Gij, dummy1, dummy2, deltai;

    /* Compute species properties mu, cp, lambda */
    for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
    {
      itSp->second.ComputeSpeciesProperties(temp);
    }

    /* Compute mu and lamda for mixture using working arrays Gij and Delta_i */
    mul = 0.0;
    lambda = 0.0;
    for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
    {
      deltai = 0.0;
      Gij = 0.0;
      for (itSp2 = mySpecies.begin(); itSp2 != mySpecies.end(); itSp2++)
      {
        Mi_o_Mj = itSp->second.GetSpeciesMolarMass()
        / itSp2->second.GetSpeciesMolarMass();
        mui_o_muj = itSp->second.GetMu() / itSp2->second.GetMu();
        dummy1 = 1.0 / sqrt(8.0 * (1.0 + Mi_o_Mj));
        dummy2 = 1.0 + sqrt(mui_o_muj * sqrt(1.0 / Mi_o_Mj));
        Gij = dummy1 * dummy2 * dummy2;
        deltai += Gij * Mi_o_Mj * itSp2->second.GetYk();
      }
      itSp->second.SetDelta(deltai);
      dummy1 = itSp->second.GetYk() / deltai;
      mul += dummy1 * itSp->second.GetMu();
      lambda += dummy1 * itSp->second.GetLambda();
    }

    return;
  }

  /***********************************************************************************************************/

  /*! \brief Compute mixture molecular mass #M and gas constant #R_over_M based on composition.
   *
   *  \warning The species mass fraction Species::Yk must first have been computed!
   */
  void ComputeMixM()
  {
    double m = 0.0;
    for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
    {
      m += itSp->second.GetYk() / itSp->second.GetSpeciesMolarMass();
    }
    M = 1.0 / m;
    R_over_M = kR / M;
    return;
  }

  /***********************************************************************************************************/

  /*! \brief Compute mixture enthalpy (chemical and sensible) based on composition and temperature.
   *
   *  Compute mixture enthalpy as a weighted sum (species mass fractions) of the species enthalpies (linear mixing rule).
   *  \param[in] temp Temperature in [K].
   *  \return Mixture enthalpy in [J/kg] for given temperature.
   *  \warning The species mass fraction Species::Yk must first have been computed!
   */
  double ComputeMixEnthalpy(double temp)
  {
    double Enthalpy;
    int i;

    Enthalpy = 0.0;
    if (temp > TPolyLinMax)
    {
      for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
      {
        Enthalpy += itSp->second.GetYk() * itSp->second.ComputeSpeciesEnthalpy(TPolyLinMax);
      }       
      Enthalpy += ComputeMixCp(TPolyLinMax) * (temp - TPolyLinMax);
    }
    else if (temp < TPolyLinMin)
    {
      for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
      {
        Enthalpy += itSp->second.GetYk() * itSp->second.ComputeSpeciesEnthalpy(TPolyLinMin);
      }       
      Enthalpy += ComputeMixCp(TPolyLinMin) * (temp - TPolyLinMin);
    }
    else
    {
      for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
      {
        Enthalpy += itSp->second.GetYk() * itSp->second.ComputeSpeciesEnthalpy(temp);
      }
    }
    return Enthalpy;
  }

  /***********************************************************************************************************/

  /*! \brief Compute mixture energy (chemical and sensible) based on composition and temperature.
   *
   *  Compute mixture energy from mixture enthalpy and temperature.
   *  \param[in] temp Temperature in [K].
   *  \return Mixture energy in [J/kg] for given temperature.
   *  \warning The species mass fraction Species::Yk must first have been computed!
   */
  double ComputeMixEnergy(double temp)
  {
    return ComputeMixEnthalpy(temp) - R_over_M * temp;
  }

  /***********************************************************************************************************/

  /*! \brief Compute mixture heat capacity based on composition and temperature.
   *
   *  Compute mixture heat capacity as a weighted sum (species mass fractions) of the species heat capacities (linear mixing rule).
   *  \param[in] temp Temperature in [K].
   *  \return Mixture heat capacity in [J/(kg*K)] for given temperature.
   *  \warning The species mass fraction Species::Yk must first have been computed!
   */
  double ComputeMixCp(double temp)
  {
    double Cpp;
    int i;

    Cpp = 0.0;
    if (temp > TPolyLinMax)
    {
      for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
      {
        Cpp += itSp->second.GetYk() * itSp->second.ComputeSpeciesCp(TPolyLinMax);
      }       
    }
    else if (temp < TPolyLinMin)
    {
      for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
      {
        Cpp += itSp->second.GetYk() * itSp->second.ComputeSpeciesCp(TPolyLinMin);
      }       
    }
    else
    {
      for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
      {
        Cpp += itSp->second.GetYk() * itSp->second.ComputeSpeciesCp(temp);
      }
    }
    return Cpp;
  }

  /***********************************************************************************************************/
  /*! \brief Compute mixture temperature based on composition and energy.
   *
   *  Compute mixture temperature using a Newton iteration. The maximum number of iteration is large enough
   *  to ensure convergence (usually takes 3 or 4 iterations).
   *  \param[in] Energy Energy of mixture in [J/kg].
   *  \param[in] Tguess Initial guess for temperature in [K].
   *  \return Mixture temperature in [K]. If iteration has not converged, the program stops with an error message.
   *  \warning The species mass fraction Species::Yk must first have been computed! 
   */
  double ComputeMixTemperature(double Energy, double Tguess = 1000.0)
  {
    int iter;
    double Tresidual, ee, cpp, error;

    iter = 0;
    if (Energy < ComputeMixEnergy(TMinClip)) return TMinClip; 
    else if (Energy > ComputeMixEnergy(TMaxClip)) return TMaxClip;
    do
    {
      ++iter; 
      //      hh = ComputeMixEnthalpy(Tguess);
      ee = ComputeMixEnergy(Tguess);
      cpp = ComputeMixCp(Tguess);

      //      Tresidual = Tguess - (hh - Tguess * R_over_M - Energy) / (cpp - R_over_M);
      Tresidual = Tguess - (ee - Energy) / (cpp - R_over_M);

      error = Tguess - Tresidual;
      if (fabs(error) <= NewtonEpsilon)
      {
        //cout << "Newton iterations: " << iter << endl;
        return Tresidual;
      }

      Tguess = Tresidual;

    }
    while (iter < NewtonIterMax);

    /* if the iteration has not converged, exit with error */
    cerr << "### Computation of temperature in Newton iteration has not converged: Tguess="
    << Tguess + error << ", Tresidual=" << Tresidual << ", energy="
    << Energy << ", #iterations=" << iter << " ###" << endl;
    throw(-1);
    return Tresidual; // clipping to ensure simulation to run
  }

  /***********************************************************************************************************/
  /*! \brief Compute mixture temperature based on composition and enthalpy.
   *
   *  Compute mixture temperature using a Newton iteration. The maximum number of iteration is large enough
   *  to ensure convergence (usually takes 3 or 4 iterations).
   *  \param[in] Enthalpy Enthalpy of mixture in [J/kg].
   *  \param[in] Tguess Initial guess for temperature in [K].
   *  \return Mixture temperature in [K]. If iteration has not converged, the program stops with an error message.
   *  \warning The species mass fraction Species::Yk must first have been computed! 
   */
  double ComputeMixTemperature_H(double Enthalpy, double Tguess = 1000.0)
  {
    int iter;
    double Tresidual, hh, cpp, error;

    iter = 0;
    if (Enthalpy < ComputeMixEnthalpy(TMinClip)) return TMinClip; 
    else if (Enthalpy > ComputeMixEnthalpy(TMaxClip)) return TMaxClip;
    do
    {
      ++iter;
      hh = ComputeMixEnthalpy(Tguess);
      cpp = ComputeMixCp(Tguess);

      Tresidual = Tguess - (hh - Enthalpy) / cpp;

      error = Tguess - Tresidual;
      if (fabs(error) <= NewtonEpsilon)
        return Tresidual;

      Tguess = Tresidual;

    }
    while (iter < NewtonIterMax);

    /* if the iteration has not converged, exit with error */
    cerr << "### Computation of temperature in Newton iteration has not converged: Tguess="
    << Tguess + error << ", Tresidual=" << Tresidual << ", enthalpy="
    << Enthalpy << ", #iterations=" << iter << " ###" << endl;
    throw(-1);
    return Tresidual;  // clipping to ensure simulation to run
  }

  /***********************************************************************************************************/

  /*! \brief Output on screen the mixture properties Species::Yk, #T, #M, #mul, #lambda, #R_over_M,
   *  #gama, #a_sound, #cp, #H, energy, #rhoD_th, #pressure.
   */
  void OutputMixProperties()
  {
    int i;
    char buffer[kNameLen];
    double sum = 0.0, Y;

    cout << "Mixture properties are:" << endl;
    for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
    {
      itSp->second.GetSpeciesName(buffer);
      Y = itSp->second.GetYk();
      sum += Y;
      cout << " Yk(" << buffer << ")=" << Y << "\t";
    }
    cout << endl;
    cout << "Sum of all species mass fraction = " << sum << endl;
    cout << "T=" << T << endl;
    cout << "M=" << M << endl;
    cout << "mu=" << mul << endl;
    cout << "lambda=" << lambda << endl;
    cout << "R/M=" << R_over_M << endl;
    cout << "gamma=" << gama << endl;
    cout << "a_sound=" << a_sound << endl;
    cout << "Cp=" << cp << endl;
    cout << "H=" << H << endl;
    cout << "E=" << ComputeMixEnergy(T) << endl;
    cout << "rho*Dk=" << rhoDk << endl;
    cout << "P=" << pressure << endl;
    return;
  }

  /***********************************************************************************************************/

  /*! \brief Compute species mass fractions to determine the mixture composition from combustion.
   *
   *  Overloaded function to compute species mass fractions from chemistry table for steady flamelet combustion model.
   *  The different chemistry tables (Cartesian, adaptive) are defined over template.
   *  \param[in] Zm Zmean (mean scalar mixture fraction from flow solver).
   *  \param[in] Zv Zvar (variance of scalar mixture fraction from flow solver).
   *  \param[in] ch Chi/progress variable (scalar dissipation rate or progress variable from flow solver).
   *  \param[in] myChemTable Chemistry table.
   */
  template <class Chemtable>
  void GetSpeciesMassFraction(double Zm, double Zv, double ch, Chemtable &myChemTable)
  {
    double Zvar, Zmean, Chi;
    double yk, sumYk = 0.0;

    // Clip mean of Z
    Zmean = max(min(Zm, 1.0), 0.0);

    // Clip variance of Z to maximum possible value of mean * (1 - mean)
    Zvar = min(Zv, Zmean * (1.0 - Zmean));
    Chi = ch;
    
#ifdef mixing_tabulated
    Zvar = 0.0;
    Chi = 0.0;
#endif

    // Interpolate for each species 
    for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
    {
      if (dominant_species.compare(itSp->first) != 0)
      {
        yk = myChemTable.Lookup(Zmean, Zvar, Chi, itSp->first);
        itSp->second.SetYk(yk);
        sumYk += yk;
      }
    }
    mySpecies[dominant_species].SetYk(1.0 - sumYk);

    ComputeMixM();

    return;
  }

  /***********************************************************************************************************/

  /*! \brief Compute species mass fractions and progress variable source term to determine the mixture composition from combustion.
   *
   *  Overloaded function to compute species mass fractions from chemistry table and source term for progress variable
   *  for FPVA combustion model.
   *  The different chemistry tables (Cartesian, adaptive) are defined over template.
   *  \param[out] CSource Source term for progress variable (interpolated from table and saved to avoid 2 interpolations).
   *  \param[in]  Zm Zmean (mean scalar mixture fraction from flow solver).
   *  \param[in]  Zv Zvar (variance of scalar mixture fraction from flow solver).
   *  \param[in]  ch Chi/progress variable (scalar dissipation rate or progress variable from flow solver).
   *  \param[in]  myChemTable Chemistry table.
   */
  template <class Chemtable>
  void GetSpeciesMassFraction(double &CSource, double Zm, double Zv, double Cm, Chemtable &myChemTable)
  {
    double Zvar, Zmean, Cmean;
    double yk, sumYk = 0.0;

    // Clip mean of Z
    Zmean = max(min(Zm, 1.0), 0.0);

    // Clip variance of Z to maximum possible value of mean * (1 - mean)
    Zvar = min(Zv, Zmean * (1.0 - Zmean));
    Cmean = Cm;
    
#ifdef mixing_tabulated
    Zvar = 0.0;
    Cmean = 0.0;
#endif

    // Interpolate for each species
    for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
    {
      if (dominant_species.compare(itSp->first) != 0)
      {
        yk = myChemTable.Lookup(Zmean, Zvar, Cmean, itSp->first);
        itSp->second.SetYk(yk);
        sumYk += yk;
      }
    }
    mySpecies[dominant_species].SetYk(1.0 - sumYk);

    ComputeMixM();
    
    /* Interpolate source term for progress variable */
    CSource = myChemTable.Lookup(Zmean, Zvar, Cm, "SRC_PROG");

    return;
  }

  /***********************************************************************************************************/

  /*! \brief Compute species mass fractions to determine the mixture composition from linear mixing.
   *
   *  Overloaded functions to compute species mass fraction from linear rules for mixing model.
   *  \param[in] Zm Zmean (mean scalar mixture fraction from flow solver).
   */
  void GetSpeciesMassFraction(double Zm)
  {
    double yk, zm;
    
    /* Interpolate for each species */
    zm = min(max(Zm, 0.0), 1.0);
//    zm = Zm;
    for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
    {
      yk = (1.0 - zm) * itSp->second.GetYk_bc_left() + zm * itSp->second.GetYk_bc_right();
      itSp->second.SetYk(yk);
    }

    ComputeMixM();
    return;
  }

  /***********************************************************************************************************/

  /*! \brief Compute species mass fractions to determine the mixture composition from input file.
   *
   *  Overloaded functions to compute species mass fraction for constant mixture for variable properties model.
   */
  void GetSpeciesMassFraction()
  {
    double yk;
    
    for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
    {
      yk = itSp->second.GetYk_bc_left();
      itSp->second.SetYk(yk);
    }

    ComputeMixM();
    return;
  }

  /***********************************************************************************************************/

  /*! \brief Compute species mass fractions to determine the mixture composition input file.
   *
   *  Function returning the value of the progress variable given Zmean, Zvar and Cmean for FPVA combustion model.
   *  Can be used with input Cmean very big in order to get the value of Cmean max given Zmean and Zvar
   *  in order to initialize/ignite the flame.
   */
  template <class Chemtable>
  void GetProgressVariable(double &Cmnew, double Zm, double Zv, double Cm, Chemtable &myChemTable)
  {
    double Zvar, Zmean, Cmean;

    // Clip mean of Z
    Zmean = max(min(Zm, 1.0), 0.0);

    // Clip variance of Z to maximum possible value of mean * (1 - mean)
    Zvar = min(Zv, Zmean * (1.0 - Zmean));
    Cmean = Cm;
    
#ifdef mixing_tabulated
    Zvar = 0.0;
    Cmean = 0.0;
#endif    

    // Interpolate progress variable
    Cmnew = myChemTable.Lookup(Zmean, Zvar, Cm, "PROG");

    return;
  }

  /***********************************************************************************************************/

  /*! \brief Set species mass fraction Yk given its name k.
   *
   *  \param[in]  Yk     the species mass fraction of species k.
   *  \param[in]  YkName name of the species k.
   */
  void SetSpeciesMassFraction(double Yk, string YkName)
  {
    int IsFound = 0;

    /* Loop over species to find the one corresponding to YkName */
    for (itSp = mySpecies.begin(); itSp != mySpecies.end(); itSp++)
    {
      if (YkName.compare(itSp->first) == 0)
      {
        itSp->second.SetYk(Yk);
        IsFound = 1;
      }
    }
    if (IsFound == 0)
    {
      cerr << "### Species " << YkName << " was not found in the list of species! ###" << endl;
      throw(-1);
    }

    return;
  }

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
