#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string>
#include <sstream>

using namespace std;


/* Define global constants: can be edited */
/***************************************************************************************************************/
#define cyl 0                 // Cylinder flag: 0 -> Cartesian or 1 -> Cylinder (not yet validated)
#define scaleX 0.001          // Scaling factor in X direction
#define scaleY 0.001          // Scaling factor in Y direction
#define scaleZ 0.001          // Scaling factor in Z direction (Thickness should be given before scaling!)

#define FILE_2D_IN  "2d.neu"  // Name of input  file (2D Gambit neutral file)
#define FILE_3D_OUT "3d.neu"  // Name of output file (3D Gambit neutral file)

#define GAMBIT_VER "2.4.6"    // Version of Gambit used to create 2D input file
/***************************************************************************************************************/


/* Define neutral file specific constants */
#define NGRPS 1               // Number of element groups
#define NBSETS 0              // Number of boundary condition sets
#define NDFCD 3               // Number of coordinate directions (2 for neutral file read, 3 for neutral file written)
#define NDFVL 3               // Number of velocity components (2 for neutral file read, 3 for neutral file written, but not used here)

#define NTYPE_TRI 3           // Element geometry type TRIANGLE
#define NDP_TRI 3             // Number of nodes that define element TRIANGLE (node order is 0 -> 1 -> 2)
#define NTYPE_PRISM 5         // Element geometry type PRISM (=WEDGE)
#define NDP_PRISM 6           // Number of nodes that define element PRISM (node order is 0 -> 1 -> 2 -> 3 -> 4 -> 5)
#define NTYPE_QUAD 2          // Element geometry type QUAD
#define NDP_QUAD 4            // Number of nodes that define element QUAD (node order is 0 -> 1 -> 2 -> 3)
#define NTYPE_BRICK 4         // ELement geometry type BRICK (=HEX)
#define NDP_BRICK 8           // Number of nodes that define element BRICK (node order is 0 -> 1 -> 3 -> 2 -> 4 -> 5 -> 7 -> 6) 

/* Define Structures for Geometric Manipulation */
#define PMAX  2000000                          /* Max # of pts in the grid */
#define NTETAMAX  52                           /* Max # of pts in the theta direction */


#ifndef M_PI
#define M_PI 3.141592654
#endif

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif


template <class T>
bool from_string(T& t, const std::string &s, std::ios_base& (*f)(std::ios_base&)) 
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}


/* Global Variables */
double ThicknessZ;        // Thickness in the Z direction (before scaling)
int    nCellsZ;           // Number of cell planes in the Z direction


/* Included functions (not used yet)*/
double Stretching(double t, double C, double D, double E)   // Stretching function
{
  double y, y1, y2;
  D = min( max(0.0,D), 1.0) * 2.0;
  E = min( max(0.0,E), 1.0) * 2.0;

  y = 1.0 + tanh(C*(2.0*t-E));
  y1 = 1.0 + tanh(C*(2.0*0.0-E));
  y2 = 1.0 + tanh(C*(2.0*1.0-E));

  y = (y-y1)/(y2-y1);

  return(y);
}

/* Grid structure */
typedef struct {
  int    nvertices, nelements;     // Total number of vertices and elements
  double *vx2d, *vy2d;             // Coordinates of vertices (X and Y)
  int    *ntype, *ndp;             // Geometry type of element and number of vertices to define it
  int    (*cell2d)[4];             // List of vertices to define element (in 2D case; 4 is maximum number of vertices required)
} grid;

void Read2DGridGambit(grid &gr);          // Read 2D Gambit mesh in Neutral format
void Write3DGridGambitCart(grid &gr);     // Write 3D Cartesian Gambit mesh in Neutral format
void Write3DGridGambitAxi(grid &gr);      // Write 3D cylindrical Gambit mesh in Neutral format (not yet completely implemented/validated)

int main(int argc, char *argv[])
{
  ThicknessZ = 0.01;
  nCellsZ = 1;

  if (argc > 1)
  {
    string str(argv[1]);
    from_string<double>(ThicknessZ, str, std::dec);
  }
  printf("Thickness in Z direction (before scaling): %.8le\n", ThicknessZ);

  if (argc > 2)
  {
    string str(argv[2]);
    from_string<int>(nCellsZ, str, std::dec);
  }
  printf("Number of cells in Z direction: %d\n", nCellsZ);
  
  printf("Scaling in X direction: %.8le\n", scaleX);
  printf("Scaling in Y direction: %.8le\n", scaleY);
  printf("Scaling in Z direction: %.8le\n", scaleZ);

  // Object gr of class grid
  grid gr;

  // Read from file and store 2D mesh
  Read2DGridGambit(gr);

  // Compute and write to file 3D mesh
  if (cyl == 0)
    Write3DGridGambitCart(gr);
  if (cyl == 1)
    Write3DGridGambitAxi(gr);

  return 0;
}

/*---------------------------------------------------------------------
 Read 2D grid input file in Gambit Neutral format
 ---------------------------------------------------------------------*/
void Read2DGridGambit(grid &gr)
{
  char line[256];
  FILE *fpin;
  int i, v, c, j;
  int ngrps, nbsets, ndfcd, ndfvl;
  int ndp, ntype;

  fpin = fopen(FILE_2D_IN, "r");
  
  /* Read header */
  fscanf(fpin, "%s\n", line);
  while (strcmp(line, "NDFVL") != 0)
  {
    fscanf(fpin, "%s\n", line);
  }
  fscanf(fpin, "%i %i %i %i %i %i \n", &gr.nvertices, &gr.nelements, &ngrps, &nbsets, &ndfcd, &ndfvl);
  if (ndfcd != 2)  // Check if it is indeed a 2D grid
  {
    printf(" The input file should be 2D but ndfcd = %i \n", ndfcd);
    printf(" Exiting! \n");
    throw(-1);
  }
  printf(" 2D Grid --> Vertices: %i        Elements: %i \n", gr.nvertices, gr.nelements);
  
  /* Allocate memory to save vertices and elements */
  gr.vx2d   = new double[gr.nvertices];
  gr.vy2d   = new double[gr.nvertices];
  gr.ntype  = new int[gr.nelements];
  gr.ndp    = new int[gr.nelements];
  gr.cell2d = new int[gr.nelements][4];

  while (strcmp(line, GAMBIT_VER) != 0)
  {
    fscanf(fpin, "%s\n", line);
  }

  /* Read nodal coordinates */
  for (i = 0; i < gr.nvertices; i++)
  {
    fscanf(fpin, "%10i%lf%lf\n", &v, &gr.vx2d[i], &gr.vy2d[i]);
  }

  fscanf(fpin, "%s\n", line);
  while (strcmp(line, GAMBIT_VER) != 0)
  {
    fscanf(fpin, "%s\n", line);
  }
  
  /* Read elements */
  for (c = 0; c < gr.nelements; c++)
  {
    fscanf(fpin, "%8i%2i%2i\n", &i, &gr.ntype[c], &gr.ndp[c]);
    
    if (gr.ntype[c] == NTYPE_QUAD)
    {
      fscanf(fpin, "%8i%8i%8i%8i\n", &gr.cell2d[c][0], &gr.cell2d[c][1], &gr.cell2d[c][3], &gr.cell2d[c][2]);   // index 2 and 3 inverted because of ordering in BRICK
      for (j = 0; j < 4; j++)
        gr.cell2d[c][j] = gr.cell2d[c][j] - 1;
    }
    else if (gr.ntype[c] == NTYPE_TRI)
    {
      fscanf(fpin, "%8i%8i%8i\n", &gr.cell2d[c][0], &gr.cell2d[c][1], &gr.cell2d[c][2]);
      for (j = 0; j < 3; j++)
        gr.cell2d[c][j] = gr.cell2d[c][j] - 1;
    }
    else
    {
      printf(" Wrong type of elements! Should be either QUAD or TRIANGLE:  ntype = %i \n", gr.ntype[c]);
      printf(" Exiting! \n");
      throw(-1);
    }
    
  }
  printf("\n 2D input file in Gambit Neutral format read! \n");

}


/*---------------------------------------------------------------------
 Write 3D grid output file in Gambit Neutral format: CARTESIAN
 ---------------------------------------------------------------------*/
void Write3DGridGambitCart(grid &gr)
{
  FILE *fpout;
  int i, j, k, v, c;
  double vx, vy, vz;
  int cellv[8];

  printf("\n Computing and writing to file 3D output file in Gambit Neutral format! \n");
  fpout = fopen(FILE_3D_OUT, "w");
  
  /* Write header */
  fprintf(fpout, "        CONTROL INFO ");
  fprintf(fpout, GAMBIT_VER);
  fprintf(fpout, "\n");
  fprintf(fpout, "** GAMBIT NEUTRAL FILE\n");
  fprintf(fpout, "2D to 3D Translator                        \n");
  fprintf(fpout, "PROGRAM:                2dto3d     VERSION:  2.0.0\n");
  fprintf(fpout, " 1 Jan 2000    00:00:00\n");          // could be improved to give correct time but not necessary
  fprintf(fpout, "     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL\n");
  fprintf(fpout, " %9i %9i %9i %9i %9i %9i \n", gr.nvertices * (nCellsZ + 1), gr.nelements * nCellsZ, NGRPS, NBSETS, NDFCD, NDFVL);
  fprintf(fpout, "ENDOFSECTION\n");
  
  printf(" 3D Grid --> Vertices: %i        Elements: %i \n", gr.nvertices * (nCellsZ + 1), gr.nelements * nCellsZ);


  /* Write nodal coordinates */
  fprintf(fpout, "   NODAL COORDINATES 1.2.1 ");
  fprintf(fpout, GAMBIT_VER);
  fprintf(fpout, "\n");

  for (k = 0; k < nCellsZ + 1; k++)
    for (i = 0; i < gr.nvertices; i++)
    {
      vx = gr.vx2d[i];
      vy = gr.vy2d[i];
      vz = (double) k / (double) nCellsZ * ThicknessZ ;
      v = i + (k) * gr.nvertices + 1;
      fprintf(fpout, "%10i%20.10e%20.10e%20.10e\n", v, vx * scaleX, vy * scaleY, vz * scaleZ);
    }
  fprintf(fpout, "ENDOFSECTION\n");

  /* Write elements */
  fprintf(fpout, "      ELEMENTS/CELLS ");
  fprintf(fpout, GAMBIT_VER);
  fprintf(fpout, "\n");
  for (k = 0; k < nCellsZ; k++)
    for (i = 0; i < gr.nelements; i++)
    {
      if (gr.ntype[i] == NTYPE_QUAD)
      {
        cellv[0] = gr.cell2d[i][0] + (k + 1) * gr.nvertices + 1;
        cellv[1] = gr.cell2d[i][1] + (k + 1) * gr.nvertices + 1;
        cellv[2] = gr.cell2d[i][2] + (k + 1) * gr.nvertices + 1;
        cellv[3] = gr.cell2d[i][3] + (k + 1) * gr.nvertices + 1;
        cellv[4] = gr.cell2d[i][0] + (k) * gr.nvertices + 1;
        cellv[5] = gr.cell2d[i][1] + (k) * gr.nvertices + 1;
        cellv[6] = gr.cell2d[i][2] + (k) * gr.nvertices + 1;
        cellv[7] = gr.cell2d[i][3] + (k) * gr.nvertices + 1;
        c = i + (k) * gr.nelements + 1;
  
        fprintf(fpout, "%8i %2i %2i %8i%8i%8i%8i%8i%8i%8i\n               %8i\n", 
            c, NTYPE_BRICK, NDP_BRICK, cellv[0], cellv[1], cellv[2], cellv[3],cellv[4], cellv[5], cellv[6], cellv[7]);
      }
      else if (gr.ntype[i] == NTYPE_TRI)
      {
        cellv[0] = gr.cell2d[i][0] + (k + 1) * gr.nvertices + 1;
        cellv[1] = gr.cell2d[i][1] + (k + 1) * gr.nvertices + 1;
        cellv[2] = gr.cell2d[i][2] + (k + 1) * gr.nvertices + 1;
        cellv[4] = gr.cell2d[i][0] + (k) * gr.nvertices + 1;
        cellv[5] = gr.cell2d[i][1] + (k) * gr.nvertices + 1;
        cellv[6] = gr.cell2d[i][2] + (k) * gr.nvertices + 1;
        c = i + (k) * gr.nelements + 1;

        fprintf(fpout, "%8i %2i %2i %8i%8i%8i%8i%8i%8i\n",
            c, NTYPE_PRISM, NDP_PRISM, cellv[0], cellv[1], cellv[2], cellv[4], cellv[5], cellv[6]);
      }
      else
      {
        printf(" Wrong type of elements! Should be either QUAD or TRIANGLE:  ntype = %i \n", gr.ntype[i]);
        printf(" Exiting! \n");
        throw(-1);        
      }
    }
  fprintf(fpout, "ENDOFSECTION\n");

  /* Write elements group */
  fprintf(fpout, "       ELEMENT GROUP ");
  fprintf(fpout, GAMBIT_VER);
  fprintf(fpout, "\n");
  fprintf(fpout, "GROUP:          1 ELEMENTS: %10i MATERIAL:          2 NFLAGS:          1\n", gr.nelements * nCellsZ);
  fprintf(fpout, "                           fluid\n");
  j = 0;
  fprintf(fpout, "%8i\n", j);
  for (k = 0; k < nCellsZ; k++)
    for (i = 0; i < gr.nelements; i++)
    {
      c = i + (k) * gr.nelements + 1;
      j++;
      fprintf(fpout, "%8i", c);
      if (j == 10)
      {
        fprintf(fpout, "\n");
        j = 0;
      }
    }
  if (j != 0)
    fprintf(fpout, "\n");
  fprintf(fpout, "ENDOFSECTION\n");

  printf("\n 3D output file written in Gambit Neutral format! \n");
}

/*---------------------------------------------------------------------
 Write 3D grid output file in Gambit Neutral format: AXI-SYMMETRIC
 ---------------------------------------------------------------------*/
void Write3DGridGambitAxi(grid &gr)
{
#if 0
  FILE *fpout;
  int i, j, k, v, c;
  int ngrps = 1, nbsets = 0, ndfcd = 3, ndfvl = 3;
  int ntype = 4, ndp = 8;
  double vx, vy, vz, theta;
  int cellv[8];
  int nonaxis, vonaxis[PMAX], ntotvertices, vteta[PMAX][NTETAMAX];

  nonaxis = 0;
  /* Check Point on the Axis */
  for (i = 0; i < nvertices; i++)
  {
    vonaxis[i] = 0;
    if (vy2d[i] == 0)
      vonaxis[i] = 1;
    nonaxis += vonaxis[i];
  }
  ntotvertices = (nvertices - nonaxis) * (nCellsZ - 1) + nvertices;

  /* Write Header */
  fpout = fopen("3d.neu", "w");
  fprintf(fpout, "CONTROL INFO 1.2.1\n");
  fprintf(fpout, "** GAMBIT NEUTRAL FILE\n");
  fprintf(fpout, "2D to 3D Translator                        \n");
  fprintf(fpout, "PROGRAM:                2dto3d     VERSION:  1.0.0\n");
  fprintf(fpout, " 1 Jan 2000    01:07:49\n");
  fprintf(fpout, "     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL\n");
  fprintf(fpout, " %9i %9i %9i %9i %9i %9i \n", ntotvertices, nelements * nCellsZ, ngrps, nbsets, ndfcd, ndfvl);
  fprintf(fpout, "ENDOFSECTION\n");
  printf(" 3D Grid --> Vertices %i Elements %i \n", ntotvertices, nelements * nCellsZ);

  /* Write Nodal Coordinates */
  v = 0;
  fprintf(fpout, "   NODAL COORDINATES 1.2.1\n");
  k = 0;
  for (i = 0; i < nvertices; i++)
  {
    theta = 2 * M_PI * (double) k / (double) nCellsZ;
    vx = vx2d[i], vy = vy2d[i] * cos(theta), vz = vy2d[i] * sin(theta);
    fprintf(fpout, "%10i%20.10e%20.10e%20.10e\n", v + 1, vx, vy, vz);
    vteta[i][k] = v + 1;
    v++;
  }
  for (k = 1; k < nCellsZ; k++)
    for (i = 0; i < nvertices; i++)
    {
      if (vonaxis[i] != 1)
      {
        theta = 2 * M_PI * (double) k / (double) nCellsZ;
        vx = vx2d[i], vy = vy2d[i] * cos(theta), vz = vy2d[i] * sin(theta);
        fprintf(fpout, "%10i%20.10e%20.10e%20.10e\n", v + 1, vx, vy, vz);
        vteta[i][k] = v + 1;
        v++;
      } else
      {
        vteta[i][k] = vteta[i][k - 1];
      }
    }
  k = nCellsZ;
  for (i = 0; i < nvertices; i++)
    vteta[i][k] = vteta[i][0];
  fprintf(fpout, "ENDOFSECTION\n");

  /* Write Elements */
  fprintf(fpout, "      ELEMENTS/CELLS 1.2.1\n");
  for (k = 0; k < nCellsZ; k++)
    for (i = 0; i < nelements; i++)
    {
      ntype = 4, ndp = 8;

      cellv[0] = vteta[cell2d[i][0]][k + 1];
      cellv[1] = vteta[cell2d[i][1]][k + 1];
      cellv[2] = vteta[cell2d[i][2]][k + 1];
      cellv[3] = vteta[cell2d[i][3]][k + 1];
      cellv[4] = vteta[cell2d[i][0]][k];
      cellv[5] = vteta[cell2d[i][1]][k];
      cellv[6] = vteta[cell2d[i][2]][k];
      cellv[7] = vteta[cell2d[i][3]][k];

      c = i + (k) * nelements + 1;

      if (cellv[0] == cellv[4] && cellv[2] == cellv[6])
      {
        ntype = 5, ndp = 6;
        fprintf(fpout, "%8i %2i %2i %8i%8i%8i%8i%8i%8i\n", c, ntype, ndp, cellv[0], cellv[1], cellv[5], cellv[2], cellv[3], cellv[7]);
      } else if (cellv[0] == cellv[4] && cellv[1] == cellv[5])
      {
        ntype = 5, ndp = 6;
        fprintf(fpout, "%8i %2i %2i %8i%8i%8i%8i%8i%8i\n", c, ntype, ndp, cellv[0], cellv[6], cellv[2], cellv[1], cellv[7], cellv[3]);
      } else if (cellv[2] == cellv[6] && cellv[3] == cellv[7])
      {
        ntype = 5, ndp = 6;
        fprintf(fpout, "%8i %2i %2i %8i%8i%8i%8i%8i%8i\n", c, ntype, ndp, cellv[0], cellv[4], cellv[2], cellv[1], cellv[5], cellv[3]);
      } else if (cellv[1] == cellv[5] && cellv[3] == cellv[7])
      {
        ntype = 5, ndp = 6;
        fprintf(fpout, "%8i %2i %2i %8i%8i%8i%8i%8i%8i\n", c, ntype, ndp, cellv[0], cellv[1], cellv[4], cellv[2], cellv[3], cellv[6]);
      } else
      {
        ntype = 4, ndp = 8;
        fprintf(fpout, "%8i %2i %2i %8i%8i%8i%8i%8i%8i%8i\n               %8i\n", c, ntype, ndp, cellv[0], cellv[1], cellv[2], cellv[3],
            cellv[4], cellv[5], cellv[6], cellv[7]);
      }

    }
  fprintf(fpout, "ENDOFSECTION\n");

  /* Write Elements Group */
  fprintf(fpout, "       ELEMENT GROUP 1.2.1\n");
  fprintf(fpout, "GROUP:          1 ELEMENTS: %10i MATERIAL:          2 NFLAGS:          1\n", nelements * nCellsZ);
  fprintf(fpout, "                           fluid\n");
  j = 0;
  fprintf(fpout, "%8i\n", j);
  for (k = 0; k < nCellsZ; k++)
    for (i = 0; i < nelements; i++)
    {
      c = i + (k) * nelements + 1;
      j++;
      fprintf(fpout, "%8i", c);
      if (j == 10)
      {
        fprintf(fpout, "\n");
        j = 0;
      }
    }
  if (j != 0)
    fprintf(fpout, "\n");
  fprintf(fpout, "ENDOFSECTION\n");
#endif
}

