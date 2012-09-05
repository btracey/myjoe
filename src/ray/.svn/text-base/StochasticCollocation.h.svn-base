#ifndef STOCHASTICCOLLOCATION_H_
#define STOCHASTICCOLLOCATION_H_


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "quad_functions.h"

#define EPS 3.0e-14 //EPS is the relative precision.



#if 0
/*
 * Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
 * arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss-
 * Legendre n-point quadrature formula.
 */
int gauss_1d_pts_wts(double *x, double *w, int n)
{
  int m, j, i;
  double z1, z, xm, xl, pp, p3, p2, p1; //High precision is a good idea for this routine.

  m = (n + 1) / 2; //The roots are symmetric in the interval, so we only have to find half of them.
  xm = 0.0;
  xl = 1.0;

  //Loop over the desired roots.
  for (i = 0; i < m; i++)
  {
    z = cos(M_PI * (i + 0.75) / (n + 0.5));

    /*Starting with the above approximation to the ith root, we enter the main loop of
     refinement by Newton's method.*/
    do
    {
      p1 = 1.0;
      p2 = 0.0;

      // Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
      for (j = 0; j < n; j++)
      {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * (j + 1.0) - 1.0) * z * p2 - j * p3) / (j + 1.0);
      }
      /*p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
       by a standard relation involving also p2, the polynomial of one lower order.*/
      pp = n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp; //Newton's method.
    } while (fabs(z - z1) > EPS);

    //Scale the root to the desired interval, and put in its symmetric counterpart.
    x[i] = xm - xl * z;
    x[n - 1 - i] = xm + xl * z;

    //Compute the weight and its symmetric counterpart.
    w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    w[n - 1 - i] = w[i];
  }

  return 0;
}

int clenshaw_curtis_1d_pts_wts(double *x, double *w, int n)
{
  int i, j;
  double s, c;

  if (n == 1)
  {
    x[0] = 0.0;
    w[0] = 2.0;
    return 0;
  }

  for (i = 0; i < n; i++)
  {
    c = (double) i / ((double) n - 1);
    if (2 * c < 1.0 + EPS && 2 * c > 1.0 - EPS)
    {
      x[i] = 0.0;
    } else
    {
      x[i] = -cos(M_PI * c);
    }
  }

  for (i = 0; i < n; i++)
  {
    if (i == 0 || i == n - 1)
    {
      w[i] = 1.0 / ((double) n * ((double) n - 2));
    } else
    {
      s = 0.0;
      for (j = 0; j < (n - 3) / 2; j++)
      {
        s = s + (1 / (4 * ((double) j + 1) * ((double) j + 1) - 1)) * cos(2 * M_PI * i * (j + 1) / ((double) n - 1));
      }
      w[i] = (2.0 / ((double) n - 1)) * (1.0 - (cos(M_PI * (double) i) / ((double) n * ((double) n - 2))) - 2 * s);
    }
  }

  return 0;

}

void update_index(int *index, int *multi_index, int dim)
{
  int i;

  for (i = 0; i < dim; i++)
  {
    if (index[i] < multi_index[i] - 1)
    {
      index[i] = index[i] + 1;
      break;
    } else
    {
      index[i] = 0;
    }
  }
}

long nchoosek(long n, long k)
{
  long r1 = 1, r2 = 1;
  long i;

  for (i = n; i > k; i--)
  {
    r1 = r1 * i;
  }
  for (i = n - k; i > 0; i--)
  {
    r2 = r2 * i;
  }
  return r1 / r2;
}

void enumerate_compositions(int d, int q, int *comps)
/*
 * Enumerate d-length compositions of k.
 */
{
  int i, j, k;
  int len, start, len_s;
  int loc_j, loc_k;
  int *loc_comps;

  if (q < d)
  {
    printf("ERROR: d-length compositions of q: q<d\n");
    return;
  }

  if (d == 1)
  {
    comps[0] = q;
  } else
  {
    len = nchoosek(q - 1, q - d);
    start = 0;
    for (i = 0; i < q - d + 1; i++)
    {
      len_s = nchoosek(q - i - 2, q - i - d);
      for (j = start; j < start + len_s; j++)
      {
        comps[j] = i + 1;
      }

      loc_comps = (int *) calloc(len_s * (d - 1), sizeof(int));

      enumerate_compositions(d - 1, q - i - 1, loc_comps);

      loc_j = 0;
      for (j = start; j < start + len_s; j++)
      {
        loc_k = 0;
        for (k = 1; k < d; k++)
        {
          comps[j + len * k] = loc_comps[loc_j + len_s * loc_k];
          loc_k++;
        }
        loc_j++;
      }
      free(loc_comps);
      start = start + len_s;
    }
  }

}
#endif

int get1DPointsAndWeights(const char *rule_name, int *level_ptr, double *p, double *w)
{

  char gauss_str[] = "GAUSS";
  int res;
  int level = *level_ptr;
  if (strcmp(gauss_str, rule_name) == 0)
  {
    res = gauss_1d_pts_wts(p, w, level);
  } else
  {
    res = clenshaw_curtis_1d_pts_wts(p, w, level);
  }
  for (int i = 0; i < level; i++)
    w[i] *= 0.5;
  return level;

}

int define_1d_points_and_weights(const char *rule_name, int *level_ptr)
{
  char gauss_str[] = "GAUSS";
  double *p, *w;
  int i, level, res;
  char p_filename[100];
  char w_filename[100];
  FILE *p_file;
  FILE *w_file;

  level = *level_ptr;

  p = (double *) calloc(level, sizeof(double));
  w = (double *) calloc(level, sizeof(double));

  if (strcmp(gauss_str, rule_name) == 0)
  {
    res = gauss_1d_pts_wts(p, w, level);
  } else
  {
    res = clenshaw_curtis_1d_pts_wts(p, w, level);
  }

  sprintf(p_filename, "%s%s%s%d%s", "points_", rule_name, "_l", level, "_d1.txt");
  sprintf(w_filename, "%s%s%s%d%s", "weights_", rule_name, "_l", level, "_d1.txt");
  p_file = fopen(p_filename, "w");
  w_file = fopen(w_filename, "w");
  for (i = 0; i < level; i++)
  {
    fprintf(p_file, "%18.16e\n", p[i]);
    fprintf(w_file, "%18.16e\n", 0.5 * w[i]);
  }
  fclose(p_file);
  fclose(w_file);

  free(p);
  free(w);
  return level;
}

int define_tensor_points_and_weights(char *rule_name, int *dim_ptr, int *level_ptr)
{
  char p_filename[100];
  char w_filename[100];
  char gauss_str[] = "gauss";
  double wt;
  double **one_d_pts, **one_d_wts;
  int dim, level, i, j, res, ngridpoints;
  int *multi_index, *index;
  FILE *p_file;
  FILE *w_file;

  dim = *dim_ptr;
  level = *level_ptr;

  multi_index = (int *) calloc(dim, sizeof(int));
  index = (int *) calloc(dim, sizeof(int));

  for (i = 0; i < dim; i++)
  {
    multi_index[i] = level;
  }
  ngridpoints = 1;
  for (i = 0; i < dim; i++)
  {
    ngridpoints = ngridpoints * multi_index[i];
  }

  //construct 1-D rules
  one_d_pts = (double **) calloc(dim, sizeof(double *));
  for (i = 0; i < dim; i++)
  {
    one_d_pts[i] = (double *) calloc(multi_index[i], sizeof(double));
  }
  one_d_wts = (double **) calloc(dim, sizeof(double *));
  for (i = 0; i < dim; i++)
  {
    one_d_wts[i] = (double *) calloc(multi_index[i], sizeof(double));
  }

  if (strcmp(gauss_str, rule_name) == 0)
  {
    for (i = 0; i < dim; i++)
    {
      res = gauss_1d_pts_wts(one_d_pts[i], one_d_wts[i], multi_index[i]);
    }
  } else
  {
    for (i = 0; i < dim; i++)
    {
      res = clenshaw_curtis_1d_pts_wts(one_d_pts[i], one_d_wts[i], multi_index[i]);
    }
  }

  for (i = 0; i < dim; i++)
  {
    index[i] = 0;
  }

  sprintf(p_filename, "%s%s%s%d%s%d%s", "points_", rule_name, "_tensor_l", level, "_d", dim, ".txt");
  sprintf(w_filename, "%s%s%s%d%s%d%s", "weights_", rule_name, "_tensor_l", level, "_d", dim, ".txt");
  p_file = fopen(p_filename, "w");
  w_file = fopen(w_filename, "w");
  for (i = 0; i < ngridpoints; i++)
  {

    wt = 1.0;
    for (j = 0; j < dim; j++)
    {
      fprintf(p_file, "%18.16e\t", one_d_pts[j][index[j]]);
      wt = 0.5 * wt * one_d_wts[j][index[j]];
    }
    update_index(index, multi_index, dim);
    fprintf(p_file, "\n");
    fprintf(w_file, "%18.16e\n", wt);
  }
  fclose(p_file);
  fclose(w_file);

  //FREE THE MEMORY!!!
  for (i = 0; i < dim; i++)
  {
    free(one_d_pts[i]);
  }
  free(one_d_pts);
  for (i = 0; i < dim; i++)
  {
    free(one_d_wts[i]);
  }
  free(one_d_wts);
  return ngridpoints;
}

int define_sparse_points_and_weights(char *rule_name, int *dim_ptr, int *level_ptr)
{
  char p_filename[100];
  char w_filename[100];
  char gauss_str[] = "gauss";
  double wt;
  double *sparse_wts, *pt;
  double **one_d_pts, **one_d_wts, **sparse_grid;
  int dim, level, i, j, k, ii, jj, q, res, len_tens, len_sparse = 1, len_comps;
  int n, one_norm_multi_index, plus_minus_one, count, same_point, found_same, same_index;
  int *multi_index, *index, *comps;
  long smolyak_factor, smolyak_size;
  FILE *p_file;
  FILE *w_file;

  dim = *dim_ptr;
  level = *level_ptr;

  multi_index = (int *) calloc(dim, sizeof(int));
  index = (int *) calloc(dim, sizeof(int));
  pt = (double *) calloc(dim, sizeof(double));

  smolyak_size = (int) pow(pow(2, level) + 1, dim);
  sparse_grid = (double **) calloc(smolyak_size, sizeof(double *));
  for (i = 0; i < smolyak_size; i++)
  {
    sparse_grid[i] = (double *) calloc(dim, sizeof(double));
  }
  sparse_wts = (double *) calloc(smolyak_size, sizeof(double));

  for (q = level + 1; q < level + dim + 1; q++)
  {

    len_comps = nchoosek(q - 1, q - dim);
    comps = (int *) calloc(len_comps * dim, sizeof(int));
    if (len_comps > 0)
    {
      enumerate_compositions(dim, q, comps);
    }

    for (i = 0; i < len_comps; i++)
    {

      len_tens = 1;
      one_norm_multi_index = 0;
      for (j = 0; j < dim; j++)
      {
        n = comps[i + len_comps * j];
        one_norm_multi_index = one_norm_multi_index + n;
        if (n == 1)
        {
          multi_index[j] = 1;
        } else
        {
          multi_index[j] = (int) pow(2, n - 1) + 1;
        }
        len_tens = len_tens * multi_index[j];
      }

      one_d_pts = (double **) calloc(dim, sizeof(double *));
      for (j = 0; j < dim; j++)
      {
        one_d_pts[j] = (double *) calloc(multi_index[j], sizeof(double));
      }
      one_d_wts = (double **) calloc(dim, sizeof(double *));
      for (j = 0; j < dim; j++)
      {
        one_d_wts[j] = (double *) calloc(multi_index[j], sizeof(double));
      }

      if (strcmp(gauss_str, rule_name) == 0)
      {
        for (j = 0; j < dim; j++)
        {
          res = gauss_1d_pts_wts(one_d_pts[j], one_d_wts[j], multi_index[j]);
        }
      } else
      {
        for (j = 0; j < dim; j++)
        {
          res = clenshaw_curtis_1d_pts_wts(one_d_pts[j], one_d_wts[j], multi_index[j]);
        }
      }

      // check each tensor point to see if it's in the grid.
      for (j = 0; j < dim; j++)
      {
        index[j] = 0;
      }
      for (j = 0; j < len_tens; j++)
      {
        wt = 1.0;
        for (k = 0; k < dim; k++)
        {
          pt[k] = one_d_pts[k][index[k]];
          wt = wt * 0.5 * one_d_wts[k][index[k]];
        }

        // build the smolyak factor
        plus_minus_one = 1;
        for (k = 0; k < level + dim - one_norm_multi_index; k++)
        {
          plus_minus_one = plus_minus_one * (-1);
        }
        smolyak_factor = plus_minus_one * nchoosek(dim - 1, level + dim - one_norm_multi_index);
        wt = wt * smolyak_factor;

        count = 1;
        found_same = 0;
        for (ii = 0; ii < len_sparse; ii++)
        {
          same_point = 1;
          for (jj = 0; jj < dim; jj++)
          {
            if (sparse_grid[ii][jj] != pt[jj])
            {
              same_point = 0;
            }
          }
          if (same_point)
          {
            found_same = 1;
            same_index = ii;
            break;
          }
        }
        if (found_same)
        {
          sparse_wts[same_index] = sparse_wts[same_index] + wt;
        } else
        {
          for (jj = 0; jj < dim; jj++)
          {
            sparse_grid[len_sparse][jj] = pt[jj];
          }
          sparse_wts[len_sparse] = wt;
          len_sparse++;
        }

        update_index(index, multi_index, dim);
      }

      // Free some memory
      for (j = 0; j < dim; j++)
      {
        free(one_d_pts[j]);
      }
      free(one_d_pts);
      for (j = 0; j < dim; j++)
      {
        free(one_d_wts[j]);
      }
      free(one_d_wts);

    }
    free(comps);

  }

  sprintf(p_filename, "%s%s%s%d%s%d%s", "points_", rule_name, "_sparse_l", level, "_d", dim, ".txt");
  sprintf(w_filename, "%s%s%s%d%s%d%s", "weights_", rule_name, "_sparse_l", level, "_d", dim, ".txt");
  p_file = fopen(p_filename, "w");
  w_file = fopen(w_filename, "w");
  for (i = 0; i < len_sparse; i++)
  {

    for (j = 0; j < dim; j++)
    {
      fprintf(p_file, "%18.16e\t", sparse_grid[i][j]);
    }
    fprintf(p_file, "\n");
    fprintf(w_file, "%18.16e\n", sparse_wts[i]);
  }
  fclose(p_file);
  fclose(w_file);

  for (i = 0; i < smolyak_size; i++)
  {
    free(sparse_grid[i]);
  }
  free(sparse_grid);
  free(sparse_wts);

  free(index);
  free(multi_index);
  free(pt);

  return len_sparse;
}




#endif /*STOCHASTICCOLLOCATION_H_*/

