#include "window.h"
#include "algorithm.h"

#define EPS 1e-15
#define MAX_IT 5000
/*
void get_my_rows (int n, int k, int p, int *i1, int *i2)
{
  *i1 = n * k;
  *i1 /= p;
  *i2 = n * (k + 1);
  *i2 = *i2 / p - 1;
}

int get_nz_matrix (int n, int m)
{
  return   (n - 1) * (m - 1) * 6 // middle points
         + (n - 1) * 4 * 2 // vertical side points
         + (m - 1) * 4 * 2 // horizontal side points
         + 2 * 3 + 2 * 2; // vertex points
}

void matr_mult (double *a, int *jnz, int n, double *x, double *b, int k, int p)
{
  int i1, i2, i, start, len, j, jj;
  double s = 0;

  get_my_rows (n, k, p, &i1, &i2);
  for (i = i1; i <= i2; i++)
    {
      s = a[i] * x[i];
      start = jnz[i];
      len = jnz[i + 1] - jnz[i];

      for (j = 0; j < len; j++)
        {
          jj = start + j;
          s += a[jj] * x[jnz[jj]];
        }
      b[i] = s;
    }

  reduce_sum <int>(p);
}

// u = u - tau * v
void sub_vect (int n, double *u, double tau, double *v, int k, int p)
{
  int i1, i2, i;

  get_my_rows (n, k, p, &i1, &i2);
  for (i = i1; i <= i2; i++)
    u[i] -= tau * v[i];

  reduce_sum <int> (p);
}

double scalp (double *x, double *y, int n, int k, int p)
{
  int i1, i2, i;
  double s = 0;
  get_my_rows (n, k, p, &i1, &i2);
  for (i = i1; i <= i2; i++)
    s += x[i] * y[i];
  reduce_sum <double> (p, &s, 1);
  return s;
}

// D * x = y => x = D^(-1) * y
void preconditioner (double *d, int n, double *y, double *x, int k, int p)
{
  int i1, i2, i;
  get_my_rows (n, k, p, &i1, &i2);
  for (i = i1; i <= i2; i++)
    x[i] = y[i] / d[i];
  reduce_sum <int> (p);
}


int solve (double *a, int *jnz, int n, double *b, double *x, double eps,
           int maxit, double *r, double *u, double *v, int k, int p)
{
  double c1, c2, tau;
  int it;

  // r = Ax - b => r = Ax, r -= 1 * b
  matr_mult  (a, jnz, n, x, r, k, p);
  sub_vect   (n, r, 1, b, k, p);
  c1 = scalp (r, r, n, k, p);

  if (fabs (c1) < eps * eps)
    {
      return 0;
    }

  for (it = 1; it < maxit; it++)
    {
      // u = D^(-1) * r
      preconditioner (a, n, r, u, k, p);

      // v = A * u
      matr_mult (a, jnz, n, u, v, k, p);

      c1 = scalp (v, r, n, k, p);
      c2 = scalp (v, v, n, k, p);

      if (fabs (c1) < eps * eps || fabs (c2) < eps * eps)
        {
          if (k == 0)
            printf ("iterations = %d, residual = %e\n", it, sqrt (c2));
          return it;
        }
      tau = c1 / c2;
      sub_vect (n, r, tau, v, k, p);
      sub_vect (n, x, tau, u, k, p);
    }

  printf ("iterations = %d, GOT MAX ITERATIONS\n", maxit);
  return maxit;
}

void *use_algorithm (void *arg)
{
  args *ar = (args *) arg;
  static double *ce;

  int err = 0;
  static double *r, *u, *v;
  if (ar->k == 0)
    {
      r = new double  [ar->P];
      u = new double  [ar->P];
      v = new double  [ar->P];
      ce = new double [ar->P];
      memset (r, 0, ar->P * sizeof(double));
      memset (u, 0, ar->P * sizeof(double));
      memset (v, 0, ar->P * sizeof(double));
      memset (ce, 0, ar->P * sizeof(double));
    }

  int N = 0;
  err = assemble_matrix (ar->N2, ar->N2, &N, &ar->jnz, &ar->a, ar->k, ar->p);
  reduce_sum<int> (ar->p, &err, 1);
  if (err)
    {
      if (ar->k == 0)
        printf ("Cannot assemble matrix: error %d\n", err);
      return 0;
    }

  matr_mult (ar->a, ar->jnz, N, ar->b, ce, ar->k, ar->p);
  reduce_sum <int> (ar->p);

  solve (ar->a, ar->jnz, N, ce, ar->x, EPS, MAX_IT, r, u, v, ar->k, ar->p);
  reduce_sum<int> (ar->p, &err, 1);

  if (ar->k == ar->p - 1)
    {
      delete [] r;
      delete [] u;
      delete [] v;
      delete [] ce;
    }
  return 0;
}
*/


                                ///MA MAN
///A*x
void matr_mult (int n, double *A, int *I, double *x,
                double *b, int k, int t) {

    int i1 = n*k/t, i2 = n*(k+1)/t, start, len;
    double s = 0;

    for (int i = i1; i < i2; ++i) {
        s = A[i]*x[i];
        start = I[i];
        len = I[i + 1] - I[i];

        for (int j = 0; j < len; ++j)
            s += A[start + j]*x[I[start + j]];

        b[i] = s;
    }

    //reduce_sum <int> (t);
}

///u = u - tau*v
void sub_vect (int n, double *u, double tau,
               double *v, int k, int t) {

    int i1 = n*k/t, i2 = n*(k+1)/t;

    for (int i = i1; i < i2; ++i)
        u[i] -= tau * v[i];

    //reduce_sum <int> (t);
}

///(x,y)
double scalar_p (int n, double *x, double *y,
                 int k, int t) {

    int i1 = n*k/t, i2 = n*(k+1)/t;
    double s = 0;

    for (int i = i1; i < i2; i++)
        s += x[i]*y[i];

    //reduce_sum <double> (t, &s, 1);

    return s;
}

///D*y = x => y = D^(-1)*x
void preconditioner (int n, double *x, double *d,
                     double *y, int k, int t) {

    int i1 = n*k/t, i2 = n*(k+1)/t;

    for (int i = i1; i < i2; i++)
        y[i] = x[i]/d[i];

    //reduce_sum <int> (t);
}

///Iterations
int solve (int n, double *A, int *I, double *b, double *x, double eps,
           int maxit, double *r, double *u, double *v, int k, int t) {

    double c1, c2, tau;

    // r = Ax - b => r = Ax, r -= 1 * b
    matr_mult  (n, A, I, x, r, k, t);

    //for (int i = 0; i < n; ++i) printf(" r[%d]=%f ", i, r[i]);

    sub_vect   (n, r, -1, b, k, t);

    //for (int i = 0; i < n; ++i) printf(" r[%d]=%f ", i, r[i]);

    c1 = scalar_p (n, r, r, k, t);

    if (fabs (c1) < eps*eps)
        return 0;

    for (int it = 1; it < maxit; it++) {
        // u = D^(-1) * r
        preconditioner (n, r, A, u, k, t);

        // v = A * u
        matr_mult (n, A, I, u, v, k, t);

        c1 = scalar_p (n, v, r, k, t);
        c2 = scalar_p (n, v, v, k, t);

        if (fabs(c1) < eps*eps || fabs(c2) < eps*eps) {
            if (k == 0)
                printf ("iterations = %d, residual = %e\n", it, sqrt (c2));

            //for (int i = 0; i < n; ++i) printf(" x[%d]=%f ", i, x[i]);

            return it;
        }

        tau = c1/c2;
        sub_vect (n, r, tau, v, k, t);
        sub_vect (n, x, -tau, u, k, t);
    }

    printf ("iterations = %d, GOT MAX ITERATIONS\n", maxit);

    return maxit;
}

///Multithread function
void *use_algorithm (void *arg) {

    args *ar = (args *) arg;

    static double *r, *u, *v;

    int n = ar->n, m = ar->m, p = ar->p,
        q = ar->q, d_x = ar->d_x, d_y = ar->d_y,
        k = ar->k, t = ar->t, err = 0;

    double x1 = ar->x1, y1 = ar->y1,
           x2 = ar->x2, y2 = ar->y2;

    double hx = (x2-x1)/(n-1), hy = (y2-y1)/(m-1), s = hx*hy;

    if (k == 0) {
        r = new double  [n*m - (d_x-1)*(d_y-1)];
        u = new double  [n*m - (d_x-1)*(d_y-1)];
        v = new double  [n*m - (d_x-1)*(d_y-1)];
//        ce = new double [n*m];
        memset (r, 0, (n*m - (d_x-1)*(d_y-1))*sizeof(double));
        memset (u, 0, (n*m - (d_x-1)*(d_y-1))*sizeof(double));
        memset (v, 0, (n*m - (d_x-1)*(d_y-1))*sizeof(double));
//        memset (ce, 0, n*m*sizeof(double));
        memset (ar->x, 0, (n*m - (d_x-1)*(d_y-1))*sizeof(double));
        memset (ar->b, 0, (n*m - (d_x-1)*(d_y-1))*sizeof(double));
    }

    err = assemble_matrix (n, m, p, q, d_x, d_y, s,
                           ar->I, ar->A, k, t);
    //reduce_sum <int> (t, &err, 1);

//    for (int i = n*m - (d_x-1)*(d_y-1) + 1; i < n*m - (d_x-1)*(d_y-1) + 1 +
//         ((n-2)*(m-2) - (d_x+1)*(d_y+1))*6
//         + ((n-2) + (d_x-1) + (m-2) + (d_y-1))*4*2
//         + 2*3 + 2*2 + 2*6 + 2*5; ++i) printf(" A[%d]=%.2f ", i, *(ar->A + i));

//    printf("\n");

//    for (int i = n*m - (d_x-1)*(d_y-1) + 1; i < n*m - (d_x-1)*(d_y-1) + 1 +
//         ((n-2)*(m-2) - (d_x+1)*(d_y+1))*6
//         + ((n-2) + (d_x-1) + (m-2) + (d_y-1))*4*2
//         + 2*3 + 2*2 + 2*6 + 2*5; ++i) {

//        if (*(ar->I + i)/10 < 1) printf(" I[%d]=%d    ", i, *(ar->I + i));

//        else printf(" I[%d]=%d   ", i, *(ar->I + i));
//    }

    if (err) {
        if (k == 0)
            printf ("Cannot assemble matrix: error %d\n", err);

        return 0;
    }
/*
    matr_mult (ar->a, ar->jnz, N, ar->b, ce, ar->k, ar->p);
    reduce_sum <int> (ar->p);
*/

    fill_vector_b (n, m, p, q, d_x, d_y, s, ar->b, ar->func);

    solve (n*m - (d_x-1)*(d_y-1), ar->A, ar->I, ar->b, ar->x,
           EPS, MAX_IT, r, u, v, k, t);

    //reduce_sum <int> (t, &err, 1);

    if (k == t - 1) {
        delete [] r;
        delete [] u;
        delete [] v;
//        delete [] ce;
    }

    return 0;
}









