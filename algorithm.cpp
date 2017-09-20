#include "window.h"
#include "algorithm.h"

#define EPS 1e-15
#define MAX_IT 5000
static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

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

    reduce_sum (t);
}

///u = u - tau*v
void sub_vect (int n, double *u, double tau,
               double *v, int k, int t) {

    int i1 = n*k/t, i2 = n*(k+1)/t;

    for (int i = i1; i < i2; ++i)
        u[i] -= tau * v[i];

    reduce_sum (t);
}

///(x,y)
double scalar_p (int n, double *x, double *y,
                 int k, int t) {

    int i1 = n*k/t, i2 = n*(k+1)/t;
    double s = 0;

    for (int i = i1; i < i2; i++)
        s += x[i]*y[i];

    return s;
}

///D*y = x => y = D^(-1)*x
void preconditioner (int n, double *x, double *d,
                     double *y, int k, int t) {

    int i1 = n*k/t, i2 = n*(k+1)/t;

    for (int i = i1; i < i2; i++)
        y[i] = x[i]/d[i];

    reduce_sum (t);
}

///Iterations
int solve (int n, double *A, int *I, double *b, double *x, double eps,
           int maxit, double *r, double *u, double *v, int k, int t) {

    double c1, c2, tau;

    static double C1 = 0, C2 = 0;

    /// r = Ax - b => r = Ax, r -= 1 * b
    matr_mult  (n, A, I, x, r, k, t);
    sub_vect   (n, r, -1, b, k, t);
    c1 = scalar_p (n, r, r, k, t);

    pthread_mutex_lock (&mutex);
    C1 += c1;
    pthread_mutex_unlock (&mutex);
    reduce_sum (t);

    if (fabs (C1) < eps*eps && k == 0)
        //return 0;
        printf ("\nError solve\n");

    for (int it = 1; it < maxit; it++) {
        // u = D^(-1) * r
        preconditioner (n, r, A, u, k, t);

        if (k == 0) {
            C1 = 0;
            C2 = 0;
        }

        // v = A * u
        matr_mult (n, A, I, u, v, k, t);

        c1 = scalar_p (n, v, r, k, t);
        c2 = scalar_p (n, v, v, k, t);

        pthread_mutex_lock (&mutex);
        C1 += c1;
        C2 += c2;
        pthread_mutex_unlock (&mutex);
        reduce_sum (t);

        if (fabs(C1) < eps*eps || fabs(C2) < eps*eps) {
            if (k == 0)
                printf ("iterations = %d, residual = %e\n", it, sqrt(C2));

            return it;
        }

        tau = C1/C2;
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
        memset (r, 0, (n*m - (d_x-1)*(d_y-1))*sizeof(double));
        memset (u, 0, (n*m - (d_x-1)*(d_y-1))*sizeof(double));
        memset (v, 0, (n*m - (d_x-1)*(d_y-1))*sizeof(double));
        memset (ar->x, 0, (n*m - (d_x-1)*(d_y-1))*sizeof(double));
        memset (ar->b, 0, (n*m - (d_x-1)*(d_y-1))*sizeof(double));
    }

    err = assemble_matrix (n, m, p, q, d_x, d_y, s,
                           ar->I, ar->A, k, t);

    if (err) {
        printf ("Cannot assemble matrix: error %d\n", err);

        return 0;
    }

    if (k == 0)
        fill_vector_b (n, m, p, q, d_x, d_y, s, ar->b, ar->func);

    reduce_sum (t);

    solve (n*m - (d_x-1)*(d_y-1), ar->A, ar->I, ar->b, ar->x,
           EPS, MAX_IT, r, u, v, k, t);

    reduce_sum (t);

    if (k == 0) {
        delete [] r;
        delete [] u;
        delete [] v;
    }

    return 0;
}









