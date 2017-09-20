#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <vector>
#include <pthread.h>

//struct args
//{
//  double *a, *b, *x;
//  std::vector<double> points;
//  double vertices[8/*vertex_pos_COUNT*/];
//  int *jnz;
//  int N2, P, p, k;
//  double (*f)(double, double);
//};
///ma man
struct args {

    double *A, *b, *x, *func;
    int *I;
    int x1, y1, x2, y2, n, m, p, q, d_x, d_y, t, k;
    double (*f)(double, double);
};

void matr_mult (int n, double *A, int *I, double *x,
                double *b, int k, int t);
void sub_vect (int n, double *u, double tau,
               double *v, int k, int t);
void preconditioner (int n, double *x, double *d,
                     double *y, int k, int t);
double scalar_p (int n, double *x, double *y,
                 int k, int t);
int solve (int n, double *A, int *I, double *b, double *x, double eps,
           int maxit, double *r, double *u, double *v, int k, int t);
void *use_algorithm (void *arg);

int get_num (int i, int j, int m,
             int p, int q, int d_x, int d_y);
int get_values (int i, int j, int n, int m,
                int p, int q, int d_x, int d_y,
                double s, double *a);
int fill_matrix (int n, int m, int p, int q,
                 int d_x, int d_y, double s, int *I,
                 double *A, int k, int t);
int fill_matrix_structure (int n, int m, int p, int q,
                           int d_x, int d_y, int *I);
void fill_vector_b (int n, int m, int p, int q,
                    int d_x, int d_y, double s,
                    double *b, double *func);
int get_links (int i, int j, int n, int m,
               int p, int q, int d_x, int d_y,
               int *x, int *y);
void get_ij (int &i, int &j, int m,
             int p, int q, int d_x, int d_y, int l);
int assemble_matrix (int n, int m, int p, int q, int d_x, int d_y,
                     double s, int *I, double *A, int k, int t);
/*
template <typename T>
void reduce_sum (int p, T *a = 0, int n = 0) {

    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;

    static int t_in = 0;
    static int t_out = 0;
    static T *p_a = 0;

    int i;

    if (p <= 1)
        return;

    pthread_mutex_lock (&m);

    if (!p_a)
        p_a = a;

    else if (a)
        for (i = 0; i < n; i++)
            p_a[i] += a[i];

    t_in++;

    if (t_in < p)
        while (t_in < p)
        pthread_cond_wait (&c_in, &m);

    else {
        t_out = 0;
        pthread_cond_broadcast (&c_in);
    }

    if (p_a != a)
        for (i = 0; i < n; i++)
            a[i] = p_a[i];

    t_out++;

    if (t_out < p)
        while (t_out < p)
            pthread_cond_wait (&c_out, &m);

    else {
        p_a = 0, t_in = 0;
        pthread_cond_broadcast (&c_out);
    }

    pthread_mutex_unlock (&m);
}
*/
#endif // ALGORITHM_H

