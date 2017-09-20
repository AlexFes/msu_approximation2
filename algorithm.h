#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <vector>
#include <pthread.h>

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
void reduce_sum (int p);


#endif // ALGORITHM_H

