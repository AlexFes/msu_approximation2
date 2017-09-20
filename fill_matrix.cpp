#include "algorithm.h"

///Close points on the plane
int get_links (int i, int j, int n, int m,
               int p, int q, int d_x, int d_y,
               int *x, int *y) {


    if ((i > 0 && i < n-1 && j > 0 && j < m-1) &&
        !(i >= p && i <= p+d_x && j >= q && j <= q+d_y)) {
        x[3] = i;
        x[5] = i + 1;
        x[4] = i + 1;
        x[2] = i;
        x[0] = i - 1;
        x[1] = i - 1;

        y[3] = j + 1;
        y[5] = j + 1;
        y[4] = j;
        y[2] = j - 1;
        y[0] = j - 1;
        y[1] = j;

        return 6;
    }

    if ((i == p+d_x && j == q) || (i == p && j == q+d_y)) {
        x[3] = i;
        x[5] = i + 1;
        x[4] = i + 1;
        x[2] = i;
        x[0] = i - 1;
        x[1] = i - 1;

        y[3] = j + 1;
        y[5] = j + 1;
        y[4] = j;
        y[2] = j - 1;
        y[0] = j - 1;
        y[1] = j;

        return 6;
    }

    if (i == 0 && j == 0) {
        x[0] = i;
        y[0] = j + 1;

        x[1] = i + 1;
        y[1] = j;

        x[2] = i + 1;
        y[2] = j + 1;

        return 3;
    }

    if (i == n-1 && j == m-1) {
        x[0] = i - 1;
        y[0] = j - 1;

        x[1] = i - 1;
        y[1] = j;

        x[2] = i;
        y[2] = j - 1;

        return 3;
    }

    if (i == n-1 && j == 0) {
        x[0] = i - 1;
        y[0] = j;

        x[1] = i;
        y[1] = j + 1;

        return 2;
    }

    if (i == 0 && j == m-1) {
        x[0] = i;
        y[0] = j - 1;

        x[1] = i + 1;
        y[1] = j;

        return 2;
    }

    if ((i == 0 && j > 0 && j < m-1) ||
       (i == p+d_x && j > q && j < q+d_y)) {
        x[1] = i;
        y[1] = j + 1;

        x[3] = i + 1;
        y[3] = j + 1;

        x[2] = i + 1;
        y[2] = j;

        x[0] = i;
        y[0] = j - 1;

        return 4;
    }

    if ((j == 0 && i > 0 && i < n-1) ||
       (j == q+d_y && i > p && i < p+d_x)) {
        x[0] = i - 1;
        y[0] = j;

        x[1] = i;
        y[1] = j + 1;

        x[3] = i + 1;
        y[3] = j + 1;

        x[2] = i + 1;
        y[2] = j;

        return 4;
    }

    if ((i == n-1 && j > 0 && j < m-1) ||
       (i == p && j > q && j < q+d_y)) {
        x[3] = i;
        y[3] = j + 1;

        x[2] = i;
        y[2] = j - 1;

        x[0] = i - 1;
        y[0] = j - 1;

        x[1] = i - 1;
        y[1] = j;

        return 4;
    }

    if ((j == m-1 && i > 0 && i < n-1) ||
       (j == q && i > p && i < p+d_x)) {
        x[3] = i + 1;
        y[3] = j;

        x[2] = i;
        y[2] = j - 1;

        x[0] = i - 1;
        y[0] = j - 1;

        x[1] = i - 1;
        y[1] = j;

        return 4;
    }

    if (i == p && j == q) {
        x[3] = i;
        y[3] = j + 1;

        x[4] = i + 1;
        y[4] = j;

        x[2] = i;
        y[2] = j - 1;

        x[0] = i - 1;
        y[0] = j - 1;

        x[1] = i - 1;
        y[1] = j;

        return 5;
    }

    if (i == p+d_x && j == q+d_y) {
        x[2] = i;
        y[2] = j + 1;

        x[4] = i + 1;
        y[4] = j + 1;

        x[3] = i + 1;
        y[3] = j;

        x[1] = i;
        y[1] = j - 1;

        x[0] = i - 1;
        y[0] = j;

        return 5;
    }

    return -1;
}

///Scalar product
int get_values (int i, int j, int n, int m,
                int p, int q, int d_x, int d_y,
                double s, double *a) {

    if ((i > 0 && i < n-1 && j > 0 && j < m-1) &&
        !(i >= p && i <= p+d_x && j >= q && j <= q+d_y)) {
        a[0] = s/2;
        a[1] = a[2] = a[3] = a[4] = a[5] = a[6] = s/12;

        return 6;
    }

    if (i == p+d_x && j == q) {
        a[0] = s*5/12;
        a[6] = a[5] = a[3] = a[1] = s/12;
        a[4] = a[2] = s/24;

        return 6;
    }

    if (i == p && j == q+d_y) {
        a[0] = s*5/12;
        a[4] = a[6] = a[1] = a[2] = s/12;
        a[5] = a[3] = s/24;

        return 6;
    }

    if (i == 0 && j == 0) {
        a[0] = s / 6;
        a[1] = s / 24;
        a[3] = s / 12;
        a[2] = s / 24;

        return 3;
    }

    if (i == n-1 && j == m-1) {
        a[0] = s / 6;
        a[3] = s / 24;
        a[1] = s / 12;
        a[2] = s / 24;

        return 3;
    }

    if (i == n-1 && j == 0) {
        a[0] = s / 12;
        a[1] = a[2] = s / 24;

        return 2;
    }

    if (i == 0 && j == m-1) {
        a[0] = s / 12;
        a[1] = a[2] = s / 24;

        return 2;
    }

    if ((i == 0 && j > 0 && j < m-1) ||
       (i == p+d_x && j > q && j < q+d_y)) {
        a[0] = s/4;
        a[2] = s/24;
        a[4] = a[3] = s/12;
        a[1] = s/24;

        return 4;
    }

    if ((i == n-1 && j > 0 && j < m-1) ||
       (i == p && j > q && j < q+d_y)) {
        a[0] = s/4;
        a[4] = a[3] = s/24;
        a[1] = a[2] = s/12;

        return 4;
    }

    if ((j == 0 && i > 0 && i < n-1) ||
       (j == q+d_y && i > p && i < p+d_x)) {
        a[0] = s / 4;
        a[2] = a[4] = s / 12;
        a[3] = a[1] = s / 24;

        return 4;
    }

    if ((j == m-1 && i > 0 && i < n-1) ||
       (j == q && i > p && i < p+d_x)) {
        a[0] = s / 4;
        a[4] = s / 24;
        a[3] = a[1] = s / 12;
        a[2] = s / 24;

        return 4;
    }

    if (i == p && j == q) {
        a[0] = s / 3;
        a[4] = a[5] = s / 24;
        a[3] = a[1] = a[2] = s / 12;

        return 5;
    }

    if (i == p+d_x && j == q+d_y) {
        a[0] = s / 3;
        a[2] = a[1] = s / 24;
        a[3] = a[5] = a[4] = s / 12;

        return 5;
    }

    return -1;
}

///Index in the matrix
int get_num (int i, int j, int m,
             int p, int q, int d_x, int d_y) {

    if (i > p && i < p+d_x)
        return (j > q ? i*m + j - (i-p)*(d_y-1) :
               i*m + j - (i-p-1)*(d_y-1));

    if (i >= p+d_x)
        return i*m + j - (d_x-1)*(d_y-1);

    return i*m + j;
}

///Index on the plane
void get_ij (int &i, int &j, int m,
             int p, int q, int d_x, int d_y, int l) {

    int first_approach, step, max_shift, definitive_shift;

    first_approach = (p + 1)*m + q + 1;     //first encounter of spaces

    if (l < first_approach) {               //the cut not encountered yet
        i = l/m;
        j = l - i*m;
    }

    else {
        step = m - d_y + 1;                 //steps before next encounter
        max_shift = (d_y - 1)*(d_x - 1);    //all points in the cut
        definitive_shift = (d_y - 1)*(1 + (l - first_approach)/step);

        if (definitive_shift > max_shift) definitive_shift = max_shift;

        l += definitive_shift;              //counting the definitive shift
        i = l/m;
        j = l - i*m;
    }
}

///Index MSR matrix
int fill_matrix_structure (int n, int m, int p, int q,
                           int d_x, int d_y, int *I) {

    int num, k, l, rows = 0, pos = (n*m - (d_x-1)*(d_y-1)) + 1;
    int x[6], y[6];

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            if (!(i > p && i < p + d_x && j > q && j < q + d_y)) {
                I[rows] = pos;
                k = get_links (i, j, n, m, p, q, d_x, d_y, x, y);

                for (l = 0; l < k; l++) {
                    num = get_num (x[l], y[l], m, p, q, d_x, d_y);
                    I[pos + l] = num;
                }

                pos += k;
                ++rows;
            }

    I[rows] = pos;

    return 0;
}

///Value MSR matrix
int fill_matrix (int n, int m, int p, int q,
                 int d_x, int d_y, double s, int *I,
                 double *A, int k, int t) {

    int i, j, nz, err = 0;
    double a[7];
    int i1 = k*(n*m - (d_x-1)*(d_y-1))/t;
    int i2 = (k + 1)*(n*m - (d_x-1)*(d_y-1))/t;

    for (int l = i1; l < i2; l++) {
        get_ij (i, j, m, p, q, d_x, d_y, l);
        nz = get_values (i, j, n, m, p, q, d_x, d_y, s, a);
        A[l] = a[0];

        if (I[l + 1] - I[l] != nz) {
            err = 1;
            break;
        }

        for (int w = 0; w < nz; w++) {
            A[I[l] + w] = a[w + 1];}
    }

    if (I[n*m - (d_x-1)*(d_y-1)] != n*m - (d_x-1)*(d_y-1) + 1 +
            ((n-2)*(m-2) - (d_x+1)*(d_y+1))*6
            + ((n-2) + (d_x-1) + (m-2) + (d_y-1))*4*2
            + 2*3 + 2*2 + 2*6 + 2*5)
        err++;

    return err;
}

///Vector b = (f, phi)
void fill_vector_b (int n, int m, int p, int q,
                    int d_x, int d_y, double s,
                    double *b, double *func) {

    int A, B, C;

    for (int i = 0; i < n-1; ++i)
        for (int j = 0; j < m-1; ++j)
            if (!(i > (p-1) && i < p + d_x
               && j > (q-1) && j < q + d_y)) {
                ///The bottom triangle
                B = get_num (i, j, m, p, q, d_x, d_y);
                A = get_num (i+1, j, m, p, q, d_x, d_y);
                C = get_num (i+1, j+1, m, p, q, d_x, d_y);

                b[A] += (2*func[A] +   func[B] +   func[C])*s/24;
                b[B] += (  func[A] + 2*func[B] +   func[C])*s/24;
                b[C] += (  func[A] +   func[B] + 2*func[C])*s/24;

                ///The top triangle
                C = get_num (i, j, m, p, q, d_x, d_y);
                A = get_num (i, j+1, m, p, q, d_x, d_y);
                B = get_num (i+1, j+1, m, p, q, d_x, d_y);

                b[A] += (2*func[A] +   func[B] +   func[C])*s/24;
                b[B] += (  func[A] + 2*func[B] +   func[C])*s/24;
                b[C] += (  func[A] +   func[B] + 2*func[C])*s/24;
            }
}

///Fill MSR matrix
int assemble_matrix (int n, int m, int p, int q, int d_x, int d_y,
                     double s, int *I, double *A, int k, int t) {

    int NZ = ((n-2)*(m-2) - (d_x+1)*(d_y+1))*6
            + ((n-2) + (d_x-1) + (m-2) + (d_y-1))*4*2
            + 2*3 + 2*2 + 2*6 + 2*5;
    int err = 0;

    if (NZ < 0)
        return -1;

    if (k == 0)
        err = fill_matrix_structure(n, m, p, q, d_x, d_y, I);

    reduce_sum (t);

    if (err)
        return -3;

    err = fill_matrix (n, m, p, q, d_x, d_y, s, I, A, k, t);

    reduce_sum (t);

    if (err)
        return -4;

    return 0;
}












