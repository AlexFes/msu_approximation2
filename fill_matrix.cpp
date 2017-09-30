#include "algorithm.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

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

        for (int w = 0; w < nz; w++)
            A[I[l] + w] = a[w + 1];
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
                    int d_x, int d_y, int k, int t, double s,
                    double *b, double (*f) (double, double),
                    double *points) {

    int i, j, v;

    int i1 = k*(n*m - (d_x-1)*(d_y-1))/t;
    int i2 = (k + 1)*(n*m - (d_x-1)*(d_y-1))/t;

    int x[6], y[6];

    ///Four triangles
    for (int l = i1; l < i2; ++l) {
        get_ij (i, j, m, p, q, d_x, d_y, l);
        v = get_links (i, j, n, m, p, q, d_x, d_y, x, y);

        ///Locate all triangles
        for (int vertex = 0; vertex < v-1; ++vertex) {
            if (y[vertex+1] - y[vertex] == 1 ||
                (y[vertex+1] - j == 1 && x[vertex+1]!=x[vertex]))
                integral (m, p, q, d_x, d_y, s, b, f, points,
                          i, j, x[vertex], y[vertex], x[vertex+1], y[vertex+1], 1);

            if (vertex!=v-2 && (y[vertex+2] - y[vertex] == 1 ||
                abs(y[vertex+2] - j) == 1))
                integral (m, p, q, d_x, d_y, s, b, f, points,
                          i, j, x[vertex], y[vertex], x[vertex+2], y[vertex+2], 0);
        }
    }


//    ///debug
//    printf("\n");
//    for (int i = 0; i < n*m - (d_x-1)*(d_y-1); ++i) printf(" %.2e ", b[i]);

}

///Calculate four integrals
void integral (int m, int p, int q,
               int d_x, int d_y, double s,
               double *b, double (*f) (double, double),
               double *points,
               int i1, int j1, int i2, int j2, int i3, int j3, int flag) {

    ///Main points
    double x1 = 0, x2 = 0, x3 = 0,
           y1 = 0, y2 = 0, y3 = 0;

    ///Middle points
    double x4, x5, x6,
           y4, y5, y6;

    ///Values
    double f1 = 0, f2 = 0, f3 = 0, f4, f5, f6,
           g1 = 0, g2 = 0, g3 = 0, g4 = 0.5, g5 = 0.5, g6 = 0.5;

    int num = 0, num1, num2;
                            ///Recognising the triangle
    ///First type
    if (flag && i2==i3) {
        num1 = get_num (i2, j2, m, p, q, d_x, d_y);
        num2 = get_num (i3, j3, m, p, q, d_x, d_y);

        if (i1 > i2) {
            x1 = points[num2*2];
            y1 = points[num2*2 + 1];

            x3 = points[num1*2];
            y3 = points[num1*2 + 1];
        }

        else {
            x1 = points[num1*2];
            y1 = points[num1*2 + 1];

            x3 = points[num2*2];
            y3 = points[num2*2 + 1];
        }

        num = get_num (i1, j1, m, p, q, d_x, d_y);
        x2 = points[num*2];
        y2 = points[num*2 + 1];

        f1 = f(x1, y1);
        f2 = f(x2, y2);
        f3 = f(x3, y3);
        g2 = 1;
        g5 = 0;
    }

    ///Second type
    if (!flag && j2==j3) {
        num1 = get_num (i2, j2, m, p, q, d_x, d_y);
        num2 = get_num (i3, j3, m, p, q, d_x, d_y);

        if (j1 > j2) {
            x1 = points[num2*2];
            y1 = points[num2*2 + 1];

            x2 = points[num1*2];
            y2 = points[num1*2 + 1];
        }

        else {
            x1 = points[num1*2];
            y1 = points[num1*2 + 1];

            x2 = points[num2*2];
            y2 = points[num2*2 + 1];
        }

        num = get_num (i1, j1, m, p, q, d_x, d_y);
        x3 = points[num*2];
        y3 = points[num*2 + 1];

        f1 = f(x1, y1);
        f2 = f(x2, y2);
        f3 = f(x3, y3);
        g3 = 1;
        g4 = 0;
    }

    ///Third type
    if ((!flag && j2!=j3) || (flag && i2!=i3)) {
        num1 = get_num (i2, j2, m, p, q, d_x, d_y);
        num2 = get_num (i3, j3, m, p, q, d_x, d_y);

        if (j1==j2) {
            x3 = points[num2*2];
            y3 = points[num2*2 + 1];

            x2 = points[num1*2];
            y2 = points[num1*2 + 1];
        }

        else {
            x3 = points[num1*2];
            y3 = points[num1*2 + 1];

            x2 = points[num2*2];
            y2 = points[num2*2 + 1];
        }

        num = get_num (i1, j1, m, p, q, d_x, d_y);
        x1 = points[num*2];
        y1 = points[num*2 + 1];

        f1 = f(x1, y1);
        f2 = f(x2, y2);
        f3 = f(x3, y3);
        g1 = 1;
        g6 = 0;
    }

    if (!flag && j2!=j3 && ((i1 == p+d_x && j1 == q && j3 > j1) ||
                            (i1 == p && j1 == q+d_y && j2 < j1)))
        return;

    x6 = x4 = (x1 + x2)/2;
    y4 = y1;

    x5 = x1;
    y6 = y5 = (y1 + y3)/2;

    f4 = f(x4, y4);
    f5 = f(x5, y5);
    f6 = f(x6, y6);

    b[num] += s*(f5*g5 + f6*g6 + f3*g3)/48 + s*(f5*(g6+g3) + f6*(g5+g3) + f3*(g5+g6))/96;
    b[num] += s*(f4*g4 + f2*g2 + f6*g6)/48 + s*(f4*(g2+g6) + f2*(g4+g6) + f6*(g4+g2))/96;
    b[num] += s*(f4*g4 + f1*g1 + f6*g6)/48 + s*(f4*(g1+g6) + f1*(g4+g6) + f6*(g4+g1))/96;
    b[num] += s*(f5*g5 + f6*g6 + f1*g1)/48 + s*(f5*(g6+g1) + f6*(g5+g1) + f1*(g5+g6))/96;

//    ///debug
//    printf ("\nnum = %d\n", num);

    return;
}

/*
///Vector b = (f, phi)
void fill_vector_b (int n, int m, int p, int q,
                    int d_x, int d_y, double s,
                    double *b, double *func) {

    int A, B, C;

   ///Single triangle
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
*/
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

//    ///debug
//    for (int i = 0; i < n*m - (d_x-1)*(d_y-1) + 1 +
//             ((n-2)*(m-2) - (d_x+1)*(d_y+1))*6
//             + ((n-2) + (d_x-1) + (m-2) + (d_y-1))*4*2
//             + 2*3 + 2*2 + 2*6 + 2*5; ++i) printf(" %d", I[i]);

    if (err)
        return -3;

    err = fill_matrix (n, m, p, q, d_x, d_y, s, I, A, k, t);

//    ///debug
//    printf("\n");
//    for (int i = 0; i < n*m - (d_x-1)*(d_y-1) + 1 +
//                 ((n-2)*(m-2) - (d_x+1)*(d_y+1))*6
//                 + ((n-2) + (d_x-1) + (m-2) + (d_y-1))*4*2
//                 + 2*3 + 2*2 + 2*6 + 2*5; ++i) printf(" %.2e ", A[i]);

    reduce_sum (t);

    if (err)
        return -4;

    return 0;
}












