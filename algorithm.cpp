#include "window.h"
#define EPS 1e-16
#define MAX_IT 5000


double Scene3D::f(double x, double y)
{
    double returnvalue;
    //returnvalue = 1;
    //returnvalue = exp(x * x - y * y);
    //returnvalue = sqrt(x * x + y * y);
    //returnvalue = x + y;
    returnvalue = sin(M_PI * x) * sin(M_PI * y);

    if(deltaEnabled)
    {
        if(fabs((x - x1) - ((x2 - x1) / 3)) < 0.05 && fabs((y - y1) - ((y2 - y1) / 3)) < 0.05)
        {
            returnvalue += deltaNums;
        }
        if(fabs((x - x1) - (2 * (x2 - x1) / 3)) < 0.05 && fabs((y - y1) - (2 * (y2 - y1) / 3)) < 0.05)
        {
            returnvalue -= deltaNums;
        }
    }

    return returnvalue;
}

void Scene3D::use_algorithm (void)
{
    static double *r, *u, *v;

    r = new double[size];
    u = new double[size];
    v = new double[size];
    memset(r, 0, size*sizeof(double));
    memset(u, 0, size*sizeof(double));
    memset(v, 0, size*sizeof(double));
    memset(x, 0, size*sizeof(double));
    memset(b, 0, size*sizeof(double));

    if(!assemble_matrix())
    {
        printf("Cannot assemble matrix\n");
        return;
    }

    fill_vector_b();
    solve(r, u, v);
    delete [] r;
    delete [] u;
    delete [] v;

    return;
}

///Iterations
int Scene3D::solve(double *r, double *u, double *v)
{
    double c1, c2;
    int start, len;
    double s;

    c1 = 0.;
    c2 = 0.;

    /// r = Ax - b => r = Ax, r -= 1 * b
    for (int i = 0; i < size; i ++)
    {
        r[i] += b[i];
    }

    for(int it = 1; it < MAX_IT; it ++)
    {
        // u = D^(-1) * r
        for(int i = 0; i < size; i ++)
        {
            u[i] = r[i] / A[i];
        }
        c1 = 0.;
        c2 = 0.;

        // v = A * u
        s = 0.;
        for(int i = 0; i < size; ++i)
        {
            s = A[i]*u[i];
            start = I[i];
            len = I[i + 1] - I[i];

            for(int j = 0; j < len; ++j)
            {
                s += A[start + j]*u[I[start + j]];
            }

            v[i] = s;
        }

        for(int i = 0; i < size; i ++)
        {
            c1 += v[i] * r[i];
            c2 += v[i] * v[i];
        }

        if(fabs(c1) < EPS*EPS || fabs(c2) < EPS*EPS)
        {
            return it;
        }

        for (int i = 0; i < size; i ++)
        {
            r[i] -= (c1 / c2) * v[i];
            x[i] += (c1 / c2) * u[i];
        }
    }

    return MAX_IT;
}

