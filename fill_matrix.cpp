#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "window.h"

bool Scene3D::assemble_matrix(void)
{
    int num, k, l, rows, pos, w;
    int x[6], y[6];
    int i, j, nz;
    double a[7];

    rows = 0;
    pos = size + 1;

    for(i = 0; i < n; i ++)
    {
        for(j = 0; j < m; j ++)
        {
            I[rows] = pos;
            k = get_links(i, j, x, y);
            for(l = 0; l < k; l ++)
            {
                num = get_num(x[l], y[l]);
                I[pos + l] = num;
            }
            pos += k;
            rows ++;
        }
    }

    I[rows] = pos;

    for(l = 0; l < size; l++)
    {
        get_ij(i, j, l);
        nz = get_values(i, j, a);
        A[l] = a[0];

        if(I[l + 1] - I[l] != nz)
        {
            return false;
        }

        for(w = 0; w < nz; w++)
        {
            A[I[l] + w] = a[w + 1];
        }
    }

    if(I[size] != size + 1 + numPoints)
    {
        return false;
    }

    return true;
}

int Scene3D::get_links(int i, int j, int *x, int *y)
{
    if((i > 0 && i < n-1 && j > 0 && j < m-1))
    {
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

    if(i == 0 && j == 0)
    {
        x[0] = i;
        y[0] = j + 1;

        x[1] = i + 1;
        y[1] = j;

        x[2] = i + 1;
        y[2] = j + 1;

        return 3;
    }

    if(i == n-1 && j == m-1)
    {
        x[0] = i - 1;
        y[0] = j - 1;

        x[1] = i - 1;
        y[1] = j;

        x[2] = i;
        y[2] = j - 1;

        return 3;
    }

    if(i == n-1 && j == 0)
    {
        x[0] = i - 1;
        y[0] = j;

        x[1] = i;
        y[1] = j + 1;

        return 2;
    }

    if(i == 0 && j == m-1)
    {
        x[0] = i;
        y[0] = j - 1;

        x[1] = i + 1;
        y[1] = j;

        return 2;
    }

    if((i == 0 && j > 0 && j < m-1))
    {
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

    if((j == 0 && i > 0 && i < n-1))
    {
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

    if((i == n-1 && j > 0 && j < m-1))
    {
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

    if((j == m-1 && i > 0 && i < n-1))
    {
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

    return -1;
}

int Scene3D::get_num(int i, int j)
{
    return i*m + j;
}

int Scene3D::get_values(int i, int j, double *a)
{
    if((i > 0 && i < n-1 && j > 0 && j < m-1))
    {
        a[0] = s/2;
        a[1] = a[2] = a[3] = a[4] = a[5] = a[6] = s/12;

        return 6;
    }

    if(i == 0 && j == 0)
    {
        a[0] = s / 6;
        a[1] = s / 24;
        a[3] = s / 12;
        a[2] = s / 24;

        return 3;
    }

    if(i == n-1 && j == m-1)
    {
        a[0] = s / 6;
        a[3] = s / 24;
        a[1] = s / 12;
        a[2] = s / 24;

        return 3;
    }

    if(i == n-1 && j == 0)
    {
        a[0] = s / 12;
        a[1] = a[2] = s / 24;

        return 2;
    }

    if(i == 0 && j == m-1)
    {
        a[0] = s / 12;
        a[1] = a[2] = s / 24;

        return 2;
    }

    if(i == 0 && j > 0 && j < m-1)
    {
        a[0] = s/4;
        a[2] = s/24;
        a[4] = a[3] = s/12;
        a[1] = s/24;

        return 4;
    }

    if(i == n-1 && j > 0 && j < m-1)
    {
        a[0] = s/4;
        a[4] = a[3] = s/24;
        a[1] = a[2] = s/12;

        return 4;
    }

    if(j == 0 && i > 0 && i < n-1)
    {
        a[0] = s / 4;
        a[2] = a[4] = s / 12;
        a[3] = a[1] = s / 24;

        return 4;
    }

    if(j == m-1 && i > 0 && i < n-1)
    {
        a[0] = s / 4;
        a[4] = s / 24;
        a[3] = a[1] = s / 12;
        a[2] = s / 24;

        return 4;
    }

    return -1;
}

void Scene3D::get_ij(int &i, int &j, int l)
{
    i = l/m;
    j = l - i*m;
}

void Scene3D::fill_vector_b()
{
    int i, j, v;
    int x[6], y[6];

    ///Four triangles
    for(int l = 0; l < size; ++l)
    {
        get_ij(i, j, l);
        v = get_links(i, j, x, y);

        ///Locate all triangles
        for(int vertex = 0; vertex < v-1; ++vertex)
        {
            if(y[vertex+1] - y[vertex] == 1 || (y[vertex+1] - j == 1 && x[vertex+1]!=x[vertex]))
            {
                integral(i, j, x[vertex], y[vertex], x[vertex+1], y[vertex+1], 1);
            }

            if(vertex!=v-2 && (y[vertex+2] - y[vertex] == 1 || abs(y[vertex+2] - j) == 1))
            {
                integral(i, j, x[vertex], y[vertex], x[vertex+2], y[vertex+2], 0);
            }
        }
    }
}

void Scene3D::integral(int i1, int j1, int i2, int j2, int i3, int j3, int flag)
{
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
    if(flag && i2==i3)
    {
        num1 = get_num(i2, j2);
        num2 = get_num(i3, j3);

        if(i1 > i2)
        {
            x1 = points[num2*2];
            y1 = points[num2*2 + 1];

            x3 = points[num1*2];
            y3 = points[num1*2 + 1];
        }
        else
        {
            x1 = points[num1*2];
            y1 = points[num1*2 + 1];

            x3 = points[num2*2];
            y3 = points[num2*2 + 1];
        }

        num = get_num(i1, j1);
        x2 = points[num*2];
        y2 = points[num*2 + 1];

        f1 = f(x1, y1);
        f2 = f(x2, y2);
        f3 = f(x3, y3);
        g2 = 1;
        g5 = 0;
    }

    ///Second type
    if (!flag && j2==j3)
    {
        num1 = get_num(i2, j2);
        num2 = get_num(i3, j3);

        if (j1 > j2)
        {
            x1 = points[num2*2];
            y1 = points[num2*2 + 1];

            x2 = points[num1*2];
            y2 = points[num1*2 + 1];
        }
        else
        {
            x1 = points[num1*2];
            y1 = points[num1*2 + 1];

            x2 = points[num2*2];
            y2 = points[num2*2 + 1];
        }

        num = get_num(i1, j1);
        x3 = points[num*2];
        y3 = points[num*2 + 1];

        f1 = f(x1, y1);
        f2 = f(x2, y2);
        f3 = f(x3, y3);
        g3 = 1;
        g4 = 0;
    }

    ///Third type
    if((!flag && j2!=j3) || (flag && i2!=i3))
    {
        num1 = get_num(i2, j2);
        num2 = get_num(i3, j3);

        if(j1==j2)
        {
            x3 = points[num2*2];
            y3 = points[num2*2 + 1];

            x2 = points[num1*2];
            y2 = points[num1*2 + 1];
        }

        else
        {
            x3 = points[num1*2];
            y3 = points[num1*2 + 1];

            x2 = points[num2*2];
            y2 = points[num2*2 + 1];
        }

        num = get_num(i1, j1);
        x1 = points[num*2];
        y1 = points[num*2 + 1];

        f1 = f(x1, y1);
        f2 = f(x2, y2);
        f3 = f(x3, y3);
        g1 = 1;
        g6 = 0;
    }

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

    return;
}
