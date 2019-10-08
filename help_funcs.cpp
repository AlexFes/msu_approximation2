#include "window.h"

void Scene3D::recount_algorithm()
{
    double hx, hy;
    int numPointsMiddle, numPointsEdge, numPointsCorners;

    numPointsMiddle = (n - 2) * (m - 2) * 6;
    numPointsEdge = ((n - 2) + (m - 2)) * 4 * 2;
    numPointsCorners = 2 * 3 + 2 * 2;
    hx = (x2 - x1) / (n - 1);
    hy = (y2 - y1) / (n - 1);

    size = n * m;
    s = hx * hy;
    numPoints = numPointsMiddle + numPointsEdge + numPointsCorners;

    update_arrays();
    get_points();
    use_algorithm();

    getVertexArray(VertexArray_appr, x, 0, 1);
    getVertexArray(VertexArray_residual, func, x, 2);

    updateGL();
}

///Get points of the rectangle
void Scene3D::get_points()
{
    int num;
    double hx = (x2-x1)/(n-1), hy = (y2-y1)/(m-1);

    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < m; ++j)
        {
            if(true)
            {
                num = get_num(i, j);
                points[num*2] = x1 + hx*i;
                points[num*2 + 1] = y1 + hy*j;
            }
        }
    }

    for(int i = 0; i < size; i ++)
    {
        func[i] = f(points[i*2], points[i*2+1]);
    }
}

///Push back a point
void Scene3D::fill_vertex_array(std::vector<GLfloat> &VertexArray, double *x, double *func, int index1, int index2, int index3)
{
    VertexArray.push_back(points[2*index1]);
    VertexArray.push_back(points[2*index1 + 1]);

    if(x)
    {
        VertexArray.push_back(fabs(func[index1] - x[index1]));
    }
    else
    {
        VertexArray.push_back(func[index1]);
    }

    VertexArray.push_back(points[2*index2]);
    VertexArray.push_back(points[2*index2 + 1]);

    if(x)
    {
        VertexArray.push_back(fabs(func[index2] - x[index2]));
    }
    else
    {
        VertexArray.push_back(func[index2]);
    }

    VertexArray.push_back(points[2*index3]);
    VertexArray.push_back(points[2*index3 + 1]);

    if(x)
    {
        VertexArray.push_back(fabs(func[index3] - x[index3]));
    }
    else
    {
        VertexArray.push_back(func[index3]);
    }
}

void Scene3D::fill_vertex_array_half_points(std::vector<GLfloat> &VertexArray, int status, int index1, int index2, int index3)
{
    double x1, x2;

    x1 = (points[2 * index1] + points[2 * index2]) / 2;
    x2 = (points[2 * index1 + 1] + points[2 * index2 + 1]) / 2;

    VertexArray.push_back(x1);
    VertexArray.push_back(x2);

    if(status == 0)
    {
        VertexArray.push_back(f(x1, x2));
    }
    else if(status == 1)
    {
        VertexArray.push_back((x[index1] + x[index2]) / 2);
    }
    else
    {
        VertexArray.push_back(fabs(f(x1, x2) - (x[index1] + x[index2]) / 2));
    }

    x1 = (points[2 * index1] + points[2 * index3]) / 2;
    x2 = (points[2 * index1 + 1] + points[2 * index3 + 1]) / 2;

    VertexArray.push_back(x1);
    VertexArray.push_back(x2);

    if(status == 0)
    {
        VertexArray.push_back(f(x1, x2));
    }
    else if(status == 1)
    {
        VertexArray.push_back((x[index1] + x[index3]) / 2);
    }
    else
    {
        VertexArray.push_back(fabs(f(x1, x2) - (x[index1] + x[index3]) / 2));
    }

    x1 = (points[2 * index2] + points[2 * index3]) / 2;
    x2 = (points[2 * index2 + 1] + points[2 * index3 + 1]) / 2;
    VertexArray.push_back(x1);
    VertexArray.push_back(x2);

    if(status == 0)
    {
        VertexArray.push_back(f (x1, x2));
    }
    else if(status == 1)
    {
        VertexArray.push_back((x[index2] + x[index3]) / 2);
    }
    else
    {
        VertexArray.push_back(fabs(f (x1, x2) - (x[index2] + x[index3]) / 2));
    }
}

///Get an array of points to draw
void Scene3D::getVertexArray(std::vector<GLfloat> &VertexArray, double *func, double *x, int status)
{
    VertexArray.clear();
    int index1, index2, index3;

    for(int i = 0; i < n-1; ++i)
    {
        for(int j = 0; j < m-1; ++j)
        {
            if(true)
            {
                index1 = get_num(i, j);
                index2 = get_num(i+1, j);
                index3 = get_num(i+1, j+1);
                fill_vertex_array(VertexArray, x, func,index1, index2, index3);
                fill_vertex_array_half_points(VertexArray, status, index1, index2, index3);

                index1 = get_num(i, j);
                index2 = get_num(i, j+1);
                index3 = get_num(i+1, j+1);
                fill_vertex_array(VertexArray, x, func, index1, index2, index3);
                fill_vertex_array_half_points(VertexArray, status, index1, index2, index3);
            }
        }
    }

    double max = 0;

    for(int i = 0; i < (int)VertexArray.size () / 9; i++)
    {
        for(int j = 0; j < 9; j += 3)
        {
            if(fabs(VertexArray[9 * i + 2 + j]) > max)
            {
                max = fabs(VertexArray[9 * i + 2 + j]);
            }
        }
    }

    if(status == 2)
    {
        for(int i = 0; i < (int)VertexArray.size () / 9; i++)
        {
            for(int j = 0; j < 9; j += 3)
            {
                VertexArray[9 * i + 2 + j]/=max;
            }
        }
        max_residual = max;
    }
}

///Free memory
template <typename T>
static inline void free_array (T *arr)
{
    if(arr)
    {
        delete[] arr;
    }
    arr = 0;
}

///Allocate memory
void Scene3D::update_arrays()
{
    free_array <double> (x);
    free_array <double> (func);
    free_array <double> (points);
    free_array <double> (A);
    free_array <int> (I);
    free_array <double> (b);

    x = new double[size];
    func = new double[size];
    points = new double[2 * size];
    A = new double[size + 1 + numPoints];
    I = new int[size + 1 + numPoints];
    b = new double[size];
}
