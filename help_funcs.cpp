#include "window.h"
#include "algorithm.h"
//#include <pthread.h>
//#include <sys/time.h>
//#include <sys/resource.h>

///Multithread
void Scene3D::recount_algorithm () {

    //get_cut(0.4, 0.3, 0.4, 0.5);

    pthread_t tid;

    update_arrays();

    get_points ();

    args *arg = new args [t];

    for (int i = 0; i < t; i++) {
        arg[i].A = A;
        arg[i].I = I;
        arg[i].b = b;
        arg[i].x = x;
        arg[i].func = func;

        arg[i].x1 = x1;
        arg[i].x2 = x2;
        arg[i].y1 = y1;
        arg[i].y2 = y2;

        arg[i].n = n;
        arg[i].m = m;
        arg[i].p = p;
        arg[i].q = q;
        arg[i].d_x = d_x;
        arg[i].d_y = d_y;

        arg[i].k = i;
        arg[i].t = t;
    }

    for (int i = 1; i < t; i++)
        pthread_create (&tid, 0, use_algorithm, arg + i);

    use_algorithm (arg + 0);

//  getVertexArray (VertexArray_real, func, 0, 0);
    getVertexArray (VertexArray_appr, x, 0, 1);
    getVertexArray (VertexArray_residual, func, x, 2);

    delete[] arg;

    updateGL();
}

///Get points of the rectangle
void Scene3D::get_points () {

    int num;
    double hx = (x2-x1)/(n-1), hy = (y2-y1)/(m-1);

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            if (!(i > p && i < p + d_x && j > q && j < q + d_y)) {
                num = get_num (i, j, m, p, q, d_x, d_y);

                points[num*2] = x1 + hx*i;
                points[num*2 + 1] = y1 + hy*j;
            }

    for (int i = 0; i < n*m - (d_x-1)*(d_y-1); i++)
        func[i] = f(points[i*2], points[i*2+1]);
}

///Push back a point
void Scene3D::fill_vertex_array (std::vector<GLfloat> &VertexArray,
     double *x, double *func, int index1, int index2, int index3) {

    VertexArray.push_back (points[2*index1]);
    VertexArray.push_back (points[2*index1 + 1]);

    if (x)
        VertexArray.push_back (1e+14*fabs(func[index1] - x[index1]));

    else
        VertexArray.push_back (func[index1]);

    VertexArray.push_back (points[2*index2]);
    VertexArray.push_back (points[2*index2 + 1]);

    if (x)
        VertexArray.push_back (1e+14*fabs(func[index2] - x[index2]));

    else
        VertexArray.push_back (func[index2]);

    VertexArray.push_back (points[2*index3]);
    VertexArray.push_back (points[2*index3 + 1]);

    if (x)
        VertexArray.push_back (1e+14*fabs(func[index3] - x[index3]));

    else
        VertexArray.push_back (func[index3]);
}

///Get an array of points to draw
void Scene3D::getVertexArray (std::vector<GLfloat> &VertexArray,
                            double *func, double *x, int status) {

    VertexArray.clear ();
    int index1, index2, index3;

    for (int i = 0; i < n-1; ++i)
        for (int j = 0; j < m-1; ++j)
            if (!(i > (p-1) && i < p + d_x
               && j > (q-1) && j < q + d_y)) {
                index1 = get_num (i, j, m, p, q, d_x, d_y);
                index2 = get_num (i+1, j, m, p, q, d_x, d_y);
                index3 = get_num (i+1, j+1, m, p, q, d_x, d_y);
                fill_vertex_array (VertexArray, x, func,
                                   index1, index2, index3);

                index1 = get_num (i, j, m, p, q, d_x, d_y);
                index2 = get_num (i, j+1, m, p, q, d_x, d_y);
                index3 = get_num (i+1, j+1, m, p, q, d_x, d_y);
                fill_vertex_array (VertexArray, x, func,
                                   index1, index2, index3);
            }

    double max = 0;

    for (int i = 0; i < (int)VertexArray.size () / 9; i++)
        for (int j = 0; j < 9; j += 3)
            if (fabs (VertexArray[9 * i + 2 + j]) > max)
                max = fabs (VertexArray[9 * i + 2]);

    if (status == 2) {
        max_residual = max;
        printf ("Max residual = %e\n", 1e-14*max_residual);
        //qDebug() << "Max residual = " << max_residual;
    }
}

///Get the cut
void Scene3D::get_cut(double X, double Y, double dX, double dY) {

    double hx = (x2-x1)/(n-1), hy = (y2-y1)/(m-1);

    p = (int) trunc ((X - x1)/hx);
    if (p == 0)
        ++p;
    else if (p >= n-2)
        p = n-3;

    q = (int) trunc ((Y - y1)/hy);
    if (q == 0)
        q++;
    else if (q >= m-2)
        q = m-3;

    d_x = (int) trunc (dX/hx);
    if (p + d_x >= n - 1)
        d_x = n - 2 - p;
    else if (d_x == 0)
        d_x = 1;

    d_y = (int) trunc (dY/hy);
    if (q + d_y >= m - 1)
        d_y = m - 2 - q;
    else if (d_y == 0)
        d_y = 1;
}

///Free memory
template <typename T>
static inline void free_array (T *arr) {

    if (arr)
        delete[] arr;

    arr = 0;
}

///Allocate memory
void Scene3D::update_arrays () {

    free_array <double> (x);
    free_array <double> (func);
    free_array <double> (points);
    free_array <double> (A);
    free_array <int> (I);
    free_array <double> (b);

    x = new double [n*m - (d_x-1)*(d_y-1)];
    func = new double [n*m - (d_x-1)*(d_y-1)];
    points = new double [2*(n*m - (d_x-1)*(d_y-1))];
    A = new double [n*m - (d_x-1)*(d_y-1) + 1 +
            ((n-2)*(m-2) - (d_x+1)*(d_y+1))*6
            + ((n-2) + (d_x-1) + (m-2) + (d_y-1))*4*2
            + 2*3 + 2*2 + 2*6 + 2*5];
    I = new int [n*m - (d_x-1)*(d_y-1) + 1 +
            ((n-2)*(m-2) - (d_x+1)*(d_y+1))*6
            + ((n-2) + (d_x-1) + (m-2) + (d_y-1))*4*2
            + 2*3 + 2*2 + 2*6 + 2*5];
    b = new double [n*m - (d_x-1)*(d_y-1)];
}

void reduce_sum (int p) {

    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;

    static int t_in = 0;
    static int t_out = 0;

    if (p <= 1)
        return;

    pthread_mutex_lock (&m);

    t_in++;

    if (t_in < p)
        while (t_in < p)
            pthread_cond_wait (&c_in, &m);

    else {
        t_out = 0;
        pthread_cond_broadcast (&c_in);
    }

    t_out++;

    if (t_out < p)
        while (t_out < p)
            pthread_cond_wait (&c_out, &m);

    else {
        t_in = 0;
        pthread_cond_broadcast (&c_out);
    }

    pthread_mutex_unlock (&m);
}









