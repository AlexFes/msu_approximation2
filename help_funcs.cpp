#include "window.h"
#include "algorithm.h"

/*
//
void Scene3D::fill_vertex_array (std::vector<GLfloat> &VertexArray,
     double *x, double *func, int index1, int index2, int index3)
{
  VertexArray.push_back (points[2 * index1]);
  VertexArray.push_back (points[2 * index1 + 1]);
  if (x)
    VertexArray.push_back (func[index1] - x[index2]);
  else
    VertexArray.push_back (func[index1]);

  VertexArray.push_back (points[2 * index2]);
  VertexArray.push_back (points[2 * index2 + 1]);
  if (x)
    VertexArray.push_back (func[index2] - x[index2]);
  else
    VertexArray.push_back (func[index2]);

  VertexArray.push_back (points[2 * index3]);
  VertexArray.push_back (points[2 * index3 + 1]);
  if (x)
    VertexArray.push_back (func[index3] - x[index3]);
  else
    VertexArray.push_back (func[index3]);
}
//
void Scene3D::fill_vertex_array_half_points (std::vector<GLfloat> &VertexArray, int status, int index1, int index2, int index3)
{
  double x1, x2;

  x1 = (points[2 * index1] + points[2 * index2]) / 2;
  x2 = (points[2 * index1 + 1] + points[2 * index2 + 1]) / 2;
  VertexArray.push_back (x1);
  VertexArray.push_back (x2);
  if (status == 0)
    VertexArray.push_back (f (x1, x2));
  else if (status == 1)
    VertexArray.push_back ((x[index1] + x[index2]) / 2);
  else
    VertexArray.push_back (f (x1, x2) - (x[index1] + x[index2]) / 2);

  x1 = (points[2 * index1] + points[2 * index3]) / 2;
  x2 = (points[2 * index1 + 1] + points[2 * index3 + 1]) / 2;
  VertexArray.push_back (x1);
  VertexArray.push_back (x2);

  if (status == 0)
    VertexArray.push_back (f (x1, x2));
  else if (status == 1)
    VertexArray.push_back ((x[index1] + x[index3]) / 2);
  else
    VertexArray.push_back (f (x1, x2) - (x[index1] + x[index3]) / 2);

  x1 = (points[2 * index2] + points[2 * index3]) / 2;
  x2 = (points[2 * index2 + 1] + points[2 * index3 + 1]) / 2;
  VertexArray.push_back (x1);
  VertexArray.push_back (x2);

  if (status == 0)
    VertexArray.push_back (f (x1, x2));
  else if (status == 1)
    VertexArray.push_back ((x[index2] + x[index3]) / 2);
  else
    VertexArray.push_back (f (x1, x2) - (x[index2] + x[index3]) / 2);
}


//
void Scene3D::getVertexArray (std::vector<GLfloat> &VertexArray,
                        double *func, double *x, int status)
{
  VertexArray.clear ();
  int row, col, index1, index2, index3;
  for (row = 0; row < N_2; row++)
    {
      for (col = 0; col < N_2; col++)
        {
          index1 = get_index (row, col, N_2);
          index2 = get_index (row + 1, col, N_2);
          index3 = get_index (row + 1, col + 1, N_2);
          fill_vertex_array (VertexArray, x, func, index1, index2, index3);
          fill_vertex_array_half_points (VertexArray, status, index1, index2, index3);

          index1 = get_index (row, col, N_2);
          index2 = get_index (row, col + 1, N_2);
          index3 = get_index ( row + 1, col + 1, N_2);
          fill_vertex_array (VertexArray, x, func, index1, index2, index3);
          fill_vertex_array_half_points (VertexArray, status, index1, index2, index3);
        }
    }

  double max = 0;
  for (int i = 0; i < (int)VertexArray.size () / 9; i++)
    {
      for (int j = 0; j < 9; j += 3)
        {
          if (fabs (VertexArray[9 * i + 2 + j]) > max)
            max = fabs (VertexArray[9 * i + 2]);
        }
    }
  if (status == 2)
    {
      max_residual = max;
      printf ("Max residual = %e\n", max_residual);
    }
}
//
int Scene3D::update_arrays ()
{
  P = (N_2 + 1) * (N_2 + 1);

  all_triangles = 2 * P;

  points.clear ();
  func.clear ();
  x.clear ();

  points.resize (2 * P);
  func.resize (P);
  x.resize (P);
  return 0;
}
//
void Scene3D::draw_points_in_nodes (std::vector<GLfloat> &VertexArray)
{
  int i;
  GLfloat x1, y1, z1, x2, y2, z2, x3, y3, z3;
  for (i = 0; i < (int)VertexArray.size () / 9; i++)
    {
      x1 = VertexArray[9 * i + 0] + ellipse_params.x_center;
      x2 = VertexArray[9 * i + 3] + ellipse_params.x_center;
      x3 = VertexArray[9 * i + 6] + ellipse_params.x_center;

      y1 = VertexArray[9 * i + 1] + ellipse_params.y_center;
      y2 = VertexArray[9 * i + 4] + ellipse_params.y_center;
      y3 = VertexArray[9 * i + 7] + ellipse_params.y_center;

      z1 = VertexArray[9 * i + 2];
      z2 = VertexArray[9 * i + 5];
      z3 = VertexArray[9 * i + 8];

      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      glBegin(GL_TRIANGLES);
      glVertex3f (x1, y1, z1);
      glVertex3f (x2, y2, z2);
      glVertex3f (x3, y3, z3);
      glEnd();
    }
}
*/

                ///MA MAN

///Multithread
void Scene3D::recount_algorithm () {

    //pthread_t tid;

    update_arrays();

    get_points ();

//    for (int i = 0; i < n*m - (d_x-1)*(d_y-1); ++i){

//        printf("  x=%f,y=%f",
//               points[2*i], points[2*i+1]);
//    }

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

//    for (int i = 1; i < t; i++)
//        pthread_create (&tid, 0, use_algorithm, arg + i);

    use_algorithm (arg + 0);

//        for (int i = 0; i < n*m - (d_x-1)*(d_y-1); ++i){

//            printf("  f[%d]=%f  ",
//                   i, x[i]);
//        }


    //getVertexArray (VertexArray_real, func.data (), 0, 0);
    getVertexArray (VertexArray_appr,
                    /*func*/x/*.data ()*/, 0, /*0*/1);
    //getVertexArray (VertexArray_residual, func.data (), x.data (), 2);

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

//                if (i==2 && j==4)
//                    printf("  x=%f,y=%f num=%d", x1 + hx*i, y1 + hy*j,
//                           num);

            }

//    for (int i = 0; i < n*m - (d_x-1)*(d_y-1); ++i){

//        printf("  x=%f,y=%f",
//               points[2*i], points[2*i+1]);
//    }

    for (int i = 0; i < n*m - (d_x-1)*(d_y-1); i++) {
        func[i] = f(points[i*2], points[i*2+1]);
        b[i] = func[i];
    }
//    for (int i = 0; i < n*m - (d_x-1)*(d_y-1); ++i){

//        printf("  x=%f,y=%f",
//               points[2*i], points[2*i+1]);
//    }

}

///Push back a point
void Scene3D::fill_vertex_array (std::vector<GLfloat> &VertexArray,
     double *x, double *func, int index1, int index2, int index3) {

    VertexArray.push_back (points[2*index1]);
    VertexArray.push_back (points[2*index1 + 1]);

    if (x)
        VertexArray.push_back (func[index1] - x[index2]);

    else
        VertexArray.push_back (func[index1]);

    VertexArray.push_back (points[2*index2]);
    VertexArray.push_back (points[2*index2 + 1]);

    if (x)
        VertexArray.push_back (func[index2] - x[index2]);

    else
        VertexArray.push_back (func[index2]);

    VertexArray.push_back (points[2*index3]);
    VertexArray.push_back (points[2*index3 + 1]);

    if (x)
        VertexArray.push_back (func[index3] - x[index3]);

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
///????????
    double max = 0;

    for (int i = 0; i < (int)VertexArray.size () / 9; i++)
        for (int j = 0; j < 9; j += 3)
            if (fabs (VertexArray[9 * i + 2 + j]) > max)
                max = fabs (VertexArray[9 * i + 2]);

    if (status == 2) {
        max_residual = max;
        printf ("Max residual = %e\n", max_residual);
    }
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
    A = new double [(n*m - (d_x-1)*(d_y-1))*
                    (n*m - (d_x-1)*(d_y-1))];
    I = new int [(n*m - (d_x-1)*(d_y-1))*
                 (n*m - (d_x-1)*(d_y-1))];
    b = new double [n*m - (d_x-1)*(d_y-1)];
}








