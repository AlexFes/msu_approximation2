#include <QtGui>
#include "window.h"
#include "algorithm.h"
#include <math.h>
#include <string.h>
#include <QDebug>
/*
static double f_0 (double x, double y)
{
  (void)x;
  (void)y;
  return sqrt (x * x + y * y);
}

Scene3D::Scene3D (QWidget* parent) : QGLWidget(parent)
{
  xRot = -90; yRot = 0; zRot = 0; zTra = 0; nSca = 1;
  N = 1;
  N_2 = 2;
  N_2 = 2;
  f = f_0;
  what_to_draw = 0;
  all_triangles = 1;
  max_function = 1.;
  max_residual = 1.;
}

Scene3D::~Scene3D ()
{
}

void Scene3D::recount_algorithm ()
{
  int i;
  pthread_t tid;

  //fill arrays of triangulation

  P = 2 * (N_2 + 1) * (N_2 + 1);
  if (update_arrays () < 0)
    return;

  get_points ();

  double *a, *b, *ans;
  int *jnz;

  int len = P + 1 + get_nz_matrix (N_2, N_2);

  if (   !(a   = new double [len])
      || !(jnz = new int [len])
      || !(b   = new double [P])
      || !(ans = new double [P]))
    return; /// \todo delete all initialized memory

  memset (a, 0, len * sizeof (double));
  memset (jnz, 0, len * sizeof (int));
  memset (b, 0, P * sizeof (double));
  memset (ans, 0, P * sizeof (double));

  init_matrix_b (points, b, P, f);

  args *arg = new args [p];
  for (i = 0; i < p; i++)
    {
      arg[i].a = a;
      arg[i].p = p;
      arg[i].b = b;
      arg[i].k = i;
      arg[i].x = ans;
      arg[i].P = P;
      arg[i].N2 = N_2;
      arg[i].jnz = jnz;
      arg[i].points = points;
      arg[i].f = f;
      memcpy (arg[i].vertices, vertices, vertex_pos_COUNT * sizeof (double));
    }

  for (i = 1; i < p; i++)
    pthread_create (&tid, 0, use_algorithm, arg + i);
  use_algorithm (arg + 0);
  //x = x_new; // copy vector

  for (int idx = 0; idx < P; idx++)
    x[idx] = ans[idx];

  updateGL();
  printf ("Ended\n");
  delete [] arg;
  delete [] a;
  delete [] jnz;
  delete [] b;

  getVertexArray (VertexArray_real, func.data (), 0, 0);
  getVertexArray (VertexArray_appr, x.data (), 0, 1);
  getVertexArray (VertexArray_residual, func.data (), x.data (), 2);
}

void Scene3D::drawFigure ()
{
  glLineWidth (1.0f);

  if (what_to_draw % 2 == 0)
    {
      glColor3f (0, 0, 255);
      draw_points_in_nodes (VertexArray_real);
      glColor3f (255, 0, 0);
      draw_points_in_nodes (VertexArray_appr);
    }
  if (what_to_draw % 2 == 1)
    {
      glColor3f (0, 1, 0);
      draw_points_in_nodes (VertexArray_residual);
    }
}

void Scene3D::paintGL ()
{
  max_function = func[0];
  for (int i = 1; i < P; i++)
    {
      if (func[i] > max_function)
        max_function = func[i];
    }
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  drawAxis ();
  drawFigure ();

  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity();

  glScalef (nSca, nSca, nSca);
  glTranslatef (0.0f, zTra, 0.0f);
  glRotatef (xRot, 1.0f, 0.0f, 0.0f);
  glRotatef (yRot, 0.0f, 1.0f, 0.0f);
  glRotatef (zRot, 0.0f, 0.0f, 1.0f);
}

void Scene3D::get_points ()
{
  double x1 = vertices[bottom_left_x],
         y1 = vertices[bottom_left_y],
         x2 = vertices[top_right_x],
         y2 = vertices[top_right_y];

  double len_x = x2 - x1,
         len_y = y2 - y1,
         step_x = len_x / N_2,
         step_y = len_y / N_2;

  double curr_x = 0.,
         curr_y = 0.;

  /// Centers are always 0, and then in drawing parallel offset is added
  double x_center = 0.,
         y_center = 0.;

  double a = ellipse_params.a,
         b = ellipse_params.b;

  double r = sqrt (  (x1 - x_center) * (x1 - x_center)
                   + (y1 - y_center) * (y1 - y_center));

  for (int row = 0; row <= N_2; row++)
    {
      curr_y = y1 + row * step_y;
      for (int col = 0; col <= N_2; col++)
        {
          curr_x = x1 + col * step_x;
          double current_r = sqrt (  (curr_x - x_center) * (curr_x - x_center)
                                   + (curr_y - y_center) * (curr_y - y_center));

          if (fabs (current_r) < 1e-12)
            {
              int ind = get_index (row, col, N_2);
              points[2 * ind] = curr_x;
              points[2 * ind + 1] = curr_y;
            }
          else
            {
              double r_ = r;
              double cos_phi = (curr_x - x_center) / current_r,
                     sin_phi = (curr_y - y_center) / current_r;

              /// down triangle
              ///
              if (row <= N_2 / 2 && (col >= row && col <= N_2 - row))
                {
                  r_ = r - r * 2 * row / N_2;
                }
              /// up triangle
              ///
              else if (row >= N_2 / 2 && (col >= N_2 - row && col <= row))
                {
                  r_ -= r * 2 * (N_2 - row) / N_2;
                }
              /// right triangle
              ///
              else if (col >= N_2 / 2 && (row >= N_2 - col && row <= col))
                {
                  r_ = r - r * 2 * (N_2 - col) / N_2;
                }
              /// left triangle
              ///
              else if (col <= N_2 / 2 && (row >= col && row <= N_2 - col))
                {
                  r_ -= r * 2 * col / N_2;
                }

              int ind = get_index (row, col, N_2);
              points[2 * ind] = a * r_ * cos_phi;
              points[2 * ind + 1] = b * r_ * sin_phi;
            }
        }
    }

  for (int i = 0; i < P; i++)
    func[i] = f (points[2 * i], points[2 * i + 1]);
}
*/

                        ///MA MAN
static double f_0 (double x, double y) {

    (void)x;
    (void)y;

    return sqrt (x*x + y*y);
    //return x + y;
}

void Scene3D::paintGL () {

    max_function = func[0];
    for (int i = 1; i < n*m - (d_x-1)*(d_y-1); i++)
        if (func[i] > max_function)
            max_function = func[i];

    //printf("  %f   ", max_function);

    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    drawAxis ();
    drawFigure ();

    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();

    glScalef (nSca/max_function, nSca/max_function, nSca/max_function);
    glTranslatef (0.0f, zTra, 0.0f);
    glRotatef (xRot, 1.0f, 0.0f, 0.0f);
    glRotatef (yRot, 0.0f, 1.0f, 0.0f);
    glRotatef (zRot, 0.0f, 0.0f, 1.0f);
}

void Scene3D::drawFigure () {

    glLineWidth (1.0f);

    if (what_to_draw % 2 == 0) {
//        glColor3f (0, 0, 255);
//        draw_points_in_nodes (VertexArray_real);
        glColor3f (255, 0, 0);
        draw_points_in_nodes (VertexArray_appr);
    }

    if (what_to_draw % 2 == 1) {
        glColor3f (0, 1, 0);
        draw_points_in_nodes (VertexArray_residual);
    }
}

void Scene3D::drawAxis() {

    glLineWidth (3.0f);

    glColor4f (1.00f, 0.00f, 0.00f, 1.0f);

    glBegin (GL_LINES);
    glVertex3f ( 20.0f,  0.0f,  0.0f);
    glVertex3f (-20.0f,  0.0f,  0.0f);
    glEnd ();

    QColor halfGreen (0, 128, 0, 255);
    qglColor (halfGreen);

    glBegin (GL_LINES);
    glVertex3f (0.0f,   20.0f,  0.0f);
    glVertex3f (0.0f, -20.0f,  0.0f);

    glColor4f (0.00f, 0.00f, 1.00f, 1.0f);
    glVertex3f ( 0.0f,  0.0f,  20.0f);
    glVertex3f ( 0.0f,  0.0f, -20.0f);
    glEnd ();
}

void Scene3D::draw_points_in_nodes (std::vector<GLfloat>
                                    &VertexArray) {

    GLfloat x1, y1, z1, x2, y2, z2, x3, y3, z3;

    for (int i = 0; i < (int) VertexArray.size()/9; i++) {
        x1 = VertexArray[9 * i + 0];
        x2 = VertexArray[9 * i + 3];
        x3 = VertexArray[9 * i + 6];

        y1 = VertexArray[9 * i + 1];
        y2 = VertexArray[9 * i + 4];
        y3 = VertexArray[9 * i + 7];

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

Scene3D::Scene3D (QWidget* parent) : QGLWidget(parent) {

    xRot = -90; yRot = 0; zRot = 0; zTra = 0; nSca = 1;
    f = f_0;
    what_to_draw = 0;
    max_function = 1.;
    max_residual = 1.;
}

Scene3D::~Scene3D () {

    delete[] x;
    delete[] func;
    delete[] points;
    delete[] A;
    delete[] I;
    delete[] b;
}







