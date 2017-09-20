#include <QtGui>
#include "window.h"
#include "algorithm.h"
#include <math.h>
#include <string.h>
#include <QDebug>

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







