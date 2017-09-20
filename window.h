#ifndef WINDOW_H
#define WINDOW_H

#include <math.h>
#include <QWidget>
#include <QPushButton>
#include <QPainter>
#include <QLabel>
#include <QRadioButton>
#include <string.h>
#include <iostream>
#include <QGLWidget>
#include <QtOpenGL>

#include <QMainWindow>
#include <vector>

struct args;

class Scene3D : public QGLWidget {

public:
    Scene3D (QWidget* parent = 0);
    ~Scene3D ();

    int input_values (int argc, char *argv[]);
    void recount_algorithm ();
    void update_arrays ();

private:
    void scale_plus();
    void scale_minus();
    void rotate_up();
    void rotate_down();
    void rotate_left();
    void rotate_right();
    void translate_down();
    void translate_up();
    void defaultScene();
    void drawAxis();

    void fill_vertex_array (std::vector<GLfloat> &VertexArray,
                            double *x, double *func,
                            int index1, int index2, int index3);

    void getVertexArray (std::vector<GLfloat> &VertexArray, double *func,
                         double *x, int status);
    void drawFigure ();
    void draw_points_in_nodes (std::vector<GLfloat> &VertexArray);
    void get_points ();

protected:
    void initializeGL();
    void resizeGL(int nWidth, int nHeight);
    void paintGL();
    void mousePressEvent(QMouseEvent* pe);
    void mouseMoveEvent(QMouseEvent* pe);
    void mouseReleaseEvent(QMouseEvent*);
    void wheelEvent(QWheelEvent* pe);

private:
    int t;      //threads
    int n;      //number of "x" points
    int m;      //number of "y" points

    int p;      //the cut
    int q;
    int d_x;
    int d_y;

    double x1;     //the rectangle
    double y1;
    double x2;
    double y2;

    double *A;     //the MSR matrix (phi,phi)
    int *I;

    double *b;     // (f,phi)
    double *x;

    double *points;
    double *func;

    double (*f) (double, double);
    int what_to_draw;
    double max_function;
    double max_residual;

    std::vector<GLfloat> VertexArray_real;
    std::vector<GLfloat> VertexArray_residual;
    std::vector<GLfloat> VertexArray_appr;

    GLfloat xRot;
    GLfloat yRot;
    GLfloat zRot;
    GLfloat zTra;
    GLfloat nSca;

    QPoint ptrMousePosition;
};

#endif // WINDOW_H
