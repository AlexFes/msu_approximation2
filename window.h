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

class Scene3D : public QGLWidget
{
public:
    Scene3D(QWidget* parent = 0);                   ///Constructor
    ~Scene3D();                                     ///Destructor

    int     input_values(int argc, char *argv[]);   ///Initialization of the parameters based on the input values
    void    recount_algorithm();                    ///The trigger to recalculate the MSR matrix
    void    update_arrays();                        ///Memory reallocation

private:
    bool    deltaEnabled;                           //whether or not the function values are shifted in several points
    int     deltaNums;                              //defines how much are values shifted
    int     n;                                      //the length of the rectangle
    int     m;                                      //the width of the rectangle
    int     numPoints;                              //amount of triangulation points on the plane.
    int     size;                                   //size of MSR matrix
    int     what_to_draw;                           //the switch between the approximation and residual
    double  x1;                                     //the 1st coordinate of the 1st point
    double  y1;                                     //the 2nd coordinate of the 1st point
    double  x2;                                     //the 1st coordinate of the 2nd point
    double  y2;                                     //the 2nd coordinate of the 2nd point
    double  s;                                      //size of rectangle between neighbouring points
    double  max_function;                           //maximal function value, is required for transformation
    double  max_residual;                           //maximal residual value, is required for transformation
    double  *A;                                     //the MSR matrix
    double  *b;                                     //Additional arrays required to calculate the matrix according to the whitepaper.
    double  *x;
    double  *points;
    double  *func;
    int     *I;                                     //companion matrix for the MSR matrix which stores relations
    std::vector<GLfloat> VertexArray_residual;      //vertexes to draw residual
    std::vector<GLfloat> VertexArray_appr;          //vertexes to draw approximation
    GLfloat xRot;                                   //the angle between Ox axis and current screen position
    GLfloat yRot;                                   //the angle between Oy axis and current screen position
    GLfloat zRot;                                   //the angle between Oz axis and current screen position
    GLfloat zTra;                                   //the variable for translation of coordinates
    GLfloat nSca;                                   //the scale of the drawn approximation/residual, changed by scrolling
    QPoint ptrMousePosition;                        //current position of the mouse, allows to rotate by tracking its change

    ///Methods for button binding:
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

    ///Methods for filling vertex arrays used to draw approximation/residual:
    void fill_vertex_array(std::vector<GLfloat> &VertexArray, double *x, double *func, int index1, int index2, int index3);
    void fill_vertex_array_half_points(std::vector<GLfloat> &VertexArray, int status, int index1, int index2, int index3);
    void getVertexArray(std::vector<GLfloat> &VertexArray, double *func, double *x, int status);
    void drawFigure();
    void draw_points_in_nodes(std::vector<GLfloat> &VertexArray);
    void get_points();

    ///Methods for solving the MSR matrix required for the approximation:
    void use_algorithm(void);                                                   ///Input function for the calculation
    bool assemble_matrix(void);                                                 ///Fills MSR matrix
    void fill_vector_b(void);                                                   ///Fills vector b
    void integral(int i1, int j1, int i2, int j2, int i3, int j3, int flag);    ///Calculate 4 integrals
    int get_links(int i, int j, int *x, int *y);                                ///Maps the vertices according to the coordinate positioning on the approximation rectangle
    int get_num(int i, int j);                                                  ///Obtains the position in the MSR matrix of the (i, j)th element of the plane
    void get_ij(int &i, int &j, int l);                                         ///Obtains the coordinates (i, j) on the plane for the l'th element in the MSR matrix
    int get_values(int i, int j, double *a);                                    ///Calculates koefficients of each triangluation point of the (i, j)th element of the plane
    int solve(double *r, double *u, double *v);                                 ///Solve the system with iterations.
    double f(double x, double y);                                               ///Stores the approximated function and returns its values with possible shifts made by user

protected:
    ///Overwritten methods of parent class:
    void initializeGL();
    void resizeGL(int nWidth, int nHeight);
    void paintGL();
    void mousePressEvent(QMouseEvent* pe);
    void mouseMoveEvent(QMouseEvent* pe);
    void mouseReleaseEvent(QMouseEvent*);
    void wheelEvent(QWheelEvent* pe);
    void keyPressEvent(QKeyEvent* pe);

private:

};

#endif // WINDOW_H
