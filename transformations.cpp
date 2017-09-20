#include "window.h"

void Scene3D::initializeGL() {

    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    qglClearColor (Qt::white);
    glEnable (GL_DEPTH_TEST);
    glShadeModel (GL_FLAT);
}

///DANGER
void Scene3D::resizeGL (int nWidth, int nHeight) {

    glMatrixMode (GL_PROJECTION);
    glLoadIdentity( );

    if (fabs (max_function) < 1e-12)
        max_function = 1.;

    GLfloat ratio = (GLfloat)nHeight / (GLfloat)nWidth;

    if (what_to_draw == 0)
        glOrtho (-2./*left*/, 2./*right*/,
             -2*ratio*max_function/*bottom*/,
              2*ratio*max_function/*top*/,
             -5./*near*/, 5./*far*/);

    else
        glOrtho (-2./*left*/, 2./*right*/,
             -2 * ratio * max_residual/*bottom*/,
              2 * ratio * max_residual/*top*/,
             -5./*near*/, 5./*far*/);
  /*
    glOrtho (-2.0 * max_function,         2.0 * max_function,
             -2.0 * max_function * ratio, 2.0 * max_function * ratio,
             -2.0 * max_function,         2.0 * max_function);
*/

    glViewport ((x1+x2)/1, (y1+y2)/1,
               (GLint)nWidth, (GLint)nHeight);
}

void Scene3D::mousePressEvent(QMouseEvent* pe)
{
  ptrMousePosition = pe->pos();
}

void Scene3D::mouseReleaseEvent(QMouseEvent* /*pe*/)
{

}

void Scene3D::mouseMoveEvent(QMouseEvent* pe) {

    xRot += 180/nSca*(GLfloat)(pe->y() - ptrMousePosition.y())/height();
    zRot += 180/nSca*(GLfloat)(pe->x() - ptrMousePosition.x())/width();

    ptrMousePosition = pe->pos();

    updateGL();
}

void Scene3D::wheelEvent(QWheelEvent* pe) {

    if ((pe->delta()) > 0)
        scale_plus ();

    else if ((pe->delta()) < 0)
        scale_minus ();

    updateGL();
}
/*
void Scene3D::keyPressEvent(QKeyEvent* pe)
{
  switch (pe->key ())
    {
      case Qt::Key_Equal:
        scale_plus();
        break;

      case Qt::Key_Minus:
        scale_minus ();
        break;

      case Qt::Key_W:
        rotate_up ();
        break;

      case Qt::Key_S:
        rotate_down ();
        break;

      case Qt::Key_A:
        rotate_left ();
        break;

      case Qt::Key_D:
          rotate_right ();
        break;

      case Qt::Key_Q:
        translate_down ();
        break;

      case Qt::Key_E:
        translate_up ();
        break;

      case Qt::Key_Space:
        defaultScene ();
        break;

      case Qt::Key_Escape:
        this->close ();
        return;

      case Qt::Key_1:
        what_to_draw++;
        if (what_to_draw == 2)
          what_to_draw = 0;
        break;

      case Qt::Key_2:
        N++;
        N_2 *= 2;
        N_2 *= 2;
        recount_algorithm ();
        break;
      case Qt::Key_3:
        if (N > 1)
          {
            N--;
            N_2 /= 2;
            N_2 /= 2;
            recount_algorithm ();
          }
        break;
      default:
        return;
    }

  updateGL ();
}
*/
void Scene3D::scale_plus()
{
  nSca = nSca * 1.1;
}

void Scene3D::scale_minus()
{
  nSca = nSca / 1.1;
}

void Scene3D::rotate_up()
{
  xRot += 1.0;
}

void Scene3D::rotate_down()
{
  xRot -= 1.0;
}

void Scene3D::rotate_left()
{
  zRot += 1.0;
}

void Scene3D::rotate_right()
{
  zRot -= 1.0;
}

void Scene3D::translate_down()
{
  zTra -= 0.05;
}

void Scene3D::translate_up()
{
  zTra += 0.05;
}

void Scene3D::defaultScene()
{
  xRot = -90;
  yRot = 0;
  zRot = 0;
  zTra = 0;
  nSca = 1;
}


