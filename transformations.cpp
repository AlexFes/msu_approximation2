#include "window.h"

void Scene3D::initializeGL()
{
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    qglClearColor(Qt::white);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_FLAT);
}

void Scene3D::resizeGL(int nWidth, int nHeight)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    GLfloat ratio=(GLfloat)nHeight/(GLfloat)nWidth;
    if (nWidth>=nHeight)
    {
       glOrtho(-1.0 / ratio, 1.0 / ratio, -1.0, 1.0, -10.0, 10.0);
    }
    else
    {
       glOrtho(-1.0, 1.0, -1.0 * ratio, 1.0 * ratio, -10.0, 10.0);
    }

    glViewport(0, 0, (GLint)nWidth, (GLint)nHeight);
}

void Scene3D::mousePressEvent(QMouseEvent* pe)
{
    ptrMousePosition = pe->pos();
}

void Scene3D::mouseReleaseEvent(QMouseEvent* /*pe*/)
{

}

void Scene3D::mouseMoveEvent(QMouseEvent* pe)
{
    xRot += 180 / nSca * (GLfloat)(pe->y() - ptrMousePosition.y()) / height();
    zRot += 180 / nSca * (GLfloat)(pe->x() - ptrMousePosition.x()) / width();

    ptrMousePosition = pe->pos();
    updateGL();
}

void Scene3D::wheelEvent(QWheelEvent* pe)
{

    if ((pe->delta()) > 0)
    {
        scale_plus();
    }
    else if ((pe->delta()) < 0)
    {
        scale_minus();
    }
    updateGL();
}

void Scene3D::keyPressEvent(QKeyEvent* pe)
{
    switch(pe->key())
    {
    case Qt::Key_Space:
        deltaEnabled = false;
        deltaNums = 0;
        defaultScene();
        break;

    case Qt::Key_Escape:
        this->close();
        return;

    case Qt::Key_1:
        what_to_draw ++;
        if(what_to_draw == 2)
        {
            what_to_draw = 0;
        }
        recount_algorithm();
        break;

    case Qt::Key_2:
        deltaEnabled = false;
        deltaNums = 0;
        n += n - 1;
        recount_algorithm();
        break;

    case Qt::Key_3:
        deltaEnabled = false;
        deltaNums = 0;
        n += 1;
        n /= 2;
        if(n < 5)
        {
            n = 5;
        }
        recount_algorithm();
        break;

    case Qt::Key_4:
        deltaEnabled = false;
        deltaNums = 0;
        m += m - 1;
        recount_algorithm();
        break;

    case Qt::Key_5:
        deltaEnabled = false;
        deltaNums = 0;
        m += 1;
        m /= 2;
        if(m < 5)
        {
            m = 5;
        }
        recount_algorithm();
        break;
    case Qt::Key_6:
        deltaEnabled = true;
        deltaNums ++;
        recount_algorithm();
        break;
    case Qt::Key_W:
        deltaEnabled = false;
        deltaNums = 0;
        y2 += 0.1;
        y1 += 0.1;
        recount_algorithm();
        break;
    case Qt::Key_S:
        deltaEnabled = false;
        deltaNums = 0;
        y2 -= 0.1;
        y1 -= 0.1;
        recount_algorithm();
        break;
    case Qt::Key_A:
        deltaEnabled = false;
        deltaNums = 0;
        x1 += 0.1;
        x2 += 0.1;
        recount_algorithm();
        break;
    case Qt::Key_D:
        deltaEnabled = false;
        deltaNums = 0;
        x1 -= 0.1;
        x2 -= 0.1;
        recount_algorithm();
        break;
    }

    updateGL();
}

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


