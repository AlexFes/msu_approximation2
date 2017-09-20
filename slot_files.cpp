#include <QObject>
#include <QPainter>
#include <stdio.h>

#include "window.h"
#include <QLineEdit>
#include <QLabel>
#include <QHBoxLayout>
#include <QSpacerItem>
#include <QPushButton>
#include <QDebug>
#include <QMessageBox>
#include <QRadioButton>


void Window::double_N()
{
    N++;
    N2 *= 2;
    update();
}

void Window::reduce_N()
{
    N--;
    N2 /= 2;
    update();
}
