#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include <QAction>
#include <QMenuBar>
#include <QMessageBox>
#include "window.h"
#include <stdio.h>

int main (int argc, char *argv[])
{

    QApplication app(argc, argv);
    QMainWindow *window = new QMainWindow;
    window->setWindowTitle("Graph");
    Scene3D draw_area;

    if (draw_area.input_values(argc, argv) != 1)
    {
        printf("Execute the program with arguments: X1, Y1, X2, Y2.\nThe function will be drawn on the rectange from (X1, Y1) to (X2, Y2).\n");
        return -1;
    }

    draw_area.resize(1920, 1080);
    draw_area.show();
    app.exec();

    delete window;
    return 0;
}
