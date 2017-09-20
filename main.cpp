#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include <QAction>
#include <QMenuBar>
#include <QMessageBox>
#include "window.h"
#include <stdio.h>

int main (int argc, char *argv[]) {

    QApplication app (argc, argv);

    QMainWindow *window = new QMainWindow;

    window->setWindowTitle ("Graph");

    Scene3D draw_area;

    if (draw_area.input_values (argc, argv) != 1)
        return -1;

    draw_area.resize (500, 500);
    draw_area.show ();
    app.exec ();

    return 0;
}
