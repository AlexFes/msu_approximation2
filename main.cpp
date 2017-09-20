#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include <QAction>
#include <QMenuBar>
#include <QMessageBox>
#include <fenv.h>
#include "window.h"
#include <stdio.h>

int main (int argc, char *argv[]) {

    QApplication app (argc, argv);

 // feenableexcept (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

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
