#include "window.h"

int Scene3D::input_values (int argc, char *argv[]) {

    if (argc != 2) {
        printf ("\nUsage: %s t\n", argv[0]);
        return -1;
    }

    if (!(t = atoi (argv[1]))) {
        printf ("\nWrong t\n");
        return -2;
    }

    if (t <= 0) {
        printf ("\nWrong t\n");
        return -3;
    }

    x1 = y1 = 0;
    x2 = y2 = 7;

    n = 50;
    m = 55;

    p = 10;
    d_x = 11;

    q = 7;
    d_y = 8;

    recount_algorithm ();

    return 1;
}
