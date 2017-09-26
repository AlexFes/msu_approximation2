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
    x2 = y2 = 4;

    n = 5;
    m = 5;

    p = 1;
    q = 1;

    d_x = 2;
    d_y = 2;

    what_to_draw = 0;

    recount_algorithm ();

    return 1;
}
