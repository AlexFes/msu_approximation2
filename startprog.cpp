#include "window.h"

int Scene3D::input_values(int argc, char *argv[])
{
    if(argc != 5)
    {
        return -1;
    }
    x1 = atoi(argv[1]);
    y1 = atoi(argv[2]);
    x2 = atoi(argv[3]);
    y2 = atoi(argv[4]);
    if((x2 - x1 < 1e-6) || (y2 - y1 < 1e-6))
    {
        return -1;
    }

    n = 5;
    m = 5;
    deltaEnabled = false;
    deltaNums = 0;
    what_to_draw = 0;

    recount_algorithm();

    return 1;
}
