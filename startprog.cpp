#include "window.h"
/*
int Scene3D::input_values (int argc,char *argv[])
{
  if (argc != 2)
    {
      printf ("Usage: %s p\n", argv[0]);
      return -1;
    }

  if (!(p = atoi (argv[1])))
    {
      printf ("Wrong p\n");
      return -2;
    }

  if (p <= 0)
    {
      printf ("p needs to be >= 1\n");
      return -3;
    }

  printf ("Input parameters of ellipse:\n");
  printf ("Input center:\nx=");
  if (!scanf ("%lf", &ellipse_params.x_center))
    {
      printf ("wrong input\n");
      return -1;
    }

  printf ("y=");
  if (!scanf ("%lf", &ellipse_params.y_center))
    {
      printf ("wrong input\n");
      return -1;
    }

  printf ("Input a:");
  if (!scanf ("%lf", &ellipse_params.a))
    {
      printf ("wrong input\n");
      return -1;
    }

  printf ("Input b:");
  if (!scanf ("%lf", &ellipse_params.b))
    {
      printf ("wrong input\n");
      return -1;
    }

  //init vertices

  vertices[bottom_left_x] = -ellipse_params.a / sqrt (2);
  vertices[bottom_left_y] = -ellipse_params.b / sqrt (2);

  vertices[bottom_right_x] = ellipse_params.a / sqrt (2);
  vertices[bottom_right_y] = -ellipse_params.b / sqrt (2);

  vertices[top_left_x] = -ellipse_params.a / sqrt (2);
  vertices[top_left_y] = ellipse_params.b / sqrt (2);

  vertices[top_right_x] = ellipse_params.a / sqrt (2);
  vertices[top_right_y] = ellipse_params.b / sqrt (2);


//  double len_cols = vertices[top_right_x] - vertices[top_left_x],
//         len_rows = vertices[top_right_y] - vertices[bottom_right_y];
  N_2 = 2;

//  printf ("N_rows = %d, N_cols = %d\n", N_2, N_2);
  if (update_arrays () < 0)
    return -5;
  recount_algorithm ();
  return 1;
}
*/
                    ///MA MAN

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
/*
    printf ("\n\nInput parameters of rectangle:\n");

    printf ("Input the bottom left point:\nx1 = ");
    if (!scanf ("%f", &x1)) {
        printf ("\nWrong input\n");
        return -1;
    }

    printf ("\ny1 = ");
    if (!scanf ("%f", &y1)) {
        printf ("\nWrong input\n");
        return -1;
    }

    printf ("\n\nInput the top right point:\nx2 = ");
    if (!scanf ("%f", &x2)) {
        printf ("\nWrong input\n");
        return -1;
    }

    printf ("\ny2 = ");
    if (!scanf ("%f", &y2)) {
        printf ("\nWrong input\n");
        return -1;
    }

    printf ("\n\nInput number of points:\nn = ");
    if (!scanf ("%d", &n)) {
        printf ("\nWrong input\n");
        return -1;
    }

    printf ("\nm = ");
    if (!scanf ("%d", &m)) {
        printf ("\nWrong input\n");
        return -1;
    }
*/
    ///init the cut

    x1 = y1 = 0;
    x2 = y2 = 7;

    n = 50;
    m = 55;
//    p = n/4;
//    d_x = n/2;

//    q = m/4;
//    d_y = m/2;

    p = 10;
    d_x = 11;

    q = 7;
    d_y = 8;

    recount_algorithm ();

    return 1;
}
