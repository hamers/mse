#include "types.h"
extern "C"
{

int test_tools();
int test_a_h_conversion();
int test_orbital_element_conversion();
int test_kepler_equation_solver();
int test_orbital_vectors_cartesian_conversion();

int test_collisions();
int test_collision_MS_MS();
int test_collision_giant_MS();

int test_collision_stars(double m1, int kw1, double m2, int kw2);

}
