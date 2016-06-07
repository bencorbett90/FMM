#include "boost_1_61_0/boost/math/special_functions/spherical_harmonic.hpp"
#include <complex>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <time.h>

typedef std::complex<double> dcomplex;

typedef unsigned uint;

typedef struct {
    double x;
    double y;
    double z;
} point_t;

typedef struct {
    double r;
    double theta;
    double phi;
} point_spher_t;

typedef struct {
    double x;
    double y;
    double z;
    double q;
} pCharge_t;


/* Return the magnitude squared of point P. */
double magnitude2(point_t *p) {
    return (p->x * p->x) + (p->y * p->y) + (p->z * p->z);
}

/* Return the magnitude of point P. */
double magnitude(point_t *p) {
    return pow(magnitude2(p), 0.5);
}

/* Convert point P in cartesian coordinates to spherical coordinates
   and store the result in P_SPHER. */
void cartesian_to_spherical(point_t *p, point_spher_t *p_spher) {
    p_spher->r = magnitude(p);
    p_spher->theta = atan2(p->y, p->x);

    double temp = pow((p->x * p->x) + (p->y * p->y), 0.5);
    p_spher->phi = atan2(temp, p->z);
}

/* Wrapper for boost::spherical_harmonic that renormalizes. */
dcomplex spherical_harmonic(uint n, int m, double theta, double phi) {
    dcomplex temp = boost::math::spherical_harmonic(n, -m, theta, phi);
    return temp * pow(4 * M_PI / (2 * n + 1), 0.5);
}



/* Construct a multipole expansion consisting of P terms, centered about the point ORIGIN. 
   The returned array has N as it's first dimension and M as the second. */
dcomplex* calc_multipoles(uint p, point_t *origin, uint num_particles, pCharge_t *particles) {
    uint num_terms = (p + 1) * (p + 1);
    dcomplex *multipole_terms = new dcomplex[num_terms];
    
    point_t shifted_pos;
    point_spher_t spher_pos;
    double weighted_charge;
    uint index;

    /* Iterate over all particles. */
    for (uint i = 0; i < num_particles; i++) {
        pCharge_t *pCharge = particles + i;
        weighted_charge = pCharge->q;

        /* Get the shifted cartesian coordinates. */
        shifted_pos = *((point_t *) &(pCharge->x));
        shifted_pos.x -= origin->x;
        shifted_pos.y -= origin->y;
        shifted_pos.z -= origin->z;


        /* Convert the shifted coordinates into spherical. */
        cartesian_to_spherical(&shifted_pos, &spher_pos);
        
        index = 0;
        /* Iterate over n, which goes from [0, p]. */
        for (int n = 0; n <= p; n++) {

            /* Iterate over m, which goes from [-n, n]. */
            for (int m = -n; m <= n; m++) {
                dcomplex temp = spherical_harmonic(n, -m, spher_pos.theta, spher_pos.phi);
                multipole_terms[index] += weighted_charge * temp;

                index++;
            }

            weighted_charge *= spher_pos.r;
        }

    }
    return multipole_terms;
}


dcomplex calc_potential_from_multipole(point_t *point, dcomplex *expansion, uint p, point_t *origin) {
    point_t shifted_pos;
    point_spher_t spher_pos;

    shifted_pos = *point;
    shifted_pos.x -= origin->x;
    shifted_pos.y -= origin->y;
    shifted_pos.z -= origin->z;

    /* After this call SPHER_POS holds the spherical coordinates of POINT from ORIGIN. */
    cartesian_to_spherical(&shifted_pos, &spher_pos);

    double r_weight = 1.0 / spher_pos.r;
    dcomplex V = dcomplex(0, 0);
    uint index = 0;

    for (int n = 0; n <= p; n++) {

        /* Iterate over m, which goes from [-n, n]. */
        for (int m = -n; m <= n; m++) {
            dcomplex temp = spherical_harmonic(n, m, spher_pos.theta, spher_pos.phi);
            V += expansion[index] * r_weight * temp;
            index++;
        }

        r_weight /= spher_pos.r;
    }
    return V;
}


void print_expansion(dcomplex *expansion, uint p) {
    int index = 0;

    printf("\nExpansion terms:\n");
    for (int n = 0; n <= p; n++) {
        printf("N = %d\n", n);
        printf("\t");

        /* Iterate over m, which goes from [-n, n]. */
        for (int m = -n; m <= n; m++) {
            printf("%g + %gI, ", expansion[index].real(), expansion[index].imag());
            index++;
        }
        printf("\n");
    }
}

/* Generate NUM_SOURCES point charges within at most a distance A of ORIGIN. */
pCharge_t* gen_source_points(uint num_sources, double a, point_t *origin) {
    
    pCharge_t *source_array = (pCharge_t *) malloc(num_sources * sizeof(pCharge_t));

    point_t temp;
    for (int i = 0; i < num_sources; ) {
        temp.x = 2 * a * (drand48() - 0.5);
        temp.y = 2 * a * (drand48() - 0.5);
        temp.z = 2 * a * (drand48() - 0.5);
        if (magnitude(&temp) >= a) continue;

        source_array[i].x = temp.x + origin->x;
        source_array[i].y = temp.y + origin->y;
        source_array[i].z = temp.z + origin->z;
        source_array[i].q = 2 * (drand48() - 0.5);
        i++;
    }
    return source_array;
}

/* Generate NUM_TEST points at least a distance A from ORIGIN. */
point_t* gen_test_points(uint num_tests, double a, point_t *origin) {
    
    point_t *test_array = (point_t *) malloc(num_tests * sizeof(point_t));

    point_t temp;
    for (int i = 0; i < num_tests; ) {
        temp.x = 2 * (drand48() - 0.5);
        temp.y = 2 * (drand48() - 0.5);
        temp.z = 2 * (drand48() - 0.5);
        if (temp.x == 0 || temp.y == 0 || temp.z == 0) continue;
        
        temp.x = a / temp.x;
        temp.y = a / temp.y;
        temp.z = a / temp.z;
        if (magnitude(&temp) <= a) continue;

        test_array[i].x = temp.x + origin->x;
        test_array[i].y = temp.y + origin->y;
        test_array[i].z = temp.z + origin->z;
        i++;
    }
    return test_array;
}


void multipole_test(uint p, point_t *origin, double radius, uint num_sources, uint num_tests) {
    pCharge_t *source_points = gen_source_points(num_sources, radius, origin);
    point_t *test_points = gen_test_points(num_tests, radius, origin);

    for (int i = 0; i < num_sources; i++) {
        if (magnitude((point_t *) (source_points + i)) >= radius) printf("Source error\n");
    }
    for (int i = 0; i < num_tests; i++) {
        if (magnitude(test_points + i) <= radius) printf("Test error\n");
    }

    dcomplex *expansion = calc_multipoles(p, origin, num_sources, source_points);

    double error = 0.0;
    double error_norm = 0.0;
    double max_error = 0.0;

    point_t shifted_pos;

    for (int i = 0; i < num_tests; i++) {
        point_t *test_point = test_points + i;
        double V_exact = 0.0;

        for (int j = 0; j < num_sources; j++) {
            pCharge_t *source_point = source_points + j;

            shifted_pos.x = test_point->x - source_point->x;
            shifted_pos.y = test_point->y - source_point->y;
            shifted_pos.z = test_point->z - source_point->z;

            V_exact += source_point->q / magnitude(&shifted_pos);
        }

        dcomplex V_aprox = calc_potential_from_multipole(test_point, expansion, p, origin);
        double cur_error = pow(V_aprox.real() - V_exact, 2);

        error += cur_error;
        error_norm += V_exact * V_exact;
        if (cur_error > max_error) max_error = cur_error;
    }

    error = pow(error / error_norm, 0.5);
    max_error = pow(max_error, 0.5);

    printf("Num terms: %u\n", p);
    printf("L2 error norm: %g\n", error);
    printf("Max error: %g\n", max_error);

    delete expansion;
    free(source_points);
    free(test_points);
}



int main(int argc, char** argv) {
    srand48(10);

    uint p = 5;
    if (argc == 2) {
        p = atoi(argv[1]);
    }

    point_t origin = {0, 0, 0};
    double radius = 1.0;
    uint num_sources = 10;
    uint num_tests = 10;

    multipole_test(p, &origin, radius, num_sources, num_tests);
    
    return 0;
}