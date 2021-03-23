// nonlinear theory of macroscopic flow induced in a drop of ferrofluid
// main file


#include "header.h"


int main() {
    int i;
    double delta_time, finish_time, mu_0, n, start_time, T, t_array[4];
    const char *abs_H_vs_z[4] = {"data/abs_H_vs_z_1.txt",
                                 "data/abs_H_vs_z_2.txt",
                                 "data/abs_H_vs_z_3.txt",
                                 "data/abs_H_vs_z_4.txt"};
    const char *u_z_vs_z[4] = {"data/u_z_vs_z_1.txt",
                               "data/u_z_vs_z_2.txt",
                               "data/u_z_vs_z_3.txt",
                               "data/u_z_vs_z_4.txt"};
    printf("---------------------------\n\n\n");
    start_time = clock();
    T = 2.0 * pi / omega; // the period of the field oscillation
    mu_0 = 4.0 * pi * 1e-7; // vacuum permeability
    n = (log(2.0 * r_0) - 1.0) / int_pow(r_0, 2); // aggregate demagnetizing factor
    printf("period = %.2lf sec\n", T);
    write_beta_ratio();
    write_phi_vs_z();
    worker_t_array(t_array, T);
    for(i = 0; i < 4; ++i) {
        write_abs_H_vs_z(t_array[i], abs_H_vs_z[i]);
        write_u_z_vs_z(mu_0, n, t_array[i], u_z_vs_z[i]);
    }
    write_vector_field(mu_0, n);
    part_diff_eq(mu_0, n, T, t_array);
    finish_time = clock();
    delta_time = finish_time - start_time;
    printf("C runtime = %.2lf sec = %.2lf min = %.2lf hours\n\n", delta_time * 1.0 / CLOCKS_PER_SEC,
           delta_time / (60.0 * CLOCKS_PER_SEC), delta_time / (3600.0 * CLOCKS_PER_SEC));
    return 0;
}
