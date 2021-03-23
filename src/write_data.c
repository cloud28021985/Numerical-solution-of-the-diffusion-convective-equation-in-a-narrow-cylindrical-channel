// nonlinear theory of macroscopic flow induced in a drop of ferrofluid
// recording calculation results into the text files


#include "header.h"


void write_abs_H_vs_z(double t, const char *abs_H_vs_z) {
    int i;
    double abs_H, h_plot, z;
    h_plot = 2.0 * L / (N_plot - 1.0);
    z = - L;
    FILE *myfile = NULL;
    myfile = fopen(abs_H_vs_z, "w");
    fprintf(myfile, "%10s %20s\n", "# z / L", "abs(H), kA/m");
    for(i = 0; i < N_plot; ++i) {
        worker_abs_H(&abs_H, t, z);
        fprintf(myfile, "%10.6lf %20.6lf\n", z / L, abs_H / 1000.0);
        z = z + h_plot;
    }
    fclose(myfile);
}


void write_beta_ratio() {
    int i;
    double beta_ratio, h_plot, r, r_min, r_max;
    r_min = 10.0;
    r_max = 100.0;
    h_plot = (r_max - r_min) / (N_plot - 1.0);
    r = r_min;
    FILE *myfile = NULL;
    myfile = fopen("data/beta_ratio_vs_r.txt", "w");
    fprintf(myfile, "%10s %20s\n", "# r", "beta_ratio");
    for(i = 0; i < N_plot; ++i) {
        worker_beta_ratio(&beta_ratio, r);
        fprintf(myfile, "%10.6lf %20.6lf\n", r, beta_ratio);
        r = r + h_plot;
    }
    fclose(myfile);
}


void write_phi_vs_z() {
    int i;
    double h_plot, phi, z;
    h_plot = 2.0 * L / (N_plot - 1.0);
    z = - L;
    FILE *myfile = NULL;
    myfile = fopen("data/phi_vs_z.txt", "w");
    fprintf(myfile, "%10s %20s\n", "# z / L", "phi");
    for(i = 0; i < N_plot; ++i) {
        worker_phi(&phi, z);
        fprintf(myfile, "%10.6lf %20.6lf\n", z / L, phi);
        z = z + h_plot;
    }
    fclose(myfile);
}


void write_phi_num_vs_z(double **phi_matrix, double *z_array, const char *phi_num_vs_z) {
    int i;
    FILE *myfile = NULL;
    myfile = fopen(phi_num_vs_z, "w");
    fprintf(myfile, "%10s %20s\n", "# z / L", "phi_num");
    for(i = 0; i < N_z; ++i) {
        fprintf(myfile, "%10.6lf %20.6lf\n", z_array[i] / L, phi_matrix[i][0]);
    }
    fclose(myfile);
}


void write_u_z_vs_z(double mu_0, double n, double t, const char *u_z_vs_z) {
    int i;
    double h_plot, rho, u_z, z;
    h_plot = 2.0 * L / (N_plot - 1.0);
    z = - L;
    rho = 0.0;
    FILE *myfile = NULL;
    myfile = fopen(u_z_vs_z, "w");
    fprintf(myfile, "%10s %20s\n", "# z / L", "u_z, mm/s");
    for(i = 0; i < N_plot; ++i) {
        worker_u_z(mu_0, n, rho, t, &u_z, z);
        fprintf(myfile, "%10.6lf %20.6lf\n", z / L, u_z * 1000.0);
        z = z + h_plot;
    }
    fclose(myfile);
}


void write_u_z_num_vs_rho(int i, double h_rho, double h_z, double mu_0, double n, double **phi_matrix, double **Psi_matrix, double *rho_array, double t, double **u_z_matrix, double *z_array, const char *u_z_vs_rho, const char * u_z_num_vs_rho) {
    int j;
    double h_plot, rho, u_z, a_U_array[N_z], der_z_U_array[N_z];
    rho = 0.0;
    h_plot = D / (2.0 * (N_plot - 1.0));
    FILE *myfile_1 = NULL;
    myfile_1 = fopen(u_z_vs_rho, "w");
    fprintf(myfile_1, "%10s %20s\n", "# 2 rho / D", "u_z, mm/s");
    for(j = 0; j < N_plot; ++j) {
        worker_u_z(mu_0, n, rho, t, &u_z, z_array[i]);
        fprintf(myfile_1, "%10.6lf %20.6lf\n", 2.0 * rho / D, u_z * 1000.0);
        rho = rho + h_plot;
    }
    fclose(myfile_1);
    worker_a_U_array(a_U_array, mu_0, n, t, z_array);
    worker_der_z_U_array(der_z_U_array, mu_0, n, t, z_array);
    worker_u_z_Psi_matrix(a_U_array, der_z_U_array, h_rho, h_z, mu_0, n, phi_matrix, Psi_matrix, rho_array, t, u_z_matrix, z_array);
    FILE *myfile_2 = NULL;
    myfile_2 = fopen(u_z_num_vs_rho, "w");
    fprintf(myfile_2, "%10s %20s\n", "# 2 rho / L", "u_z_num, mm/s");
    for(j = 0; j < N_rho; ++j) {
        fprintf(myfile_2, "%10.6lf %20.6lf\n", 2.0 * rho_array[j] / D, u_z_matrix[i][j] * 1000.0);
    }
    fclose(myfile_2);
}


void write_u_z_num_vs_t(int i, double h_rho, double h_z, double mu_0, double n, double **phi_matrix, double **Psi_matrix, double *rho_array, double T, double **u_z_matrix, double *z_array) {
    int k;
    double h_plot, t, a_U_array[N_z], der_z_U_array[N_z];
    t = 0.0;
    h_plot = T / (2.0 * (N_plot - 1.0));
    FILE *myfile = NULL;
    myfile = fopen("data/u_z_num_vs_t.txt", "w");
    fprintf(myfile, "%10s %20s\n", "# t, sec", "u_z, mm/s");
    for(k = 0; k < N_plot; ++k) {
        worker_a_U_array(a_U_array, mu_0, n, t, z_array);
        worker_der_z_U_array(der_z_U_array, mu_0, n, t, z_array);
        worker_u_z_Psi_matrix(a_U_array, der_z_U_array, h_rho, h_z, mu_0, n, phi_matrix, Psi_matrix, rho_array, t, u_z_matrix, z_array);
        fprintf(myfile, "%10.6lf %20.6lf\n", t, u_z_matrix[i][0] * 1000.0);
        t = t + h_plot;
    }
    fclose(myfile);
}


void write_u_z_num_vs_z(double h_rho, double h_z, double mu_0, double n, double **phi_matrix, double **Psi_matrix, double *rho_array, double t, double **u_z_matrix, double *z_array, const char *u_z_num_vs_z) {
    int i;
    double a_U_array[N_z], der_z_U_array[N_z];
    worker_a_U_array(a_U_array, mu_0, n, t, z_array);
    worker_der_z_U_array(der_z_U_array, mu_0, n, t, z_array);
    worker_u_z_Psi_matrix(a_U_array, der_z_U_array, h_rho, h_z, mu_0, n, phi_matrix, Psi_matrix, rho_array, t, u_z_matrix, z_array);
    FILE *myfile = NULL;
    myfile = fopen(u_z_num_vs_z, "w");
    fprintf(myfile, "%10s %20s\n", "# z / L", "u_z_num, mm/s");
    for(i = 0; i < N_z; ++i) {
        fprintf(myfile, "%10.6lf %20.6lf\n", z_array[i] / L, u_z_matrix[i][0] * 1000.0);
    }
    fclose(myfile);
}


void write_vector_field(double mu_0, double n) {
    int i, j;
    double h_vf_rho, h_vf_z, t, z_max, z_array[N_vf_z], rho_matrix[N_vf_z][N_vf_rho], u_rho_matrix[N_vf_z][N_vf_rho], u_z_matrix[N_vf_z][N_vf_rho], z_matrix[N_vf_z][N_vf_rho];
    z_max = 0.1 * L;
    h_vf_rho = D / (N_vf_rho - 1.0);
    h_vf_z = 2.0 * z_max / (N_vf_z - 1.0);
    t = 0.0;
    FILE *myfile_1 = NULL;
    FILE *myfile_2 = NULL;
    FILE *myfile_3 = NULL;
    FILE *myfile_4 = NULL;
    myfile_1 = fopen("data/vector_field/z_matrix.txt", "w");
    myfile_2 = fopen("data/vector_field/rho_matrix.txt", "w");
    myfile_3 = fopen("data/vector_field/u_z_matrix.txt", "w");
    myfile_4 = fopen("data/vector_field/u_rho_matrix.txt", "w");
    for(i = 0; i < N_vf_z; ++i) {
        z_array[i] = - z_max + i * h_vf_z;
    }
    for(i = 0; i < N_vf_z; ++i) {
        for(j = 0; j < N_vf_rho; ++j) {
            z_matrix[i][j] = z_array[i];
            rho_matrix[i][j] = - D / 2.0 + j * h_vf_rho;
            worker_u_z(mu_0, n, rho_matrix[i][j], t, &u_z_matrix[i][j], z_matrix[i][j]);
            worker_u_rho(mu_0, n, rho_matrix[i][j], t, &u_rho_matrix[i][j], z_matrix[i][j]);
        }
    }
    for(i = 0; i < N_vf_z; ++i) {
        for(j = 0; j < N_vf_rho; ++j) {
            fprintf(myfile_1, "%15.10lf", z_matrix[i][j] / L);
            fprintf(myfile_2, "%15.10lf", 2.0 * rho_matrix[i][j] / D);
            fprintf(myfile_3, "%15.10lf", u_z_matrix[i][j]);
            fprintf(myfile_4, "%15.10lf", u_rho_matrix[i][j]);

        }
        fprintf(myfile_1, "\n");
        fprintf(myfile_2, "\n");
        fprintf(myfile_3, "\n");
        fprintf(myfile_4, "\n");
    }
    fclose(myfile_1);
    fclose(myfile_2);
    fclose(myfile_3);
    fclose(myfile_4);
}
