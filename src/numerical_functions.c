// nonlinear theory of macroscopic flow induced in a drop of ferrofluid
// numerical functions


#include "header.h"


// organising of the numerical algorithms
void part_diff_eq(double mu_0, double n, double T, double *t_array) {
    int i, i_1, i_2, j, k;
    double h_rho, h_z, int_phi, int_phi_0, t[1], t_max_total, rho_array[N_rho], t_max[2], z_array[N_z],
        **p = (double**)malloc(N_z * N_rho * sizeof(double)),
        **phi_matrix = (double**)malloc(N_z * N_rho * sizeof(double)),
        **phi_up_matrix = (double**)malloc(N_z * N_rho * sizeof(double)),
        **Psi_matrix = (double**)malloc(N_z * N_rho * sizeof(double)),
        **u_rho_matrix = (double**)malloc(N_z * N_rho * sizeof(double)),
        **u_z_matrix = (double**)malloc(N_z * N_rho * sizeof(double));
    const char *phi_num_vs_z[2] = {"data/phi_num_1_vs_z.txt",
                                   "data/phi_num_2_vs_z.txt"};
    const char *u_z_vs_rho[4] = {"data/u_z_vs_rho_1.txt",
                                 "data/u_z_vs_rho_2.txt",
                                 "data/u_z_vs_rho_3.txt",
                                 "data/u_z_vs_rho_4.txt"};
    const char *u_z_num_vs_rho[4] = {"data/u_z_num_vs_rho_1.txt",
                                     "data/u_z_num_vs_rho_2.txt",
                                     "data/u_z_num_vs_rho_3.txt",
                                     "data/u_z_num_vs_rho_4.txt"};
    const char *u_z_num_vs_z[4] = {"data/u_z_num_vs_z_1.txt",
                                   "data/u_z_num_vs_z_2.txt",
                                   "data/u_z_num_vs_z_3.txt",
                                   "data/u_z_num_vs_z_4.txt"};
    for(i = 0; i < N_z; ++i) {
        p[i] = (double*)malloc(N_rho * sizeof(double));
        phi_matrix[i] = (double*)malloc(N_rho * sizeof(double));
        phi_up_matrix[i] = (double*)malloc(N_rho * sizeof(double));
        Psi_matrix[i] = (double*)malloc(N_rho * sizeof(double));
        u_rho_matrix[i] = (double*)malloc(N_rho * sizeof(double));
        u_z_matrix[i] = (double*)malloc(N_rho * sizeof(double));
    }
    t[0] = 0.0;
    t_max[0] = 20.0; // time = 20 sec
    //t_max[0] = 0.5;
    t_max[1] = 60.0; // time = 60 sec
    //t_max[1] = 1.0;
    t_max_total = t_max[1];
    worker_h_z(&h_z);
    worker_h_rho(&h_rho);
    worker_z_array(h_z, z_array);
    worker_rho_array(h_rho, rho_array);
    worker_phi_matrix(phi_matrix, z_array);
    worker_int_phi(h_rho, h_z, &int_phi, phi_matrix, rho_array);
    int_phi_0 = int_phi;
    for(i = 1; i < N_z; ++i) {
        for(j = 1; j < N_rho; ++j) {
            phi_up_matrix[i][j] = phi_matrix[i][j];
        }
    }
    printf("sigma = %.6lf\n", sigma);
    i_1 = 220;
    //i_1 = 55;
    printf("z_array[i_1] = %.6lf\n", z_array[i_1]);
    i_2 = 200;
    //i_2 = 50;
    printf("z_array[i_2] = %.6lf\n\n", z_array[i_2]);
    printf("numerical solution of the partial differential equation ...\n");
    for(k = 0; k < 2; ++k) {
        time_steps(h_rho, h_z, int_phi, int_phi_0, mu_0, n, p, phi_matrix, phi_up_matrix, Psi_matrix, rho_array, t, t_max[k], t_max_total, u_rho_matrix, u_z_matrix, z_array);
        write_phi_num_vs_z(phi_matrix, z_array, phi_num_vs_z[k]);
    }
    printf("numerical solution of the partial differential equation complete\n\n");
    for(k = 0; k < 4; ++k) {
        write_u_z_num_vs_z(h_rho, h_z, mu_0, n, phi_matrix, Psi_matrix, rho_array, t_array[k], u_z_matrix, z_array, u_z_num_vs_z[k]);
        write_u_z_num_vs_rho(i_1, h_rho, h_z, mu_0, n, phi_matrix, Psi_matrix, rho_array, t_array[k], u_z_matrix, z_array, u_z_vs_rho[k], u_z_num_vs_rho[k]);
        write_u_z_num_vs_t(i_2, h_rho, h_z, mu_0, n, phi_matrix, Psi_matrix, rho_array, T, u_z_matrix, z_array);
    }
    free(p);
    free(phi_matrix);
    free(phi_up_matrix);
    free(Psi_matrix);
    free(u_rho_matrix);
    free(u_z_matrix);
}


void time_steps(double h_rho, double h_z, double int_phi, double int_phi_0, double mu_0, double n, double **p, double **phi_matrix, double **phi_up_matrix, double **Psi_matrix, double *rho_array, double *t, double t_max, double t_max_total, double **u_rho_matrix, double **u_z_matrix, double *z_array) {
    int i, j, k;
    double a_U_array[N_z], der_z_U_array[N_z];
    k = 0;
    while(t[0] < t_max) {
        worker_a_U_array(a_U_array, mu_0, n, t[0], z_array);
        worker_der_z_U_array(der_z_U_array, mu_0, n, t[0], z_array);
        nonlinear_system(a_U_array, der_z_U_array, h_rho, h_z, mu_0, n, p, phi_matrix, phi_up_matrix, Psi_matrix, rho_array, t[0], u_rho_matrix, u_z_matrix, z_array);
        t[0] = t[0] + tau;
        for(i = 0; i < N_z; ++i) {
            for(j = 0; j < N_rho; ++j) {
                phi_matrix[i][j] = phi_up_matrix[i][j];
            }
        }
        k = k + 1;
        if(k % 1000 == 0) {
            worker_int_phi(h_rho, h_z, &int_phi, phi_matrix, rho_array);
            printf("t = %4.2f sec, int_phi/int_phi_0 = %8.6f, progress = %4.2f %%\n", t[0], int_phi / int_phi_0, t[0] * 100.0 / t_max_total);
        }
    }
}


void worker_a_U_array(double *a_U_array, double mu_0, double n, double t, double *z_array) {
    int i;
    for(i = 0; i < N_z; ++i) {
        worker_a_U(&a_U_array[i], mu_0, n, t, z_array[i]);
    }
}


void worker_alpha_Psi_array(int i, double *alpha_Psi_array, double *der_z_phi_array, double h_rho, double *rho_array) {
    int j;
    double f[N_rho];
    for(j = 1; j < N_rho; ++j) {
        f[j] = der_z_phi_array[j] * rho_array[j];
    }
    numerical_integration(i, f, h_rho, alpha_Psi_array);
}


void worker_B_Psi(int i, double *a_U_array, double *B_Psi, double *delta_Psi_array, double *der_z_U_array, double *epsilon_Psi_array, double *eta_Psi_array, double *zeta_Psi_array) {
    *B_Psi = 32.0 / (eta * int_pow(D, 2)) * (8.0 / int_pow(D, 2) * der_z_U_array[i] * zeta_Psi_array[N_rho - 1] + 8.0 * a_U_array[i] / int_pow(D, 2) * eta_Psi_array[N_rho - 1] - der_z_U_array[i] * delta_Psi_array[N_rho - 1] - a_U_array[i] * epsilon_Psi_array[N_rho - 1]);
}


void worker_beta_Psi_array(int i, double *beta_Psi_array, double h_rho, double **phi_matrix, double *rho_array) {
    int j;
    double f[N_rho];
    for(j = 1; j < N_rho; ++j) {
        f[j] = phi_matrix[i][j] * rho_array[j];
    }
    numerical_integration(i, f, h_rho, beta_Psi_array);
}


void worker_D_Psi(int i, double *a_U_array, double *D_Psi, double *delta_Psi_array, double *der_z_U_array, double *epsilon_Psi_array, double *eta_Psi_array, double *zeta_Psi_array) {
    *D_Psi = 1.0 / eta * (der_z_U_array[i] * delta_Psi_array[N_rho - 1] + a_U_array[i] * epsilon_Psi_array[N_rho - 1] - 16.0 / int_pow(D, 2) * der_z_U_array[i] * zeta_Psi_array[N_rho - 1] - 16.0 * a_U_array[i] / int_pow(D, 2) * eta_Psi_array[N_rho - 1]);
}


void worker_delta_Psi_array(int i, double *beta_Psi_array, double *delta_Psi_array, double h_rho, double *rho_array) {
    int j;
    double f[N_rho];
    for(j = 1; j < N_rho; ++j) {
        f[j] = beta_Psi_array[j] / rho_array[j];
    }
    numerical_integration(i, f, h_rho, delta_Psi_array);
}


void worker_der_z_phi_array(int i, double *der_z_phi_array, double h_z, double **phi_matrix) {
    int j;
    for(j = 0; j < N_rho; ++j) {
        der_z_phi_array[j] = (phi_matrix[i + 1][j] - phi_matrix[i - 1][j]) / (2.0 * h_z);
    }
}


// z-derivative of the potential energy of elongated particles
void worker_der_z_U_array(double *der_z_U_array, double mu_0, double n, double t, double *z_array) {
    int i;
    for(i = 0; i < N_z; ++i) {
        worker_der_z_U(&der_z_U_array[i], mu_0, n, t, z_array[i]);
    }
}


void worker_gamma_psi_array(int i, double *alpha_Psi_array, double *gamma_Psi_array, double h_rho, double *rho_array) {
    int j;
    double f[N_rho];
    for(j = 1; j < N_rho; ++j) {
        f[j] = alpha_Psi_array[j] * rho_array[j];
    }
    numerical_integration(i, f, h_rho, gamma_Psi_array);
}


void worker_epsilon_Psi_array(int i, double *gamma_Psi_array, double *epsilon_Psi_array, double h_rho, double *rho_array) {
    int j;
    double f[N_rho];
    for(j = 1; j < N_rho; ++j) {
        f[j] = gamma_Psi_array[j] / rho_array[j];
    }
    numerical_integration(i, f, h_rho, epsilon_Psi_array);
}


void worker_eta_Psi_array(int i, double *epsilon_Psi_array, double *eta_Psi_array, double h_rho, double *rho_array) {
    int j;
    double f[N_rho];
    for(j = 1; j < N_rho; ++j) {
        f[j] = epsilon_Psi_array[j] * rho_array[j];
    }
    numerical_integration(i, f, h_rho, eta_Psi_array);
}


// step along the coordinate rho
void worker_h_rho(double *h_rho) {
    *h_rho = D / (2.0 * (N_rho - 1.0));
}


// step along the coordinate z
void worker_h_z(double *h_z) {
    *h_z = 2.0 * L / (N_z - 1.0);
}


// total amount of ferrofluid
void worker_int_phi(double h_rho, double h_z, double *int_phi, double **phi_matrix, double *rho_array) {
    int i;
    double Sum, int_phi_rho[N_z];
    worker_int_phi_rho(h_rho, int_phi_rho, phi_matrix, rho_array);
    Sum = int_phi_rho[0] / 2.0;
    for(i = 1; i < N_z - 1; ++i) {
        Sum = Sum + int_phi_rho[i];
    }
    Sum = Sum + int_phi_rho[N_z - 1] / 2.0;
    *int_phi = h_z * Sum;
}


// average value for the coordinate rho of the ferrofluid concentration
void worker_int_phi_rho(double h_rho, double *int_phi_rho, double **phi_matrix, double *rho_array) {
    int i, j;
    double Sum_1[N_z] = {0.0}, Sum_2[N_z] = {0.0}, Sum_3[N_z] = {0.0};
    for(i = 0; i < N_z; ++i) {
        Sum_1[i] = Sum_1[i] + h_rho * phi_matrix[i][0];
        for(j = 1; j < N_rho - 1; ++j) {
            Sum_2[i] = Sum_2[i] + rho_array[j] * phi_matrix[i][j];
        }
        Sum_3[i] = Sum_3[i] + (rho_array[N_rho - 2] + 2.0 * rho_array[N_rho - 1]) * phi_matrix[i][N_rho - 1];
        int_phi_rho[i] = pi * h_rho * (Sum_1[i] / 3.0 + 2.0 * Sum_2[i] + Sum_3[i] / 3.0);
    }
}


// the ferrofluid concentration
void worker_phi_matrix(double **phi_matrix, double *z_array) {
    int i, j;
    for(i = 0; i < N_z; ++i) {
        for(j = 0; j < N_rho; ++j) {
            worker_phi(&phi_matrix[i][j], z_array[i]);
        }
    }
}


// nodes on the coordinate rho
void worker_rho_array(double h_rho, double *rho_array) {
    int i;
    for(i = 0; i < N_rho; ++i) {
        rho_array[i] = i * h_rho;
    }
}


// z-component of the suspension flow velocity and the standard stream function psi
void worker_u_z_Psi_matrix(double *a_U_array, double *der_z_U_array, double h_rho, double h_z, double mu_0, double n, double **phi_matrix, double **Psi_matrix, double *rho_array, double t, double **u_z_matrix, double *z_array) {
    int i, j;
    double B_Psi, D_Psi, alpha_Psi_array[N_rho], beta_Psi_array[N_rho], delta_Psi_array[N_rho], der_z_phi_array[N_rho], gamma_Psi_array[N_rho], epsilon_Psi_array[N_rho], eta_Psi_array[N_rho], zeta_Psi_array[N_rho];
    for(j = 0; j < N_rho; ++j) {
        u_z_matrix[0][j] = 0.0;
        u_z_matrix[N_z - 1][j] = 0.0;
        Psi_matrix[0][j] = 0.0;
        Psi_matrix[N_z - 1][j] = 0.0;
    }
    for(i = 1; i < N_z - 1; ++i) {
        worker_der_z_phi_array(i, der_z_phi_array, h_z, phi_matrix);
        worker_alpha_Psi_array(i, alpha_Psi_array, der_z_phi_array, h_rho, rho_array);
        worker_beta_Psi_array(i, beta_Psi_array, h_rho, phi_matrix, rho_array);
        worker_gamma_psi_array(i, alpha_Psi_array, gamma_Psi_array, h_rho, rho_array);
        worker_delta_Psi_array(i, beta_Psi_array, delta_Psi_array, h_rho, rho_array);
        worker_epsilon_Psi_array(i, gamma_Psi_array, epsilon_Psi_array, h_rho, rho_array);
        worker_zeta_Psi_array(i, delta_Psi_array, h_rho, rho_array, zeta_Psi_array);
        worker_eta_Psi_array(i, epsilon_Psi_array, eta_Psi_array, h_rho, rho_array);
        worker_B_Psi(i, a_U_array, &B_Psi, delta_Psi_array, der_z_U_array, epsilon_Psi_array, eta_Psi_array, zeta_Psi_array);
        worker_D_Psi(i, a_U_array, &D_Psi, delta_Psi_array, der_z_U_array, epsilon_Psi_array, eta_Psi_array, zeta_Psi_array);
        for(j = 0; j < N_rho - 1; ++j) {
            u_z_matrix[i][j] = 1.0 / eta * der_z_U_array[i] * delta_Psi_array[j] + a_U_array[i] / eta * epsilon_Psi_array[j] + B_Psi * int_pow(rho_array[j], 2) / 4.0 + D_Psi;
        }
        u_z_matrix[i][N_rho - 1] = 0.0;
        Psi_matrix[i][0] = 0.0;
        for(j = 2; j < N_rho - 1; ++j) {
            Psi_matrix[i][j] = 1.0 / eta * der_z_U_array[i] * zeta_Psi_array[j] / rho_array[j] + a_U_array[i] / eta * eta_Psi_array[j] / rho_array[j] + B_Psi * int_pow(rho_array[j], 3) / 16.0 + D_Psi * rho_array[j] / 2.0;
        }
        Psi_matrix[i][1] = Psi_matrix[i][2] / 2.0;
        Psi_matrix[i][N_rho - 1] = 0.0;
    }
}


// rho-component of the suspension flow velocity
void worker_u_rho_matrix(double h_z, double **Psi_matrix, double **u_rho_matrix) {
    int i, j;
    for(j = 0; j < N_rho; ++j) {
        u_rho_matrix[0][j] = 0.0;
        u_rho_matrix[N_z - 1][j] = 0.0;
    }
    for(i = 1; i < N_z - 1; ++ i) {
        for(j = 0; j < N_rho; ++j) {
            u_rho_matrix[i][j] = - (Psi_matrix[i + 1][j] - Psi_matrix[i - 1][j]) / (2.0 * h_z);
        }
    }
}


// nodes on the coordinate z
void worker_z_array(double h_z, double *z_array) {
    int i;
    for(i = 0; i < N_z; ++i) {
        z_array[i] = - L + i * h_z;
    }
}


void worker_zeta_Psi_array(int i, double *delta_Psi_array, double h_rho, double *rho_array, double *zeta_Psi_array) {
    int j;
    double f[N_rho];
    for(j = 1; j < N_rho; ++j) {
        f[j] = delta_Psi_array[j] * rho_array[j];
    }
    numerical_integration(i, f, h_rho, zeta_Psi_array);
}
