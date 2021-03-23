// nonlinear theory of macroscopic flow induced in a drop of ferrofluids
// numerical algorithms


#include "header.h"


// convergence condition
_Bool converge(double **x, double **p) {
    int i, j;
    double norm;
    norm = 0.0;
    for (i = 0; i < N_z; ++i) {
        for (j = 0; j < N_rho; ++j) {
            norm += int_pow(x[i][j] - p[i][j], 2);
        }
    }
    return (sqrt(norm) > error);
}


// exponentiation by squaring
double int_pow(double x, int y) {
    double res = 1.0;
    while (y) {
        if (y & 1) {
            res *= x;
        }
        y >>= 1;
        x *= x;
    }
    return res;
}


void gauss_seidel_method(double h_rho, double h_z, double **p, double **phi_matrix, double **phi_up_matrix, double *rho_array, double **u_rho_matrix, double **u_z_matrix) {
    int i, iter, j;
    iter = 0;
    do {
        for(i = 0; i < N_z; ++i) {
            for(j = 0; j < N_rho; ++j) {
                p[i][j] = phi_up_matrix[i][j];
            }
        }
        phi_up_matrix[0][0] = 1.0 / (1.0 + 2.0 * k_D * tau / int_pow(h_rho, 2) + 2.0 * k_D * tau / int_pow(h_z, 2)) * phi_matrix[0][0] + 1.0 / (int_pow(h_rho, 2) / (2.0 * k_D * tau) + 1.0 + int_pow(h_rho, 2) / int_pow(h_z, 2)) * phi_up_matrix[0][1] + 1.0 / (int_pow(h_z, 2) / (2.0 * k_D * tau) + 1.0 * int_pow(h_z, 2) / int_pow(h_rho, 2) + 1.0) * phi_up_matrix[1][0];
        //-----------------------------
        phi_up_matrix[N_z - 1][0] = 1.0 / (1.0 + 2.0 * k_D * tau / int_pow(h_rho, 2) + 2.0 * k_D * tau / int_pow(h_z, 2)) * phi_matrix[N_z - 1][0] + 1.0 / (int_pow(h_rho, 2) / (2.0 * k_D * tau) + 1.0 + int_pow(h_rho, 2) / int_pow(h_z, 2)) * phi_up_matrix[N_z - 1][1] + 1.0 / (int_pow(h_z, 2) / (2.0 * k_D * tau) + 1.0 * int_pow(h_z, 2) / (int_pow(h_rho, 2)) + 1.0) * phi_up_matrix[N_z - 2][0];
        for(j = 1; j < N_rho - 1; ++j) {
            phi_up_matrix[0][j] = 1.0 / (1.0 + 2.0 * k_D * tau / int_pow(h_rho, 2) + 2.0 * k_D * tau / int_pow(h_z, 2)) * phi_matrix[0][j] + 1.0 / (int_pow(h_rho, 2) / (k_D * tau) + 2.0 + 2.0 * int_pow(h_rho, 2) / int_pow(h_z, 2)) * (phi_up_matrix[0][j - 1] + phi_up_matrix[0][j + 1]) + 1.0 / (int_pow(h_z, 2) / (2.0 * k_D * tau) + int_pow(h_z, 2) / int_pow(h_rho, 2) + 1.0) * phi_up_matrix[1][j];
            //--------------------------------------------------------------------
            phi_up_matrix[N_z - 1][j] = 1.0 / (1.0 + 2.0 * k_D * tau / int_pow(h_rho, 2) + 2.0 * k_D * tau / int_pow(h_z, 2)) * phi_matrix[N_z - 1][j] + 1.0 / (int_pow(h_rho, 2) / (k_D * tau) + 2.0 + 2.0 * int_pow(h_rho, 2) / int_pow(h_z, 2)) * (phi_up_matrix[N_z - 1][j - 1] + phi_up_matrix[N_z - 1][j + 1]) + 1.0 / (int_pow(h_z, 2) / (2.0 * k_D * tau) + int_pow(h_z, 2) / int_pow(h_rho, 2) + 1.0) * phi_up_matrix[N_z - 2][j];
        }
        phi_up_matrix[0][N_rho - 1] = 1.0 / (1.0 + 2.0 * k_D * tau / int_pow(h_rho, 2) + 2.0 * k_D * tau / int_pow(h_z, 2)) * phi_matrix[0][N_rho - 1] + 1.0 / (int_pow(h_rho, 2) / (2.0 * k_D * tau) + 1.0 + int_pow(h_rho, 2) / int_pow(h_z, 2)) * phi_up_matrix[0][N_rho - 2] + 1.0 / (int_pow(h_z, 2) / (2.0 * k_D * tau) + 1.0 * int_pow(h_z, 2) / int_pow(h_rho, 2) + 1.0) * phi_up_matrix[1][N_rho - 1];
        //---------------------------------
        phi_up_matrix[N_z - 1][N_rho - 1] = 1.0 / (1.0 + 2.0 * k_D * tau / int_pow(h_rho, 2) + 2.0 * k_D * tau / int_pow(h_z, 2)) * phi_matrix[N_z - 1][N_rho - 1] + 1.0 / (int_pow(h_rho, 2) / (2.0 * k_D * tau) + 1.0 + int_pow(h_rho, 2) / int_pow(h_z, 2)) * phi_up_matrix[N_z - 1][N_rho - 2] + 1.0 / (int_pow(h_z, 2) / (2.0 * k_D * tau) + 1.0 * int_pow(h_z, 2) / int_pow(h_rho, 2) + 1.0) * phi_up_matrix[N_z - 2][N_rho - 1];
        for(i = 1; i < N_z - 1; ++i) {
            phi_up_matrix[i][0] = 1.0 / (1.0 + 2.0 * k_D * tau / int_pow(h_rho, 2) + 2.0 * k_D * tau / int_pow(h_z, 2)) * phi_matrix[i][0] + 1.0 / (int_pow(h_rho, 2) / (2.0 * k_D * tau) + 1.0 + int_pow(h_rho, 2) / int_pow(h_z, 2)) * phi_up_matrix[i][1] + 1.0 / (h_z / tau + 2.0 * k_D * h_z / int_pow(h_rho, 2) + 2.0 * k_D / h_z) * ((k_D / h_z + u_z_matrix[i][0] / 2.0) * phi_up_matrix[i - 1][0] + (k_D / h_z - u_z_matrix[i][0] / 2.0) * phi_up_matrix[i + 1][0]);
            for(j = 1; j < N_rho - 1; ++j) {
                phi_up_matrix[i][j] = 1.0 / (1.0 + 2.0 * k_D * tau / int_pow(h_rho, 2) + 2.0 * k_D * tau / int_pow(h_z, 2)) * phi_matrix[i][j] + 1.0 / (h_rho / tau + 2.0 * k_D / h_rho + 2.0 * k_D * h_rho / int_pow(h_z, 2)) * ((k_D / h_rho - k_D / (2.0 * rho_array[j]) + u_rho_matrix[i][j] / 2.0) * phi_up_matrix[i][j - 1] + (k_D / h_rho + k_D / (2.0 * rho_array[j]) - u_rho_matrix[i][j] / 2.0) * phi_up_matrix[i][j + 1]) + 1.0 / (h_z / tau + 2.0 * k_D * h_z / int_pow(h_rho, 2) + 2.0 * k_D / h_z) * ((k_D / h_z + u_z_matrix[i][j] / 2.0) * phi_up_matrix[i - 1][j] + (k_D / h_z - u_z_matrix[i][j] / 2.0) * phi_up_matrix[i + 1][j]);
            }
            phi_up_matrix[i][N_rho - 1] = 1.0 / (1.0 + 2.0 * k_D * tau / int_pow(h_rho, 2) + 2.0 * k_D * tau / int_pow(h_z, 2)) * phi_matrix[i][N_rho - 1] + 1.0 / (int_pow(h_rho, 2) / (2.0 * k_D * tau) + 1.0 + int_pow(h_rho, 2) / int_pow(h_z, 2)) * phi_up_matrix[i][N_rho - 2] + 1.0 / (int_pow(h_z, 2) / (k_D * tau) + 2.0 * int_pow(h_z, 2) / int_pow(h_rho, 2) + 2.0) * (phi_up_matrix[i - 1][N_rho - 1] + phi_up_matrix[i + 1][N_rho - 1]);
        }
        iter = iter + 1;
        if(iter > 1000) {
            printf("stop linear: iter > 1000\n\n");
            exit(1);
        }
    }
    while (converge(phi_up_matrix, p));
}


// numerical solution of the nonlinear system
void nonlinear_system(double *a_U_array, double *der_z_U_array, double h_rho, double h_z, double mu_0, double n, double **p, double **phi_matrix, double **phi_up_matrix, double **Psi_matrix, double *rho_array, double t, double **u_rho_matrix, double **u_z_matrix, double *z_array) {
    int i, iter, j;
    iter = 0;
    do {
        for(i = 0; i < N_z; ++i) {
            for(j = 0; j < N_rho; ++j) {
                p[i][j] = phi_up_matrix[i][j];
            }
        }
        worker_u_z_Psi_matrix(a_U_array, der_z_U_array, h_rho, h_z, mu_0, n, phi_matrix, Psi_matrix, rho_array, t, u_z_matrix, z_array);
        worker_u_rho_matrix(h_z, Psi_matrix, u_rho_matrix);
        gauss_seidel_method(h_rho, h_z, p, phi_matrix, phi_up_matrix, rho_array, u_rho_matrix, u_z_matrix);
        iter = iter + 1;
        if(iter > 1000) {
            printf("stop nonlinear: iter > 1000\n\n");
            exit(1);
        }
    }
    while (converge(phi_up_matrix, p));
}


// numerical integration
void numerical_integration(int i, double *f, double h_rho, double *int_array) {
    int j;
    int_array[0] = 0.0;
    int_array[1] = h_rho * f[1] / 2.0;
    int_array[2] = h_rho / 3.0 * (4.0 * f[1] + f[2]);
    for(j = 3; j < N_rho; ++j) {
        int_array[j] = int_array[j - 2] + h_rho / 3.0 * (f[j - 2] + 4.0 * f[j - 1] + f[j]);
    }
}
