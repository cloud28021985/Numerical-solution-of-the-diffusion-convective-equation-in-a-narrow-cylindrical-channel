// nonlinear theory of macroscopic flow induced in a drop of ferrofluid
// analytical functions


#include "header.h"


void worker_a_1z(double *a_1z, double z) {
    *a_1z = I_0 / (2.0 * l) * ((L + l + z) / sqrt(int_pow(L + l + z, 2) + int_pow(R, 2)) - (L + z) / sqrt(int_pow(L + z, 2) + int_pow(R, 2)));
}


void worker_a_1rho(double *a_1rho, double z) {
    *a_1rho = I_0 * int_pow(R, 2) / (4.0 * l) * (1.0 / pow(int_pow(z + L, 2) + int_pow(R, 2), 1.5) - 1.0 / pow(int_pow(z + L + l, 2) + int_pow(R, 2), 1.5));
}


void worker_a_2z(double *a_2z, double z) {
    *a_2z = I_0 / (2.0 * l) * ((L + l - z) / sqrt(int_pow(L + l - z, 2) + int_pow(R, 2)) - (L - z) / sqrt(int_pow(L - z, 2) + int_pow(R, 2)));
}


void worker_a_2rho(double *a_2rho, double z) {
    *a_2rho = I_0 * int_pow(R, 2) / (4.0 * l) * (1.0 / pow(int_pow(L + l - z, 2) + int_pow(R, 2), 1.5) - 1.0 / pow(int_pow(L - z, 2) + int_pow(R, 2), 1.5));
}


void worker_a_U(double *a_U, double mu_0, double n, double t, double z) {
    double a_1rho, a_1z, a_2z, a_2rho, b_1z, b_2z;
    worker_a_1z(&a_1z, z);
    worker_a_2z(&a_2z, z);
    worker_a_1rho(&a_1rho, z);
    worker_a_2rho(&a_2rho, z);
    worker_b_1z(&b_1z, z);
    worker_b_2z(&b_2z, z);
    *a_U = mu_0 * chi_p / (1 + chi_p * n) * ((a_1z * b_1z + int_pow(a_1rho, 2)) * int_pow(cos(omega * t), 2) + ((a_1z * b_2z + a_2z * b_1z) / 2.0 + a_1rho * a_2rho) * sin(2.0 * omega * t) + (a_2z * b_2z + int_pow(a_2rho, 2)) * int_pow(sin(omega * t), 2));
}


void worker_A_psi(double *A_psi, double mu_0, double n, double t, double z) {
    double a_U;
    worker_a_U(&a_U, mu_0, n, t, z);
    *A_psi = a_U * phi_0 * z / (eta * int_pow(sigma, 2)) * exp(- int_pow(z, 2) / int_pow(sigma, 2));
}


// absolute value of magnetiс field
void worker_abs_H(double *abs_H, double t, double z) {
    double a_1z, a_2z, h_1z, H_1z, h_2z, H_2z;
    worker_a_1z(&a_1z, z);
    worker_a_2z(&a_2z, z);
    h_1z = a_1z;
    h_2z = a_2z;
    H_1z = h_1z * cos(omega * t);
    H_2z = h_2z * sin(omega * t);
    *abs_H = fabs(H_1z + H_2z);
}


void worker_b_1z(double *b_1z, double z) {
    *b_1z = 3.0 * I_0 * int_pow(R, 2) / (4.0 * l) * ((L + l + z) / pow(int_pow(L + l + z, 2) + int_pow(R, 2), 2.5) - (L + z) / pow(int_pow(L + z, 2) + int_pow(R, 2), 2.5));
}


void worker_b_2z(double *b_2z, double z) {
    *b_2z = 3.0 * I_0 * int_pow(R, 2) / (4.0 * l) * ((L + l - z) / pow(int_pow((L + l - z), 2) + int_pow(R, 2), 2.5) - (L - z) / pow(int_pow((L - z), 2) + int_pow(R, 2), 2.5));
}


// ratio of the components of the mobility tensor
void worker_beta_ratio(double *beta_ratio, double r) {
    *beta_ratio = 2.0 * (2.0 * log(2.0 * r) - 1.0) / (2.0 * log(2.0 * r) + 1.0);
}


void worker_der_z_a_1z(double *der_z_a_1z, double z) {
    *der_z_a_1z = I_0 * int_pow(R, 2) / (2.0 * l) * (1.0 / pow(int_pow(L + l + z, 2) + int_pow(R, 2), 1.5) - 1.0 / pow(int_pow(L + z, 2) + int_pow(R, 2), 1.5));
}


void worker_der_z_a_2z(double *der_z_a_2z, double z) {
    *der_z_a_2z = I_0 * int_pow(R, 2) / (2.0 * l) * (1.0 / pow(int_pow(L - z, 2) + int_pow(R, 2), 1.5) - 1.0 / pow(int_pow(L + l - z, 2) + int_pow(R, 2), 1.5));
}


void worker_der_z_A_psi(double *der_z_A_psi, double mu_0, double n, double t, double z) {
    double a_U;
    worker_a_U(&a_U, mu_0, n, t, z);
    *der_z_A_psi = a_U * phi_0 / (eta * int_pow(sigma, 2)) * (1.0 - 2.0 * int_pow(z, 2) / int_pow(sigma, 2) * exp(- int_pow(z, 2) / int_pow(sigma, 2)));
}


// z-derivative of the potential energy of elongated particles
void worker_der_z_U(double *der_z_U, double mu_0, double n, double t, double z) {
    double a_1z, a_2z, der_z_a_1z, der_z_a_2z;
    worker_a_1z(&a_1z, z);
    worker_a_2z(&a_2z, z);
    worker_der_z_a_1z(&der_z_a_1z, z);
    worker_der_z_a_2z(&der_z_a_2z, z);
    *der_z_U = - mu_0 * chi_p * (a_1z * cos(omega * t) + a_2z * sin(omega * t)) / (1.0 + chi_p * n) * (der_z_a_1z * cos(omega * t) + der_z_a_2z * sin(omega * t));
}


// volume concentration of the particles
void worker_phi(double *phi, double z) {
    *phi = phi_0 * exp(- int_pow(z, 2) / int_pow(sigma, 2));
}


void worker_t_array(double *t_array, double T) {
    int i;
    for(i = 0; i < 4; ++i) {
        t_array[i] = i * T / 8.0;
    }
}


// z-component of the suspension flow velocity at the initial time
void worker_u_z(double mu_0, double n, double rho, double t, double *u_z, double z) {
    double A_psi;
    worker_A_psi(&A_psi, mu_0, n, t, z);
    *u_z = - A_psi / 16.0 * (int_pow(rho, 4) - int_pow(D, 2) * int_pow(rho, 2) / 3.0 + int_pow(D, 4) / 48.0);
}


// rho-component of the suspension flow velocity at the initial time
void worker_u_rho(double mu_0, double n, double rho, double t, double *u_rho, double z) {
    double der_z_A_psi;
    worker_der_z_A_psi(&der_z_A_psi, mu_0, n, t, z);
    *u_rho = rho / 96.0 * der_z_A_psi * (int_pow(rho, 4) - int_pow(D, 2) * int_pow(rho, 2) / 2.0 + int_pow(D, 4) / 16.0);
}
