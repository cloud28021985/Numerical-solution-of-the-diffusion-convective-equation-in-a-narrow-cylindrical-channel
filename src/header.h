// nonlinear theory of macroscopic flow induced in a drop of ferrofluid
// header file


#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define chi_p 25.0 // magnetic susceptibility of these particles
#define eta 0.001 // the viscosity of water
#define D 0.00175 // the the diameter of the cylinder
#define error 1e-10 // error of iterative algorithms
#define I_0 18000.0 // current amplitude
#define k_D 1e-7 // diffusion coefficient
//#define k_D 1e-6 // diffusion coefficient
#define l 0.05 // length of the solenoid
#define L 0.1 // the distance origin and the solenoid
#define N_plot 400 // the number of points for analytical graphs
#define N_rho 100 // the number of nodes along the coordinate rho
//#define N_rho 50 // the number of nodes along the coordinate rho
#define N_z 401 // the number of nodes along the coordinate z
//#define N_z 101 // the number of nodes along the coordinate z
#define N_vf_rho 15 // the parameter for the vector field of fluid velocity
#define N_vf_z 21 // the parameter for the vector field of fluid velocity
#define omega 15.0 // the angular frequency of the alternating magnetic field
#define phi_0 0.05 // the initial volume concentration
#define pi M_PI // the number pi
#define R 0.1 // radius of the solenoid
#define r_0 10.0 // aspect ratio of a particle
#define sigma 0.01 // typical length of the ferroparticles cloud
#define tau 0.0001 // time step
//#define tau 0.0004 // time step


double int_pow(double x, int y);
void nonlinear_system(double *a_U_array, double *der_z_U_array, double h_rho, double h_z, double mu_0, double n, double **p, double **phi_matrix, double **phi_up_matrix, double **Psi_matrix, double *rho_array, double t, double **u_rho_matrix, double **u_z_matrix, double *z_array);
void numerical_integration(int i, double *f, double h_rho, double *int_array);
void part_diff_eq(double mu_0, double n, double T, double *t_array);
void time_steps(double h_rho, double h_z, double int_phi, double int_phi_0, double mu_0, double n, double **p, double **phi_matrix, double **phi_up_matrix, double **Psi_matrix, double *rho_array, double *t, double t_max, double t_max_total, double **u_rho_matrix, double **u_z_matrix, double *z_array);
void worker_a_U(double *a_U, double mu_0, double n, double t, double z);
void worker_a_U_array(double *a_U_array, double mu_0, double n, double t, double *z_array);
void worker_abs_H(double *abs_H, double t, double z);
void worker_b_1z(double *b_1z, double z);
void worker_b_2z(double *b_2z, double z);
void worker_beta_ratio(double *beta_ratio, double r);
void worker_der_z_U(double *der_z_U, double mu_0, double n, double t, double z);
void worker_der_z_U_array(double *der_z_U_array, double mu_0, double n, double t, double *z_array);
void worker_h_rho(double *h_rho);
void worker_h_z(double *h_z);
void worker_int_phi(double h_rho, double h_z, double *int_phi, double **phi_matrix, double *rho_array);
void worker_int_phi_rho(double h_rho, double *int_phi_rho, double **phi_matrix, double *rho_array);
void worker_phi(double *Phi, double z);
void worker_phi_matrix(double **phi_matrix, double *z_array);
void worker_rho_array(double h_rho, double *rho_array);
void worker_t_array(double *t_array, double T);
void worker_u_z(double mu_0, double n, double rho, double t, double *u_z, double z);
void worker_u_z_Psi_matrix(double *a_U_array, double *der_z_U_array, double h_rho, double h_z, double mu_0, double n, double **phi_matrix, double **Psi_matrix, double *rho_array, double t, double **u_z_matrix, double *z_array);
void worker_u_rho(double mu_0, double n, double rho, double t, double *u_rho, double z);
void worker_u_rho_matrix(double h_z, double **Psi_matrix, double **u_rho_matrix);
void worker_z_array(double h_z, double *z_array);
    void worker_zeta_Psi_array(int i, double *delta_Psi_array, double h_rho, double *rho_array, double *zeta_Psi_array);
void write_abs_H_vs_z(double t, const char *abs_H_vs_z);
void write_beta_ratio();
void write_phi_vs_z();
void write_phi_num_vs_z(double **phi_matrix, double *z_array, const char *phi_num_vs_z);
void write_u_z_vs_z(double mu_0, double n, double t, const char *u_z_vs_z);
void write_u_z_num_vs_rho(int i, double h_rho, double h_z, double mu_0, double n, double **phi_matrix, double **Psi_matrix, double *rho_array, double t, double **u_z_matrix, double *z_array, const char *u_z_vs_rho, const char * u_z_num_vs_rho);
void write_u_z_num_vs_t(int i, double h_rho, double h_z, double mu_0, double n, double **phi_matrix, double **Psi_matrix, double *rho_array, double T, double **u_z_matrix, double *z_array);
void write_u_z_num_vs_z(double h_rho, double h_z, double mu_0, double n, double **phi_matrix, double **Psi_matrix, double *rho_array, double t, double **u_z_matrix, double *z_array, const char *u_z_num_vs_z);
void write_vector_field(double mu_0, double n);
