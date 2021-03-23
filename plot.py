# nonlinear theory of macroscopic flow induced in a drop of ferrofluid
# function graphing


import matplotlib.pyplot as plt
import numpy as np
import timeit
start = timeit.default_timer()


abs_H_vs_z_txt = ['data/abs_H_vs_z_1.txt', 'data/abs_H_vs_z_2.txt', 'data/abs_H_vs_z_3.txt', 'data/abs_H_vs_z_4.txt']
abs_H_vs_z_pdf = ['figs/abs_H_vs_z_1.pdf', 'figs/abs_H_vs_z_2.pdf', 'figs/abs_H_vs_z_3.pdf', 'figs/abs_H_vs_z_4.pdf']
for i in range(4):
    x, y = np.loadtxt(abs_H_vs_z_txt[i], comments = '#', unpack = True)
    fig = plt.figure(figsize = (6.0, 4.5))
    plt.plot(x, y, color = 'black', linewidth = 2.0)
    plt.xlabel('$z$ / $L$', fontsize = 16)
    plt.ylabel('$|H|$, kA/m', fontsize = 16)
    #plt.ylabel('$|H|$, кА/м', fontsize = 16)
    plt.tight_layout()
    plt.grid(True)
    plt.savefig(abs_H_vs_z_pdf[i])
    plt.close('all')


x, y = np.loadtxt('data/beta_ratio_vs_r.txt', comments = '#', unpack = True)
fig = plt.figure(figsize = (6.0, 4.5))
plt.plot(x, y, color = 'black', linewidth = 2.0)
plt.xlabel('$r$', fontsize = 16)
plt.ylabel('$β_{||}/β_⊥$', fontsize = 16)
plt.tight_layout()
plt.grid(True)
plt.savefig('figs/beta_ratio_vs_r.pdf')
plt.close('all')


x, y = np.loadtxt('data/phi_vs_z.txt', comments = '#', unpack = True)
x_num_1, y_num_1  = np.loadtxt('data/phi_num_1_vs_z.txt', comments = '#', unpack = True)
x_num_2, y_num_2  = np.loadtxt('data/phi_num_2_vs_z.txt', comments = '#', unpack = True)
fig = plt.figure(figsize = (6.0, 4.5))
plt.plot(x, y, color = 'black', linewidth = 3.0)
#plt.plot(x_num_1, y_num_1, '--', color = 'black', linewidth = 1.5)
#plt.plot(x_num_2, y_num_2, '-.', color = 'black', linewidth = 1.5)
plt.plot(x_num_1, y_num_1, color = 'blue', linewidth = 1.5)
plt.plot(x_num_2, y_num_2, color = 'red', linewidth = 1.5)
plt.xlabel('$z$ / $L$', fontsize = 16)
plt.ylabel('$φ$', fontsize = 16)
plt.tight_layout()
plt.grid(True)
plt.savefig('figs/phi_vs_z.pdf')
plt.close('all')


u_z_vs_z_txt = ['data/u_z_vs_z_1.txt', 'data/u_z_vs_z_2.txt', 'data/u_z_vs_z_3.txt', 'data/u_z_vs_z_4.txt']
u_z_num_vs_z_txt = ['data/u_z_num_vs_z_1.txt', 'data/u_z_num_vs_z_2.txt', 'data/u_z_num_vs_z_3.txt', 'data/u_z_num_vs_z_4.txt']
u_z_vs_z_pdf = ['figs/u_z_vs_z_1.pdf', 'figs/u_z_vs_z_2.pdf', 'figs/u_z_vs_z_3.pdf', 'figs/u_z_vs_z_4.pdf']
for i in range(4):
    x, y = np.loadtxt(u_z_vs_z_txt[i], comments = '#', unpack = True)
    x_num, u_z_num  = np.loadtxt(u_z_num_vs_z_txt[i], comments = '#', unpack = True)
    fig = plt.figure(figsize = (6.0, 4.5))
    plt.plot(x, y, color = 'black', linewidth = 3.0)
    #plt.plot(x_num, u_z_num, '-.', color = 'black', linewidth = 1.5)
    plt.plot(x_num, u_z_num, color = 'red', linewidth = 1.5)
    plt.xlabel('$z$ / $L$', fontsize = 16)
    #plt.ylabel('$u_z$, mm/s', fontsize = 16)
    plt.ylabel('$u_z$, мм/с', fontsize = 16)
    plt.tight_layout()
    plt.grid(True)
    plt.savefig(u_z_vs_z_pdf[i])
    plt.close('all')


u_z_vs_rho_txt = ['data/u_z_vs_rho_1.txt', 'data/u_z_vs_rho_2.txt', 'data/u_z_vs_rho_3.txt', 'data/u_z_vs_rho_4.txt']
u_z_num_vs_rho_txt = ['data/u_z_num_vs_rho_1.txt', 'data/u_z_num_vs_rho_2.txt', 'data/u_z_num_vs_rho_3.txt', 'data/u_z_num_vs_rho_4.txt']
u_z_vs_rho_pdf = ['figs/u_z_vs_rho_1.pdf', 'figs/u_z_vs_rho_2.pdf', 'figs/u_z_vs_rho_3.pdf', 'figs/u_z_vs_rho_4.pdf']
for i in range(4):
    x, y = np.loadtxt(u_z_vs_rho_txt[i], comments = '#', unpack = True)
    x_num, u_z_num  = np.loadtxt(u_z_num_vs_rho_txt[i], comments = '#', unpack = True)
    fig = plt.figure(figsize = (6.0, 4.5))
    plt.plot(x, y, color = 'black', linewidth = 3.0)
    #plt.plot(x_num, u_z_num, '-.', color = 'black', linewidth = 1.5)
    plt.plot(x_num, u_z_num, color = 'red', linewidth = 1.5)
    plt.xlabel('$2 ρ$ / $D$', fontsize = 16)
    #plt.ylabel('$u_z$, mm/s', fontsize = 16)
    plt.ylabel('$u_z$, мм/с', fontsize = 16)
    plt.tight_layout()
    plt.grid(True)
    plt.savefig(u_z_vs_rho_pdf[i])
    plt.close('all')


z = np.loadtxt('data/vector_field/z_matrix.txt', comments = '#', unpack = True)
rho = np.loadtxt('data/vector_field/rho_matrix.txt', comments = '#', unpack = True)
u_z = np.loadtxt('data/vector_field/u_z_matrix.txt', comments = '#', unpack = True)
u_rho = np.loadtxt('data/vector_field/u_rho_matrix.txt', comments = '#', unpack = True)
fig = plt.figure(figsize = (6.0, 4.5))
plt.quiver(z, rho, u_z, u_rho)
plt.xlabel('$z$ / $L$', fontsize = 16)
plt.ylabel('$2ρ$ / $D$', fontsize = 16)
plt.tight_layout()
plt.savefig('figs/vector_field.pdf')
plt.close('all')


x, y = np.loadtxt('data/u_z_num_vs_t.txt', comments = '#', unpack = True)
fig = plt.figure(figsize = (6.0, 4.5))
plt.plot(x, y, color = 'black', linewidth = 2.0)
plt.axhline(0, color = 'black', linewidth = 1.0)
plt.xlabel('$t$, sec', fontsize = 16)
plt.ylabel('$u_z$, mm/s', fontsize = 16)
plt.tight_layout()
plt.grid(True)
plt.savefig('figs/u_z_num_vs_t.pdf')
plt.close('all')


stop = timeit.default_timer()
runtime = stop - start
print('plot runtime = {:4.2f} sec\n' .format(runtime))
