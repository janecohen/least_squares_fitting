#%% imports and parameters

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 13,
    'axes.titlesize': 15,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 11,
})

#%% plotting first spectrum to show weight

data_theor1 = np.loadtxt("S1.txt")
data_theor_weighted1 = np.loadtxt("S1_weighted.txt")
 
x1 = data_theor1[:, 0]
w1 = data_theor1[:, 1]
S1_theor = data_theor1[:, 2]
S1_theor_weighted = data_theor_weighted1[:, 2]
S1_exp = data_theor_weighted1[:, 3]

plt.figure(figsize=(6,5))

plt.scatter(w1, S1_exp, s=12, facecolors='none', edgecolors="#3764B8", linewidths=0.7, label='Experimental data', zorder=3)
plt.plot(w1, S1_theor, color='#E69F00', linewidth=2.2, label=r'Fit, $W=1$', zorder=2)
plt.plot(w1, S1_theor_weighted, color='#009E73', linewidth=2.2, label=r'Fit, $W=F(\omega)^2$', zorder=1)
plt.xlabel(r'$\omega$ [cm$^{-1}$]', fontsize=13)
plt.ylabel(r'Intensity F($\omega$) [a.u. cm$^{-1}$]', fontsize=13)
plt.title('Least-Squares Fit of Spectrum 1', fontsize=14)

plt.grid(alpha=0.25)
plt.legend(frameon=False, fontsize=11)
plt.tick_params(labelsize=11)

plt.tight_layout()
#plt.savefig("spectra1_plot.pdf", bbox_inches="tight")
plt.show()

# %% plotting all 5 spectrum

# file contents: x, w, spectrum_theor, spectrum_exp

plt.figure(figsize=(6,5))

# spectrum 1
data1 = np.loadtxt("S1_weighted.txt")
x1 = data1[:, 0]
w1 = data1[:, 1]
S1_theor_weighted = data1[:, 2]
S1_exp = data1[:, 3]

plt.scatter(w1, S1_exp, s=4, color='b', label="Spectrum 1 - 1 atm", zorder=1)
plt.plot(w1, S1_theor_weighted, color='k', zorder=2)

# spectrum 2
data2 = np.loadtxt("S2_weighted.txt")
x2 = data2[:, 0]
w2 = data2[:, 1]
S2_theor_weighted = data2[:, 2]
S2_exp = data2[:, 3]

plt.scatter(w2, S2_exp, s=4, color='purple', label="Spectrum 2 - 3 atm", zorder=3)
plt.plot(w2, S2_theor_weighted, color='k', zorder=4)

# spectrum 3
data3 = np.loadtxt("S3_weighted.txt")
x3 = data3[:, 0]
w3 = data3[:, 1]
S3_theor_weighted = data3[:, 2]
S3_exp = data3[:, 3]

plt.scatter(w3, S3_exp, s=4, color='m', label="Spectrum 3 - 6 atm", zorder=5)
plt.plot(w3, S3_theor_weighted, color='k', zorder=6)

# spectrum 4
data4 = np.loadtxt("S4_weighted.txt")
x4 = data4[:, 0]
w4 = data4[:, 1]
S4_theor_weighted = data4[:, 2]
S4_exp = data4[:, 3]

plt.scatter(w4, S4_exp, s=4, color='orange', label="Spectrum 4 - 10 atm", zorder=7)
plt.plot(w4, S4_theor_weighted, color='k', zorder=8)

# spectrum 5
data5 = np.loadtxt("S5_weighted.txt")
x5 = data5[:, 0]
w5 = data5[:, 1]
S5_theor_weighted = data5[:, 2]
S5_exp = data5[:, 3]

plt.scatter(w5, S5_exp, s=4, color='r', label="Spectrum 5 - 15 atm", zorder=9)
plt.plot(w5, S5_theor_weighted, color='k', linewidth=1, zorder=10)

plt.xlabel(r'$\omega$ [cm$^{-1}$]', fontsize=13)
plt.ylabel(r'$F(\omega)$ [a.u. cm$^{-1}$]', fontsize=13)
plt.title('Least-Squares Fit of Spectrum 1-5', fontsize=14)

plt.grid(alpha=0.25)
plt.legend(frameon=False, fontsize=11)
plt.tick_params(labelsize=11)

plt.tight_layout()
plt.savefig("all_spectra_plot.pdf", bbox_inches="tight")
plt.show()


# %% pressure vs max position linear regression

p = np.array([1, 3, 6, 10, 15])  # atm
spectrum = [S1_theor_weighted, S2_theor_weighted,S3_theor_weighted,S4_theor_weighted, S5_theor_weighted]
w_list = [w1,w2,w3,w4,w5]
max_vals = np.zeros(5)
w_max = np.zeros(5)

for i in range(len(spectrum)):
    index = np.argmax(spectrum[i])
    max_vals[i] = spectrum[i][index]
    w_max[i] = w_list[i][index]
    
slope, intercept = np.polyfit(p, w_max, 1)
print("slope:", slope)
print("intercept:", intercept)

lin_fit = p*slope+intercept

plt.figure(figsize=(6,5))

plt.plot(p,lin_fit, color='#E69F00',zorder=0, label=r'$\omega_m$=$\omega_0$+p$\Delta$$\omega$')
plt.scatter(p,w_max, facecolors='none', edgecolors="#3764B8", label="Least-square fit data")
plt.xlabel(r'Pressure (atm)', fontsize=13)
plt.ylabel(r'Maximum position $\omega_m$ [cm$^{-1}$]', fontsize=13)
plt.title('Linear Regression of Maximum Position', fontsize=14)

plt.grid(alpha=0.25)
plt.legend(frameon=False, fontsize=11)
plt.tick_params(labelsize=11)

plt.tight_layout()
plt.savefig("linear_regression_max.pdf", bbox_inches="tight")
plt.show()

# %% pressure vs half width linear regression
gamma = [0.20037941502663698, 0.59804514065912961, 1.1985108530819473, 2.0003349003299271, 3.0094367574964158]
    
slope, intercept = np.polyfit(p, gamma, 1)
print("slope:", slope)
print("intercept:", intercept)

lin_fit = p*slope+intercept

plt.figure(figsize=(6,5))

plt.plot(p,lin_fit, color='#E69F00',zorder=0, label=r'$\gamma$=p$\gamma_1$')
plt.scatter(p,gamma, facecolors='none', edgecolors="#3764B8", label="Least-square fit data")
plt.xlabel(r'Pressure (atm)', fontsize=13)
plt.ylabel(r'Line half-width $\gamma$ [cm$^{-1}$]', fontsize=13)
plt.title('Linear Regression of Line Width', fontsize=14)

plt.grid(alpha=0.25)
plt.legend(frameon=False, fontsize=11)
plt.tick_params(labelsize=11)

plt.tight_layout()
plt.savefig("linear_regression_width.pdf", bbox_inches="tight")
plt.show()
# %%