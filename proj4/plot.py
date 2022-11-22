import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats


def equilibration_time(filename1, filename2, filename3, filename4):
    data1 = pd.read_csv(filename1, header = None, skiprows = 1)
    data2 = pd.read_csv(filename2, header = None, skiprows = 1)
    data3 = pd.read_csv(filename3, header = None, skiprows = 1)
    data4 = pd.read_csv(filename4, header = None, skiprows = 1)
    cycles= []
    e1, e2, e3, e4 = [], [], [], []
    m1, m2, m3, m4 = [], [], [], []

    for i in range(0, len(data1[0])):
        a = data1[0][i].split()
        b = data2[0][i].split()
        c = data3[0][i].split()
        d = data4[0][i].split()
        cycles.append(a[0])
        e1.append(float(a[1]))
        e2.append(float(b[1]))
        e3.append(float(c[1]))
        e4.append(float(d[1]))
        m1.append(float(a[2]))
        m2.append(float(b[2]))
        m3.append(float(c[2]))
        m4.append(float(d[2]))

    eps1, eps2, eps3, eps4 = np.array(e1), np.array(e2), np.array(e3), np.array(e4)
    mag1, mag2, mag3, mag4 = np.array(m1), np.array(m2), np.array(m3), np.array(m4)
    x = np.array(cycles)

    fig, ax = plt.subplots(2, 2, figsize=(11,6))
    ax[0, 0].plot(x, eps1, label="<ϵ> for T=1.0, not ordered")
    ax[0, 0].plot(x, eps2, label="<ϵ> for T=1.0, ordered")
    ax[1, 0].plot(x, eps3, label="<ϵ> for T=2.4, not ordered")
    ax[1, 0].plot(x, eps4, label="<ϵ> for T=2.4, ordered")

    ax[0, 0].set_xscale('linear')
    ax[1, 0].set_xscale('linear')
    ax[0, 0].legend()
    ax[1, 0].legend()
    #ax[0].set_title('Burn in time for <ϵ> measured in MC cycles')
    #ax[1].set_title('Burn in time for <ϵ> measured in MC cycles')
    ax[0, 0].set_xlabel('Cycles')
    ax[1, 0].set_xlabel('Cycles')
    ax[0, 0].grid()
    ax[1, 0].grid()
    ax[0, 0].set_ylabel('<ϵ> [$J$]')
    ax[1, 0].set_ylabel('<ϵ> [$J$]')

    ax[0, 1].plot(x, mag1, label="<|m|> for T=1.0, not ordered")
    ax[0, 1].plot(x, mag2, label="<|m|> for T=1.0, ordered")
    ax[1, 1].plot(x, mag3, label="<|m|> for T=2.4, not ordered")
    ax[1, 1].plot(x, mag4, label="<|m|> for T=2.4, ordered")

    ax[0, 1].set_xscale('linear')
    ax[1, 1].set_xscale('linear')
    ax[0, 1].legend()
    ax[1, 1].legend()
    #ax[0].set_title('Burn in time for <|m|> measured in MC cycles')
    #ax[1].set_title('Burn in time for <|m|> measured in MC cycles')
    ax[0, 1].set_xlabel('Cycles')
    ax[1, 1].set_xlabel('Cycles')
    ax[0, 1].grid()
    ax[1, 1].grid()
    ax[0, 1].set_ylabel('<|m|> [$1$]')
    ax[1, 1].set_ylabel('<|m|> [$1$]')
    fig.savefig("Equilibration_time.pdf")
    #plt.show()



def histograms(filename1, filename2):
    E_lowT = pd.read_csv("probability_density_T1_L20.csv", header = None)
    E_highT = pd.read_csv("probability_density_T2.4_L20.csv", header = None)


    e1, e2 = [], []

    for i in range(0, len(E_lowT[0])):
        a = E_lowT[0][i].split()
        b = E_highT[0][i].split()
        e1.append(float(a[1]))
        e2.append(float(b[1]))

    E_lowT, E_highT = np.array(e1), np.array(e2)


    hist_lowT, bin_edges_lowT = np.histogram(E_lowT, bins=np.arange(np.min(E_lowT)-4, np.max(E_lowT)-4 + 8, 8))
    hist_highT, bin_edges_highT = np.histogram(E_highT, bins=np.arange(np.min(E_highT)-4, np.max(E_highT)-4 + 8, 8))


    prob_lowT = hist_lowT/float(np.sum(hist_lowT))
    prob_highT = hist_highT/float(np.sum(hist_highT))

    bin_midpoints_lowT = (bin_edges_lowT[1:]+bin_edges_lowT[:-1])/2.0
    bin_midpoints_highT = (bin_edges_highT[1:]+bin_edges_highT[:-1])/2.0

    bin_width_lowT = bin_edges_lowT[1] - bin_edges_lowT[0]
    bin_width_highT = bin_edges_highT[1] - bin_edges_highT[0]


    plt.figure(1)
    plt.bar(bin_midpoints_lowT, prob_lowT, width=bin_width_lowT, label='$T$ = 1.0 $J/k$')
    plt.bar(bin_midpoints_highT, prob_highT, width=bin_width_highT, label='$T$ = 2.4 $J/k$')
    #plt.title("Normalized histogram of ϵ for L=20 systems with T=1.0 and T=2.4")
    plt.xlabel(r'System energy $<ϵ>$ [$J$]')
    plt.ylabel(r'Probability p(ϵ;T)')
    plt.yticks(np.arange(0, 1, 0.1))
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig("Probability_histogram.pdf")
    #plt.show()


def phase_transitions(f1, f2, f3, f4):
    data1 = pd.read_csv(f1, header = None)
    data2 = pd.read_csv(f2, header = None)
    data3 = pd.read_csv(f3, header = None)
    data4 = pd.read_csv(f4, header = None)

    e1, e2, e3, e4 = [], [], [], []
    m1, m2, m3, m4 = [], [], [], []
    cv1, cv2, cv3, cv4 = [], [], [], []
    chi1, chi2, chi3, chi4 = [], [], [], []
    t = []


    for i in range(0,len(data1[0])):
        a = data1[0][i].split()
        b = data2[0][i].split()
        c = data3[0][i].split()
        d = data4[0][i].split()

        t.append(float(a[0]))
        e1.append(float(a[1]))
        e2.append(float(b[1]))
        e3.append(float(c[1]))
        e4.append(float(d[1]))

        m1.append(float(a[2]))
        m2.append(float(b[2]))
        m3.append(float(c[2]))
        m4.append(float(d[2]))

        cv1.append(float(a[3]))
        cv2.append(float(b[3]))
        cv3.append(float(c[3]))
        cv4.append(float(d[3]))

        chi1.append(float(a[4]))
        chi2.append(float(b[4]))
        chi3.append(float(c[4]))
        chi4.append(float(d[4]))

    temp = np.array(t)
    eps1, eps2, eps3, eps4 = np.array(e1), np.array(e2), np.array(e3), np.array(e4)
    mag1, mag2, mag3, mag4 = np.array(m1), np.array(m2), np.array(m3), np.array(m4)
    c_v1, c_v2, c_v3, c_v4 = np.array(cv1), np.array(cv2), np.array(cv3), np.array(cv4)
    chi_1, chi_2, chi_3, chi_4 = np.array(chi1), np.array(chi2), np.array(chi3), np.array(chi4)

    fig, ax = plt.subplots(2, 2, figsize=(11,6))

    ax[0, 0].plot(temp, eps1, label="L=40")
    ax[0, 0].plot(temp, eps2, label="L=60")
    ax[0, 0].plot(temp, eps3, label="L=80")
    ax[0, 0].plot(temp, eps4, label="L=100")
    ax[0, 0].legend()
    ax[0, 0].grid()
    ax[0, 0].set_ylabel("Energy <ϵ> [$J$]")
    ax[0, 0].set_xlabel("Temperature []$J/k_B$]")

    ax[1, 0].plot(temp, mag1, label="L=40")
    ax[1, 0].plot(temp, mag2, label="L=60")
    ax[1, 0].plot(temp, mag3, label="L=80")
    ax[1, 0].plot(temp, mag4, label="L=100")
    ax[1, 0].legend()
    ax[1, 0].grid()
    ax[1, 0].set_ylabel("Magnetization <|m|> [$1$]")
    ax[1, 0].set_xlabel("Temperature [$J/k_B$]")

    ax[0, 1].plot(temp, c_v1, label="L=40")
    ax[0, 1].plot(temp, c_v2, label="L=60")
    ax[0, 1].plot(temp, c_v3, label="L=80")
    ax[0, 1].plot(temp, c_v4, label="L=100")
    ax[0, 1].legend()
    ax[0, 1].grid()
    ax[0, 1].set_ylabel("C_v [$J^2/K_B$]")
    ax[0, 1].set_xlabel("Temperature [$J/k_B$]")

    ax[1, 1].plot(temp, chi_1, label="L=40")
    ax[1, 1].plot(temp, chi_2, label="L=60")
    ax[1, 1].plot(temp, chi_3, label="L=80")
    ax[1, 1].plot(temp, chi_4, label="L=100")
    ax[1, 1].legend()
    ax[1, 1].grid()
    ax[1, 1].set_ylabel("χ [$1/K_B$]")
    ax[1, 1].set_xlabel("Temperature [$J/k_B$]")

    fig.savefig("phase_transitions.pdf")


histograms("probability_density_T1_L20.csv", "probability_density_T2.4_L20.csv")
equilibration_time("Equilibrium_time_study_L20_T1.0_not_ordered.csv","Equilibrium_time_study_L20_T1.0_ordered.csv", "Equilibrium_time_study_L20_T2.4_not_ordered.csv","Equilibrium_time_study_L20_T2.4_ordered.csv")
phase_transitions("Paralellization_L40.csv", "Paralellization_L60.csv", "Paralellization_L80.csv", "Paralellization_L100.csv")
