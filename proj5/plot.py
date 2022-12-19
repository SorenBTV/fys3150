import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pyarma as pa


def prob(f1, f2):
    files = [f1, f2]
    type = ["No slit", "Double slit"]
    t = np.linspace(0,0.008, int(0.008/(2.5e-5))+1)

    fig, ax = plt.subplots()

    for i, j in zip(files, type):
        file = pa.cx_cube()
        file.load(i)
        p = np.array(file)
        p = np.swapaxes(p,0,1)
        z = (p*np.conjugate(p)).real
        p = np.sum(z, axis=(1,2))
        ax.plot(t, np.abs(p-1), label=j)

    ax.set_ylabel(r"$|\sum_{i,j} p_{ij}^{t} - 1$)|")
    ax.set_xlabel(r"$time [1]$")
    ax.legend(), ax.grid(), plt.savefig("figures/probability_dist.pdf")


def two_dimensional(f1):
    file = pa.cx_cube()
    file.load("data/"+f1+".dat")
    p = np.array(file)
    p = np.swapaxes(p, 0, 1)
    z = np.sqrt((p*np.conjugate(p)).real)
    #print(z_d_l.shape)


    t = np.linspace(0, 0.002, int(0.002/2.5e-5)+1)
    dt = 2.5e-5

    t0 = 0
    t1 = np.where(t==0.001)[0][0]
    t2 = np.where(t==0.002)[0][0]
    t_list = [t0, t1, t2]

    fig, ax = plt.subplots(ncols=len(t_list), constrained_layout=True, sharey=False, figsize = (10,3))
    j = 0
    ax[0].set_ylabel(r"$y$", fontsize=10)
    for i in t_list:
        norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z[i]))
        img = ax[j].imshow(z[i], extent=[0,1,0,1])
        ax[j].set_xlabel(r"$x$", fontsize=10)
        cbar = fig.colorbar(img, ax=ax[j], location="top", shrink = 0.6)
        txt = ax[j].text(0.95, 0.95, "t = {:.3e}".format(t[i]), color="white",
                        horizontalalignment="right", verticalalignment="top", fontsize=6)
        j += 1

    #plt.show()
    plt.savefig("figures/" + f1 + ".pdf")


    #Real and imaginary plots
    fig, ax = plt.subplots(ncols=len(t_list), constrained_layout=True, sharey=False, figsize = (10,6))
    j = 0
    ax[0].set_ylabel(r"$y$", fontsize=10)
    for i in t_list:
        norm = matplotlib.cm.colors.Normalize(vmin=np.min(p[i].real), vmax=np.max(p[i].real))
        img = ax[j].imshow(p[i].real, extent=[0,1,0,1], norm=norm)
        ax[j].set_xlabel(r"$x$", fontsize=10)
        cbar = fig.colorbar(img, ax=ax[j], location="top", shrink = 0.6)
        txt = ax[j].text(0.95, 0.95, "t = {:.3e}".format(t[i]), color="white",
                        horizontalalignment="right", verticalalignment="top", fontsize=6)
        j += 1
    #plt.show()
    plt.savefig("figures/" + f1 + "_real.pdf")


    fig, ax = plt.subplots(ncols=len(t_list), constrained_layout=True, sharey=False, figsize = (10,6))
    j = 0
    ax[0].set_ylabel(r"$y$", fontsize=10)
    for i in t_list:
        norm = matplotlib.cm.colors.Normalize(vmin=np.min(p[i].imag), vmax=np.max(p[i].imag))
        img = ax[j].imshow(p[i].imag, extent=[0,1,0,1], norm=norm)
        ax[j].set_xlabel("x", fontsize=10)
        cbar = fig.colorbar(img, ax=ax[j], location="top", shrink = 0.6)
        txt = ax[j].text(0.95, 0.95, "t = {:.3e}".format(t[i]), color="white",
                        horizontalalignment="right", verticalalignment="top", fontsize=6)
        j += 1
    #plt.show()
    plt.savefig("figures/" + f1 + "_imag.pdf")


def one_dimensional(f1, f2, f3):
    files = ["data/"+f1+".dat", "data/"+f2+".dat", "data/"+f3+".dat"]
    types = [f1, f2, f3]
    t = np.linspace(0, 0.002, int(0.002/2.5e-5)+1)
    t_point = np.where(t==0.002)[0][0]

    for file, type in zip(files, types):
        p = pa.cx_cube()
        p.load(file)
        p = np.array(p)
        p = np.swapaxes(p, 0, 1)
        z = (p*np.conjugate(p)).real

        pic = z[-1]
        x = np.linspace(0, 1, len(pic))
        x_point = np.where(x == 0.8)[0][0]
        screen = pic[:, x_point]
        screen /= np.sum(screen)

        plt.plot(x, screen, label=type)
    plt.grid()
    plt.legend()
    plt.xlabel(r"$y$")
    plt.ylabel(r"$p(y|x=0.8;t=0.002)$")
    plt.savefig("figures/detection_probabilities.pdf")
    #plt.show()


prob("data/test7.1.dat", "data/test7.2.dat")
two_dimensional("single_slit")
#two_dimensional("test7.1")
two_dimensional("double_slit")
two_dimensional("triple_slit")
one_dimensional("single_slit", "double_slit", "triple_slit")
