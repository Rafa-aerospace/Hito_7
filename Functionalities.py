from Geometry import Surface_selection
from Aerodynamics import LIM_Application
import plotly.graph_objects as go
import plotly.io as pio

from numpy import zeros

import matplotlib.pyplot as plt
from matplotlib import rc # LaTeX tipography
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True); plt.rc('font', family='serif')

import matplotlib
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)


def Execute(name, app, nx, nphi, alpha, beta, Mach, gamma):

    if app == "Evolución de los Coeficientes con el Mach":
        Coefs_fM(name, alpha, beta, nx, nphi)
    elif app == "Cáculo y representación 3-D":
        Basic(name, alpha, beta, Mach, nx, nphi,gamma)
    elif app == "Convergencia con el número de elementos":
        Convergence(name, Mach, alpha, beta)
    else:
        print("Please introduce one application between the default")



def Coefs_fM(name, alpha, beta, nx, nphi, gamma=1.4): # Para ver cómo varían los coeficientes en función del Mach

    Mach_list = [2, 3, 4, 5, 6, 7, 8, 9, 10]

    g = Surface_selection(name, nx, nphi)

    Coefs = {}

    for Mach in Mach_list:
        Cps,  Coefs[Mach] = LIM_Application(g, Mach, gamma, alpha, beta)

    fig, ax = plt.subplots(1,1, figsize=(9,6), constrained_layout='true')
    ax.set_xlim( Mach_list[0]*0.9 , Mach_list[-1]*1.1 )
    ax.set_ylim( -0.2 , 0.75 )
    ax.set_title("Aerodynamic Coefficients variation with $M$", fontsize=20)
    ax.grid()
    ax.set_xlabel("$M$ [-]" ,fontsize=20)
    ax.set_ylabel("C [-]" ,fontsize=20)

    for key in Coefs:

        ax.plot( key, Coefs[key][0], 'bo-')
        ax.plot( key, Coefs[key][1], 'ks-')
        ax.plot( key, Coefs[key][2], 'rd--')

    # ax.legend(loc=0, fancybox=False ,edgecolor="black", ncol = 1, fontsize=20)
    plt.show()


def Convergence(name, Mach, alpha, beta, gamma=1.4): # Para ver cómo varían los coeficientes en función del Mach

    nx_list =  [200, 400, 600, 800, 1000]
    nphi_list =  [6, 10, 50]

    Coefs_nx = {}

    for i, nx in enumerate(nx_list): # Bucle completo en nx

        Coefs_nx[str(nx)] = []

        t = 0

        for j, nphi in enumerate(nphi_list[-3:]):

            Coefs_nx[str(nx)].append([])

            g = Surface_selection(name, nx, nphi)
            Cps, Coefs_nx[str(nx)][t][:] = LIM_Application(g, Mach, gamma, alpha, beta)

            t = t+1

    fig_convergence_NX, ax = plt.subplots(1,1, figsize=(9,6), constrained_layout='true')
    ax.set_xlim( nx_list[0]*0.95, nx_list[-1]*1.05 )
    ax.set_title("Aerodynamic Coefficients variation with $N_x$", fontsize=20)
    ax.grid()
    ax.set_xlabel("$N_x$ [-]" ,fontsize=20)
    ax.set_ylabel("$C_D$ [-]" ,fontsize=20)

    for key in Coefs_nx:

        ax.plot( int(key), Coefs_nx[key][0][0], 'bo--' )
        ax.plot( int(key), Coefs_nx[key][1][0], 'ks--' )
        ax.plot( int(key), Coefs_nx[key][2][0], 'rd--' )

    plt.show()

    fig_convergence_NX, ax = plt.subplots(1,1, figsize=(9,6), constrained_layout='true')
    ax.set_xlim( nx_list[0]*0.95, nx_list[-1]*1.05 )
    ax.set_title("Aerodynamic Coefficients variation with $N_x$", fontsize=20)
    ax.grid()
    ax.set_xlabel("$N_x$ [-]" ,fontsize=20)
    ax.set_ylabel("$C_L$ [-]" ,fontsize=20)

    for key in Coefs_nx:

        ax.plot( int(key), Coefs_nx[key][0][2], 'bo--' )
        ax.plot( int(key), Coefs_nx[key][1][2], 'ks--' )
        ax.plot( int(key), Coefs_nx[key][2][2], 'rd--' )

    plt.show()


def Basic(name, alpha, beta, Mach, nx, nphi, gamma=1.4):

    g = Surface_selection(name, nx, nphi)

    pio.renderers.default='browser'

    Cps,  Coefs = LIM_Application(g, Mach, gamma, alpha, beta)

    X, Y, Z = [ g.coords[i] for i in range(3) ]

    Cps_ = zeros([len(X), len(X)])
    Cps_[0:-1,0:-1] = Cps[:,:]
    Cps_[-1,0:-1] = Cps[-2,:]
    Cps_[0:-1,-1] = Cps[:,-2]
    Cps_[-1,-1] = 1

    fig = go.Figure(go.Surface(x=X, y=Y, z=Z, surfacecolor=Cps_.transpose()))
    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=-1.5, y=-1.5, z=0.5)
    )
    fig.update_layout(scene_camera=camera, title_text= g.surf_name)

    fig.show()
    return Cps, Coefs




