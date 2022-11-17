# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 17:03:33 2022

@author: Rafael Rivero de Nicolás
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc # LaTeX tipography
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.rc('text', usetex=True); plt.rc('font', family='serif')

import matplotlib
matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)

# %%

#### Funciones para ploteo ####
def Plot_geometry_1D(x, y, title=r'2D Geometry', axis_labels=[r'$x$ [-]', r'$z(x)$ [-]']):

    fig, ax = plt.subplots(1,1, figsize=(9,6), constrained_layout='true')
    ax.set_xlim( x[0] , x[-1] );
    ax.set_ylim( -y[0] - x[-1] + x[0], y[0] + x[-1] - x[0] );
    ax.set_title(title, fontsize=20)
    ax.grid();
    ax.set_xlabel(axis_labels[0] ,fontsize=20)
    ax.set_ylabel(axis_labels[1] ,fontsize=20)
    ax.plot( x, y, c="blue" )
    ax.plot( x, -y, c="blue" )
    plt.plot()



def Simple_plot_1D(x, y, title=None, axis_labels=[None, None]):

    fig, ax = plt.subplots(1,1, figsize=(9,6), constrained_layout='true')
    ax.set_xlim( x[0] , x[-1] );
    ax.set_ylim( np.min(y)*0.95, np.max(y)*1.05 );
    ax.set_title(title, fontsize=20)
    ax.grid();
    ax.set_xlabel(axis_labels[0] ,fontsize=20)
    ax.set_ylabel(axis_labels[1] ,fontsize=20)
    ax.plot( x, y, c="blue" )
    plt.plot()
################################



#### Funciones para definir la geometría ####
def Horizontal_Parabola(x): # Libro Anderson, Cap 3, page 63

    return np.sqrt( (x+1)/0.769 )


def Quarter_of_a_Circle(x):

    return np.sqrt( x )

def Flechita(x): # Se la ha inventado rafa

    f1 = x[x<=0.5]
    f2 = 2*x[np.logical_and(x>0.5, x<=0.75)]-0.5
    f3 = 1-4*(x[x>0.75]-0.75)

    return np.hstack((f1, f2, f3))
########################





#### Funciones serias ####
def Derivative(x, f):

    return ( f(x+1E-5)-f(x) ) / (1E-5)


def Betta_calculation(x, f):

    betta = np.arctan( Derivative( x, f ) )

    betta[betta<0] = 0 # betta<0 implies Cp=0 based on Newton Theory

    return betta



def Cp_Dist_1D_by_Newton(x, f, M):

    Plot_geometry_1D(x, f(x))

    Cp = 2*(np.sin(Betta_calculation(x, f)))**2

    Simple_plot_1D(x, Cp)

    return Cp



def Cp_Dist_1D_by_Newton_Modified(x, f, M, gamma=1.4, Mach=10):

    Plot_geometry_1D(x, f(x), title=r'2D Geometry', axis_labels=[r'$x$ [-]', r'$z(x)$ [-]'])

    Cp_max = 2/(g*M**2) * ( (((g+1)**2 * M**2)/(4*g*M**2-2*(g-1)))**(g/(g-1)) * (1-g+2*g*M**2)/(g+1) - 1 ) # From page 62 Andeson Hypersonic...

    Cp_modified = Cp_max*( np.sin(Betta_calculation(x, f)) )**2

    Simple_plot_1D(x, Cp_modified)

    return Cp_modified
########################



# %% Inputs

M_inf = 8; g = 1.4
x = np.linspace(-1, 1, 70)
x_1 = np.linspace(0, 1, 50)

# Cp = Cp_Dist_1D_by_Newton(x, Horizontal_Parabola, M_inf)

# Cp = Cp_Dist_1D_by_Newton(x_1, Quarter_of_a_Circle, M_inf)

Cpmod = Cp_Dist_1D_by_Newton_Modified(x_1, Flechita, M_inf)

