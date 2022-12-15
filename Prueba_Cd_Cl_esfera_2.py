# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 16:08:17 2022

@author:  Rafael Rivero de Nicolás, Guillermo García del Río, Inés Arauzo Andrés
"""

from numpy import deg2rad, array, size, dot, zeros, min, max, logical_and, sqrt, hstack, sin, cos, arctan, linspace, geomspace, reshape, pi, floor, ceil, sum
from numpy.random import rand

import Aerodynamics as aero
import Math_Utilities as mth
import Normals_3D as nor
import Geometry as gm
import OurDecorators as pl


n_x = 4; dx = 1/n_x; n_phi = 10; dphi = 2*pi/(n_phi+1); # radianes

x = linspace(0, 1, n_x+1)
alpha = 0; beta = 0
# z = gm.Quarter_of_a_Circle(x)
# dz = mth.Analytical_Derivative( x, gm.Quarter_of_a_Circle )

z = gm.conito(x)
dz = mth.Analytical_Derivative( x, gm.conito )

pl.Plot_geometry_1D(x, z)


N_panels = (n_x)*(n_phi) # Número de paneles totales en los que se ha mallado la superficie

# %% Función de los coeficientes

# for i in range(n_x):

#     S_rev = S_rev + 2*pi * z[i] * sqrt( 1 + (dz[i])**2 ) * dx

# %%
# Cps = rand(n_x, n_phi)
# Normals = rand(n_x, n_phi, 3)

Normals = reshape(linspace(1, N_panels*3, N_panels*3), [n_x, n_phi, 3])
Cps = reshape(linspace(1,N_panels,N_panels), [n_x, n_phi])

N_m = reshape(Normals, [N_panels, 3])
Cps_v = reshape(Cps, [N_panels])

# %%
d_Surface = zeros(N_panels)

for n in range(N_panels):

    i = int(floor( n/ (n_phi) ))
    # print(i)

    hipotenusa = sqrt( ( x[i+1] - x[i] )**2 + ( z[i+1] - z[i] )**2 )

    d_Surface[n] = dphi * ( z[i+1] + z[i] )/2 * hipotenusa

S_rev = sum(d_Surface)

Coefs_XYZ = ( (-1) * dot( (Cps_v* d_Surface), N_m)  ) / S_rev


AA = d_Surface * Cps_v

BB = dot(AA, N_m)

CD_CY_CL = dot( mth.TWB_matrix(alpha, beta), Coefs_XYZ )
CD_CY_CL[0] = -CD_CY_CL[0]; CD_CY_CL[1] = CD_CY_CL[1]; CD_CY_CL[2] = -CD_CY_CL[2];



# S_an = 4 * pi * 1**2 / 2 # Solo para la esfera
# S_an_conito = pi*1*sqrt(2)

















