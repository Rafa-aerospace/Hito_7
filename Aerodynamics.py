"""
Created on Mon Nov 14 17:03:33 2022
@authors: Rafael Rivero de Nicolás, Guillermo García del Río, Inés Arauzo Andrés
"""
from numpy import deg2rad, array, size, dot, zeros, min, max, logical_and, sqrt, hstack, sin, cos, arccos, arctan, linspace, geomspace, reshape, pi, floor
from numpy.linalg import norm
from time import process_time
from Math_Utilities import Theta_calculation, TWB_matrix
from OurDecorators import Simple_plot_1D, Plot_geometry_1D
from Normals_3D import Rotated_Normals_Loops, Rotated_Normals_Giant
from numba import njit
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

def Cp_Dist_1D_by_Newton(x, f, M):
    """Calculates the pressure coefficient over a 1D shape (defined by the evaluation of the function f over the parrtition x) at Mach M

    Args:
        x (array): linspace in the coordinate x
        f (array): result of evaluating a function in x
        M (float): Mach number

    Returns:
        _type_: _description_
    """

    Plot_geometry_1D(x, f(x))

    Cp = 2*(sin(Theta_calculation(x, f)))**2

    #Simple_plot_1D(x, Cp)

    return Cp

def Cp_Dist_1D_by_Newton_Modified(x, f, M, gamma=1.4, Mach=10):

    Plot_geometry_1D(x, f(x), title=r'2D Geometry', axis_labels=[r'$x$ [-]', r'$z(x)$ [-]'])

    Cp_max = 2/(gamma*M**2) * ( (((gamma+1)**2 * M**2)/(4*gamma*M**2-2*(gamma-1)))**(gamma/(gamma-1)) * (1-gamma+2*gamma*M**2)/(gamma+1) - 1 ) # From page 62 Andeson Hypersonic...

    Cp_modified = Cp_max*( sin(Theta_calculation(x, f)) )**2

    Simple_plot_1D(x, Cp_modified)

    return Cp_modified

def Cp_3D(alpha: float, beta: float, nx:int, nphi: int, f):
    """This function calculates the Cp for each normal of the surface. Remember that the Normals_Tensor already contemplates both the angle of atack and the sideslip.

    Args:
        alpha (float): Angle of attack.       [Deg]
        beta (float): Sideslip.               [Deg]
        nx(int): number of x partitions
        nphi(int): number of phi partitions
        f(function): f(x) which will be revoluted

    Returns:

    Cp_matrix: Matrix (nx,nphi). Each row is a "x" partition and has has "nphi" Cps, this for every "x" partition (nx)
    """
    # Normals_tensor =  Rotated_Normals_Loops(alpha, beta, nx, nphi, f)
    Normals_tensor =  Rotated_Normals_Giant(alpha, beta, nx, nphi, f) # Option 2 for the rotation

    Cp_matrix = zeros([nx,nphi])

    V0 = array([-1,0,0])

    for i in range(nx):
        for j in range(nphi):

            normal = (Normals_tensor[i,j])
            psi = arccos(dot(V0,normal)/norm(V0))

            #print(psi)

            if (pi/2 - psi) < 0:
                theta = 0
            else:
                theta = pi/2 - psi

            # print(pi/2 - psi)

            Cp_matrix[i,j] = 2*(sin(theta))**2

    return Cp_matrix

def Cp_3D_for_coefs(alpha: float, beta: float, nx:int, nphi: int, f):
    """This function calculates the Cp for each normal of the surface. Remember that the Normals_Tensor already contemplates both the angle of atack and the sideslip.

    Args:
        alpha (float): Angle of attack.       [Deg]
        beta (float): Sideslip.               [Deg]
        nx(int): number of x partitions
        nphi(int): number of phi partitions
        f(function): f(x) which will be revoluted

    Returns:

    Cp_matrix: Matrix (nx,nphi). Each row is a "x" partition and has has "nphi" Cps, this for every "x" partition (nx)
    """
    Normals_tensor =  Rotated_Normals_Loops(alpha, beta, nx, nphi, f)
    # Normals_tensor =  Rotated_Normals_Giant(alpha, beta, nx, nphi, f) # Option 2 for the rotation

    Cp_matrix = zeros([nx,nphi])

    V0 = array([-1,0,0])

    for i in range(nx):
        for j in range(nphi):

            normal = (Normals_tensor[i,j])
            psi = arccos(dot(V0,normal)/norm(V0))

            #print(psi)

            if (pi/2 - psi) < 0:
                theta = 0
            else:
                theta = pi/2 - psi

            # print(pi/2 - psi)

            Cp_matrix[i,j] = 2*(sin(theta))**2

    return Cp_matrix, Normals_tensor




def Cp_3D_Modified(alpha: float, beta: float, nx:int, nphi: int, f, M: float, gamma: float):
    """This function calculates the modified Cp for each normal of the surface. Remember that the Normals_Tensor already contemplates both the angle of atack and the sideslip.

    Args:
        alpha (float): Angle of attack.       [Deg]
        beta (float): Sideslip.               [Deg]
        nx(int): number of x partitions
        nphi(int): number of phi partitions
        f(function): f(x) which will be revoluted
        M (float): : Mach number
        gamma (float): fluid's gamma

    Output:
        Cp_matrix: Matrix (nx,nphi). Each row is a "x" partition and has has "nphi" Cps, this for every "x" partition (nx)

    """
    Normals_tensor =  Rotated_Normals_Loops(alpha, beta, nx, nphi, f)

    Cp_matrix = zeros([nx,nphi])

    Cp_max = 2/(gamma*M**2) * ( (((gamma+1)**2 * M**2)/(4*gamma*M**2-2*(gamma-1)))**(gamma/(gamma-1)) * (1-gamma+2*gamma*M**2)/(gamma+1) - 1 )

    V0 = array([1,0,0])

    for i in range(nx):
        for j in range(nphi):

            normal = (Normals_tensor[i,j])
            psi = arccos(dot(V0,normal)/norm(V0))

            if (pi/2 - psi) < 0:
                theta = 0
            else:
                theta = pi/2 - psi

            Cp_matrix[i,j] = Cp_max*(sin(theta))**2

    return  Cp_matrix


def aero_coeffs(Cp, N, x, z, dz, n_x, n_phi, dx, dphi, alpha, beta):

    N_panels = n_x*n_phi

    N_matrix = reshape(N, [N_panels, 3])
    Cp_v = reshape(Cp, [N_panels])

    d_Surface = zeros(N_panels)

    for n in range(N_panels):

        i = int(floor( n/ (n_phi) ))
        # print(i)

        hipotenusa = sqrt( ( x[i+1] - x[i] )**2 + ( z[i+1] - z[i] )**2 )

        d_Surface[n] = dphi * ( z[i+1] + z[i] )/2 * hipotenusa

    S_rev = sum(d_Surface)

    Coefs_XYZ = ( (-1) * dot( (Cp_v* d_Surface), N_matrix)  ) / S_rev

    CD_CY_CL = dot( TWB_matrix(alpha, beta), Coefs_XYZ )
    CD_CY_CL[0] = -CD_CY_CL[0]; CD_CY_CL[1] = CD_CY_CL[1]; CD_CY_CL[2] = -CD_CY_CL[2];

    return CD_CY_CL

#%% Tests

# from numpy import zeros, array
# T = zeros([2,3,3])         ## El 2 son las particiones en X , el primer 3 las particiones en Phi y el último 3
#                            ## es el número de componentes del vector normal (i j k)

# A = array([[1, 2, 3],
#            [4, 5, 6],
#            [7, 8, 9]])

# C = array([[4,0,0],
#            [5,0,0],
#            [6,0,0]])

# T[0] = A
# T[1] = C

# for i in range(size(T,0)):
#     for j in range(size(T,1)):
#         print(T[i,j])
