"""
Created on Mon Nov 14 17:03:33 2022
@authors: Rafael Rivero de Nicolás, Guillermo García del Río, Inés Arauzo Andrés
"""

from numpy import array, size, dot, zeros, linspace, reshape
from Math_Utilities import rotation_rev, rotation_WB, TBW_matrix, Normal

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

#@njit(parallel=True)
def Revoluting_Normals(nx:int, nphi: int, f):
    """Given the normal vector field of a 1D shape in an array, returns the  normal vector field of the 3D revoluted shape

    Args:
        nx(int): number of x partitions
        nphi(int): number of phi partitions
        f(function): f(x) which will be revoluted

    Returns:
        Normals_revoluted(array): tensor with all the normals for every phi partition in each x partition.
            nx     = size(Normals_revoluted,0)     # Number of partitions in X = Number of matrixes in the tensor
            nphi   = size(Normals_revoluted,1)     # Number of partitions in the revolution angle = Number of files in the matrix
            ncomps = size(Normals_revoluted,2)     # Number of components (3 = i j k)
    """

    phi_list = linspace(0,360, nphi)
    # print('changethis')
    x = linspace(0, 1, nx + 1)
    Normals_array = Normal(x, f)
    Normals_revoluted = zeros([nx, nphi, 3])
    

    for i, n in enumerate(Normals_array):

        for j, phi in enumerate(phi_list):

            Normals_revoluted[i, j, :] = rotation_rev(phi,n)


    return Normals_revoluted

# print(Revoluting_Normals(5, 4, conito))

# Giro de las normales para un alpha o beta
def Rotated_Normals_Loops(alpha: float, beta: float, nx:int, nphi: int, f):
    """ Rotation of all the normals a certain angle of attack and a certain sideslip.

    Args:
        alpha (float): Angle of attack.       [Deg]
        beta (float): Sideslip.               [Deg]
        nx(int): number of x partitions
        nphi(int): number of phi partitions
        f(function): f(x) which will be revoluted

    Returns:
        Normals_tensor_rotated: Tensor in which each matrix represents all the normals along a section of the revolution surface
                                taking into acount the angle of attack and the sideslip.
    """
    Normals_tensor = Revoluting_Normals(nx, nphi, f)

    ncomps = size(Normals_tensor,2)     # Number of components (3 = i j k)

    Normals_tensor_rotated = zeros([nx, nphi, ncomps])

    for j in range(nx):
        for i in range(nphi):
            #Normals_tensor_rotated[j,i,:] = rotation_BW(alpha,beta,Normals_tensor[j,i,:])
            Normals_tensor_rotated[j,i,:] = rotation_WB(alpha,beta,Normals_tensor[j,i,:])

    return Normals_tensor_rotated

def Rotated_Normals_Giant(alpha: float, beta: float, nx:int, nphi: int, f):
    """ Rotation of all the normals a certain angle of attack and a certain sideslip.

    Args:
        alpha (float): Angle of attack.       [Deg]
        beta (float): Sideslip.               [Deg]
        nx(int): number of x partitions
        nphi(int): number of phi partitions
        f(function): f(x) which will be revoluted

    Returns:
        Normals_tensor_rotated: Tensor in which each matrix represents all the normals along a section of the revolution surface
                                taking into acount the angle of attack and the sideslip.
    """
    Normals_tensor = Revoluting_Normals(nx, nphi, f)

    nx     = size(Normals_tensor,0)     # Number of partitions in X = Number of matrixes in the tensor
    nphi   = size(Normals_tensor,1)     # Number of partitions in the revolution angle = Number of files in the matrix
    ncomps = size(Normals_tensor,2)     # Number of components (3 = i j k)

    Output = zeros( [nx, nphi, ncomps] )

    Tbw = TBW_matrix(alpha, beta)

    # Giant matrix building
    Tbw_giant = tbw_giant_loop (zeros([3*nphi, 3*nphi]), Tbw, nphi)

    for i in range( size(Normals_tensor,0) ):

        p_O = reshape(Output[i,:,:], [3*nphi])         # Pointer to the Output tensor. It will be a column vector
        p_T = reshape(Normals_tensor[i,:,:], [3*nphi])  # Pointer to the Input tensor. It will be a column  vector
        p_O[:] = p_O[:] + dot(Tbw_giant, p_T)          # Vector corresponding to all normals rotated for a x value

    return Output

@njit(parallel=True)
def tbw_giant_loop (A:array, Tbw:array, nphi:int):
    """_summary_

    Args:
        A (array): init tensor full of zeros
        Tbw (array): body-wind rotation matrix
        nphi(int): number of phi partitions

    Returns:
        A: Rotation Matrix for all the normals
    """

    for i in range( nphi ):

        A [3*i:3*i+3, 3*i:3*i+3] = Tbw
    return A


# %%

# def fx(x):
#     return x


# nx = 5; nphi = 4; f = fx

# """Given the normal vector field of a 1D shape in an array, returns the  normal vector field of the 3D revoluted shape

# Args:
#     nx(int): number of x partitions
#     nphi(int): number of phi partitions
#     f(function): f(x) which will be revoluted

# Returns:
#     Normals_revoluted(array): tensor with all the normals for every phi partition in each x partition.
#         nx     = size(Normals_revoluted,0)     # Number of partitions in X = Number of matrixes in the tensor
#         nphi   = size(Normals_revoluted,1)     # Number of partitions in the revolution angle = Number of files in the matrix
#         ncomps = size(Normals_revoluted,2)     # Number of components (3 = i j k)
# """

# phi_list = linspace(0,360, nphi)
# print('changethis')
# x = linspace(0, 1, nx + 1)
# Normals_array = Normal(x, f)
# Normals_revoluted = zeros([nx, nphi, 3])

# for i, n in enumerate(Normals_array):

#     for j, phi in enumerate(phi_list):

#         Normals_revoluted[i, j, :] = rotation_rev(phi,n)
#         # print(i, j)




















