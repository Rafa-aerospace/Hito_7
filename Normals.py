from numpy import array, size, dot, zeros, linspace, reshape
from numpy.linalg import norm
from Rotations import rotation_rev, rotation_WB, TBW_matrix
from numba import njit


def Normal(x, f):
    """Given a function f evaluated on an x partition, returns the normal vector of every partition

    Args:
        x (array): linspace in the coordinate x
        f (array): result of evaluating a function in x

    Returns:
        Normal (array): Array of normal vectors
    """

    Normals = zeros([len(x)-1, 3])
    for i in range (len(x)-1):
        dx = x[i+1] - x[i]
        dy = 0
        dz = f[i+1] - f[i]
        Normals[i, :] = (-dz, dy, dx)/norm((-dz, dy, dx))

    return Normals

def Revoluting_Normals(nx:int, nphi: int, f):
    """Given the normal vector field of a 1D shape in an array, returns the  normal vector field of the 3D revoluted shape

    Args:
        nx(int): number of x partitions
        nphi(int): number of phi partitions
        f(1D-ndarray): f(x) which will be revoluted

    Returns:
        Normals_revoluted(array): tensor with all the normals for every phi partition in each x partition.
            nx     = size(Normals_revoluted,0)     # Number of partitions in X = Number of matrixes in the tensor
            nphi   = size(Normals_revoluted,1)     # Number of partitions in the revolution angle = Number of files in the matrix
            ncomps = size(Normals_revoluted,2)     # Number of components (3 = i j k)
    """

    phi_list = linspace(0, 360, nphi+1)
    dphi = phi_list[1] - phi_list[0]
    x = linspace(0, 1, nx + 1)
    dx = x[1]-x[0]
    Normals_array = Normal(x, f)

    for i in range(nx):
        Normals_array[i,:] = rotation_rev(dphi/2, Normals_array[i,:])

    Normals_revoluted = zeros([nx, nphi, 3])
    
    for i in range(len(Normals_array)):
        Normals_revoluted[i, 0, :] = Normals_array[i, :]

        for j, phi in enumerate(phi_list[1:]):

            Normals_revoluted[i, j, :] = rotation_rev(phi, Normals_revoluted[i, 0, :])

    return Normals_revoluted, dx, dphi

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