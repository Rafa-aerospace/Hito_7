from numpy import array, zeros, linspace
from Rotations import rotation_rev
import surf2stl
import os

def coords_tensor_generator(nx:int  ,nphi: int, coords_matrix:array):

    phi_list = linspace(0, 360, nphi + 1)

    coords_revoluted = zeros([nx+1, nphi+1, 3])

    for i in range(nx+1):

        for j in range(nphi+1):

            coords_revoluted[i, j, :] = rotation_rev( phi_list[j], coords_matrix[i])

    return coords_revoluted

def Mesh_gen(x: array, z: array, nx: int, nphi: int):

    y = zeros(len(z))
    
    coords = zeros( [nx+1, 3] )
    for i in range(nx+1):
        coords[i,0] = x[i]
        coords[i,1] = y[i]
        coords[i,2] = z[i]

    coords_revoluted = coords_tensor_generator(nx, nphi, coords)

    X = zeros([nphi+1,nx+1])
    Y = zeros([nphi+1,nx+1])
    Z = zeros([nphi+1,nx+1])

    for k in range(nx+1):
        X[:,k] = coords_revoluted[k,:,0]
        Y[:,k] = coords_revoluted[k,:,1]
        Z[:,k] = coords_revoluted[k,:,2]
    
    return [X, Y, Z]

def py2stl(z: array, name: str):

    folder = 'STL'
    filename = f'{name}.stl'
    filepath = f'{folder}\{filename}'

    if not os.path.exists(filepath):

        nx = 100
        nphi = 500
        x = linspace(0, 1, nx + 1)
        znum = int( len(z) / (nx + 1) )

        z_rep = [z[i*znum] for i in range(nx+1)] # z for representation
        z_rep[-1] = z[-1]

        [X, Y, Z] = Mesh_gen(x, z, nx, nphi)

        surf2stl.write(filepath, X, Y, Z)

