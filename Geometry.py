
from Shapes import Shape_definition
from Normals import Revoluting_Normals
from Body_Representation import py2stl, Mesh_gen

class Revolution_Surface(object):

    def __init__(self, x, z, nx, nphi, dx, dphi, N, L, surf_name = "Unknown"):
        '''
        Arguments
        ----------
        x : 1-D ndarray
            DESCRIPTION.
        z : 1-D ndarray
            DESCRIPTION.
        N : 3-D ndarray
            DESCRIPTION.
        dx : Integer
            DESCRIPTION.
        dphi : Integer
            DESCRIPTION.
        surf_name : String
            Name of the surface. The default is "Unknown".
        Coords: List
            List with the matrix coordinates of the surface.
        '''
        self.x = x
        self.z = z
        self.nx = nx
        self.nphi = nphi
        self.N = N
        self.dx = dx
        self.dphi = dphi
        self.surf_name = surf_name
        self.coords = [ L[0], L[1], L[2] ] # List with X, Y, Z mesh for 3-D representation


def Surface_selection(name, nx, nphi):

    z, x = Shape_definition(name, nx)
    Normals_tensor, dx, dphi = Revoluting_Normals(nx, nphi, z)

    [X, Y, Z] = Mesh_gen(x, z, nx, nphi)

    Surface = Revolution_Surface(x, z, nx, nphi, dx, dphi, Normals_tensor, [X, Y, Z], name)

    py2stl(z, name)

    return Surface

# z = shape_definition(name, nx)
















