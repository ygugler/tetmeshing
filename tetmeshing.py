"""
Test
"""

import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# meshing
import pygalmesh
import pygmsh  # is this really necessary?
import gmsh  # background mesh specifying element sizes! --> very interesting / open, merge, write, addvolume , addsurfaceloop
# https://github.com/nschloe/pygmsh/issues/279
# import tetgen  #  pip install tetgen --> look up how that works
import meshio

# marching cubes
from skimage import measure
import mcubes

def sitk2np(img_sitk, pixel_type=np.float32):
    """

    :param img_sitk:
    :param pixel_type:
    :return:
    """
    img_np = sitk.GetArrayFromImage(img_sitk)  # transform to numpy array
    img_np = img_np.transpose(2, 1, 0)  # switch ordering of indexes: SITK is in x, y, z / np switches to z, y, x
    img_np = img_np.astype(pixel_type)
    return img_np


img_file = "/home/ygugler/Desktop/FEAMUR2020/Gerico_Images/PA08A02A-vol-output.mhd"
img_sitk = sitk.ReadImage(img_file)

img_np = sitk2np(img_sitk)
img_np_bin = np.copy(img_np)
img_np_bin[img_np_bin != 0] = 1

verts, faces, normals, values = measure.marching_cubes(img_np_bin, step_size=3)

fig = plt.figure(figsize=(100, 100))
ax = fig.add_subplot(111, projection='3d')


mesh = Poly3DCollection(verts[faces])
mesh.set_edgecolor('k')
ax.add_collection3d(mesh)


ax.set_xlabel("x-axis: a = 6 per ellipsoid")
ax.set_ylabel("y-axis: b = 10")
ax.set_zlabel("z-axis: c = 16")

ax.set_xlim(0, 160)  # a = 6 (times two for 2nd ellipsoid)
ax.set_ylim(0, 160)  # b = 10
ax.set_zlim(0, 160)  # c = 16

plt.tight_layout()
plt.show()

points = np.array(verts)
cells = [("triangle", np.array(faces))]
meshio.write_points_cells('out.vtu', points, cells)
meshio.write_points_cells('out.inp', points, cells)
meshio.write_points_cells('out.stl', points, cells)

# gmsh
gmsh.initialize()
gmsh.merge("/home/ygugler/Code/tetmeshing/out.stl")
n = gmsh.model.getDimension()
s = gmsh.model.getEntities(n)
l = gmsh.model.geo.addSurfaceLoop([s[i][1] for i in range(len(s))])
gmsh.model.geo.add_surface_loop([faces[i][1] for i in range(np.size(faces[0]))])
gmsh.model.geo.addVolume([l])
gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(3)
gmsh.write("/home/ygugler/Code/tetmeshing/out.msh")
gmsh.finalize()


mesh_new = meshio.read("/home/ygugler/Code/tetmeshing/out.msh")
mesh_new_abq = meshio.write("/home/ygugler/Code/tetmeshing/out2.inp", mesh_new)


# pygalmesh
# pymesh = pygalmesh.generate_volume_mesh_from_surface_mesh('out.vtu')
"""
pymesh = pygalmesh.generate_volume_mesh_from_surface_mesh('out.vtu', min_facet_angle=50.0,
    max_radius_surface_delaunay_ball=3.0,
    max_facet_distance=0.5,
    max_circumradius_edge_ratio=4.0,
    max_cell_circumradius=5.0)
"""

pymesh.write("femur.inp")