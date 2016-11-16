import yt
import numpy as np
import linecache
from math import ceil

import time
import sys
from IPython.core.display import Image
from yt.visualization.volume_rendering.api import Scene, VolumeSource
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
# usage: python script_name data_file_path frame_num

# the data file should be a gaussian cube file format
# http://paulbourke.net/dataformats/cube/

# IMPORTANT: PLEASE DO READ ALL THE COMMENTS IF YOU ARE NOT FAMILIAR WITH YT OR THIS PROJECT
# THE COMMENTS ARE INTENDED TO HELP RATHER THAN TO CONFUSE

# original authors contact: zwen6@illinois.edu tcai4@illinois.edu

# a purely mathematical helper function
# checks if the given point is within a 3-D sphere
# c is the center of the sphere, which should be a tuple of 3 elements
# r is the radius which should be a float
# p is the point to be tested which should be a tuple of 3 elements
# returns a boolean
def within_sphere(c, r, p):
	dist_sq = (p[0] - c[0]) * (p[0] - c[0]) + (p[1] - c[1]) * (p[1] - c[1]) + (p[2] - c[2]) * (p[2] - c[2])
	r_sq = r * r
	if dist_sq <= r_sq: # yes
		return True
	else:
		return False

# for each of the atom centers
# iterates each cell in the cube surrounding it
# the cube: (2 * r) by (2 * r) by (2 * r)
# and sees if the cell is within the atom sphere
# if yes assign HUGE density
# if not assign original value
# this function forms crude 'spheres' using arbitrary density values
# there are HARD-CODED values here subject to change
# THIS FUNCTION HAS SPACE FOR IMPROVEMENT
# please note here the length unit should be different for three axis
# i used the unit for z-axis here for all three
# since yt's unit system only support one uniform unit for the 3-D space
def add_spheres(dataset, atoms, r):
	for a in atoms:
		for i in range(max(int((a[0] - r) / .20800511), 0), min(int((a[0] + r) / .20800511), 180)):
			for j in range(max(int((a[1] - r) / .20800511), 0), min(int((a[1] + r) / .20800511), 180)):
				for k in range(max(int((a[2] - r) / .20800511), 0), min(int((a[2] + r) / .20800511), 720)):

					ii = i * .20800511
					jj = j * .20800511
					kk = k * .20800511

					dist_sq = (ii - a[0]) * (ii - a[0]) + (jj - a[1]) * (jj - a[1]) + (kk - a[2]) * (kk - a[2])
					r_sq = r * r
					if dist_sq <= r_sq: # yes
						dataset[i][j][k] = 1

print str(sys.argv)
# print sys.argv[1]

# a timer to check the runtime of the script
t0 = time.time()

# a domain of data
# this is a spatial representation of a 3-D uniform grid
# within each entry (cell) is the electronic density at this location
dom = np.ndarray(shape=(180, 180, 720), dtype=float, order='F')

# counter variables
# keeps track of the spatial location
x = 0
y = 0
z = 0
count = 0

# a list of tuples that records the locations of centers of atoms
atoms = []

# reading the atom locations from the data file
for i in range(7, 207): # excluding the hydrogen atom on line 207
	curr_line = linecache.getline(sys.argv[1], i, None)
	ll = curr_line.split(" ")
	atoms.append((float(ll[2]), float(ll[3]), float(ll[4]))) # record the centers of the atoms

# parsing the file
# filling in the spatial data to the domain
for i in range(208, 3888208):
	curr_line = linecache.getline(sys.argv[1], i, None)
	ll = curr_line.split(" ")

	dom[x][y][z] = float(ll[0])

	dom[x][y][z + 1] = float(ll[1])

	dom[x][y][z + 2] = float(ll[2])

	dom[x][y][z + 3] = float(ll[3])

	dom[x][y][z + 4] = float(ll[4])

	dom[x][y][z + 5] = float(ll[5])

	# update the coordiniates by the gaussian file format
	z = z + 6
	if z == 720:
		z = 0
		y = y + 1
		if y == 180:
			y = 0
			x = x + 1
	# print "at line " + str(i)


# add the density "spheres" to cells
# but also check if the cell is within any atom
# if yes assign fixed electronic density to the cell
# this is a work-around for adding 3-D spheres to the rendering
# add_spheres(dom, atoms, 1)
	
# extract the hydrogen atom location from the file
curr_line = linecache.getline(sys.argv[1], 207, None)
ll = curr_line.split(" ")
hx = float(ll[2])
hy = float(ll[3])
hz = float(ll[4])

# yt-specific data structures and parameters setup
field = 'density'
use_log = True
data = dict(density = (dom, "g/cm**3"))
bbox = np.array([[0, 180], [0, 180], [0, 720]])

# initialize the yt-specific data structure (ds)
# for more info read: 
# http://yt-project.org/docs/dev/reference/api/generated/yt.frontends.stream.data_structures.load_uniform_grid.html?highlight=load_uniform_grid

# please note here the length unit should be different for three axis
# i used the unit for z-axis here for all three
# since yt's unit system only support one uniform unit for the 3-D space
ds = yt.load_uniform_grid(data, dom.shape, length_unit=0.20800511, bbox=bbox, nprocs=64)

# initialize the yt scene
sc = Scene()
vol = VolumeSource(ds, field=field)

# camera position
#hc = ds.arr([hx, hy, hz], 'cm') # hydrogen location
hc = ds.arr([19.13987520, 19.13987520, 82.13300000], 'cm') # hydrogen location

# Find the bounds in log space of for your field
dd = ds.all_data()
mi, ma = dd.quantities.extrema(field)
if use_log:
    mi,ma = np.log10(mi), np.log10(ma)

# instantiating the ColorTransferfunction
tf = yt.ColorTransferFunction((mi, ma))
'''
tfh = TransferFunctionHelper(ds)
tfh.set_field('density')
tfh.set_log(True)
tfh.set_bounds([mi, ma])
tfh.build_transfer_function()
'''
# tf.add_layers(30, colormap = 'RdBu')
# tf.add_gaussian(-2, 0.0, [0.0, 0.0, 0.9, 0.9])
tf.add_gaussian(-0, 1e-2, [0.0, 0.0, 0.8, 0.8])
tf.add_gaussian(-2, 1e-3, [0.1, 0.05, 0.3, 0.6])
tf.add_gaussian(-4, 1e-4, [0.3, 0.05, 0.05, 0.6])
tf.add_gaussian(-6, 1e-4, [0.1, 0.2, 0.05, 0.6])
tf.add_gaussian(-8, 1e-4, [0.0, 0.1, 0.1, 0.6])
# tf.add_gaussian(1e-7, 1e-8, [0.0, 1.0, 0.0, 1.0])
# tf.add_gaussian(1e-6, 1
# tfh.plot("vol.png", "density")
'''
tfh.tf.map_to_colormap(-8., -3., colormap = 'Blues', scale=1.0)
tfh.tf.map_to_colormap(-3., -1., colormap = 'Greens', scale=2.0)
tfh.tf.map_to_colormap(-1., 1., colormap = 'Reds', scale=2.0)
'''
# set the tf for the volume rendering
vol.set_transfer_function(tf)

# plot the transfer function for debugging and demo
tf.plot("func_org.png")


# setting up the camera parameters: center, looking direction, width, resolution
# c = (ds.domain_right_edge + ds.domain_left_edge)/2.0
# c = ds.arr([-100, 90, 180], 'code_length')
# L = np.array([1.0, 1.0, 1.0])
L = [0., 0., 1.]
W = ds.quan(0.5, 'unitary')
N = 512
north = [0., 1., 0.]

# Create a camera object
#cam = sc.add_camera(ds, lens_type='stereo-spherical')
cam = sc.add_camera(ds, lens_type='spherical')

cam.resolution = [4096, 2048]
#cam.resolution = [2048, 1024]

# cam.resolution = [1600, 1600]
# cam.positon = c
cam.position = hc
cam.switch_orientation(normal_vector=L, north_vector=north)
# cam.lens.disparity = ds.domain_width[0] * 1.e-3
sc.add_source(vol)
sc.render()

# arbitrary output file name
# subject to change
sc.save('spherical_rendering_frame' + str(sys.argv[2]) + '.png', sigma_clip=0)

# the following lines are legacy from earlier development

# snapshots = list(enumerate(cam.zoomin(3.0, 1, clip_ratio=0.0)))
# snapshots[0][1].write_png('test_rendering_zoomed.png')
#frame = 0
#for i, snapshot in enumerate(cam.zoomin(10.0, 5)):
#    snapshot.write_png('camera_movement_%04i.png' % frame)
#    frame += 1

# Now let's add some isocontours, and take a snapshot, saving the image
# to a file.

# im = cam.snapshot('test_rendering_frame' + str(sys.argv[2]) + '.png')

# To add the domain box to the image:
# inim = cam.draw_domain(im)
# nim.write_png('test_rendering_with_domain.png')

# To add the grid outlines to the image:
# nim = cam.draw_grids(im)
# nim.write_png('test_rendering_with_grids.png')


# timer
t1 = time.time()
total = t1-t0
print total

