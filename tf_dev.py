import yt
import numpy as np
import linecache

import time
import sys
from IPython.core.display import Image
from yt.visualization.volume_rendering.api import Scene, VolumeSource
from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
# usage: python script_name data_file_path frame_num

print str(sys.argv)
# print sys.argv[1]

t0 = time.time()


use_log = True

# experiment with tf functions
mi = YTQuantity(7.11301e-09, 'g/cm**3')
ma = YTQuantity(3.33129, 'g/cm**3')

if use_log:
    mi,ma = np.log10(mi), np.log10(ma)

# Instantiate the ColorTransferfunction.
print mi, ma
tf = yt.ColorTransferFunction((mi, ma))
'''
tfh = TransferFunctionHelper(ds)
tfh.set_field('density')
tfh.set_log(True)
tfh.set_bounds([mi, ma])
tfh.build_transfer_function()
'''
# tf.add_layers(30, colormap = 'RdBu')
tf.add_gaussian(-8, 0.001, [0.0, 0.0, 1.0, 0.05])
tf.add_gaussian(-7, 0.001, [0.0, 1.0, 0.0, 0.05])
tf.add_gaussian(-4, 0.001, [1.0, 0.0, 0.0, 0.05])
tf.plot("tf.png")

t1 = time.time()

total = t1-t0

print total
