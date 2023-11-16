#!/usr/bin/env python

import numpy as np
import pandas as pd
import warmth
import io
import random

#
# Load model
#

model_str = """
top          topage  k_cond    rhp           phi       decay     solidus      liquidus   
 360.000000     0.0  2.564000  1.498000e-06  0.584000  0.622000  2708.000000  2708.000000
 523.190000     5.0  2.277900  1.585208e-06  0.606586  0.662183  2707.024001  2707.024001
 845.880000    23.0  1.856198  1.608295e-06  0.637742  0.700175  2707.734916  2707.734916
1113.070000    34.0  0.769366  9.523252e-07  0.328388  0.389374  1266.639495  1266.639495
1364.940000    56.0  2.685191  1.428223e-06  0.568786  0.594719  2709.049274  2709.049274
1739.870000    66.0  1.959962  1.634753e-06  0.617657  0.637082  2709.410640  2709.410640
1991.920000   145.0  2.389656  1.524629e-06  0.592339  0.623368  2708.640057  2708.640057
2142.170000   164.0  2.139735  1.593953e-06  0.604479  0.662501  2638.327714  2638.327714
2572.900000   174.0  3.129789  1.008131e-06  0.493581  0.449826  2716.135856  2716.135856
3094.860000   201.0  3.481859  8.653857e-07  0.454647  0.385675  2708.997194  2708.997194
5416.485239   250.0  2.564000  1.498000e-06  0.584000  0.622000  2708.000000  2708.000000
"""
model_df = pd.read_fwf(io.StringIO(model_str.strip()))

#
# Set parameters
#   hc          initial crust thickness (eg, 30000 meters)
#   qbase       heat flow at base of the crust/moho, (eg, 0.030 W/m^2)
#   time_start  start simulation (eg, 0 250 Ma)
#   time_end    end of simulation (eg, 0 Ma)
#   rift_start  time of rift starting (eg, 250 Ma)
#   rift_end    time of rift ending (eg, 240 Ma)
#

params = {
    'hc':           random.uniform(25000, 45000),  
    'qbase':        random.uniform(0.020, 0.060), 
    'time_start':   250,     
    'time_end':       0,     
    'rift_start':   250,     
    'rift_end':     240      
}

#
# Run simulation
#

node = warmth.single_node()
node.sediments_inputs = model_df
model = warmth.Model()
model.parameters.time_end   = params['time_end']
model.parameters.time_start = params['time_start']
node.hc    = params['hc']
node.qbase = params['qbase']
node.rift  = np.array([[params['rift_start'], params['rift_end']]])
model.builder.nodes = [[node]]
model.builder.set_eustatic_sea_level(warmth.data.haq87)
model.simulator.run(parallel=False)

#
# Capture results
#

node = model.builder.nodes[0][0]
depths = node.result.temperature(0)['depth'][1:]
temps = node.result.temperature(0)['values'][1:]
depths_temps = [(d, t) for d, t in zip(depths, temps)][:50]

#
# Output results
#

with open("mywell_temp.out", "w", encoding="utf-8") as f:
    for d, t in depths_temps:
        print(t, file=f)
