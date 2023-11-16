import numpy as np
import pandas as pd
import warmth
import random

#
# Load model
#

model_df = pd.read_fwf("mymodel.txt")

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

print(f"{params['hc']=}")
print(f"{params['qbase']=}")

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

print("Depth Temp")
for d, t in depths_temps[:5]:
    print(f"{int(d):>5} {int(t):>4}")
print("...")
for d, t in depths_temps[-5:]:
    print(f"{int(d):>5} {int(t):>4}")
