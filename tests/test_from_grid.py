
import os
import multiprocessing
from pathlib import Path
import sys
import numpy as np

import temperer


maps_dir = Path("../data/mapA")

model = temperer.Model()

inputs = model.builder.input_horizons_template

inputs.loc[0]=[0,"0.gri",None,"Onlap"]
inputs.loc[1]=[66,"66.gri",None,"Onlap"]
inputs.loc[2]=[100,"100.gri",None,"Onlap"]
inputs.loc[3]=[163,"163.gri",None,"Erosive"]
inputs.loc[4]=[168,"168.gri",None,"Erosive"]
inputs.loc[5]=[170,"170.gri",None,"Onlap"]
inputs.loc[6]=[182,"182.gri",None,"Erosive"]
model.builder.input_horizons=inputs


inc = 2000
model.builder.define_geometry(maps_dir/"0.gri",xinc=inc,yinc=inc,fformat="irap_binary")

model.builder.extract_nodes(4,maps_dir)

from temperer.data import haq87
model.builder.set_eustatic_sea_level(haq87)

for i in model.builder.iter_node():
    i.rift=[[182,175]]


# simulate every 10 nodes
for index in model.builder.grid.indexing_arr:
    if (index[0] % 10 > 0):
        pass
    else:
        if isinstance(model.builder.nodes[index[0]][index[1]],bool) is False:
            model.builder.nodes[index[0]][index[1]] = False
    if (index[1] % 10 > 0):
        pass
    else:
        if isinstance(model.builder.nodes[index[0]][index[1]],bool) is False:
            model.builder.nodes[index[0]][index[1]] = False

indexer_full_sim = [i.indexer for i in model.builder.iter_node() if i._full_simulation is True]

for ni in range(len(model.builder.nodes)):
    for nj in range(len(model.builder.nodes[ni])):
        if model.builder.nodes[ni][nj] is not False:
            nn = model.builder.nodes[ni][nj]
            if nn.Y>40000:
                nn.sediments['phi'] = np.ones([len(nn.sediments['phi']),1]) * 0.5
                # nn.sediments['k_cond'] = np.ones([len(nn.sediments['k_cond']),1]) * 2.2


model.simulator.run(save=False,purge=True)

# %%
for i in model.builder.iter_node():
    if i is not False:
        print(i.result.heatflow(0))

# interpolate and extrapolate the missing nodes
# find nearby nodes from the array indexer_full_sim, which is sorted by x-index
import itertools
for ni in range(len(model.builder.nodes)):
    for nj in range(len(model.builder.nodes[ni])):
        if model.builder.nodes[ni][nj] is False:
            closest_x_up = []
            for j in range(ni,len(model.builder.nodes)):
                matching_x = [ i[0] for i in indexer_full_sim if i[0]==j ]
                closest_x_up = closest_x_up + list(set(matching_x))
                if len(matching_x)>0:
                    break
            closest_x_down = []
            for j in range(ni-1,-1,-1):
                matching_x = [ i[0] for i in indexer_full_sim if i[0]==j ]
                closest_x_down = closest_x_down + list(set(matching_x))
                if len(matching_x)>0:
                    break
            closest_y_up = []
            for j in range(nj,len(model.builder.nodes[ni])):
                matching_y = [ i[1] for i in indexer_full_sim if (i[1]==j and ((i[0] in closest_x_up) or i[0] in closest_x_down)) ]
                closest_y_up = closest_y_up + list(set(matching_y))
                if len(matching_y)>0:
                    break
            closest_y_down = []
            for j in range(nj-1,-1,-1):
                matching_y = [ i[1] for i in indexer_full_sim if (i[1]==j and (i[0] in closest_x_up or i[0] in closest_x_down) ) ]
                closest_y_down = closest_y_down + list(set(matching_y))
                if len(matching_y)>0:
                    break

            interpolationNodes = [  model.builder.nodes[i[0]][i[1]] for i in itertools.product(closest_x_up+closest_x_down, closest_y_up+closest_y_down)  ]
            interpolationNodes = [nn for nn in interpolationNodes if nn is not False]
            node = temperer.build.interpolateNode(interpolationNodes)
            node.X, node.Y = model.builder.grid.location_grid[ni,nj,:]
            model.builder.nodes[ni][nj] = node
        else:
            node = temperer.build.interpolateNode([model.builder.nodes[ni][nj]])  # "interpolate" the node from itself to make sure the same member variables exist at the end
            model.builder.nodes[ni][nj] = node

max_age = model.builder.input_horizons['Age'].to_list()[-1]

#
# assert that the grid of nodes is fully defined;  
# assert that the .crust_ls property has been defined for every node
#
for index in model.builder.grid.indexing_arr:
    nn = model.builder.nodes[index[0]][index[1]]
    assert nn is not False
    assert type(nn)==temperer.build.single_node
    assert nn.crust_ls.shape[0] > max_age-1
    if nn.Y>40000:
        nn.sediments['phi'] = np.ones([len(nn.sediments['phi']),1]) * 0.50
        # nn.sediments['k_cond'] = np.ones([len(nn.sediments['k_cond']),1]) * 2.2

#
# run the 3D simulation
#
from subsheat3D.fixed_mesh_model import UniformNodeGridFixedSizeMeshModel
from subsheat3D.Helpers import NodeGrid

try:
    os.mkdir('out')
except FileExistsError:
    pass
try:
    os.mkdir('mesh')
except FileExistsError:
    pass
try:
    os.mkdir('temp')
except FileExistsError:
    pass

# convert to 1D array of nodes and add padding!

pad = 1
nodes = []
for j in range(-pad, model.builder.grid.num_nodes_y+pad):
    for i in range(-pad, model.builder.grid.num_nodes_x+pad):
        node_i = max(min(i, model.builder.grid.num_nodes_x-1), 0)
        node_j = max(min(j, model.builder.grid.num_nodes_y-1), 0)
        if (node_i != i) or (node_j != j):
            # this is a padded node: create it from the nearest node in the unpadded grid
            node_new = temperer.build.interpolateNode([model.builder.nodes[node_j][node_i]])
            node_new.X = model.builder.grid.origin_x + i * model.builder.grid.step_x 
            node_new.Y = model.builder.grid.origin_y + j * model.builder.grid.step_y
            nodes.append(node_new)
        else:
            nodes.append(model.builder.nodes[node_j][node_i])

nodeGrid = NodeGrid(0, 0, model.builder.grid.num_nodes_x+2*pad, model.builder.grid.num_nodes_y+2*pad, 100, 100, 5, 100, 100)
nodeGrid.num_nodes_x = model.builder.grid.num_nodes_x + 2*pad
nodeGrid.num_nodes_y = model.builder.grid.num_nodes_y + 2*pad


start_time, end_time = max_age, 0
modelNamePrefix = 'test1_'
mms2, mms_tti = [], []
out_dir = './out/'

no_steps = 4
dt = 314712e8 / no_steps

for tti in range(start_time, end_time-1,-1): 
    rebuild_mesh = (tti==start_time)
    
    # build_tti = 0     # use the full stack at every time step
    build_tti = tti   # do sedimentation: update stack positions at time steps

    if rebuild_mesh:
        print("Rebuild/reload mesh at tti=", tti)          
        mm2 = UniformNodeGridFixedSizeMeshModel(nodeGrid, nodes, modelName=modelNamePrefix+str(tti), sedimentsOnly=False)
        mm2.buildMesh(build_tti)
    else:
        print("Re-generating mesh vertices at tti=", tti)
        mm2.updateMesh(build_tti)

    mm2.baseFluxMagnitude = 0.04 if (tti>50) else 0.04
    print("===",tti,"  ",mm2.baseFluxMagnitude,"=========== ")
    if ( len(mms2) == 0):
        mm2.setupSolverAndSolve(no_steps=no_steps, time_step = dt, skip_setup=False)   
    else:    
        mm2.setupSolverAndSolve( no_steps=no_steps, time_step=dt, skip_setup=(not rebuild_mesh))
    if (True):
        mm2.writeLayerIDFunction(out_dir+"LayerID-"+str(tti)+".xdmf", tti=tti)
        mm2.writeTemperatureFunction(out_dir+"Temperature-"+str(tti)+".xdmf", tti=tti)
        mm2.writePoroFunction(out_dir+"Poro0-"+str(tti)+".xdmf", tti=tti)
        if (tti==0):
            fn = out_dir+"mesh-pos-"+str(tti)+".npy"
            np.save(fn, mm2.mesh.geometry.x)
            print("np save", fn, mm2.mesh.geometry.x.shape)
            fn = out_dir+"T-value-"+str(tti)+".npy"
            np.save(fn, mm2.uh.x.array[:])            
            fn = out_dir+"cell-pos-"+str(tti)+".npy"
            print("cell pos shape", mm2.getCellMidpoints().shape )
            np.save(fn, mm2.getCellMidpoints())
            print("thermal cond shape", mm2.thermalCond.x.array.shape)
            fn = out_dir+"k-value-"+str(tti)+".npy"
            np.save(fn, mm2.thermalCond.x.array[:])
            fn = out_dir+"phi-value-"+str(tti)+".npy"
            np.save(fn, mm2.mean_porosity.x.array[:])

    mms2.append(mm2)
    mms_tti.append(tti)

EPCfilename = mm2.write_hexa_mesh_resqml("temp/")
print("RESQML model written to: " , EPCfilename)
# read_mesh_resqml(EPCfilename, meshTitle="hexamesh")  # test reading of the .epc file


# hx = model.builder.grid.num_nodes_x // 2
# hy = model.builder.grid.num_nodes_y // 2
hx = model.builder.grid.num_nodes_x - 1 - pad
hy = model.builder.grid.num_nodes_y - 1 - pad
# hx = 1
# hy = 1

nn = model.builder.nodes[hy][hx]
dd = nn.depth_out[:,0]

node_ind = hy*model.builder.grid.num_nodes_x + hx
v_per_n = int(mm2.mesh_vertices.shape[0]/(model.builder.grid.num_nodes_y*model.builder.grid.num_nodes_x))

temp_1d = np.nan_to_num(nn.temperature_out[:,0], nan=5.0)
temp_3d_ind = np.array([ np.where([mm2.mesh_reindex==i])[1][0] for i in range(node_ind*v_per_n, (node_ind+1)*v_per_n) ] )
dd_mesh = mm2.mesh.geometry.x[temp_3d_ind,2]
temp_3d_mesh = mm2.u_n.x.array[temp_3d_ind]

temp_1d_at_mesh_pos = np.interp(dd_mesh, dd, temp_1d)
dd_subset = np.where(dd_mesh<5000)
print(f'Max. absolute error in temperature at 3D mesh vertex positions: {np.amax(np.abs(temp_1d_at_mesh_pos-temp_3d_mesh))}')
print(f'Max. absolute error at depths < 5000m: {np.amax(np.abs(temp_1d_at_mesh_pos[dd_subset]-temp_3d_mesh[dd_subset]))}')


# when not in ipython:
# import matplotlib as plt
# plt.use('qtagg')
#
# %matplotlib
# import matplotlib as plt
# plt.use('TkAgg')

# plt.pyplot.plot( dd, temp_1d, label=f'1D simulation'); 
# # plt.pyplot.plot( dd, temp_3d, label=f'3D simulation'); 
# plt.pyplot.plot( dd_mesh, temp_3d_mesh, 'o', label=f'3D simulation (nodes)'); 
# plt.pyplot.legend(loc="lower right", prop={'size': 7})
# plt.pyplot.show()



