import sys, os
from bmtk.simulator import bionet
import numpy as np
import pandas as pd
import h5py
from neuron import h
from scipy.stats import skew
import synapses
from bmtk.simulator.bionet.pyfunction_cache import add_weight_function

synapses.load()

pc = h.ParallelContext()  # object to access MPI methods
MPI_size = int(pc.nhost())
MPI_rank = int(pc.id())

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        inp = sys.argv[-1]
    else:
        raise Exception("no work" + str(sys.argv[-1]))

fname = str(inp)

config_file = 'simulation_config.json'

def gaussianBL(edge_props, source, target):
    w0 = edge_props["syn_weight"]
    sigma = edge_props["weight_sigma"]
    try:
        maximum = edge_props["weight_max"]
        return min(maximum, np.random.normal(w0, sigma, 1))
    except:
        return np.random.normal(w0, sigma, 1)

def lognormal(edge_props, source, target):
    m = edge_props["syn_weight"]
    s = edge_props["weight_sigma"]
    mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
    std = np.sqrt(np.log((s/m)**2 + 1))

    try:
        maximum = edge_props["weight_max"]
        return max(min(maximum, np.random.lognormal(mean, std, 1)), 0)
    except:
        return max(0, np.random.lognormal(mean, std, 1))

add_weight_function(lognormal)
add_weight_function(gaussianBL)


conf = bionet.Config.from_json(config_file, validate=True)
conf.build_env()

graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)

# cells = graph.get_local_cells()

# exc_strengths = {}
# inh_strengths = {}

# for gid, cell in cells.items():
#     exc_strens = []
#     inh_strens = []
#     for con in cell._connections:
#         if con._edge_prop.source_population == 'exc_stim':
#             #exc_strens.append(con.syn_weight)
#             exc_strens.append(con._syn.initW)
#         elif con._edge_prop.source_population == 'inh_stim':
#             #inh_strens.append(con.syn_weight)
#             inh_strens.append(con._syn.initW)
#         else:
#             raise Exception("Source pop is: " + str(con._edge_prop.source_population))

    
#     exc_strengths[gid] = exc_strens
#     inh_strengths[gid] = inh_strens

sim.run()
pc.barrier()

raster_file = './output/spikes.h5'

frs = {}
local_gids = list(graph.local_gids)

for key in local_gids:
    frs[key] = 0

try:
    f = h5py.File(raster_file,'r')
    
    spike_keys = list(f['spikes'].keys())
    if len(spike_keys) > 1:
        raise Exception("Spike keys: " + str(spike_keys))

    spike_key = list(f['spikes'].keys())[0]
    timestamps = f['spikes'][spike_key]['timestamps'].value
    gids = f['spikes'][spike_key]['node_ids'].value

    for i in range(len(gids)):
        if gids[i] in local_gids and timestamps[i] >= 200:
            frs[gids[i]] += 1
except:
    print("No spikes.")

dicts = [{"gid": gid, "FR": frs[gid] / 5} for gid in local_gids]
# dicts = [{"gid": gid, "FR": frs[gid] / 5, "num_exc": len(exc_strengths[gid]), "num_inh": len(inh_strengths[gid]),
#             "avg_exc": np.mean(exc_strengths[gid]), "avg_inh": np.mean(inh_strengths[gid]), 
#             "max_exc": np.max(exc_strengths[gid]), "max_inh": np.max(inh_strengths[gid]),
#             "std_exc": np.std(exc_strengths[gid]), "std_inh": np.std(inh_strengths[gid]),
#             "skew_exc": skew(exc_strengths[gid]), "skew_inh": skew(inh_strengths[gid])} for gid in local_gids]
df = pd.DataFrame(dicts)
#df.set_index("gid")
df.to_csv(fname+str(MPI_rank)+'.csv', index=False)

#import pdb; pdb.set_trace()

pc.barrier()

if MPI_rank == 0:
    base_df = pd.read_csv(fname+"0.csv", index_col="gid")
    res_df = pd.concat([base_df] + [pd.read_csv(fname+str(rank)+".csv", index_col="gid") for rank in range(1, MPI_size)])
    frs_df = pd.read_csv('frs_temp.csv', index_col="gid")
    res_df = res_df.join(frs_df)
    [os.remove(fname+str(rank)+".csv") for rank in range(MPI_size)]
    os.remove("frs_temp.csv")
    res_df.to_csv(fname+".csv")
