import sys
from bmtk.simulator import bionet
import numpy as np
import pandas as pd
import h5py

from bmtk.simulator.bionet.pyfunction_cache import add_weight_function

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
    return np.random.normal(w0, sigma, 1)

def lognormal(edge_props, source, target):
    m = edge_props["syn_weight"]
    s = edge_props["weight_sigma"]
    mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
    std = np.sqrt(np.log((s/m)**2 + 1))
    return np.random.lognormal(mean, std, 1)

add_weight_function(lognormal)
add_weight_function(gaussianBL)


conf = bionet.Config.from_json(config_file, validate=True)
conf.build_env()

graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)

cells = graph.get_local_cells()

exc_strengths = {}
inh_strengths = {}

for gid, cell in cells.items():
    exc_strens = []
    inh_strens = []
    for con in cell._connections:
        if con._edge_prop.source_population == 'exc_stim':
            exc_strens.append(con.syn_weight)
        elif con._edge_prop.source_population == 'inh_stim':
            inh_strens.append(con.syn_weight)
        else:
            raise Exception("Source pop is: " + str(con._edge_prop.source_population))

    
    exc_strengths[gid] = exc_strens
    inh_strengths[gid] = inh_strens

sim.run()

import pdb; pdb.set_trace()

raster_file = './output/spikes.h5'

frs = {}
all_gids = [i for i in range(len(exc_strengths))]

for key in all_gids:
    frs[key] = 0

try:
    f = h5py.File(raster_file,'r')
    timestamps = f['spikes']['inh_stim']['timestamps'].value
    gids = f['spikes']['inh_stim']['node_ids'].value

    for i in range(len(gids)):
        if timestamps[i] >= 200:
            frs[gids[i]] += 1
except:
    print("No spikes.")

df = pd.DataFrame()
dicts = [{"FR": frs[gid], "num_exc": len(exc_strengths[gid]), "num_inh": len(inh_strengths[gid]),
            "avg_exc": np.mean(exc_strengths[gid]), "avg_inh": np.mean(inh_strengths[gid]), 
            "max_exc": np.max(exc_strengths[gid]), "max_inh": np.max(inh_strengths[gid])} for gid in all_gids]
df = pd.DataFrame(dicts)
df.to_csv(fname, index=False)