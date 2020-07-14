import os, sys
from bmtk.simulator import bionet
import numpy as np
import pandas as pd
import h5py

def run(config_file):


    
    from bmtk.simulator.bionet.pyfunction_cache import add_weight_function

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

    # exc_strengths = []
    # for edge in graph.find_edges(source_nodes='exc_stim')[0].get_edges():
    #     exc_strengths.append(edge.syn_weight())

    # inh_strengths = []
    # for edge in graph.find_edges(source_nodes='inh_stim')[0].get_edges():
    #     inh_strengths.append(edge.syn_weight())

    exc_strengths = []
    inh_strengths = []

    pn = graph.get_cell_gid(0)
    for con in pn._connections:
        if con._edge_prop.source_population == 'exc_stim':
            exc_strengths.append(con.syn_weight)
        elif con._edge_prop.source_population == 'inh_stim':
            inh_strengths.append(con.syn_weight)
        else:
            raise Exception("Source pop is: " + str(con._edge_prop.source_population))

    # sim = bionet.BioSimulator.from_config(conf, network=graph)
    sim.run()
    #bionet.nrn.quit_execution()

    raster_file = './output/spikes.h5'
    try:
        fr = 0
        f = h5py.File(raster_file,'r')
        timestamps = f['spikes']['inh_stim']['timestamps'].value
        for t in timestamps:
            if t >= 200:
                fr += 1
    except:
        fr = 0

    res_file = "results.csv"
    df = pd.read_csv(res_file)
    new_line = {"FR": fr, "num_exc": len(exc_strengths), "num_inh": len(inh_strengths),
                "avg_exc": np.mean(exc_strengths), "avg_inh": np.mean(inh_strengths), 
                "max_exc": np.max(exc_strengths), "max_inh": np.max(inh_strengths)}
    print(new_line)
    df = df.append(new_line, ignore_index=True)
    df.to_csv(res_file, index=False)


if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        run(sys.argv[-1])
    else:
        run('simulation_config.json')
