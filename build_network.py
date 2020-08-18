from bmtk.builder import NetworkBuilder
import numpy as np
import sys
import synapses             
import pandas as pd                                                                                                       
synapses.load()
syn = synapses.syn_params_dicts()

np.random.seed(3546)

if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        inp = sys.argv[-1]
    else:
        raise Exception("no work" + str(sys.argv[-1]))

N = int(inp)

# Initialize our network

net = NetworkBuilder("biophysical")

def lognormal(m, s):
        mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
        std = np.sqrt(np.log((s/m)**2 + 1))
        return max(np.random.lognormal(mean, std, 1), 0)

#num_inh = [int(lognormal(56, 7.5)) for i in range(N)]
num_inh = [56 for i in range(N)]
#print(num_inh)
inh_bounds = []
sum = 0
for num in num_inh:
        sum += num
        inh_bounds.append(sum)

#num_exc = [int(lognormal(1000, 80)) for i in range(N)]
num_exc = [1000 for i in range(N)]
#print(num_exc)
exc_bounds = []
sum = 0
for num in num_exc:
        sum += num
        exc_bounds.append(sum)


exc_frs = np.random.uniform(0, 2, N)
inh_frs = np.random.uniform(8, 12, N)
print("FRS:",exc_frs,inh_frs)

fr_df = pd.DataFrame()
fr_df['gid'] = [i for i in range(N)]
fr_df.set_index('gid')
fr_df['exc_fr'] = exc_frs
fr_df['inh_fr'] = inh_frs

fr_df.to_csv('frs_temp.csv')

# exc_fr = 1.09
# inh_fr = 10

##################################################################################
###################################Pyr Type C#####################################

net.add_nodes(N=N, pop_name='Pyrc',
    potental='exc',
    model_type='biophysical',
    model_template='ctdb:Biophys1.hoc',
    model_processing='aibs_allactive',
    dynamics_params='491766131_fit.json',
    morphology='Rbp4-Cre_KL100_Ai14-203503.04.01.01_527109145_m.swc')


##################################################################################
###################################External Networks##############################

#print("Internal nodes built")

# External excitatory inputs
exc_stim = NetworkBuilder('exc_stim')
exc_stim.add_nodes(N=np.sum(num_exc),
                pop_name='exc_stim',
                potential='exc',
                model_type='virtual')

# External inhibitory inputs
inh_stim = NetworkBuilder('inh_stim')
inh_stim.add_nodes(N=np.sum(num_inh),
                pop_name='inh_stim',
                potential='exc',
                model_type='virtual')

##################################################################################
###################################Edges##########################################

def correct_cell(source, target, bounds):
        sid = source.node_id
        tid = target.node_id

        lower_bound = 0
        if tid > 0:
                lower_bound = bounds[tid - 1]

        upper_bound = bounds[tid]

        if sid < upper_bound and sid >= lower_bound:
                #print("connecting cell {} to {}".format(sid,tid))
                return 1
        else:
                return None

#Create connections between Inh --> Pyr cells
# net.add_edges(source=inh_stim.nodes(), target=net.nodes(),
#        connection_rule=correct_cell,
#        connection_params={'bounds': inh_bounds},
#        syn_weight=1e-3,
#        #weight_function='lognormal',
#        weight_sigma=1e-3,
#        weight_max=20e-3,
#        dynamics_params='GABA_InhToExc.json',
#        model_template='Exp2Syn',
#        distance_range=[0.0, 300.0],
#        target_sections=['somatic'],
#        delay=2.0)
# Create connections between Exc --> Pyr cells
net.add_edges(source=inh_stim.nodes(), target=net.nodes(),
                connection_rule=correct_cell,
                connection_params={'bounds': inh_bounds},
                syn_weight=1,
                delay=0.1,
                dynamics_params='INT2PN.json',
                model_template=syn['INT2PN.json']['level_of_detail'],
                distance_range=[0.0, 300.0],
                target_sections=['somatic'])

# Create connections between Exc --> Pyr cells
net.add_edges(source=exc_stim.nodes(), target=net.nodes(),
                connection_rule=correct_cell,
                connection_params={'bounds': exc_bounds},
                syn_weight=1,
                target_sections=['dend'],
                delay=0.1,
                distance_range=[0.0, 300.0],
                dynamics_params='PN2PN.json',
                model_template=syn['PN2PN.json']['level_of_detail'])


# Build and save our networks
net.build()
net.save_nodes(output_dir='network')
net.save_edges(output_dir='network')

exc_stim.build()
exc_stim.save_nodes(output_dir='network')

inh_stim.build()
inh_stim.save_nodes(output_dir='network')


from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
from bmtk.utils.reports.spike_trains.spikes_file_writers import write_csv

exc_psg = PoissonSpikeGenerator(population='exc_stim')
# exc_psg.add(node_ids=range(np.sum(num_exc)),  
#         firing_rate=exc_fr / 1000,    
#         times=(200.0, 5200.0))    
for i in range(N):
        exc_psg.add(node_ids=range(exc_bounds[i] - num_exc[i], exc_bounds[i]),
                firing_rate=exc_frs[i] / 1000,
                times=(200.0, 5200.0))
exc_psg.to_sonata('exc_stim_spikes.h5')

inh_psg = PoissonSpikeGenerator(population='inh_stim')
# inh_psg.add(node_ids=range(np.sum(num_inh)), 
#         firing_rate=inh_fr / 1000,  
#         times=(200.0, 5200.0))   
for i in range(N):
        inh_psg.add(node_ids=range(inh_bounds[i] - num_inh[i], inh_bounds[i]),
                firing_rate=inh_frs[i] / 1000,
                times=(200.0, 5200.0))
inh_psg.to_sonata('inh_stim_spikes.h5')


from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
                network_dir='./network',
                tstop=5200.0, dt = 0.1,
                report_vars=['v'],
                spikes_inputs=[('exc_stim', 'exc_stim_spikes.h5'), ('inh_stim', 'inh_stim_spikes.h5')],
                components_dir='biophys_components',
                compile_mechanisms=True)
