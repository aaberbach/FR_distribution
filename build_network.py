from bmtk.builder import NetworkBuilder
import os, sys
import numpy as np
from bmtk.builder.auxi.node_params import positions_cuboid, positions_list

#np.random.seed(123412)

# Initialize our network

net = NetworkBuilder("SPWR_biophysical")

def lognormal(m, s):
        mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
        std = np.sqrt(np.log((s/m)**2 + 1))
        return np.random.lognormal(mean, std, 1)

# num_inh = int(np.random.lognormal(43, 13, 1))
# num_exc = int(np.random.lognormal(25, 10, 1))
num_inh = int(lognormal(43, 13))
num_exc = int(lognormal(25, 10))
exc_fr = 2
inh_fr = 10

##################################################################################
###################################Pyr Type C#####################################

net.add_nodes(N=1, pop_name='PyrC',
              mem_potential='e',
              model_type='biophysical',
              model_template='hoc:stylized_typeC',
              morphology=None)


##################################################################################
###################################External Networks##############################

#print("Internal nodes built")

# External excitatory inputs
exc_stim = NetworkBuilder('exc_stim')
exc_stim.add_nodes(N=num_exc,
                   pop_name='exc_stim',
                   potential='exc',
                   model_type='virtual')

# External inhibitory inputs
inh_stim = NetworkBuilder('inh_stim')
inh_stim.add_nodes(N=num_inh,
                   pop_name='inh_stim',
                   potential='exc',
                   model_type='virtual')

##################################################################################
###################################Edges##########################################

# Create connections between Inh --> Pyr cells
net.add_edges(source=inh_stim.nodes(), target=net.nodes(),
              connection_rule=1,
              syn_weight=5.0e-03,
              weight_function='lognormal',
              weight_sigma=1.0e-03,
              dynamics_params='GABA_InhToExc.json',
              model_template='Exp2Syn',
              distance_range=[0.0, 300.0],
              target_sections=['somatic'],
              delay=2.0)

# Create connections between Exc --> Pyr cells
net.add_edges(source=exc_stim.nodes(), target=net.nodes(),
                   connection_rule=1,
                   syn_weight=12.0e-03,
                   weight_function='gaussianBL',
                   weight_sigma=1.0e-03,
                   target_sections=['somatic'],
                   delay=2.0,
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='Exp2Syn')


# Build and save our networks
net.build()
net.save_nodes(output_dir='network')
net.save_edges(output_dir='network')

exc_stim.build()
exc_stim.save_nodes(output_dir='network')

inh_stim.build()
inh_stim.save_nodes(output_dir='network')


from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator

exc_psg = PoissonSpikeGenerator(population='exc_stim')
exc_psg.add(node_ids=range(num_exc),  
        firing_rate=int(exc_fr) / 1000,    
        times=(200.0, 1200.0))    
exc_psg.to_sonata('exc_stim_spikes.h5')

inh_psg = PoissonSpikeGenerator(population='inh_stim')
inh_psg.add(node_ids=range(num_inh), 
        firing_rate=int(inh_fr) / 1000,  
        times=(200.0, 1200.0))   
inh_psg.to_sonata('inh_stim_spikes.h5')


from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
		network_dir='./network',
		tstop=1200.0, dt = 0.1,
		report_vars=['v'],
        spikes_inputs=[('exc_stim', 'exc_stim_spikes.h5'), ('inh_stim', 'inh_stim_spikes.h5')],
		components_dir='biophys_components',
		compile_mechanisms=True)