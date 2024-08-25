# IMPORTS

import csv
import os
import random
from neuron import h, gui
from neuron.units import ms, mV, um, V, cm
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from tqdm.notebook import tqdm

from time import time
import os.path

from socket import gethostname

h.load_file('model-files/basal_project.hoc')

id = int(os.getenv('SLURM_ARRAY_TASK_ID'))
g_index = int(np.floor(id/15)) # Determines which g_inhib value to use

# Establishes trial number
if id == 1:
    trial = 1
else:
    trial = id % 15 + 1

# BARDIA CELL INITIATION (ALTERED)

ratio=0
loc=10
nsyn=0
sec=h.a1_111

dendrec=True
iotim = 0
tic = time()

h.synrec = 1

h.tsamp = h.tstop/h.dt+1
h.synrec = 1
tsamp = int(h.tsamp)

r = h.Random(h.luckyoffset)
r.negexp(h.meanisi)

# h.a1_111.nseg = int(np.floor(h.a1_111.L))
# h.a1_111.nseg = int(np.floor(h.a1_111.L)/2)
seclist = [h.a1_111,
           h.a9_122,
           h.a10_11,
           h.a5_1,
           h.a8_11,
           h.a3_11,
           h.a9_122,
           h.a9_121,
           h.a8_122,
           h.a8_121,
           h.a7_111]

sl2 = h.SectionList()
for sec in seclist:
    sl2.append(sec = sec)
poppedsecs = sl2.unique()
h.refreshnseg(h.makeactivelist(sl2))
h.cvode.cache_efficient(1)
h.cvode_active(0)
h.poisson = 0

######################### BIOPHYSICAL FUNCTIONS #########################

def add_spines(branches):
    """
    Adds spines to all branches within the branch list
    5 spines per segment on each branch

    Returns nested list of spines within each branch
    """

    # Intialize vectors for the spine head/neck
    spine_head = []
    spine_neck = []

    # Iterates through the number of branches of interest (list)
    num_branch = len(branches)
    for j in range(num_branch):

        # branches[j].nseg = int((np.floor(branches[j].L)) / 4)
        spine_head_temp = []
        spine_neck_temp = []

        # Number of spines: 5 per segment
        n_spines = branches[j].nseg * 5

        # Iterate through the number of spines for the particular branch
        for i in range(n_spines):

            # Create head and neck biophysical parameters
            spine_head_temp.append(h.Section())
            spine_neck_temp.append(h.Section())

            spine_head_temp[-1].L = 0.5
            spine_head_temp[-1].diam = 0.5
            spine_head_temp[-1].Ra = 100
            spine_head_temp[-1].insert('pas')
            for seg in spine_head_temp[-1]:
                seg.pas.e = -70
                seg.pas.g = 0.00005

            spine_neck_temp[-1].L = 1.5
            spine_neck_temp[-1].diam = 0.25
            spine_neck_temp[-1].Ra = 100
            spine_neck_temp[-1].insert('pas')
            for seg in spine_neck_temp[-1]:
                seg.pas.e = -70
                seg.pas.g = 0.00005

            # Connect head and neck together; then connect neck to branch
            spine_neck_temp[-1].connect(branches[j]((i + 1) / n_spines))
            spine_head_temp[-1].connect(spine_neck_temp[-1])

        spine_head.append(spine_neck_temp)
        spine_neck.append(spine_head_temp)

    return spine_neck, spine_head


def add_inhib_soma(num):
    """
    Adds GABA synapses to the soma
    num: number of synapses

    Returns list of GABA synapses on the soma

    """

    syn_gaba_soma = []
    for i in range(num):
        syn_gaba_soma.append(h.GABAa_S(h.soma(0.5)))

    return syn_gaba_soma


def add_syn(branches, inhib=False):
    """
    Adds AMPA, NMDA, and GABA synaptic channels to a list of branches
    branches: list of branches
    inhib: boolean, determines if inhibitory synapses will be added to the branches

    Returns three nested vectors containing the synapses
    """

    syn_ampa = []
    syn_nmda = []
    syn_gaba = []

    num_branch = len(branches)

    # Iterates through the number of branches and applies synapses to every spine
    for j in range(num_branch):

        syn_ampa_temp = []
        syn_nmda_temp = []
        syn_gaba_temp = []

        # Indexes spine heads of the current branch
        n_spines = len(spine_head[j])
        for i in range(n_spines):

            # adds AMPA, NMDA, and GABA synapses to each spine
            syn_ampa_temp.append(h.Exp2Syn(spine_head[j][i](0.5)))
            syn_nmda_temp.append(h.Exp2SynNMDA(spine_head[j][i](0.5)))
            syn_ampa_temp[-1].e = 0
            syn_ampa_temp[-1].tau1 = 0.05
            syn_ampa_temp[-1].tau2 = 0.5
            syn_nmda_temp[-1].e = 0
            syn_nmda_temp[-1].tau1 = 2.1
            syn_nmda_temp[-1].tau2 = 18.8

            if inhib:
                syn_gaba_temp.append(h.GABAa_S(spine_head[j][i](0.5)))

        syn_ampa.append(syn_ampa_temp)
        syn_nmda.append(syn_nmda_temp)
        syn_gaba.append(syn_gaba_temp)

    return syn_ampa, syn_nmda, syn_gaba


######################### CURRENT INJECTIONS ############################

def frequency_counter(t_vec, v_vec):
    """
    Returns the frequency of spikes measured at the soma

    t_vec: time vector
    v_vec: voltage vector

    """

    # Chops off the first 100ms, in which the model is still reaching steady-state
    t_tot = t_vec[-1]
    start_t_i = np.where(np.floor(t_vec) == 100)[0][0]

    spike_counter = 0
    length_v = len(v_vec)

    # Iterates from 100ms to end of time for the vector; calculates spikes
    for i in range(start_t_i, length_v - 1):
        if v_vec[i] < 0 and v_vec[i + 1] >= 0:
            spike_counter += 1

    # Determines the frequency of spikes
    frequency = (spike_counter / (t_tot - 100)) * 10 ** 3

    return frequency


def current_trial(n_gaba,
                  g,
                  firing,
                  dtb,
                  dendrite_gaba=False,
                  distance=0.5):
    """
    Determines the amount of current needed to inject for background conditions of ~0.5 Hz

    max_c: the maximum amount of injected current at the soma (nA)
    inc: the size step of the amount of current (nA)
    avg: the number of times each FI point will be averaged along the curve
    g: conductance of GABA synapses
    dtb: distribution of GABA synpse firing

    """

    h.tstop = 1100  # ms
    h.v_init = -70  # mV  

    # Iterates through increasing current injections
    loopFlag = True
        
    inc = 0.01               

    frequency_list = []
    amp_list = []
    
    # Uses previous current value as maximum and evaluates where the cell begins to jump off from 0.5 Hz
    
    for amp in tqdm(np.linspace(1, 3.5, int(1 / inc + 1))):

        # Performs averaged trials at a specific current value
        # Creates a background GABA stimulus
        ncstim_ = None
        ncstim_, syn_stim_gaba = background_simulation(basal_dendrites,
                                                       n_gaba,
                                                       g,
                                                       firing,
                                                       dtb,
                                                       dendrite_gaba,
                                                       distance)

        # Creates current injection
        temp = []
        temp = h.IClamp(h.soma(1))
        temp.delay = 0
        temp.dur = 1e9
        temp.amp = amp

        # Calculates voltage trace
        soma_v = h.Vector().record(h.soma(1)._ref_v)
        t_vec = h.Vector().record(h._ref_t)
        h.finitialize()
        h.run()

        temp = None
        ncstim_ = None
        syn_stim_gaba = None

        frequency_list.append(frequency_counter(t_vec, soma_v))
        amp_list.append(amp)
    
    return amp_list, frequency_list


##################### SYNAPSE BEHAVIOR ##########################

def background_simulation(branches,
                          n_gaba,
                          g,
                          firing=50,
                          dtb=3,
                          dendrite_gaba=False,
                          distance=0.5):
    """
    Simulates background synapses along basal dendrites,
    in which the soma has basket cell-like inhibitory synapses

    Firing rate is fixed to 5 Hz with a Guassian distribution of start times

    branches: list of branches
    n_gaba: number of GABA synapses at the soma
    g: conductance of GABA synapses (µS)
    firing: firing rate of GABA-gamma synapses [50 Hz default]
    dtb: std of distribution for firing rates [Hz]
    dendrite_gaba: determines if GABA synapses are added to the basal dendrites
    distance: determines the distance GABA synapses are added to basal dendrites (can be decimal or exact value)
    """

    # GABA Stimulations; all fire at 5 Hz (no noise) at Gaussian-distributed start times
    syn_stim_gaba = []

    # Determine interval between GABA synapses (dt)
    interval_gaba = 1000 / firing

    # Calculates number of times GABA fires within time period
    num_firing = int(np.floor(h.tstop / interval_gaba))

    # Iterates through GABA firing
    for i in range(num_firing + 1):

        # For each synapse, generates distribution of firing for each time point
        gaba_start_distribution = np.random.normal(0, dtb, n_gaba)

        # Iterate through all firing instances for each synapse
        for j in range(n_gaba):
            syn_stim_gaba.append(create_stim(interval=interval_gaba,
                                             num=1,
                                             start=gaba_start_distribution[j] + interval_gaba * i,
                                             noise=0,
                                             s=1))

    num_branch = len(branches)
    ncstim = []

    # Iterate through GABA firing periods
    for q in range(num_firing + 1):

        n_s = 0  # Number of synapses
        c = 0  # Spread of repeated branch synapses
        i = 1  # Index of branches

        # Iterate through the number of synapses, then shift the firing period
        while n_s < n_gaba:

            # Place stimulus at x distance along dendrite of the way along the basal dendrite
            branch_i = branches[i]
            branch_nseg = branch_i.nseg

            loc = int(np.floor(distance * branch_nseg))

            # Append firing instances to SOMA GABA synapses
            ncstim.append(h.NetCon(syn_stim_gaba[n_s + q * n_gaba], syn_gaba_soma[n_s], 0, 0, g))
            n_s += 1

            # Reinitialize branch index to first branch, increase branch spread by one
            if i == num_branch - 1:
                i = 1
                c = c + 1


    return ncstim, syn_stim_gaba


def create_stim(interval=10, num=3, start=5, noise=0, s=1):
    """
    Creates a stimulator object used as the pre-synaptic stimulus
    interval: amount of time between stimulation events (ms)
    num: the number of stimulation events
    start: delay of stimulation (ms)
    s: seed number
    """

    stim = h.NetStim()
    stim.interval = interval
    stim.number = num
    stim.start = start
    stim.noise = noise

    if s == 1:
        stim.seed(np.random.randint(1000000))
    else:
        stim.seed(s)

    return stim

basal_dendrites = [h.a1_111,
                   h.a9_122,
                   h.a10_11,
                   h.a5_1,
                   h.a8_11,
                   h.a3_11,
                   h.a9_122,
                   h.a9_121,
                   h.a8_122,
                   h.a8_121,
                   h.a7_111]

n_gaba = 30

spine_neck, spine_head = add_spines(basal_dendrites)
syn_ampa, syn_nmda, syn_gaba = add_syn(basal_dendrites, inhib=True)
syn_gaba_soma = add_inhib_soma(n_gaba)

inhib_g = [0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015]
x = inhib_g[g_index]

a_list, f_list = current_trial(n_gaba=30,
                               g=x,
                               firing=50,
                               dtb=3,
                               dendrite_gaba=False,
                               distance=0.5)   

filename = "current-collection/group" + str(g_index) + "-trial" + str(trial) + ".csv"

with open(filename, 'w', newline='') as f:
     
    # using csv.writer method from CSV package
    write = csv.writer(f, delimiter=',')
    write.writerow(a_list)
    write.writerow(f_list)