'''
This file lists out the each of the cells of the model (Pyr, Inh, RE, TC).
'''

from netpyne import specs, sim 
from neuron import h 
import numpy as np 
import config
import syn_conns

netParams = specs.NetParams()
simConfig = specs.SimConfig()

simConfig.oneSynPerNetcon = False # Allows for mutliple synapses to connect to one NetCon 
netParams.defaultThreshold = 0 # mV, setting early as NetPyNE's deafault is -10 mV, whereas NEURON in 0 mV
raster_filename = 'raster_test_again.pkl' # raster_filename located here to change raster names quickly as the rest of the code should remain constant
# config.NP_fac = 0.8 # used for the time being to alter AMPA_d2 synaptic weights for each sleep state (lines 169, 170, and 172) as well as pyr2pyr_nmda_d1_str (line 111)
# config.NP_fac=0.8 should be the default

NPyr = config.Npyr # number of pyramidal cells in network 
NInh = config.Ninh # number of inhibitory cells in network 
NRE = config.Nre # number of reticular cells in nework 
NTC = config.Ntc # number of thalamocortical cells in nework 

##############################################################################
# Population Parameters
##############################################################################
netParams.popParams['PyrCell'] = {'cellType': 'Pyr', 'numCells': NPyr}
netParams.popParams['InhCell'] = {'cellType': 'Inh', 'numCells': NInh}
netParams.popParams['RECell'] = {'cellType': 'RE', 'numCells': NRE}
netParams.popParams['TCCell'] = {'cellType': 'TC', 'numCells': NTC}

##############################################################################
# Cell Parameters
##############################################################################
'''
Pyr Cell 
'''
soma_L = soma_diam = np.sqrt(100/np.pi)
geom_factor = (soma_L/2)/(np.pi * (soma_diam/2)**2) + ((165*soma_L)/2)/(np.pi * (soma_diam/2)**2)
cellRule = {'conds': {'cellType': 'Pyr'}, 'secs': {}, 'globals': {}}
cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}, 'ions':{}} # template for soma characteristics 
cellRule['secs']['dend']= {'geom': {}, 'topol': {}, 'mechs': {}, 'ions': {}, 'synMechParams': {}} # template for the dendrite characterisitcs 
cellRule['secs']['soma']['geom'] = { # geometry of the soma 
        'diam': soma_diam, # microns;gives an area of 100 microns^2 
        'L': soma_L, # microns
        'cm': 0.75, # uF/cm2
        'Ra': 1e7/geom_factor * 1e-4}# this post is extremely helpful for setting the axial resistivity to match the absolute values described in the line above: https://www.neuron.yale.edu/phpBB/viewtopic.php?f=8&t=3904&p=16841&hilit=axial+resistance#p16841
        #axial resistance from middle of soma to end where it meets dendrite is soma.L/2; then divide by area and multiple by axial resistivity to get resistance. Then do the same to calculate resistance from soma/dendrite junction to halfway along dendrite
        # 1e-4 converts the ohm * microns into ohm * cm 
cellRule['secs']['dend']['geom'] = { # geometry of the dend 
        'diam': soma_diam, # microns;gives an area of 100 microns^2  }
        'L': 165*soma_L, # microns
        'cm': 0.75, # uF/cm2
        'Ra': 1e7/geom_factor * 1e-4}# this post is extremely helpful for setting the axial resistivity to match the absolute values described in the line above: https://www.neuron.yale.edu/phpBB/viewtopic.php?f=8&t=3904&p=16841&hilit=axial+resistance#p16841
        #axial resistance from middle of soma to end where it meets dendrite is soma.L/2; then divide by area and multiple by axial resistivity to get resistance. Then do the same to calculate resistance from soma/dendrite junction to halfway along dendrite
        # 1e-4 converts the ohm * microns into ohm * cm 
cellRule['secs']['dend']['topol'] = {'parentSec': 'soma', 'parentX': 1.0, 'childX': 0.0} # connect the dend to the soma

# Soma mechanisms
cellRule['secs']['soma']['mechs']['kdr'] = { # delayed-rectifying potassium current
        'gkbar': 0.200} # S/cm2 
cellRule['secs']['soma']['mechs']['naf'] = { # regular and fast-Na+ current
        'gnabar': 3.0} # S/cm2 
cellRule['secs']['soma']['mechs']['nap'] = { # persistant Na+ current
        'gnabar': 0.0003} # S/cm2 
cellRule['secs']['soma']['ions']['na'] = { # allows for the change in reversal potential, internal concentration, and external concentration of sodium
        'e': 50} # mV
cellRule['secs']['soma']['ions']['k'] = { # allows for the change in reversal potential, internal concentration, and external concentration of potassium
        'e': -90} # mV

# Dend mechanisms 
cellRule['secs']['dend']['mechs']['naf'] = { # regular and fast sodium current
        'gnabar': 0.0008} # S/cm2 
cellRule['secs']['dend']['mechs']['nap'] = { # persistant sodium current
        'gnabar': 0.000042} # S/cm2
cellRule['secs']['dend']['mechs']['hva'] = { # HVA calcium current
        'gcabar': 0.000012} #, S/cm2 
        #'eca': 140} # mV; must be commented out until forum post answered 
cellRule['secs']['dend']['mechs']['km'] = { # slow non-inactivating potassium current
        'gkbar': 0.00002} # S/cm2
cellRule['secs']['dend']['mechs']['cad'] = { # calcium accumulation mechanism
        'taur': 165, # 165 
        'depth': 1.0, # um
        'cainf': 2.4e-4} # mM
cellRule['secs']['dend']['mechs']['kca'] = { # clacium-dependent potassium current
        'gkcabar': 0.00005} # S/cm2 
cellRule['secs']['dend']['mechs']['kL'] = { # potassium leak current
        'gkL': config.gkl_pyr_baseline * config.init_gkl} # S/cm2
cellRule['secs']['dend']['mechs']['pas'] = { # passive leak current
        'g': 0.000011, # S/cm2
        'e': -67} # mV
cellRule['secs']['dend']['ions']['ca'] = { # allows for change in reversal potential, internal concentration, and external concentration of calcium 
        'e': 140} # mV
cellRule['secs']['dend']['ions']['k'] = { # allows for change in reversal potential, internal concentration, and external concentration of potassium
        'e': -90} # mV 
cellRule['secs']['dend']['ions']['na'] = { # allows for change in reversal potential, internal concentration, and external concentration of sodium 
        'e': 50} # mV 

cellRule['globals']['krev'] = -95 # mV

netParams.cellParams['Pyr'] = cellRule 

'''
Inh Cell
'''
soma_L = soma_diam = np.sqrt(100/np.pi)
geom_factor = (soma_L/2)/(np.pi * (soma_diam/2)**2) + ((50*soma_L)/2)/(np.pi * (soma_diam/2)**2)
cellRule = {'conds': {'cellType': 'Inh'}, 'secs': {}, 'globals': {}}
cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}, 'ions': {}} # template for soma characteristics 
cellRule['secs']['dend'] = {'geom': {}, 'topol': {}, 'mechs': {}, 'ions': {}, 'globals': {}} # template for dend characteristics
cellRule['secs']['soma']['geom'] = { # geometrty of soma 
        'diam': soma_diam, # microns;gives an area of 100 microns^2 
        'L': soma_L, # microns
        'cm': 0.75, # uF/cm2
        'Ra': 1e7/geom_factor * 1e-4}# this post is extremely helpful for setting the axial resistivity to match the absolute values described in the line above: https://www.neuron.yale.edu/phpBB/viewtopic.php?f=8&t=3904&p=16841&hilit=axial+resistance#p16841
        #axial resistance from middle of soma to end where it meets dendrite is soma.L/2; then divide by area and multiple by axial resistivity to get resistance. Then do the same to calculate resistance from soma/dendrite junction to halfway along dendrite
        # 1e-4 converts the ohm * microns into ohm * cm 
cellRule['secs']['dend']['geom'] = { # geometry of the dend 
        'diam': soma_diam, # microns;gives an area of 100 microns^2  }
        'L': 50*soma_L, # microns
        'cm': 0.75, # uF/cm2
        'Ra': 1e7/geom_factor * 1e-4}# this post is extremely helpful for setting the axial resistivity to match the absolute values described in the line above: https://www.neuron.yale.edu/phpBB/viewtopic.php?f=8&t=3904&p=16841&hilit=axial+resistance#p16841
        #axial resistance from middle of soma to end where it meets dendrite is soma.L/2; then divide by area and multiple by axial resistivity to get resistance. Then do the same to calculate resistance from soma/dendrite junction to halfway along dendrite
        # 1e-4 converts the ohm * microns into ohm * cm 
cellRule['secs']['dend']['topol'] = {'parentSec': 'soma', 'parentX': 1.0, 'childX': 0.0} # connect the dend to the distal end of the soma

# Soma mechanisms
cellRule['secs']['soma']['mechs']['kdr'] = { # delayed rectifying potassium current 
        'gkbar': 0.20} # S/cm2 
cellRule['secs']['soma']['mechs']['naf'] = { # regular and fast sodium current
        'gnabar': 2.5} # S/cm2
cellRule['secs']['soma']['ions']['na'] = { # allows for change in reversal potential, internal concentration, and external concentration
        'e': 50} # mV 
cellRule['secs']['soma']['ions']['k'] = { # allows for change in reversal potential, internal concentration, and external concentration 
        'e': -90} # mV

# Dend mechanisms
cellRule['secs']['dend']['mechs']['naf'] = { # regular and fast sodium currents
        'gnabar': 0.0008} # S/cm2 
cellRule['secs']['dend']['mechs']['hva'] = { # HVA calcium current 
        'gcabar': 0.000012} # S/cm2 
cellRule['secs']['dend']['mechs']['km'] = { # slow, non-inactivating potassium current
        'gkbar': 0.000015} # S/cm2
cellRule['secs']['dend']['mechs']['cad'] = { # calcium accumulation mechanism
        'taur': 165, # ms 
        'depth': 1.0, # um
        'cainf': 2.4e-4} # mM
cellRule['secs']['dend']['mechs']['kca'] = { # calcium-dependent potassium current 
        'gkcabar': 0.00005} # S/cm2
cellRule['secs']['dend']['mechs']['kL'] = { # potassium leak current
        'gkL': config.gkl_inh_baseline * config.init_gkl} # S/cm2
cellRule['secs']['dend']['mechs']['pas'] = { # passive leak current
        'g': 0.000009, # S/cm2
        'e': -70.0} # mV
cellRule['secs']['dend']['ions']['na'] = { # allows for change in reversal potential, internal concentration, and external concentration
        'e': 50} # mV 
cellRule['secs']['dend']['ions']['k'] = { # allows for change in reversal potential, internal concentration, and external concentration
        'e': -90} # mV
cellRule['secs']['dend']['ions']['ca'] = { # allows for change in reversal potential, internal concentration, and external concentration
        'e': 140} # mV '

cellRule['secs']['dend']['globals']['krev'] = -95 # mV

netParams.cellParams['Inh'] = cellRule

'''
RE Cell
'''
cellRule = {'conds': {'cellType': 'RE'}, 'secs': {}, 'globals' : {}} # creates a template for the cellRule or applications of the cell 
cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}, 'ions': {}} # creates a soma with geometry, mechansims, and ions
cellRule['secs']['soma']['geom'] = {'diam': np.sqrt(1.43e4/np.pi), # one-compartment of 1.43e4 um2
                                    'L': np.sqrt(1.43e4/np.pi),
                                    'Ra': 100, 
                                    'cm': 1} # uF/cm2
cellRule['secs']['soma']['mechs']['pas'] = {
        'g': 5e-05, 
        'e': -77} # leak current
cellRule['secs']['soma']['mechs']['kL'] = {
        'gkL': config.gkl_RE_baseline * config.init_gkl_RE}  # S/cm2, awake state of kL
cellRule['secs']['soma']['mechs']['naf_re'] = { # fast-Na+ current
        'gnabar': 0.1} # S/cm2
cellRule['secs']['soma']['mechs']['kdr_re'] = { # delayed rectifying K+ current
        'gkbar': 0.01} # S/cm2 
cellRule['secs']['soma']['mechs']['it_re'] = { # low-threshold Ca++ current  
        'gcabar': 0.0022} # S/cm2
cellRule['secs']['soma']['mechs']['cad'] = { # fast mech for submebranal Ca++ concentration 
        'depth': 1, 
        'taur': 5,
        'cainf': 2.4e-4} # mM; equlilibrium Ca++ concentration 
cellRule['secs']['soma']['ions']['na'] = { # allows for change in reversal potential, internal concentration, and external concentration
        'e': 50} # mV
cellRule['secs']['soma']['ions']['k'] = { # allows for change in reversal potential, internal concentration, and external concentration
        'e': -95} # mV
cellRule['secs']['soma']['ions']['ca'] = { # allows for chane in reversal potential, internal concentration, and external concentration
        'o': 2.0} # mM

cellRule['globals']['krev'] = -95

netParams.cellParams['RE'] = cellRule

'''
TC Cell
'''
cellRule = {'conds': {'cellType': 'TC'}, 'secs': {}, 'globals': {}} 
cellRule['secs']['soma'] = {'geom': {}, 'mechs': {}, 'ions': {}} # template for soma characteristics 
cellRule['secs']['soma']['geom'] = { # geometry of TC soma 
        'diam': np.sqrt(2.9e4/np.pi), # one-compartment of 2.9e4 um2
        'L': np.sqrt(2.9e4/np.pi),
        'Ra': 100, 
        'cm': 1} # uF/cm2}
cellRule['secs']['soma']['mechs']['pas'] = { # leak current
        'e': -70, # mV
        'g': 1e-5} # S/cm2
cellRule['secs']['soma']['mechs']['kL'] = { # K+ leak current
        #'h.krev': -95, # global reversal potential, I don't think this is important, because the reversal potential is already altered in 
        'gkL': config.gkl_TC_baseline * config.init_gkl_TC} # S/cm2, factor of 0.79 applied to simulate awake state from config
cellRule['secs']['soma']['mechs']['naf_tc'] = { # fast-Na+ current 
        'gnabar': 0.090} # S/cm2 
cellRule['secs']['soma']['mechs']['kdr_tc'] = { # delayed-rectifying K+ current
        'gkbar': 0.012} # S/cm2 
cellRule['secs']['soma']['mechs']['it_tc'] = { # low threshold Ca++ current 
        'qm': 3.55,
        'qh': 3.0,
        'gcabar': 0.0025} # S/cm2 
cellRule['secs']['soma']['mechs']['cad'] = { # calcium decay
        'depth': 2.0, # um
        'taur': 5, # ms
        'cainf': 2.4e-4} # #mM
cellRule['secs']['soma']['mechs']['iar'] = { # anomalous rectifier channel
        'ghbar': 1.6e-5, # S/cm2  
        'fac_gh_TC': config.init_gh_TC} # mV; this should -8mV, -3mV, -2mV, and 0mV for awake, N2, N3, and REM sleep, respectively (Krishnan eLife 2016,

# Note: in order to change the reversal potential and current of ions, each ion must be called
cellRule['secs']['soma']['ions']['na'] = { # allows for a change in the reversal potential. internal concentration, and external concentration of sodium
        'e': 50} # mV 
cellRule['secs']['soma']['ions']['k'] = { # allows for change in the reversal potential, internal concentration, and external concentration of potassium
        'e': -95} # mV   
cellRule['secs']['soma']['ions']['h'] = { # anomalous rectifier channel 
        'e': -40} # mV
cellRule['secs']['soma']['ions']['ca'] = { # alloows for change in reversal potential, internal concentration, and external concentration of calcium
        'e': 120, # mV
        'i': 2.4e-4, # mM
        'o': 2.0} # mM

# This block functions as h. block in cell_classes (lines 447-454)
cellRule['globals']['cac'] = 0.0015 # half-activation of Ca++ binding; Krishnan currents.h line 376
cellRule['globals']['k2'] = 0.0004 # decay of Ca++ binding on protein (Krishnan currents.cpp line 167)
cellRule['globals']['Pc'] = 0.007 # half-activation of binding on Ih channel (Krishnan CellSyn.h line 234)
cellRule['globals']['k4'] = 0.001 # decay of protein binding on Ih channel (Krishnan CellSyn.h line 235)
cellRule['globals']['nca'] = 4 # "h." rather than "soma." because this is a GLOBAL variable in the mod file; nb of binding sites for Ca++ on protein (Krishnan currents.cpp line 168)
cellRule['globals']['nexp'] = 1 # nb of binding sites on Ih channel(Krishnan currents.cpp line 168)
cellRule['globals']['ginc'] = 2.0 # augm of conductance of bound Ih (Krishnan CellSyn.h line 229)
cellRule['globals']['taum'] = 20.0 # Krishnan currents.cpp line 168

cellRule['globals']['krev'] = -95

netParams.cellParams['TC'] = cellRule

##############################################################################
# Synaptic Parameters
##############################################################################
netParams.synMechParams['PY2IN_AMPA_D2'] = {'mod': 'AMPA_D2', 'gmax': syn_conns.pyr2inh_ampa_gmax, 'psp_weight': config.pyr2inh_ampa_d2_mini_str/config.pyr2inh_ampa_d2_str} # IN->PY AMPA_D2 synapse, Note: 'gid' and 'syn_index' have no weight here but is simply a place holder, also psp_weight = 0 for the initial conditions, but it evolves over time
netParams.synMechParams['PY2IN_NMDA_D1'] = {'mod': 'NMDA_D1', 'gmax': syn_conns.pyr2inh_nmda_gmax} # PY->IN NMDA_D1 synapse
netParams.synMechParams['PY2PY_AMPA_D2'] = {'mod': 'AMPA_D2', 'gmax': syn_conns.pyr2pyr_ampa_gmax, 'psp_weight': config.pyr2pyr2_ampa_d2_mini_str/config.pyr2pyr_ampa_d2_str} # PY->PY AMPA_D2 synapse, Note: 'gid' and 'syn_index' have no weight here but is simply a place holder, also psp_weight = 0 for the initial conditions, but it evolves over time
netParams.synMechParams['PY2PY_NMDA_D1'] = {'mod': 'NMDA_D1', 'gmax': syn_conns.pyr2pyr_nmda_gmax} # PY->PY NMDA_D1 synapse
netParams.synMechParams['IN2PY_GABA_A_D2'] = {'mod': 'GABA_A_D2', 'gmax': syn_conns.inh2pyr_gaba_a_gmax, 'psp_weight':config.inh2pyr_gaba_a_d2_mini_str/config.inh2pyr_gaba_a_d2_str} # IN->PY GABA_A_D2 synapse, Note: 'gid' and 'syn_index' have no weight here but is simply a place holder, also psp_weight = 0 for the initial conditions, but it evolves over time

netParams.synMechParams['RE2TC_GABA_A'] = {'mod': 'GABA_A', 'gmax': syn_conns.re2tc_gaba_a_gmax, 'Erev': -83.0} # RE->TC GABA_A synapse
netParams.synMechParams['RE2TC_GABA_B'] = {'mod': 'GABA_B', 'gmax': syn_conns.re2tc_gaba_b_gmax, 'Erev': -95.0} # RE->TC GABA_B synapse
netParams.synMechParams['RE2RE_GABA_A'] = {'mod': 'GABA_A', 'gmax': syn_conns.re2re_gaba_gmax} # RE->RE GABA_A synapse

netParams.synMechParams['TC2RE_AMPA'] = {'mod': 'AMPA', 'gmax': syn_conns.tc2re_ampa_gmax} # TC->RE AMPA synapse
netParams.synMechParams['TC2PY_AMPA'] = {'mod': 'AMPA', 'gmax': syn_conns.tc2pyr_ampa_gmax} # TC->PY AMPA synapse
netParams.synMechParams['TC2IN_AMPA'] = {'mod': 'AMPA', 'gmax': syn_conns.tc2inh_ampa_gmax} # TC->IN AMPA synapse

netParams.synMechParams['PY2TC_AMPA'] = {'mod': 'AMPA', 'gmax': syn_conns.pyr2tc_ampa_gmax} # PY->TC AMPA synapse
netParams.synMechParams['PY2RE_AMPA'] = {'mod': 'AMPA', 'gmax': syn_conns.pyr2re_ampa_gmax} # PY->RE AMPA synapse

##############################################################################
# Stimulation Parameters
##############################################################################
'''0.0417 uA/cm2 must be injected into pyramidal cell in order for its rheobase (and entire f-I curve) to match
that of Bazhenov's reduced model (see July 10, 2019 entry of bazhenov_sleep_model_journal.docx)'''
dur_PY = 1e9 # ms; duration of the stimulation 
delay_PY = 0 # ms; delay in the stimulation 
amp_PY = 0.0417 # uA/cm2; stimulus current 
netParams.stimSourceParams['IClamp_Pyr'] = {'type': 'IClamp',
                                        'del': delay_PY, # ms
                                        'dur': dur_PY, # creates an infinite duration; ms 
                                        'amp': (amp_PY * np.pi * netParams.cellParams['Pyr']['secs']['dend']['geom']['diam'] *netParams.cellParams['Pyr']['secs']['dend']['geom']['L'] *1e-5)} # the current in nA and the location of the stimulation on a section

'''0.1091 uA/cm2 must be injected into basket cell in order for its rheobase (and entire f-I curve) to match
that of Bazhenov's reduced model (see July 10, 2019 entry of bazhenov_sleep_model_journal.docx)'''
dur_IN = 1e9 # ms; duration of the stimulation 
delay_IN = 0 # ms; delay in the stimulation 
amp_IN = 0.1091 # uA/cm2; stimulus current 
netParams.stimSourceParams['IClamp_Inh'] = {'type': 'IClamp',
                                        'del': delay_IN, # ms
                                        'dur': dur_IN, # creates an infinite duration; ms 
                                        'amp': (amp_IN * np.pi * netParams.cellParams['Inh']['secs']['dend']['geom']['diam'] *netParams.cellParams['Inh']['secs']['dend']['geom']['L'] *1e-5)} # the current in nA and the location of the stimulation on a section

netParams.stimTargetParams['IClamp->Pyr'] = { # Inserts the 
        'source': 'IClamp_Pyr', # Calls the stimSource above
        'sec' : 'dend', 
        'loc': 0.5,
        'conds': {'cellType': 'Pyr'}} # insert the tag of cell in question, if targeting a single cell use cell id

netParams.stimTargetParams['IClamp->Inh'] = {
        'source': 'IClamp_Inh', # Calls the stimSource above
        'sec' : 'dend', 
        'loc': 0.5,
        'conds': {'cellType': 'Inh'}} # insert the tag of cell in question, if targeting a single cell use cell id 

##############################################################################
# Connectivity Rules (Same order as Synaptic Parameters)
##############################################################################
netParams.connParams['PY2IN_AMPA_D2'] = {
        'preConds': {'pop':'PyrCell'}, 
        'postConds': {'pop':'InhCell'},
        'synMech': 'PY2IN_AMPA_D2', # Extremely important to cell the syanpse here, as it allows for the weight, delay, loc, sec, and other conenctions parameters to be added to the synapse
        'weight': 1.0, # scales gmax
        'delay': 0.1, # ms
        'sec': 'dend', # attaches synapse to the location on the post-synaptic cell (default is soma(0.5))
        'loc': 0.5,
        'connList': syn_conns.pyr2inh_ampa_d2_conn # defines a list of connections from [pre-synaptic cell id, post-synaptic cell id], Note: The cell ids are relative to its population
        }

netParams.connParams['PY2IN_NMDA_D1'] = {
        'preConds': {'pop': 'PyrCell'},
        'postConds': {'pop': 'InhCell'},
        'synMech': 'PY2IN_NMDA_D1', # Extremely important to list the syanpse here, as it allows for the weight, delay, loc, sec, and other conenctions parameters to be added to the synapse
        'weight': 1.0, # scales gmax 
        'delay': 0.1, # ms
        'sec': 'dend', # attaches synapse to the location on the post-synaptic cell (default is soma(0.5))
        'loc': 0.5,
        'connList': syn_conns.pyr2inh_nmda_d1_conn # defines a list of connections from [pre-synaptic cell id, post-synaptic cell id], Note: The cell ids are relative to its population
        }

netParams.connParams['PY2PY_AMPA_D2'] = { 
        'preConds': {'pop': 'PyrCell'},
        'postConds': {'pop': 'PyrCell'},
        'synMech': 'PY2PY_AMPA_D2', # Extremely important to list the syanpse here, as it allows for the weight, delay, loc, sec, and other conenctions parameters to be added to the synapse
        'weight': 1.0, # scales gmax
        'delay': 0.1, # ms
        'sec': 'dend', # attaches synapse to the location on the post-synaptic cell (default is soma(0.5))
        'loc': 0.5,
        'connList': syn_conns.pyr2pyr_ampa_d2_conn # defines a list of connections from [pre-synaptic cell id, post-synaptic cell id], Note: The cell ids are relative to its population
        } 

netParams.connParams['PY2PY_NMDA_D1'] = { 
        'preConds': {'pop': 'PyrCell'},
        'postConds': {'pop': 'PyrCell'},
        'synMech': 'PY2PY_NMDA_D1', # Extremely important to list the syanpse here, as it allows for the weight, delay, loc, sec, and other conenctions parameters to be added to the synapse
        'weight': 1.0, # scales gmax
        'delay': 0.1, # ms
        'sec': 'dend', # attaches synapse to the location on the post-synaptic cell (default is soma(0.5))
        'loc': 0.5,
        'connList': syn_conns.pyr2pyr_nmda_d1_conn # defines a list of connections from [pre-synaptic cell id, post-synaptic cell id], Note: The cell ids are relative to its population
        }

netParams.connParams['IN2PY_GABA_A_D2'] = {
        'preConds': {'pop':'InhCell'},
        'postConds': {'pop':'PyrCell'},
        'synMech': 'IN2PY_GABA_A_D2', # Extremely important to list the syanpse here, as it allows for the weight, delay, loc, sec, and other conenctions parameters to be added to the synapse
        'weight': 1.0, # scales gmax
        'delay': 0.1, # ms
        'sec': 'dend', # attaches synapse to the location on the post-synaptic cell (default is soma(0.5))
        'loc': 0.5,
        'connList': syn_conns.inh2pyr_gaba_a_d2_conn # defines a list of connections from [pre-synaptic cell id, post-synaptic cell id], Note: The cell ids are relative to its population
        }

netParams.connParams['RE2TC_GABA_A'] = {
        'preConds': {'pop': 'RECell'},
        'postConds': {'pop': 'TCCell'},
        'synMech': 'RE2TC_GABA_A', # Extremely important to list the syanpse here, as it allows for the weight, delay, loc, sec, and other conenctions parameters to be added to the synapse
        'weight': 1.0, # scales gmax
        'delay': 0.1, # ms
        'sec': 'soma', # attaches synapse to the location on the post-synaptic cell (default is soma(0.5))
        'loc': 0.5, 
        'connList': syn_conns.re2tc_gaba_a_conn # defines a list of connections from [pre-synaptic cell id, post-synaptic cell id], Note: The cell ids are relative to its population
        }

netParams.connParams['RE2TC_GABA_B'] = {
        'preConds': {'pop': 'RECell'},
        'postConds': {'pop': 'TCCell'},
        'synMech': 'RE2TC_GABA_B', # Extremely important to list the syanpse here, as it allows for the weight, delay, loc, sec, and other conenctions parameters to be added to the synapse
        'weight': 1.0, # scales gmax
        'delay': 0.1, # ms
        'sec': 'soma', # attaches synapse to the location on the post-synaptic cell (default is soma(0.5))
        'loc': 0.5, 
        'connList': syn_conns.re2tc_gaba_b_conn # defines a list of connections from [pre-synaptic cell id, post-synaptic cell id], Note: The cell ids are relative to its population
        }

netParams.connParams['RE2RE_GABA_A'] = {
        'preConds': {'cellType': 'RE'}, 
        'postConds': {'cellType': 'RE'},
        'synMech': 'RE2RE_GABA_A', # Extremely important to list the syanpse here, as it allows for the weight, delay, loc, sec, and other conenctions parameters to be added to the synapse
        'weight': 1.0, # scales the gmax
        'delay': 0.1, # ms 
        'sec': 'soma', # attaches synapse to the location on the post-synaptic cell (default is soma(0.5))
        'loc': 0.5,
        'connList': syn_conns.re2re_gaba_a_conn # defines a list of connections from [pre-synaptic cell id, post-synaptic cell id], Note: The cell ids are relative to its population
        }

netParams.connParams['TC2RE_AMPA'] = {
        'preConds': {'pop': 'TCCell'},
        'postConds': {'pop': 'RECell'},
        'synMech': 'TC2RE_AMPA', # Extremely important to cell the syanpse here, as it allows for the weight, delay, loc, sec, and other conenctions parameters to be added to the synapse
        'weight': 1.0, # scales gmax
        'delay': 0.1, # ms
        'sec': 'soma', # attaches synapse to the location on the post-synaptic cell (default is soma(0.5))
        'loc': 0.5,
        'connList': syn_conns.tc2re_ampa_conn # defines a list of connections from [pre-synaptic cell id, post-synaptic cell id], Note: The cell ids are relative to its population
        }

netParams.connParams['TC2PY_AMPA'] = {
        'preConds': {'pop': 'TCCell'},
        'postConds': {'pop': 'PyrCell'},
        'synMech': 'TC2PY_AMPA', # Extremely important to cell the syanpse here, as it allows for the weight, delay, loc, sec, and other conenctions parameters to be added to the synapse
        'weight': 1.0, # scales gmax 
        'delay': 0.1, # ms
        'sec': 'dend', # attaches synapse to the location on the post-synaptic cell (default is soma(0.5))
        'loc': 0.5,
        'connList': syn_conns.tc2pyr_ampa_conn # defines a list of connections from [pre-synaptic cell id, post-synaptic cell id], Note: The cell ids are relative to its population
        }

netParams.connParams['TC2IN_AMPA'] = {
        'preConds': {'pop':'TCCell'},
        'postConds': {'pop':'InhCell'},
        'synMech':'TC2IN_AMPA', # Extremely important to cell the syanpse here, as it allows for the weight, delay, loc, sec, and other conenctions parameters to be added to the synapse
        'weight': 1.0, # scales gmax
        'delay': 0.1, # ms
        'sec':'dend', # attaches synapse to the location on the post-synaptic cell (default is soma(0.5))
        'loc': 0.5,
        'connList': syn_conns.tc2inh_ampa_conn # defines a list of connections from [pre-synaptic cell id, post-synaptic cell id], Note: The cell ids are relative to its population
        }

netParams.connParams['PY2TC_AMPA'] = {
        'preConds': {'pop': 'PyrCell'},
        'postConds': {'pop': 'TCCell'},
        'synMech': 'PY2TC_AMPA', # Extremely important to cell the syanpse here, as it allows for the weight, delay, loc, sec, and other conenctions parameters to be added to the synapse
        'weight': 1.0, # scales gmax 
        'delay': 0.1, # ms
        'sec': 'soma',
        'loc': 0.5,
        'connList': syn_conns.pyr2tc_ampa_conn # defines a list of connections from [pre-synaptic cell id, post-synaptic cell id], Note: The cell ids are relative to its population
        }

netParams.connParams['PY2RE_AMPA'] = {
        'preConds': {'pop': 'PyrCell'},
        'postConds': {'pop': 'RECell'},
        'synMech': 'PY2RE_AMPA', # Extremely important to cell the syanpse here, as it allows for the weight, delay, loc, sec, and other conenctions parameters to be added to the synapse
        'weight': 1.0, # scales gmax
        'delay': 0.1, # ms
        'sec': 'soma', # attaches synapse to the location on the post-synaptic cell (default is soma(0.5))
        'loc': 0.5,
        'connList': syn_conns.pyr2re_ampa_conn # defines a list of connections from [pre-synaptic cell id, post-synaptic cell id], Note: The cell ids are relative to its population
        }

##############################################################################
# Simulation Configuration Parameters 
##############################################################################
simConfig.recordStep = 0.025 # ms, time interval to pickle data, lose resolution with default recordStep = 0.1 (see simConfig.py by ctrl + left click)
simConfig.verbose = False  # Helpful in NSG runs 
# Sets up the simulation parameters
simConfig.duration = config.duration # ms
simConfig.dt = 0.025 # ms
#simConfig.analysis['plotConn'] = {'feature': 'numConns', 'groupBy': 'cell', 'synMech': 'PY2TC_AMPA', 'saveData': 'C:/Users/logan/Desktop/Research2020/bazh-model/save_images/plotConn_data1.pkl', 'showFig': True} 
simConfig.analysis['plotRaster'] = {'maxSpikes': 1e10, 'marker': 'o', 'saveData': raster_filename, 'showFig': True} # Plots the raster slows down simulation

##############################################################################
# Simulation Parameters 
##############################################################################
sim.create(netParams, simConfig) # Instantiates the network without running the simulation

for cell_id in range(0,NPyr+NInh,1): # Inserts h.ion_style into all Pyr and Inh cells 
    h.ion_style("ca_ion",3,1,0,0,0,sec=sim.net.cells[cell_id].secs['dend']['hObj']) # this ensures that cai is treated as a dynamical variable, while eca is a fixed parameter (see Krishnan currents.cpp line 205)
for cell_id in range(NPyr+NInh, NPyr+NInh+NRE+NTC, 1): # Inserts h.ion_style into all RE and TC cells
    h.ion_style("ca_ion",3,1,0,0,0,sec=sim.net.cells[cell_id].secs['soma']['hObj']) # this ensures that cai is treated as a dynamical variable, while eca is a fixed parameter (see Krishnan currents.cpp line 205)

for pyr_id in range(0,NPyr,1):
    sim.net.cells[pyr_id].secs['dend']['synMechs'][0]['hObj'].setrand(pyr_id,0) # Calls the setrand function in AMPA_D2 mod file, see the following post for more information: https://groups.google.com/forum/?nomobile=true#!topic/netpyne-forum/goaTVTNvc2M
    h.SS_denom_AMPA_D2 = 250.0 # Krishnan uses 250.0 ms (need to use "h.*" bc. this is a GLOBAL variable in mod file); this parameter controls the mean inter-event interval of stochastic EPSP's (so smaller value will lead to more frequent stochastic EPSP's)
    h.mini_fre_AMPA_D2 = config.pyr2pyr2_ampa_d2_mini_f # Krishnan uses 20.0 ms
    if config.pyr2pyr_ampa_d2_str > 0: # prevent division by zero
        sim.net.cells[pyr_id].secs['dend']['synMechs'][0]['psp_weight'] = config.pyr2pyr2_ampa_d2_mini_str/config.pyr2pyr_ampa_d2_str # Chris Fink wrote the mod file so that 'psp_weight' is a normalized value that tells you how much stronger (for values greater than 1) or weaker (values less than 1) the stochastic stimulation is than the regular presynaptic weight
    else:
        sim.net.cells[pyr_id].secs['dend']['synMechs'][0]['psp_weight'] = 0 # if regular pyr2pyr_ampa_d2 connections have strength of zero, then stochastic stimulation of these synapses will also be absent    
    # NMDA_D1 PY2PY synapse usage
    h.Use_NMDA_D1 = config.pyr2pyr_nmda_d1_Use #depression is actually turned off, as far as I can tell, in Krishnan's model (which implies Use=0)
    # PY->PY AMPA_D2 gmax list
    sim.net.cells[pyr_id].secs['dend']['synMechs'][0]['gmax'] = syn_conns.pyr2pyr_ampa_d2_gmax[pyr_id]
    # PY->PY NMDA_D1 gmax list
    sim.net.cells[pyr_id].secs['dend']['synMechs'][1]['gmax'] = syn_conns.pyr2pyr_nmda_d1_gmax[pyr_id]
    # TC->PY AMPA gmax list 
    sim.net.cells[pyr_id].secs['dend']['synMechs'][3]['gmax'] = syn_conns.tc2pyr_ampa_gmax_wghts[pyr_id]
     # Randomized g_pas
    r = h.Random()
    r.Random123(pyr_id,2,0) #set stream of random numbers; first argument is gid, make second argument 2 because 0 and 1 are already taken by AMPA_D2 and GABA_D2 synapses (see cell_classes.py)
    sim.net.cells[pyr_id].secs['dend']['mechs']['pas']['g'] = 0.000011 + (r.uniform(0,1)+1) * 0.000003 #add cell-to-cell variability in this parameter, as prescribed in Krishnan CellSyn.h line 365
    
    
for pyr_id_GABA_A_D2 in range(0,NPyr,1):
    sim.net.cells[pyr_id_GABA_A_D2].secs['dend']['synMechs'][2]['hObj'].setrand(pyr_id_GABA_A_D2,1) # Calls the setrand function in GABA_A_D2 mod file, see the following post for more information: https://groups.google.com/forum/?nomobile=true#!topic/netpyne-forum/goaTVTNvc2M
    h.SS_denom_GABA_A_D2 = 250.0 # Krishnan uses 250.0 ms (need to use "h.*" bc. this is a GLOBAL variable in mod file); see note on SS_denom_AMPA_D2 above
    h.mini_fre_GABA_A_D2 = config.inh2pyr_gaba_a_d2_mini_f # Krishnan uses 20.0 ms
    if config.inh2pyr_gaba_a_d2_str > 0: # prevent division by zero
        sim.net.cells[pyr_id_GABA_A_D2].secs['dend']['synMechs'][2]['psp_weight'] = config.inh2pyr_gaba_a_d2_mini_str/config.inh2pyr_gaba_a_d2_str # see note for 'psp_weight' above; I'm (Chris Fink) not sure what value this should be...
    else: 
        sim.net.cells[pyr_id_GABA_A_D2].secs['dend']['synMechs'][2]['psp_weight'] = 0 # if regular inh2pyr_gaba_a_d2 connections have zero strength, then stochastic stimulation of these synapses will also be absent
    # IN->PY GABA_D2 gmax list
    sim.net.cells[pyr_id_GABA_A_D2].secs['dend']['synMechs'][2]['gmax'] = syn_conns.inh2pyr_gaba_a_d2_gmax[pyr_id_GABA_A_D2]
    
for inh_id in range(NPyr,NPyr+NInh,1):
    sim.net.cells[inh_id].secs['dend']['synMechs'][0]['hObj'].setrand(inh_id,0) # Calls the setrand function in ANPA_D2 mod file, see the following post for more information: https://groups.google.com/forum/?nomobile=true#!topic/netpyne-forum/goaTVTNvc2M
    h.SS_denom_AMPA_D2 = 250.0 # Krishnan uses 250.0 ms (need to use "h.*" bc. this is a GLOBAL variable in mod file)
    h.mini_fre_AMPA_D2 = config.pyr2inh_ampa_d2_mini_f # Krishnan uses 20.0 ms
    if config.pyr2inh_ampa_d2_str > 0: # prevent division by zero
       sim.net.cells[inh_id].secs['dend']['synMechs'][0]['psp_weight'] = config.pyr2inh_ampa_d2_mini_str/config.pyr2inh_ampa_d2_str # Krishnan uses 0.20 uS for stochastic stimulation strength, and 0.12 uS for PY->INH syn strength; see note on 'psp_weight' for pyramidal cell AMPA_D2 synapse
    else: 
       sim.net.cells[inh_id].secs['dend']['synMechs'][0]['psp_weight'] = 0 # if regular pyr2inh_ampa_d2 connections have strength of zero, then stochastic stimulation of these synapses will also be absent
    # NMDA_D1 PY2IN synapse usage
    h.Use_NMDA_D1 = config.pyr2inh_nmda_d1_Use
    # PY->IN AMOA_D2 gmax list 
    sim.net.cells[inh_id].secs['dend']['synMechs'][0]['gmax'] = syn_conns.pyr2inh_ampa_d2_gmax[inh_id-NPyr]
    # PY->IN NMDA_D1 gmax list 
    sim.net.cells[inh_id].secs['dend']['synMechs'][1]['gmax'] = syn_conns.pyr2inh_nmda_d1_gmax[inh_id-NPyr]
    # TC->IN AMPA gmax list 
    sim.net.cells[inh_id].secs['dend']['synMechs'][2]['gmax'] = syn_conns.tc2inh_ampa_gmax_wghts[inh_id-NPyr]
    # Randomized g_pas
    r = h.Random()
    r.Random123(inh_id,1,0) #set stream of random numbers; first argument is gid, make second argument 2 because 0 and 1 are already taken by AMPA_D2 and GABA_D2 synapses (see cell_classes.py)
    sim.net.cells[inh_id].secs['dend']['mechs']['pas']['g'] = 0.000011 + (r.uniform(0,1)+1) * 0.000003 #add cell-to-cell variability in this parameter, as prescribed in Krishnan CellSyn.h line 365
    
    
for re_id in range(NPyr+NInh,NPyr+NInh+NRE,1):
    # RE->RE GABA_A gmax list
    sim.net.cells[re_id].secs['soma']['synMechs'][0]['gmax'] = syn_conns.re2re_gaba_a_gmax[re_id-NPyr-NInh]
    # TC->RE AMPA gmax list
    sim.net.cells[re_id].secs['soma']['synMechs'][1]['gmax'] = syn_conns.tc2re_ampa_gmax_wghts[re_id-NPyr-NInh]
    # PY->RE AMPA gmax list
    sim.net.cells[re_id].secs['soma']['synMechs'][2]['gmax'] = syn_conns.pyr2re_ampa_gmax_wghts[re_id-NPyr-NInh]
    
for tc_id in range(NPyr+NInh+NRE,NPyr+NInh+NRE+NTC,1):
    # RE->TC GABA_A gmax list 
    sim.net.cells[tc_id].secs['soma']['synMechs'][0]['gmax'] = syn_conns.re2tc_gaba_a_gmax_wghts[tc_id-NPyr-NInh-NRE]
    # RE->TC GABA_B gmax list 
    sim.net.cells[tc_id].secs['soma']['synMechs'][1]['gmax'] = syn_conns.re2tc_gaba_b_gmax_wghts[tc_id-NPyr-NInh-NRE]
    # PY->TC AMPA gmax list
    sim.net.cells[tc_id].secs['soma']['synMechs'][2]['gmax'] = syn_conns.pyr2tc_ampa_gmax_wghts[tc_id-NPyr-NInh-NRE]

sim.simulate()
sim.analyze()

