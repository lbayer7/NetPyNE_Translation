# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 09:22:33 2020

This file will aim to generate the necessary synaptic connections for all the cells in the network.

@author: logan
"""

import config
import numpy as np 

##### specify PYR->INH AMPA_D2 connections 
rad=config.pyr2inh_ampa_d2_rad #number of outgoing connections from each PY cell 
pyr2inh_ampa_d2_conn = [] #generates empty list for the synaptic connections of PY->IN_AMPA_D2, Note: This will be a list within a list as per NetPyNE documentation.
pyr2inh_ampa_d2_gmax = [] #generates empty list for the weight of the synaptic conductance 
p2i_ratio=config.Npyr/config.Ninh #we will assume that Npyr > Ninh; this ratio will be important for determining the center of the set of source PYR cells for each post-synaptic INH cell
for inh_gid in range(config.Npyr, config.Npyr+config.Ninh,1): #loop over all post-synaptic INH cells on this machine (have to loop over post-syn cells bc. pc.gid_connect only works with post-synaptic gid on machine)
    i_inh = inh_gid - config.Npyr #easier to work by indexing the INH cells from 0 to Ninh-1, rather than from Npyr to Npyr+Ninh-1
    tmp = list(range( int(np.floor((i_inh-rad)*p2i_ratio)), int(np.floor((i_inh+rad)*p2i_ratio)) +1) ) #center_pyr should be roughly round(i_inh*p2i_ratio); then you need to consider the radius of connectivity surrounding that, and consider that larger p2i_ratio increases the number of PYR cells sending to a given INH cell (for a particular radius of connectivity); "+1" is due to definition of Python's "range" function 
    pyr_set = [val % config.Npyr for val in tmp] #apply modulus operation to eliminate negative numbers; this implements periodic boundary conditions; indexes PYR cells from 0 to Npyr-1
    k_PY_IN_AMPA = 0 #make sure in-degree is initialized to zero
    for i_pyr in pyr_set:
        pyr2inh_ampa_d2 = [i_pyr,i_inh] #list is read in NetPyNE as follows: [presynaptic cell id, postsynaptic cell id], Note: See http://www.netpyne.org/reference.html#network-parameters under 'Connectivity rules -- connList' for description of ids
        pyr2inh_ampa_d2_conn.append(pyr2inh_ampa_d2)
        k_PY_IN_AMPA += 1 
    #reduce gmax for this cell's synapse, so that the total synaptic strength is equal to that specified in the config file
    #(gmax was already set in cell_classes.py, which used the value specified in config.py)
    #note that this should also normalize the mini's, since their strength is set in the mod file by psp_weight, which is set in cell_classes.py
    #for fac_AMPA_D2, see Krishnan's currents.cpp line 545 and main.cpp lines 579-583, 639-776
    pyr2inh_ampa_gmax = config.init_AMPA_D2 * config.pyr2inh_ampa_d2_str / k_PY_IN_AMPA 
    pyr2inh_ampa_d2_gmax.append(pyr2inh_ampa_gmax)

##### specify PYR->INH NMDA_D1 connections
rad=config.pyr2inh_nmda_d1_rad #number of outgoing connnections from each PY cell
pyr2inh_nmda_d1_conn = [] #generates empty list for the synaptic connections of PY->IN_NMDA_D1, Note: This will be a list within a list as per NetPyNE documentation.
pyr2inh_nmda_d1_gmax = [] #generates empty list for the weight of the synaptic conductance
p2i_ratio = config.Npyr/config.Ninh #we will assume that Npyr > Ninh; this ratio will be important for determining the center of the set of source PYR cells for each post-synaptic INH cell
for inh_gid in range(config.Npyr, config.Npyr+config.Ninh,1): #loop over all post-synaptic INH cells on this machine (have to loop over post-syn cells bc. pc.gid_connect only works with post-synaptic gid on machine)
    i_inh = inh_gid - config.Npyr #easier to work by indexing the INH cells from 0 to Ninh-1, rather than from Npyr to Npyr+Ninh-1
    tmp = list(range( int(np.floor((i_inh-rad)*p2i_ratio)), int(np.floor((i_inh+rad)*p2i_ratio)) +1) ) #center_pyr should be roughly round(i_inh*p2i_ratio); then you need to consider the radius of connectivity surrounding that, and consider that larger p2i_ratio increases the number of PYR cells sending to a given INH cell (for a particular radius of connectivity); "+1" is due to definition of Python's "range" function 
    pyr_set = [val % config.Npyr for val in tmp] #apply modulus operation to eliminate negative numbers; this implements periodic boundary conditions; indexes PYR cells from 0 to Npyr-1
    k_PY_IN_NMDA = 0 #make sure in-degree is initialized to zero
    for i_pyr in pyr_set:
        k_PY_IN_NMDA += 1
        pyr2inh_nmda_d1 = [i_pyr,i_inh] #list is read in NetPyNE as follows: [presynaptic cell id, postsynaptic cell id], Note: See http://www.netpyne.org/reference.html#network-parameters under 'Connectivity rules -- connList' for description of ids
        pyr2inh_nmda_d1_conn.append(pyr2inh_nmda_d1)
    #reduce gmax for this cell's synapse, so that the total synaptic strength is equal to that specified in the config file
    #(gmax was already set in cell_classes.py, which used the value specified in config.py)
    pyr2inh_nmda_gmax = config.pyr2inh_nmda_d1_str / k_PY_IN_NMDA
    pyr2inh_nmda_d1_gmax.append(pyr2inh_nmda_gmax)

##### specify PYR->PYR AMPA_D2 connections
rad=config.pyr2pyr_ampa_d2_rad #number of outgoing connections from each PY cell
pyr2pyr_ampa_d2_conn = [] #generates empty list for the synaptic connections of PY->PY_AMPA_D2, Note: This will be a list within a list as per NetPyNE documentation.
pyr2pyr_ampa_d2_gmax = [] #generates empty list for the weight of the synaptic conductance
for pyr_gid in range(0, config.Npyr, 1): #loop over all post-synaptic PYR cells on this machine (have to loop over post-syn cells bc. pc.gid_connect only works with post-synaptic gid on machine)
    i_post_pyr = pyr_gid #no subtraction needed, because pyramidal cells are the first set of gid's
    tmp = list(range(i_post_pyr-rad,i_post_pyr+rad+1)) #generate list of PYR sources (but this will in general include negative values)
    pyr_pre_set = [val % config.Npyr for val in tmp] #apply modulus operation to eliminate negative numbers; this implements periodic boundary conditions; indexes PYR cells from 0 to Nre-1
    k_PY_PY_AMPA = 0 #make sure in-degree is initialized to zero 
    for i_pre_pyr in pyr_pre_set:
        if i_pre_pyr != i_post_pyr: #prevent self connections
            pyr2pyr_ampa_d2 = [i_pre_pyr,i_post_pyr] #list is read in NetPyNE as follows: [presynaptic cell id, postsynaptic cell id], Note: See http://www.netpyne.org/reference.html#network-parameters under 'Connectivity rules -- connList' for description of ids
            pyr2pyr_ampa_d2_conn.append(pyr2pyr_ampa_d2)
            k_PY_PY_AMPA += 1
    #reduce gmax for this cell's synapse, so that the total synaptic strength is equal to that specified in the config file
    #(gmax was already set in cell_classes.py, which used the value specified in config.py)
    #note that this should also normalize the mini's, since their strength is set in the mod file by psp_weight, which is set in cell_classes.py
    #for fac_AMPA_D2, see Krishnan's currents.cpp line 545 and main.cpp lines 579-583, 639-776
    pyr2pyr_ampa_gmax = config.init_AMPA_D2 * config.pyr2pyr_ampa_d2_str / k_PY_PY_AMPA
    pyr2pyr_ampa_d2_gmax.append(pyr2pyr_ampa_gmax)

##### specify PYR->PYR NMDA_D1 connections   
rad=config.pyr2pyr_nmda_d1_rad #number of outgoing connections from each PY cell
pyr2pyr_nmda_d1_conn = [] #generates empty list for the synaptic connections of PY->PY_NMDA_D1, Note: This will be a list within a list as per NetPyNE documentation.
pyr2pyr_nmda_d1_gmax = [] #generates empty list for the weight of the synaptic conductance
for pyr_gid in range(0, config.Npyr, 1): #loop over all post-synaptic PYR cells on this machine (have to loop over post-syn cells bc. pc.gid_connect only works with post-synaptic gid on machine)
    i_post_pyr = pyr_gid #no subtraction needed, because pyramidal cells are the first set of gid's
    tmp = list(range(i_post_pyr-rad,i_post_pyr+rad+1)) #generate list of PYR sources (but this will in general include negative values)
    pyr_pre_set = [val % config.Npyr for val in tmp] #apply modulus operation to eliminate negative numbers; this implements periodic boundary conditions; indexes PYR cells from 0 to Nre-1
    k_PY_PY_NMDA = 0 #make sure in-degree is initialized to zero 
    for i_pre_pyr in pyr_pre_set:
        if i_pre_pyr != i_post_pyr: #prevent self connections
            pyr2pyr_nmda_d1 = [i_pre_pyr,pyr_gid] #list is read in NetPyNE as follows: [presynaptic cell id, postsynaptic cell id], Note: See http://www.netpyne.org/reference.html#network-parameters under 'Connectivity rules -- connList' for description of ids
            pyr2pyr_nmda_d1_conn.append(pyr2pyr_nmda_d1)
            k_PY_PY_NMDA += 1
    #reduce gmax for this cell's synapse, so that the total synaptic strength is equal to that specified in the config file
    #(gmax was already set in cell_classes.py, which used the value specified in config.py)
    pyr2pyr_nmda_gmax = config.pyr2pyr_nmda_d1_str / k_PY_PY_NMDA
    pyr2pyr_nmda_d1_gmax.append(pyr2pyr_nmda_gmax)
    
##### specify INH->PYR GABA_A_D2 connections 
rad=config.inh2pyr_gaba_a_d2_rad #number of outgoing connections from each INH cell
inh2pyr_gaba_a_d2_conn = [] #generates empty list for the synaptic connections of IN->PY GABA_A_D2, Note: This will be a list within a list as per NetPyNE documentation.
inh2pyr_gaba_a_d2_gmax = [] #generates empty list for the weight of the synaptic conductance
p2i_ratio = config.Npyr/config.Ninh #we will assume that Npyr > Ninh; this ratio will be important for determining the center of the set of source INH cells for each post-synaptic PYR cell
for pyr_gid in range(0, config.Npyr, 1): #loop over all post-synaptic PYR cells on this machine (have to loop over post-syn cells bc. pc.gid_connect only works with post-synaptic gid on machine)
    i_pyr = pyr_gid #no subtraction needed, because pyramidal cells are the first set of gid's
    tmp = list(range( int(np.floor((i_pyr-rad)/p2i_ratio)) , int(np.floor((i_pyr+rad)/p2i_ratio))+1 ) )  #center_inh should be roughly round(i_pyr/p2i_ratio); then you need to consider the radius of connectivity surrounding that, and consider that larger p2i_ratio decreases the number of INH cells sending to a given PYR cell (for a particular radius of connectivity); "+1" is due to definition of Python's "range" function 
    inh_set = [val % config.Ninh for val in tmp] #apply modulus operation to eliminate negative numbers; this implements periodic boundary conditions; indexes INH cells from 0 to Ninh-1
    k_IN_PY = 0 #make sure in-degree is initialized to zero 
    for i_inh in inh_set:
        inh2pyr_gaba_a_d2 = [i_inh,pyr_gid]
        inh2pyr_gaba_a_d2_conn.append(inh2pyr_gaba_a_d2)   
        k_IN_PY += 1
    #reduce gmax for this cell's synapse, so that the total synaptic strength is equal to that specified in the config file
    #(gmax was already set in cell_classes.py, which used the value specified in config.py)
    #note that this should also normalize the mini's, since their strength is set in the mod file by psp_weight, which is set in cell_classes.py
    #for fac_GABA_D2, see Krishnan's currents.cpp line 755, and main.cpp
    inh2pyr_gaba_a_gmax = config.init_GABA_D2 * config.inh2pyr_gaba_a_d2_str / k_IN_PY
    inh2pyr_gaba_a_d2_gmax.append(inh2pyr_gaba_a_gmax)

##### specify RE->TC GABA_A connections (NOTE: It is assumed Nre=Ntc) 
rad=config.re2tc_gaba_a_rad #number of outgoing connections from each RE cell
re2tc_gaba_a_conn = [] #generates empty list for the synaptic connections of RE->TC GABA_A, Note: This will be a list within a list as per NetPyNE documentation.
re2tc_gaba_a_gmax_wghts = [] #generates empty list for the weight of the synaptic conductance
for tc_gid in range(config.Npyr+config.Ninh+config.Nre, config.Npyr+config.Ninh+config.Nre+config.Ntc,1): #loop over all post-synaptic TC cells on this machine (have to loop over post-syn cells bc. pc.gid_connect only works with post-synaptic gid on machine)
    i_tc = tc_gid-(config.Npyr+config.Ninh+config.Nre) #easier to work by indexing the TC cells from 0 to Ntc-1, rather than from Npyr+Ninh+Nre to Npyr+Ninh+Nre+Ntc-1
    tmp = list(range(i_tc-rad,i_tc+rad+1)) #generate list of RE sources (but this will in general include negative values); note that this assumes Nre=Ntc
    re_set=[val % config.Nre for val in tmp] #apply modulus operation to eliminate negative numbers; this implements periodic boundary conditions; indexes RE cells from 0 to Nre-1
    k_RE_TC_GABA_A = 0 #make sure in-degree is initialized to zero
    #connect each RE cell to each target TC cell
    for i_re in re_set:
        re2tc_gaba_a = [i_re,i_tc]
        re2tc_gaba_a_conn.append(re2tc_gaba_a)
        k_RE_TC_GABA_A += 1 #update in-degree for this cell
    #reduce gmax for this cell's synapse, so that the total synaptic strength is equal to that specified in the config file
    #(gmax was already set in cell_classes.py, which used the value specified in config.py)
    #for fac_GABA_TC, see currents.cpp lines 341 & 369, and main.cpp
    re2tc_gaba_a_gmax = config.init_GABA_TC * config.re2tc_gaba_a_str / k_RE_TC_GABA_A
    re2tc_gaba_a_gmax_wghts.append(re2tc_gaba_a_gmax)
    #Note: I currently have the lists commented out as they are not used in the current version due to the single stage testing 
    #of this model
    
##### specify RE->TC GABA_B connections (NOTE: It is assumed Nre=Ntc) 
rad=config.re2tc_gaba_b_rad #number of outgoing connections from each RE cell
re2tc_gaba_b_conn = [] #generates empty list for the synaptic connections of RE->TC GABA_B, Note: This will be a list within a list as per NetPyNE documentation.
re2tc_gaba_b_gmax_wghts = [] #generates empty list for the weight of the synaptic conductance
for tc_gid in range(config.Npyr+config.Ninh+config.Nre, config.Npyr+config.Ninh+config.Nre+config.Ntc,1): #loop over all post-synaptic TC cells on this machine (have to loop over post-syn cells bc. pc.gid_connect only works with post-synaptic gid on machine)
    i_tc = tc_gid-(config.Npyr+config.Ninh+config.Nre) #easier to work by indexing the TC cells from 0 to Ntc-1, rather than from Npyr+Ninh+Nre to Npyr+Ninh+Nre+Ntc-1
    tmp = list(range(i_tc-rad,i_tc+rad+1)) #generate list of RE sources (but this will in general include negative values); note that this assumes Nre=Ntc
    re_set=[val % config.Nre for val in tmp] #apply modulus operation to eliminate negative numbers; this implements periodic boundary conditions; indexes RE cells from 0 to Nre-1
    k_RE_TC_GABA_B = 0 #make sure in-degree is initialized to zero
    for i_re in re_set:  #must add Npyr+Ninh to i_re in order to get RE cell's gid
        re2tc_gaba_b = [i_re,i_tc]
        re2tc_gaba_b_conn.append(re2tc_gaba_b)
        k_RE_TC_GABA_B += 1 #update in-degree for this cell
    #see note above for GABA_synapse
    #for fac_GABA_TC, see currents.cpp lines 341 & 369, and main.cpp
    re2tc_gaba_b_gmax = config.init_GABA_TC * config.re2tc_gaba_b_str / k_RE_TC_GABA_B
    re2tc_gaba_b_gmax_wghts.append(re2tc_gaba_b_gmax)
    
##### specify RE->RE GABA_A connections
rad=config.re2re_gaba_a_rad #number of outgoing connections from each RE cell
re2re_gaba_a_conn = [] #generates empty list for the synaptic connections of RE->RE GABA_A, Note: This will be a list within a list as per NetPyNE documentation.
re2re_gaba_a_gmax = [] #generates empty list for the weight of the synaptic conductance
for re_gid in range(config.Npyr+config.Ninh, config.Npyr+config.Ninh+config.Nre,1): #loop over all post-synaptic RE cells on this machine (have to loop over post-syn cells bc. pc.gid_connect only works with post-synaptic gid on machine)
    i_post_re = re_gid - (config.Npyr+config.Ninh) #easier to work by indexing the RE cells from 0 to Nre-1, rather than from Npyr+Ninh to Npyr+Ninh+Nre-1
    tmp = list(range(i_post_re-rad,i_post_re+rad+1)) #generate list of RE sources (but this will in general include negative values)
    re_pre_set=[val % config.Nre for val in tmp] #apply modulus operation to eliminate negative numbers; this implements periodic boundary conditions; indexes RE cells from 0 to Nre-1
    k_RE_RE = 0 #make sure in-degree is initialized to zero
    for i_pre_re in re_pre_set:
        if i_pre_re != i_post_re: #prevent self-connection
            re2re_gaba_a = [i_pre_re,i_post_re]
            re2re_gaba_a_conn.append(re2re_gaba_a)
            k_RE_RE += 1 #update in-degree for this cell
    #reduce gmax for this cell's synapse, so that the total synaptic strength is equal to that specified in the config file
    #(gmax was already set in cell_classes.py, which used the value specified in config.py)
    #for fac_GABA_TC, see currents.cpp lines 341 & 369, and main.cpp
    re2re_gaba_gmax = config.init_GABA_TC * config.re2re_gaba_a_str / k_RE_RE
    re2re_gaba_a_gmax.append(re2re_gaba_gmax)
    
##### specify TC->RE AMPA connections (NOTE: It is assumed Nre=Ntc)  
rad = config.tc2re_ampa_rad    #number of outgoing connections from each TC cell
tc2re_ampa_conn = [] #generates empty list for the synaptic connections of RE->TC AMPA, Note: This will be a list within a list as per NetPyNE documentation.
tc2re_ampa_gmax_wghts = [] #generates empty list for the weight of the synaptic conductance
for re_gid in range(config.Npyr+config.Ninh, config.Npyr+config.Ninh+config.Nre,1): #loop over all post-synaptic RE cells on this machine (have to loop over post-syn cells bc. pc.gid_connect only works with post-synaptic gid on machine)
    i_re = re_gid - (config.Npyr+config.Ninh) #easier to work by indexing the RE cells from 0 to Nre-1, rather than from Npyr+Ninh to Npyr+Ninh+Nre-1    
    tmp = list(range(i_re-rad,i_re+rad+1)) #generate list of TC sources (but this will in general include negative values)
    tc_set=[val % config.Ntc for val in tmp] #apply modulus operation to eliminate negative numbers; this implements periodic boundary conditions; indexes TC cells from 0 to Ntc-1
    k_TC_RE = 0 #make sure in-degree is initialized to zero
    for i_tc in tc_set: 
        tc2re_ampa = [i_tc,i_re]
        tc2re_ampa_conn.append(tc2re_ampa)
        k_TC_RE += 1
    #reduce gmax for this cell's synapse, so that the total synaptic strength is equal to that specified in the config file
    #(gmax was already set in cell_classes.py, which used the value specified in config.py)
    #for fac_AMPA_TC, see Krishnan's currents.cpp line 429, and main.cpp 
    tc2re_ampa_gmax = config.init_AMPA_TC * config.tc2re_ampa_str / k_TC_RE
    tc2re_ampa_gmax_wghts.append(tc2re_ampa_gmax)
    
##### specify TC->PYR AMPA connections
rad=config.tc2pyr_ampa_rad #number of outgoing connections from each TC cell
tc2pyr_ampa_conn = [] #generates empty list for the synaptic connections of TC->PY AMPA, Note: This will be a list within a list as per NetPyNE documentation.
tc2pyr_ampa_gmax_wghts = [] #generates empty list for the weight of the synaptic conductance 
p2t_ratio = config.Npyr/config.Ntc #we will assume that Npyr > Ntc; this ratio will be important for determining the center of the set of source TC cells for each post-synaptic PYR cell
for pyr_gid in range(0, config.Npyr, 1): #loop over all post-synaptic PYR cells on this machine (have to loop over post-syn cells bc. pc.gid_connect only works with post-synaptic gid on machine)
    i_pyr = pyr_gid #no subtraction needed, because PYR cells are the first set of gid's
    tmp = list(range(  int(np.floor((i_pyr-rad)/p2t_ratio)), int(np.floor((i_pyr+rad)/p2t_ratio))+1))  #center_tc should be roughly round(i_pyr/p2t_ratio); then you need to consider the radius of connectivity surrounding that, and consider that larger p2t_ratio decreases the number of TC cells sending to a given PYR cell (for a particular radius of connectivity); "+1" is due to definition of Python's "range" function 
    tc_set = [val % config.Ntc for val in tmp] #apply modulus operation to eliminate negative numbers; this implements periodic boundary conditions; indexes TC cells from 0 to Ntc-1
    k_TC_PY = 0 #make sure in-degree is initialized to zero
    for i_tc in tc_set:
        tc2pyr_ampa = [i_tc,i_pyr]
        tc2pyr_ampa_conn.append(tc2pyr_ampa)
        k_TC_PY += 1
    #reduce gmax for this cell's synapse, so that the total synaptic strength is equal to that specified in the config file
    #(gmax was already set in cell_classes.py, which used the value specified in config.py)
    #for fac_AMPA_TC, see Krishnan's currents.cpp line 429, and main.cpp; note that network.cfg lists this as a D2 synapse with mini_f=0, which
    #doesn't make sense; but if it was somehow coded as a D2 synapse, then we should instead apply fac_AMPA_D2;
    #so it is unclear whether we should apply fac_AMPA_TC or fac_AMPA_D2 here (I'm guessing D2)
    tc2pyr_ampa_gmax = config.init_AMPA_D2 * config.tc2pyr_ampa_str / k_TC_PY
    tc2pyr_ampa_gmax_wghts.append(tc2pyr_ampa_gmax)
    
##### specify TC->INH AMPA connections (NOTE: this assumes that Ntc=Ninh)
rad=config.tc2inh_ampa_rad #number of outgoing connections from each TC cell
tc2inh_ampa_conn = [] #generates empty list for the synaptic connections of TC->IN AMPA, Note: This will be a list within a list as per NetPyNE documentation.
tc2inh_ampa_gmax_wghts = [] #generates empty list for the weight of the synaptic conductance
for inh_gid in range(config.Npyr, config.Npyr+config.Ninh,1): #loop over all post-synaptic INH cells on this machine (have to loop over post-syn cells bc. pc.gid_connect only works with post-synaptic gid on machine)
    i_inh = inh_gid-config.Npyr #easier to work by indexing the INH cells from 0 to Ninh-1, rather than from Npyr to Npyr+Ninh-1    
    tmp = list(range(i_inh-rad,i_inh+rad+1)) #generate list of TC sources (but this will in general include negative values)
    tc_set=[val % config.Ntc for val in tmp] #apply modulus operation to eliminate negative numbers; this implements periodic boundary conditions; indexes TC cells from 0 to Ntc-1
    k_TC_IN = 0 #make sure in-degree is initialized to zero
    for i_tc in tc_set:
       tc2inh_ampa = [i_tc,i_inh]
       tc2inh_ampa_conn.append(tc2inh_ampa)
       k_TC_IN += 1
    #reduce gmax for this cell's synapse, so that the total synaptic strength is equal to that specified in the config file
    #(gmax was already set in cell_classes.py, which used the value specified in config.py)
    #for fac_AMPA_TC, see Krishnan's currents.cpp line 429, and main.cpp; note that network.cfg lists this as a D2 synapse with mini_f=0, which
    #doesn't make sense; but if it was somehow coded as a D2 synapse, then we should instead apply fac_AMPA_D2
    #so it is unclear whether we should apply fac_AMPA_TC or fac_AMPA_D2 here (I'm guessing D2)
    tc2inh_ampa_gmax = config.init_AMPA_D2 * config.tc2inh_ampa_str / k_TC_IN
    tc2inh_ampa_gmax_wghts.append(tc2inh_ampa_gmax)

##### specify PYR->TC AMPA connections
rad=config.pyr2tc_ampa_rad #number of outgoing connections from each PYR cell
pyr2tc_ampa_conn = [] #generates empty list for the synaptic connections of PY->TC AMPA, Note: This will be a list within a list as per NetPyNE documentation.
pyr2tc_ampa_gmax_wghts = [] #generates empty list for the weight of the synaptic conductance
p2t_ratio = config.Npyr/config.Ntc #we will assume that Npyr>Ntc; this ratio will be important for determining the center of the set of source PYR cells for each post-synaptic TC cell
for tc_gid in range(config.Npyr+config.Ninh+config.Nre, config.Npyr+config.Ninh+config.Nre+config.Ntc,1): #loop over all post-synaptic TC cells on this machine (have to loop over post-syn cells bc. pc.gid_connect only works with post-synaptic gid on machine)
    i_tc = tc_gid - (config.Npyr+config.Ninh+config.Nre) #easier to work by indexing the TC cells from 0 to Ntc-1, rather than from Npyr+Ninh+Nre to from Npyr+Ninh+Nre+Ntc-1
    tmp = list(range( int(np.floor((i_tc-rad)*p2t_ratio)), int(np.floor((i_tc+rad)*p2t_ratio)) +1) ) #center_pyr should be roughly round(i_tc*p2i_ratio); then you need to consider the radius of connectivity surrounding that, and consider that larger p2t_ratio increases the number of PYR cells sending to a given TC cell (for a particular radius of connectivity); "+1" is due to definition of Python's "range" function 
    pyr_set = [val % config.Npyr for val in tmp] #apply modulus operation to eliminate negative numbers; this implements periodic boundary conditions; indexes PYR cells from 0 to Npyr-1
    k_PY_TC = 0 #make sure in-degree is initialized to zero
    for i_pyr in pyr_set:
        pyr2tc_ampa = [i_pyr,i_tc]
        pyr2tc_ampa_conn.append(pyr2tc_ampa)
        k_PY_TC += 1
    #reduce gmax for this cell's synapse, so that the total synaptic strength is equal to that specified in the config file
    #(gmax was already set in cell_classes.py, which used the value specified in config.py)
    #for fac_AMPA_TC, see currents.cpp line 429, and main.cpp
    pyr2tc_ampa_gmax = config.init_AMPA_TC * config.pyr2tc_ampa_str / k_PY_TC
    pyr2tc_ampa_gmax_wghts.append(pyr2tc_ampa_gmax)
    
##### specify PYR->RE AMPA connections
rad=config.pyr2re_ampa_rad #number of outgoing connections from each PYR cell
pyr2re_ampa_conn = [] #generates empty list for the synaptic connections of PY->RE AMPA, Note: This will be a list within a list as per NetPyNE documentation.
pyr2re_ampa_gmax_wghts = [] #generates empty list for the weight of the synaptic conductance
p2r_ratio = config.Npyr/config.Nre #we will assume that Npyr>Nre; this ratio will be important for determining the center of the set of source PYR cells for each post-synaptic RE cell
for re_gid in range(config.Npyr+config.Ninh, config.Npyr+config.Ninh+config.Nre,1): #loop over all post-synaptic RE cells on this machine (have to loop over post-syn cells bc. pc.gid_connect only works with post-synaptic gid on machine)
    i_re = re_gid - (config.Npyr+config.Ninh) #easier to work by indexing the TC cells from 0 to Nre-1, rather than from Npyr+Ninh to from Npyr+Ninh+Nre-1
    tmp = list(range( int(np.floor((i_re-rad)*p2r_ratio)), int(np.floor((i_re+rad)*p2r_ratio)) +1) ) #center_pyr should be roughly round(i_re*p2i_ratio); then you need to consider the radius of connectivity surrounding that, and consider that larger p2r_ratio increases the number of PYR cells sending to a given RE cell (for a particular radius of connectivity); "+1" is due to definition of Python's "range" function 
    pyr_set = [val % config.Npyr for val in tmp] #apply modulus operation to eliminate negative numbers; this implements periodic boundary conditions; indexes PYR cells from 0 to Npyr-1
    k_PY_RE = 0 #make sure in-degree is initialized to zero
    for i_pyr in pyr_set:
        pyr2re_ampa= [i_pyr,i_re]
        pyr2re_ampa_conn.append(pyr2re_ampa)
        k_PY_RE += 1
    #reduce gmax for this cell's synapse, so that the total synaptic strength is equal to that specified in the config file
    #(gmax was already set in cell_classes.py, which used the value specified in config.py)
    #for fac_AMPA_TC, see currents.cpp line 429, and main.cpp
    pyr2re_ampa_gmax = config.init_AMPA_TC * config.pyr2re_ampa_str / k_PY_RE
    pyr2re_ampa_gmax_wghts.append(pyr2re_ampa_gmax)
    