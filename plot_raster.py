# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 12:54:02 2020

This file is used to plot the raster data. 

@author: Logan Bayer
"""
import pickle
from matplotlib import pyplot as plt
from pathlib import Path
import numpy as np

Npyr = 500 
Ninh = 100 
Nre = 100 
Ntc = 100  
pyrList = []
inhList = [] 
reList = [] 
tcList = [] 
pyrTime = [] 
inhTime = [] 
reTime = [] 
tcTime = [] 

filename = 'raster_test.pkl'
infile = open(filename, 'rb')
test = pickle.load(infile)
infile.close()

spk_times = test['spkTimes']
spk_ids = test['spkInds']
spk_times_arr = np.array(spk_times)
spk_ids_arr = np.array(spk_ids)

for i_len in range(len(spk_times_arr)):
    if spk_ids_arr[i_len] < Npyr:
        pyrList.append(spk_ids_arr[i_len])
        pyrTime.append(spk_times_arr[i_len])
    if spk_ids_arr[i_len] >= Npyr and spk_ids_arr[i_len] < (Npyr+Ninh):
        inhList.append(spk_ids_arr[i_len])
        inhTime.append(spk_times_arr[i_len])
    if spk_ids_arr[i_len] >= (Npyr+Ninh) and spk_ids_arr[i_len] < (Npyr+Ninh+Nre):
        reList.append(spk_ids_arr[i_len])
        reTime.append(spk_times_arr[i_len])
    if spk_ids_arr[i_len] >= (Npyr+Ninh+Nre) and spk_ids_arr[i_len] < (Npyr+Ninh+Nre+Ntc):
        tcList.append(spk_ids_arr[i_len])
        tcTime.append(spk_times_arr[i_len])

plt.figure() 
plt.title('NP_v4_awake_post_summer', loc = 'center')
plt.scatter(pyrTime, pyrList, marker='o', s=5, color='red')
plt.scatter(inhTime, inhList, marker='o', s=5, color='blue')
plt.scatter(reTime, reList, marker='o', s=5, color='green')
plt.scatter(tcTime, tcList, marker='o', s=5, color='orange')