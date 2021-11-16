#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 17:50:27 2021

@author: cristiano
"""
import numpy as np
import HessCMD
from MOONcatLIGHT import McL
from matplotlib_param import plot_layout
import time

plot_layout()

#initiallize time
s0 = time.time()

#initiallize class
cat = McL()

#If you want to start with data already extracted
#cat.import_catalogue_fromline(ID, ra, dec, l, b, J, H, Ks, rmag, ejk)

#if you want to read informations by the namelist file (If you import catalogue from line, you have to set -> from_line="yes")
cat.read_input_namelist("namelist", from_line="no")

#here you make the cut in radius
cat.spherical_cut()
cutted_for_sky = cat.get_data()
cat.get_info()

#here you make the cut in mag
cat.mag_cut()
cutted_mag_and_radius = cat.get_data()

#here you get the cutted data into an ndarray
cat.mag_cut()

#here you get the de-reddened mag
J0, H0, K0 = cat.extinction_correction()


#here you get the indexes of target & acquisition stars from CMD from line
i_RC, j_RC = cat.select_stars_fromCMD(COL="J-K", MAG="K", COL_VAL=[0.7,1.2], MAG_VAL=[12,13.5], return_index=True, print_columns="no")
i_RGB, j_RGB = cat.select_stars_fromCMD(COL="J-K", MAG="K", COL_VAL=[0.3,2.0], MAG_VAL=[10,12], return_index=True)

#here you find skypoints
sky = cat.get_sky( cutted_for_sky[1], cutted_for_sky[2], rr = 1.4, n_fiber_sky=5, show=True)

#name of the outputs
name_output = "PROVA1.moons"

#Here you give flag and priority
ID_RC, Ra_RC, Dec_RC, flag_RC, priority_RC, mag_input_RC = cat.make_priority_and_flag(i_RC, 1, j_RC, 0)
ID_RGB, Ra_RGB, Dec_RGB, flag_RGB, priority_RGB, mag_input_RGB = cat.make_priority_and_flag(i_RGB, 3, j_RGB, 2)

#here you merge catalogues derived from the two colors selections
IDw, Raw, Decw, flagw, priorityw, mag_inputw = cat.merge_catalogues(ID_RC, Ra_RC, Dec_RC, flag_RC, priority_RC, mag_input_RC,ID_RGB, Ra_RGB, Dec_RGB, flag_RGB, priority_RGB, mag_input_RGB)

#here you write the input catalogue and parameter file for MOONLIGHT
#cat.write_MOONLIGHT_input(IDw, Raw, Decw, flagw, priorityw, mag_inputw, path_save="./", cat_name=name_output)
#cat.write_MOONLIGHT_param(path_save="./", cat_name=name_output, param_name=name_output+".ini")

#print the time
s1 = time.time()
print(round(s1-s0,1), "seconds to make catalogues for MOONS")



#plot the extintion map
cat.plot_extinction(Ai="Ah")


#get the data and plot CMD
ID_RC, ra_RC, dec_RC, l_RC, b_RC, Jmag_RC, Hmag_RC, Kmag_RC, J0mag_RC, H0mag_RC, K0mag_RC, rmag_RC, E_jk_RC = cat.get_data_dered(i_RC)
ID_RGB, ra_RGB, dec_RGB, l_RGB, b_RGB, Jmag_RGB, Hmag_RGB, Kmag_RGB, J0mag_RGB, H0mag_RGB, K0mag_RGB, rmag_RGB, E_jk_RGB = cat.get_data_dered(i_RGB)
ID_all_cut, ra_all_cut, dec_all_cut, l_all_cut, b_all_cut, Jmag_all_cut, Hmag_all_cut, Kmag_all_cut, J0mag_all_cut, H0mag_all_cut, K0mag_all_cut, rmag_all_cut, E_jk_all_cut = cat.get_data_dered()

K0_all  = K0mag_all_cut
JK0_all = J0mag_all_cut - K0mag_all_cut

K0_RC   = K0mag_RC
JK0_RC  = J0mag_RC - K0mag_RC

K0_RGB  = K0mag_RGB
JK0_RGB = J0mag_RGB - K0mag_RGB

HessCMD.plotHess(K0_all, JK0_all,K0_RC,JK0_RC,K0_RGB,JK0_RGB, levels=np.arange(20,300,10), cbarrtitle='Number', xlab='$J_{0}$ $-$ $K_{0} (mag)$', ylab='$Ks_{0} (mag)$', saveas='PROVA1.png', cbarr='Yes',xlims=[-1,5], ylims=[15,5],colormap='jet',ftitle= name_output+".pdf")
