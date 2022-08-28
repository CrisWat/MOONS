#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 12:02:56 2022

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
ID_cut, ra_cut, dec_cut, l_cut, b_cut, Jmag_cut, Hmag_cut, Kmag_cut, rmag_cut, E_jk_cut = cat.get_data()

#here you get the de-reddened mag
J0, H0, K0 = cat.extinction_correction()


#here you get the indexes of target & acquisition stars from CMD from line
i_RC, j_RC = cat.select_stars_fromCMD(COL="J-K", MAG="K", COL_VAL=[0.7,1.2], MAG_VAL=[12,13.5], return_index=True, print_columns="no")
i_RGB, j_RGB = cat.select_stars_fromCMD(COL="J-K", MAG="K", COL_VAL=[0.3,2.0], MAG_VAL=[10,12], return_index=True)



#get the data and plot CMD
ID_RC, ra_RC, dec_RC, l_RC, b_RC, Jmag_RC, Hmag_RC, Kmag_RC, J0mag_RC, H0mag_RC, K0mag_RC, rmag_RC, E_jk_RC = cat.get_data_dered(i_RC)
ID_RGB, ra_RGB, dec_RGB, l_RGB, b_RGB, Jmag_RGB, Hmag_RGB, Kmag_RGB, J0mag_RGB, H0mag_RGB, K0mag_RGB, rmag_RGB, E_jk_RGB = cat.get_data_dered(i_RGB)
ID_all_cut, ra_all_cut, dec_all_cut, l_all_cut, b_all_cut, Jmag_all_cut, Hmag_all_cut, Kmag_all_cut, J0mag_all_cut, H0mag_all_cut, K0mag_all_cut, rmag_all_cut, E_jk_all_cut = cat.get_data_dered()


#here you find skypoints
idd_sky, ra_sky, ra_dec, flag_sky, priority_sky, mag_sky = cat.get_sky(cutted_for_sky[1], cutted_for_sky[2], rr = 1, n_fiber_sky=1, show=False)

#name of the outputs
name_output = 'field_'+str(cat.CENTRE_CUT[0])+'_'+str(cat.CENTRE_CUT[1])+'_r'+str(cat.RADIUS_CUT)+'.moons'
#name_output = "PROVA3.moons"

#ID_RC, Ra_RC, Dec_RC, flag_RC, priority_RC, mag_input_RC = cat.make_priority_and_flag(i_RC, 1, 'T')






#plot the extintion map
cat.plot_extinction(Ai="Ah", namefig=name_output+'extinction.png')


#plot the Hess CMD
K0_all  = K0mag_all_cut
H0_all  = H0mag_all_cut
J0_all  = J0mag_all_cut
JK0_all = J0mag_all_cut - K0mag_all_cut
JH0_all = J0mag_all_cut - H0mag_all_cut

K0_RC   = K0mag_RC
JK0_RC  = J0mag_RC - K0mag_RC

K0_RGB  = K0mag_RGB
JK0_RGB = J0mag_RGB - K0mag_RGB


v = np.where((Hmag_all_cut>0) & (Jmag_all_cut>0) )
JH_all_cut = Jmag_all_cut-Hmag_all_cut
v2 = np.where((H0_all>0) & (J0_all>0) )


HessCMD.plotHess3(Hmag_all_cut[v], JH_all_cut[v], levels=np.arange(20,130, 10), cbarrtitle='Number', xlab='J-H', ylab='H', saveas=name_output+'.png', cbarr='Yes',xlims=[-0.51,5.99], ylims=[18.9,12.1],colormap='jet',ftitle = name_output+".pdf")
#HessCMD.plotHess3(H0_all[v2], J0_all[v2]-H0_all[v2], levels=np.arange(20,210, 10), cbarrtitle='Number', xlab='J-H', ylab='H', saveas=name_output+'.png', cbarr='Yes',xlims=[-3,3], ylims=[15.9,9],colormap='jet',ftitle = name_output+".pdf")












