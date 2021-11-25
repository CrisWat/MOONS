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
idd_sky, ra_sky, ra_dec, flag_sky, priority_sky, mag_sky = cat.get_sky(cutted_for_sky[1], cutted_for_sky[2], rr = 3, n_fiber_sky=30, show=False)

#name of the outputs
name_output = "PROVA3.moons"

#ID_RC, Ra_RC, Dec_RC, flag_RC, priority_RC, mag_input_RC = cat.make_priority_and_flag(i_RC, 1, 'T')


blend_RC = cat.single_malt_star(ID_all_cut, ra_all_cut, dec_all_cut, Jmag_all_cut, ID_RC, ra_RC, dec_RC, Jmag_RC, rr=1.4/3600, keyw=['a0','a1','a2','a3'])
blend_RGB = cat.single_malt_star(ID_all_cut, ra_all_cut, dec_all_cut, Jmag_all_cut, ID_RGB, ra_RGB, dec_RGB, Jmag_RGB, rr=1.4/3600, keyw=['a0','a1','a2','a3'])




#Here you give flag and priority
NSC_radius, l_centre, b_centre = 0.07, 0, 0

A0 = np.where((rmag_cut!=-999) & (rmag_cut<16) )

T0_RC_NSC = np.where( (blend_RC=='a0') & (np.sqrt((l_RC-l_centre)**2+(b_RC-b_centre)**2)<=NSC_radius))
T0_RGB_NSC = np.where( (blend_RGB=='a0') & (np.sqrt((l_RGB-l_centre)**2+(b_RGB-b_centre)**2)<=NSC_radius))

T1_RC_NSC = np.where( (blend_RC=='a1') & (np.sqrt((l_RC-l_centre)**2+(b_RC-b_centre)**2)<=NSC_radius))
T1_RGB_NSC = np.where( (blend_RGB=='a1') & (np.sqrt((l_RGB-l_centre)**2+(b_RGB-b_centre)**2)<=NSC_radius))

T2_RC = np.where( (blend_RC=='a0') & (np.sqrt((l_RC-l_centre)**2+(b_RC-b_centre)**2)>NSC_radius))
T2_RGB = np.where( (blend_RGB=='a0') & (np.sqrt((l_RGB-l_centre)**2+(b_RGB-b_centre)**2)>NSC_radius))

T3_RC = np.where( (blend_RC=='a1') & (np.sqrt((l_RC-l_centre)**2+(b_RC-b_centre)**2)>NSC_radius))
T3_RGB = np.where( (blend_RGB=='a1') & (np.sqrt((l_RGB-l_centre)**2+(b_RGB-b_centre)**2)>NSC_radius))

T4_RC_NSC = np.where( (blend_RC=='a2') & (np.sqrt((l_RC-l_centre)**2+(b_RC-b_centre)**2)<=NSC_radius))
T4_RGB_NSC = np.where( (blend_RGB=='a2') & (np.sqrt((l_RGB-l_centre)**2+(b_RGB-b_centre)**2)<=NSC_radius))

T5_RC = np.where( (blend_RC=='a2') & (np.sqrt((l_RC-l_centre)**2+(b_RC-b_centre)**2)>NSC_radius))
T5_RGB = np.where( (blend_RGB=='a2') & (np.sqrt((l_RGB-l_centre)**2+(b_RGB-b_centre)**2)>NSC_radius))



a = cat.make_priority_and_flag(ID_cut, ra_cut, dec_cut, rmag_cut, A0, 0, 'A')
b = cat.make_priority_and_flag(ID_RC, ra_RC, dec_RC, Jmag_RC, T0_RC_NSC, 1, 'T')
c = cat.make_priority_and_flag(ID_RGB, ra_RGB, dec_RGB, Jmag_RGB, T0_RGB_NSC, 1, 'T')
d = cat.make_priority_and_flag(ID_RC, ra_RC, dec_RC, Jmag_RC, T1_RC_NSC, 1, 'T')
e = cat.make_priority_and_flag(ID_RGB, ra_RGB, dec_RGB, Jmag_RGB, T1_RGB_NSC, 1, 'T')
d = cat.make_priority_and_flag(ID_RC, ra_RC, dec_RC, Jmag_RC, T2_RC, 2, 'T')
e = cat.make_priority_and_flag(ID_RGB, ra_RGB, dec_RGB, Jmag_RGB, T2_RGB, 2, 'T')
f = cat.make_priority_and_flag(ID_RC, ra_RC, dec_RC, Jmag_RC, T3_RC, 3, 'T')
g = cat.make_priority_and_flag(ID_RGB, ra_RGB, dec_RGB, Jmag_RGB, T3_RGB, 3, 'T')
h = cat.make_priority_and_flag(ID_RC, ra_RC, dec_RC, Jmag_RC, T4_RC_NSC, 4, 'T')
i = cat.make_priority_and_flag(ID_RGB, ra_RGB, dec_RGB, Jmag_RGB, T4_RGB_NSC, 4, 'T')
l = cat.make_priority_and_flag(ID_RC, ra_RC, dec_RC, Jmag_RC, T5_RC, 5, 'T')
m = cat.make_priority_and_flag(ID_RGB, ra_RGB, dec_RGB, Jmag_RGB, T5_RGB, 5, 'T')






file=open('./'+name_output, "w")
file.write("# ID          RA        DEC        flag   priority   mag  \n".format())
for val in zip(a[0], a[1], a[2], a[4], a[5]):
    file.write('  {:10}   {:=9.5f}   {:=9.5f}   A   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
for val in zip(b[0], b[1], b[2], b[4], b[5]):
    file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
for val in zip(c[0], c[1], c[2], c[4], c[5]):
    file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
for val in zip(idd_sky, ra_sky, ra_dec, priority_sky, mag_sky):
    file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
for val in zip(d[0], d[1], d[2], d[4], d[5]):
    file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
for val in zip(e[0], e[1], e[2], e[4], e[5]):
    file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
for val in zip(f[0], f[1], f[2], f[4], f[5]):
    file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
for val in zip(g[0], g[1], g[2], g[4], g[5]):
    file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
for val in zip(h[0], h[1], h[2], h[4], h[5]):
    file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
for val in zip(i[0], i[1], i[2], i[4], i[5]):
    file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
for val in zip(l[0], l[1], l[2], l[4], l[5]):
    file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
for val in zip(m[0], m[1], m[2], m[4], m[5]):
    file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
file.close()



#here you write the input catalogue and parameter file for MOONLIGHT
cat.write_MOONLIGHT_param(path_save="./", cat_name=name_output, param_name=name_output+".ini")

#print the time
s1 = time.time()
print(round(s1-s0,1), "seconds to make catalogues for MOONS")




#plot the extintion map
cat.plot_extinction(Ai="Ah", namefig="extintion_plot")


#plot the Hess CMD
K0_all  = K0mag_all_cut
JK0_all = J0mag_all_cut - K0mag_all_cut

K0_RC   = K0mag_RC
JK0_RC  = J0mag_RC - K0mag_RC

K0_RGB  = K0mag_RGB
JK0_RGB = J0mag_RGB - K0mag_RGB

HessCMD.plotHess(K0_all, JK0_all,K0_RC,JK0_RC,K0_RGB,JK0_RGB, levels=np.arange(20,300,10), cbarrtitle='Number', xlab='$J_{0}$ $-$ $K_{0} (mag)$', ylab='$Ks_{0} (mag)$', saveas='PROVA1.png', cbarr='Yes',xlims=[-1,5], ylims=[15,5],colormap='jet',ftitle = name_output+".pdf")







#here you merge catalogues derived from the two colors selections TODO
#IDw, Raw, Decw, flagw, priorityw, mag_inputw = cat.merge_catalogues(ID_RC, Ra_RC, Dec_RC, flag_RC, priority_RC, mag_input_RC,ID_RGB, Ra_RGB, Dec_RGB, flag_RGB, priority_RGB, mag_input_RGB)
#cat.write_MOONLIGHT_input(IDw, Raw, Decw, flagw, priorityw, mag_inputw, path_save="./", cat_name=name_output)




