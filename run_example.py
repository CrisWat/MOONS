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
import os


plot_layout()



path_save = "/Users/cristiano/software/MOONLIGHT/test5/"


#initiallize time
s0 = time.time()

#initiallize class
cat = McL()

#If you want to start with data already extracted
#cat.import_catalogue_fromline(ID, ra, dec, l, b, J, H, Ks, rmag, ejk)

#if you want to read informations by the namelist file (If you import catalogue from line, you have to set -> from_line="yes")
cat.read_input_namelist("namelist", from_line="no")

#name of the outputs
name_output = 'field_'+str(cat.CENTRE_CUT[0])+'_'+str(cat.CENTRE_CUT[1])+'_r'+str(cat.RADIUS_CUT)+'.moons'
#name_output = "PROVA3.moons"

#here you make the cut in radius
cat.spherical_cut()
cutted_for_sky = cat.get_data()
cat.get_info()

#here you make the cut in mag
#cat.mag_cut()
ID_cut, ra_cut, dec_cut, l_cut, b_cut, Jmag_cut, Hmag_cut, Kmag_cut, E_jk_cut = cat.get_data()

#here you get the de-reddened mag
#J0, H0, K0 = cat.extinction_correction()



#here you get the indexes of target stars from CMD
#i_target = cat.select_stars_fromCMD(COL="J-H", MAG="H", COL_VAL=[2.5,4], MAG_VAL=[14.5,16.5], return_index=True, print_columns="no")
i_target = cat.select_stars_fromCMD(return_index=True)



#get the data and plot CMD
#ID_target, ra_target, dec_target, l_target, b_target, Jmag_target, Hmag_target, Kmag_target, J0mag_target, H0mag_target, K0mag_target, E_jk_target = cat.get_data_dered(i_target)
ID_target, ra_target, dec_target, l_target, b_target, Jmag_target, Hmag_target, Kmag_target, E_jk_target = cat.get_data(i_target)
ID_all_cut, ra_all_cut, dec_all_cut, l_all_cut, b_all_cut, Jmag_all_cut, Hmag_all_cut, Kmag_all_cut, E_jk_all_cut = cat.get_data()

#select_tellurics
i_tell, j_tell = cat.select_tellurics()
ID_tell, ra_tell, dec_tell, l_tell, b_tell, Jmag_tell, Hmag_tell, Kmag_tell, E_jk_tell = cat.get_data(i_tell)
ID_tell_all, ra_tell_all, dec_tell_all, l_tell_all, b_tell_all, Jmag_tell_all, Hmag_tell_all, Kmag_tell_all, E_jk_tell_all = cat.get_data(j_tell)

file=open(path_save+name_output[:-6]+'.tellurics', "w")
file.write("# ID          RA        DEC       flag   mag  \n".format())
for val in zip(ID_tell_all, ra_tell_all, dec_tell_all, Hmag_tell_all):
      file.write('  {:10}   {:=9.5f}   {:=9.5f}   Tell   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3]))
file.close()

#here you find skypoints
#idd_sky, ra_sky, dec_sky, flag_sky, priority_sky, mag_sky = cat.get_sky(cutted_for_sky[1], cutted_for_sky[2], rr = 3, n_fiber_sky=20, show=False)
idd_sky, ra_sky, dec_sky, flag_sky, priority_sky, mag_sky, idd_sky_all, ra_sky_all, dec_sky_all, flag_sky_all, priority_sky_all, mag_sky_all = cat.get_sky_grid(cutted_for_sky[1], cutted_for_sky[2], allsky=True)

file=open(path_save+name_output[:-6]+'.sky', "w")
file.write("# ID          RA        DEC       flag   mag  \n".format())
for val in zip(idd_sky_all, ra_sky_all, dec_sky_all, mag_sky_all):
      file.write('  {:10}   {:=9.5f}   {:=9.5f}   Sky   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3]))
file.close()

#find the "blend properties"
blend_target = cat.single_malt_star(ID_all_cut, ra_all_cut, dec_all_cut, Jmag_all_cut, ID_target, ra_target, dec_target, Jmag_target, keyw=['a0','a1','a2','a3'])


NSC = 'yes'
#Here you give flag and priority
if ((cat.CENTRE_CUT[0] == 0.0) & (cat.CENTRE_CUT[1] == 0.0) & (NSC=='yes')):
    NSC_radius, l_centre, b_centre = 0.07, 0, 0


    T1_target_NSC = np.where( (blend_target=='a0') & (np.sqrt((l_target-l_centre)**2+(b_target-b_centre)**2)<=NSC_radius))
    T2_target_NSC = np.where( (blend_target=='a1') & (np.sqrt((l_target-l_centre)**2+(b_target-b_centre)**2)<=NSC_radius))
    T3_target = np.where( (blend_target=='a0') & (np.sqrt((l_target-l_centre)**2+(b_target-b_centre)**2)>NSC_radius))
    T4_target = np.where( (blend_target=='a1') & (np.sqrt((l_target-l_centre)**2+(b_target-b_centre)**2)>NSC_radius))
    a = cat.make_priority_and_flag(ID_target, ra_target, dec_target, Jmag_target, T1_target_NSC, 3, 'T')
    b = cat.make_priority_and_flag(ID_target, ra_target, dec_target, Jmag_target, T2_target_NSC, 5, 'T')
    c = cat.make_priority_and_flag(ID_target, ra_target, dec_target, Jmag_target, T3_target, 4, 'T')
    d = cat.make_priority_and_flag(ID_target, ra_target, dec_target, Jmag_target, T4_target, 6, 'T')



    file=open(path_save+name_output, "w")
    file.write("# ID          RA        DEC        flag   priority   mag  \n".format())
    for val in zip(idd_sky, ra_sky, dec_sky, mag_sky):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   1   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3]))
    for val in zip(ID_tell, ra_tell, dec_tell, Hmag_tell):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   2   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3]))
    for val in zip(a[0], a[1], a[2], a[4], a[5]):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
    for val in zip(b[0], b[1], b[2], b[4], b[5]):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
    for val in zip(c[0], c[1], c[2], c[4], c[5]):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
    for val in zip(d[0], d[1], d[2], d[4], d[5]):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
    file.close()

    #nams1 = path_save+name_output
    #cat.control_doubles(nams1)


    file=open(path_save+name_output[:-6]+'.type', "w")
    file.write("# ID          RA        DEC        flag   priority   mag  \n".format())
    for val in zip(idd_sky, ra_sky, dec_sky, mag_sky):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   Sky    1   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3]))
    for val in zip(ID_tell, ra_tell, dec_tell, Hmag_tell):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   Tell   2   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3]))
    for val in zip(a[0], a[1], a[2], a[4], a[5]):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   RGB    {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
    for val in zip(b[0], b[1], b[2], b[4], b[5]):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   RGB    {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
    for val in zip(c[0], c[1], c[2], c[4], c[5]):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   RGB    {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
    for val in zip(d[0], d[1], d[2], d[4], d[5]):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   RGB    {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
    file.close()

    #nams2 = path_save+name_output[:-6]+'.type'
    #cat.control_doubles(nams2)



else:
    T1_target = np.where((blend_target=='a0'))
    T2_target = np.where((blend_target=='a1'))
    a = cat.make_priority_and_flag(ID_target, ra_target, dec_target, Jmag_target, T1_target, 3, 'T')
    b = cat.make_priority_and_flag(ID_target, ra_target, dec_target, Jmag_target, T2_target, 4, 'T')



    file=open(path_save+name_output, "w")
    file.write("# ID          RA        DEC        flag   priority   mag  \n".format())
    for val in zip(idd_sky, ra_sky, dec_sky, mag_sky):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   1   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3]))
    for val in zip(ID_tell, ra_tell, dec_tell, Hmag_tell):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   2   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3]))
    for val in zip(a[0], a[1], a[2], a[4], a[5]):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
    for val in zip(b[0], b[1], b[2], b[4], b[5]):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   T   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
    file.close()


    #nams1 = path_save+name_output
    #cat.control_doubles(nams1)


    file=open(path_save+name_output[:-6]+'.type', "w")
    file.write("# ID          RA        DEC        flag   priority   mag  \n".format())
    for val in zip(idd_sky, ra_sky, dec_sky, mag_sky):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   Sky    1   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3]))
    for val in zip(ID_tell, ra_tell, dec_tell, Hmag_tell):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   Tell   2   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3]))
    for val in zip(a[0], a[1], a[2], a[4], a[5]):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   RGB    {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
    for val in zip(b[0], b[1], b[2], b[4], b[5]):
        file.write('  {:10}   {:=9.5f}   {:=9.5f}   RGB    {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4]))
    file.close()


    #nams2 = path_save+name_output[:-6]+'.type'
    #cat.control_doubles(nams2)




#here you write the input catalogue and parameter file for MOONLIGHT
cat.write_MOONLIGHT_param(savepath=path_save, cat_name=name_output, param_name=name_output+".ini")

#print the time
s1 = time.time()
print(round(s1-s0,1), "seconds to make catalogues for MOONS")


"""

#plot the extintion map
cat.plot_extinction(Ai="Ah", namefig=name_output+'extinction.png', savepath=path_save)


#plot the Hess CMD
K0_all  = K0mag_all_cut
JK0_all = J0mag_all_cut - K0mag_all_cut

K0_RC   = K0mag_RC
JK0_RC  = J0mag_RC - K0mag_RC

K0_RGB  = K0mag_RGB
JK0_RGB = J0mag_RGB - K0mag_RGB

HessCMD.plotHess(K0_all, JK0_all,K0_RC,JK0_RC,K0_RGB,JK0_RGB, levels=np.arange(20,300,10), cbarrtitle='Number', xlab='$J_{0}$ $-$ $K_{0} (mag)$', ylab='$Ks_{0} (mag)$', savepath=path_save,namefig=name_output+'.png', cbarr='Yes',xlims=[-1,5], ylims=[15,5],colormap='jet',ftitle = name_output+".pdf")

"""
