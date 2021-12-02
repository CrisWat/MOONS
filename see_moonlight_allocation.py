#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 10:18:21 2021

@author: cristiano
"""
import numpy as np
import matplotlib.pyplot as plt
import HessCMD
from MOONcatLIGHT import McL
from matplotlib_param import plot_layout
import time
from astroML.crossmatch import crossmatch_angular
from astropy.visualization import hist


plot_layout()


#initiallize class
cat = McL()

#If you want to start with data already extracted
#cat.import_catalogue_fromline(ID, ra, dec, l, b, J, H, Ks, rmag, ejk)

#if you want to read informations by the namelist file (If you import catalogue from line, you have to set -> from_line="yes")
cat.read_input_namelist("namelist", from_line="no")
ID, ra, dec, l, b, Jmag, Hmag, Kmag, rmag, E_jk = cat.get_data()

#here you make the cut in radius
c_cut = (0.0,0.0)
r_cut = 0.25
cat.spherical_cut(CENTRE_CUT=c_cut, RADIUS_CUT=r_cut)
info = cat.get_info(return_dict=True)

#here you make the cut in mag
cat.mag_cut()
ID_cut, ra_cut, dec_cut, l_cut, b_cut, Jmag_cut, Hmag_cut, Kmag_cut, rmag_cut, E_jk_cut = cat.get_data()
J0, H0, K0 = cat.extinction_correction()
#ID_cut, ra_cut, dec_cut, l_cut, b_cut, Jmag_cut, Hmag_cut, Kmag_cut, J0mag_cut, H0mag_cut, K0mag_cut, rmag_cut, E_jk_cut = cat.get_data_dered()

field_center = info['Field center']



path = '/Users/cristiano/software/MOONLIGHT/test3/field_'+str(c_cut[0])+'_'+str(c_cut[1])+'_r'+str(r_cut)+'/results/'
path = '/Users/cristiano/software/MOONLIGHT/test3/field_'+str(c_cut[0])+'_'+str(c_cut[1])+'_r'+str(r_cut)+'_second/results/'

RAa,DECa,prioritya,maga,obsa,infova = np.genfromtxt(path+'Catalogue.txt', usecols=(1,2,4,5,6,7), dtype='float', comments='#', unpack=True)
IDa,flaga = np.genfromtxt(path+'Catalogue.txt', usecols=(0,3), dtype='str', comments='#', unpack=True)
#ID_tar,x_tar,y_tar,size_tar,flag_tar,priority_tar,obs_tar = np.genfromtxt(path+'Catalogue_pix.txt', usecols=(0,1,2,3,4,5,6), dtype='float', comments='#', unpack=True)


#IDt,RAt,DECt,priorityt,magt,obst,infovt = np.genfromtxt(path+'field_'+str(c_cut[0])+'_'+str(c_cut[1])+'_r'+str(r_cut)+'.type', usecols=(0,1,2,4,5,6,7), dtype='float', comments='#', unpack=True)
typet = np.genfromtxt(path+'field_'+str(c_cut[0])+'_'+str(c_cut[1])+'_r'+str(r_cut)+'.type', usecols=(3), dtype='str', comments='#', unpack=True)




ia = np.where((flaga=='A') & (obsa==2))
i = np.where((flaga=='T') & (obsa==2))
i1 = np.where((flaga=='T') & (obsa==2))
i2 = np.where((flaga=='T') & (obsa==4))
i3 = np.where((flaga=='T') & (obsa!=0) & (obsa!=1) & (obsa!=2))

#j = np.where((flag_tar==1) & (priority_tar==1) & (obs_tar==2))

j = np.zeros((len(RAa[i1])))
for n in range(len(RAa[i1])):
    if (IDa[i1][n][0]) == 's':
        j[n] = 1
sky = np.where(j==1)

fig, ax = plt.subplots(figsize=(6,6), dpi=500)
plt.scatter(ra, dec, s=0.0001, color='gray')
plt.scatter(ra_cut, dec_cut, s=0.001, color='k')
#plt.scatter(x_tar[j],y_tar[j], s=4, color='r')
#plt.scatter(RAa[i], DECa[i], s=4, color='r')
plt.scatter(RAa[i1][typet[i1]=='RC'], DECa[i1][typet[i1]=='RC'], marker='o', s=8, color='r')
plt.scatter(RAa[i1][typet[i1]=='RGB'], DECa[i1][typet[i1]=='RGB'], marker='^', s=8, color='r')
plt.scatter(RAa[i2][typet[i2]=='RC'], DECa[i2][typet[i2]=='RC'], marker='o', s=8, color='b')
plt.scatter(RAa[i2][typet[i2]=='RGB'], DECa[i2][typet[i2]=='RGB'], marker='^', s=8, color='b')
plt.scatter(RAa[i1][sky], DECa[i1][sky], s=4, color='c')
#plt.scatter(RAa[ia], DECa[ia], s=4, color='r')





"""
print()
print('Allocated Fibers: ',len(RAa[i]))
print('RC NSC: ',len(RAa[i1][typet[i1]=='RC']))
print('RGB NSC: ',len(RAa[i1][typet[i1]=='RGB']))
print('RC OST: ',len(RAa[i2][typet[i2]=='RC']))
print('RGB OST: ',len(RAa[i2][typet[i2]=='RGB']))
print('Sky: ',len(RAa[i1][sky]))
print()

print()
print('Allocated Fibers: ',len(RAa[i]))
print('RC : ',len(RAa[i1][typet[i1]=='RC']))
print('RGB: ',len(RAa[i1][typet[i1]=='RGB']))
print('RC : ',len(RAa[i2][typet[i2]=='RC']))
print('RGB: ',len(RAa[i2][typet[i2]=='RGB']))
print('Sky: ',len(RAa[i1][sky]))
print()
"""

d, m = 0.3, 345
plt.xlim(field_center[0]-d, field_center[0]+d)
plt.ylim(field_center[1]-d, field_center[1]+d)
#plt.xlim((alpha[m])-d, (alpha[m])+d)
#plt.ylim((delta[m])-d, (delta[m])+d)
plt.minorticks_on()
plt.show()
plt.close()




IDa = np.genfromtxt(path+'Catalogue.txt', usecols=(0), dtype='float', comments='#', unpack=True)
 
    
k = [np.where(ID_cut==IDa[i][n])[0] for n in range(len(IDa[i]))]
k = np.concatenate(k, axis=0) 


    
HessCMD.plotHess2(K0, J0-K0, K0[k], J0[k]-K0[k], levels=np.arange(20,300,10), cbarrtitle='Number', xlab='$J_{0}$ $-$ $K_{0} (mag)$', ylab='$Ks_{0} (mag)$', saveas=None, cbarr='Yes',xlims=[-1,5], ylims=[15,5],colormap='jet',ftitle = None)

    

cat.plot_extinction_fibers(l_cut[k], b_cut[k], Ai="Ah", namefig = None)

   

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

"""
num = np.genfromtxt(path+'moonLight.paf', usecols=(1), dtype='str', delimiter='.', comments='#', skip_header=54, unpack=True)
keyw, keyv = np.genfromtxt(path+'moonLight.paf', usecols=(0,1), dtype='str', comments=';', skip_header=54, deletechars='"', unpack=True)

num = [int(num[n][3:]) for n in range(len(num))]

for n in range(len(keyw)):
    if num[n]<=9:
        keyw[n]=keyw[n][9:]
    elif (num[n]>=10) & (num[n]<=99):
        keyw[n]=keyw[n][10:]
    elif num[n]>=100:
        keyw[n]=keyw[n][11:]
    
keyword_list = keyw[0:9]


for n in range(len(keyv)):
    if (keyv[n][0]=='"'):
        keyv[n] = keyv[n][1:-1]






alpha, delta = (keyv[np.where(keyw==keyword_list[1])]), (keyv[np.where(keyw==keyword_list[2])])

alpha = [ float(alpha[n]) for n in range(len(alpha)) ]
delta = [ float(delta[n]) for n in range(len(delta)) ]
"""





