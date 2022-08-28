#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 12:35:09 2021
@author: Cristiano Fanelli
"""

import numpy as np
from astropy.io import fits
#import matplotlib
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.spatial import KDTree
from collections import Counter
from scipy.spatial import distance_matrix
from sklearn.cluster import KMeans


class McL:


    def __init__ (self):
        print()
        print("********************************")
        print("***MOONcatLIGHT Class Started***")
        print("********************************")
        print()



    def import_catalogues(self, pathfile, cols=(0,1,2,3,4,5,6,7,9)):
        """
        This function identifies the extension of the catalog file
        and the number of catalogues to merge
        --------------------------
        pathfile : the entire path of the catalogue
        cols : the number of columns used for the data extraction
        """
        extension = pathfile[-1:-5:-1][::-1]
        if extension=='fits':
            type_input='fits'
        else:
            type_input='ascii'

        comma = McL.number_of_cat(self,pathfile)


        if len(comma)==1:

            #Da Fare il check sui cataloghi multipli in txt
            if type_input=='ascii':
                self.id_obj, self.ra, self.dec, self.l, self.b, self.Jmag, self.Hmag, self.Kmag, self.E_jk = np.genfromtxt(pathfile, usecols=cols, dtype='float', unpack=True)
            elif type_input=='fits':
                f           =   fits.open(pathfile)
                f_in        =   f[1].data
                self.id_obj =   f_in['SOURCEID'].astype(int)
                self.ra     =   f_in['RA2000']
                self.dec    =   f_in['DEC2000']
                self.l      =   f_in['L']
                self.b      =   f_in['B']
                self.Jmag   =   f_in['J']
                self.Hmag   =   f_in['H']
                self.Kmag   =   f_in['K']
                self.E_jk   =   f_in['EJK']
                self.Hmag[np.isnan(self.Hmag)]   =  -999
                self.Jmag[np.isnan(self.Jmag)]   =  -999
                self.Kmag[np.isnan(self.Kmag)]   =  -999
                self.E_jk[np.isnan(self.E_jk)]   =  -999


            McL.clean_for_duplicate(self)
        elif len(comma)>1:
            paths = McL.import_paths(self,pathfile, comma)

            if type_input=='ascii':
                for n in range(len(paths)):
                    at, bt, ct, dt, et, ft, gt, ht, it, lt = np.genfromtxt(paths[n], usecols=cols, dtype='float', unpack=True)
                    self.id_obj.append(at)
                    self.ra.append(bt)
                    self.dec.append(ct)
                    self.l.append(dt)
                    self.b.append(et)
                    self.Jmag.append(ft)
                    self.Hmag.append(gt)
                    self.Kmag.append(ht)
                    self.E_jk.append(lt)

            elif type_input=='fits':
                at, bt, ct, dt, et, ft, gt, ht, it, lt = [], [], [], [], [], [], [], [], [], []
                for n in range(len(paths)):
                    f       =   fits.open(paths[n])
                    f_in    =   f[1].data
                    at.append(f_in['SOURCEID'].astype(int))
                    bt.append(f_in['RA2000'])
                    ct.append(f_in['DEC2000'])
                    dt.append(f_in['L'])
                    et.append(f_in['B'])
                    ft.append(f_in['J'])
                    gt.append(f_in['H'])
                    ht.append(f_in['K'])
                    lt.append(f_in['EJK'])
                at   =   np.asarray(at, dtype=object)
                bt   =   np.asarray(bt, dtype=object)
                ct   =   np.asarray(ct, dtype=object)
                dt   =   np.asarray(dt, dtype=object)
                et   =   np.asarray(et, dtype=object)
                ft   =   np.asarray(ft, dtype=object)
                gt   =   np.asarray(gt, dtype=object)
                ht   =   np.asarray(ht, dtype=object)
                lt   =   np.asarray(lt, dtype=object)
                if len(paths) == 2:
                    self.id_obj   =   np.concatenate((at[0],at[1]))
                    self.ra       =   np.concatenate((bt[0],bt[1]))
                    self.dec      =   np.concatenate((ct[0],ct[1]))
                    self.l        =   np.concatenate((dt[0],dt[1]))
                    self.b        =   np.concatenate((et[0],et[1]))
                    self.Jmag     =   np.concatenate((ft[0],ft[1]))
                    self.Hmag     =   np.concatenate((gt[0],gt[1]))
                    self.Kmag     =   np.concatenate((ht[0],ht[1]))
                    self.E_jk     =   np.concatenate((lt[0],lt[1]))
                elif len(paths) == 3:
                    self.id_obj   =   np.concatenate((at[0],at[1],at[2]))
                    self.ra       =   np.concatenate((bt[0],bt[1],bt[2]))
                    self.dec      =   np.concatenate((ct[0],ct[1],ct[2]))
                    self.l        =   np.concatenate((dt[0],dt[1],dt[2]))
                    self.b        =   np.concatenate((et[0],et[1],et[2]))
                    self.Jmag     =   np.concatenate((ft[0],ft[1],ft[2]))
                    self.Hmag     =   np.concatenate((gt[0],gt[1],gt[2]))
                    self.Kmag     =   np.concatenate((ht[0],ht[1],ht[2]))
                    self.E_jk     =   np.concatenate((lt[0],lt[1],lt[2]))
                self.Jmag[np.isnan(self.Jmag)]   =  -999
                self.Hmag[np.isnan(self.Hmag)]   =  -999
                self.Kmag[np.isnan(self.Kmag)]   =  -999
                self.E_jk[np.isnan(self.E_jk)]   =  -999

                McL.clean_for_duplicate(self)

        return self.id_obj, self.ra, self.dec, self.l, self.b, self.Jmag, self.Hmag, self.Kmag, self.E_jk



    def import_catalogue_fromline(self, id_obj, ra, dec, l, b, Jmag, Hmag, Kmag, E_jk):
        """
        This function import data of a catalog file directly from line.
        Here is to make variables explicit
        """

        self.id_obj    =     id_obj
        self.ra        =     ra
        self.dec       =     dec
        self.l         =     l
        self.b         =     b
        self.Jmag      =     Jmag
        self.Hmag      =     Hmag
        self.Kmag      =     Kmag
        self.E_jk      =     E_jk

        McL.clean_for_duplicate(self)


    def number_of_cat(self,pathfile):
        """
        This function identifies the number of the catalog file
        to merge
        --------------------------
        pathfile : the entire path of the catalogue
        """
        comma = []
        for n in range(len(pathfile)):
            if pathfile[n]==",":
                comma.append(n)
        comma.append(len(pathfile))

        return comma


    def import_paths(self, pathfile, comma):
        """
        This function split the catalogues path, if they are
        more than 1
        --------------------------
        pathfile : the entire path of the catalogue
        comma : dimension of the path array
        """
        paths = ['a']*len(comma)
        for m in range(len(paths)):
            if m==0:
                paths[m]= pathfile[0:comma[m]]
            elif (m>0) & (m!=len(paths)):
                paths[m]= pathfile[comma[m-1]+1:comma[m]]
            elif m==len(paths):
                paths[m]= pathfile[comma[m]+1:]

        return paths




    def import_catalogue(self, pathfile, cols=(0,1,2,3,4,5,6,7,9)):
        """
        This function import single catalog file
        --------------------------
        pathfile : the entire path of the catalogue
        cols : the number of columns used for the data extraction
        """
        extension = pathfile[-1:-5:-1][::-1]
        if extension=='fits':
            type_input='fits'
        else:
            type_input='ascii'

        if type_input=='ascii':
            self.id_obj, self.ra, self.dec, self.l, self.b, self.Jmag, self.Hmag, self.Kmag, self.E_jk = np.genfromtxt(pathfile, usecols=cols, dtype='float', unpack=True)
        elif type_input=='fits':
            f           =   fits.open(pathfile)
            f_in        =   f[1].data
            self.id_obj =   f_in['SOURCEID'].astype(int)
            self.ra     =   f_in['RA2000']
            self.dec    =   f_in['DEC2000']
            self.l      =   f_in['L']
            self.b      =   f_in['B']
            self.Jmag   =   f_in['J']
            self.Hmag   =   f_in['H']
            self.Kmag   =   f_in['K']
            self.E_jk   =   f_in['EJK']
            self.Hmag[np.isnan(self.Hmag)]   =  -999
            self.Jmag[np.isnan(self.Jmag)]   =  -999
            self.Kmag[np.isnan(self.Kmag)]   =  -999
            self.E_jk[np.isnan(self.E_jk)]   =  -999

        McL.clean_for_duplicate(self)



    def clean_for_duplicate(self):
        """
        This function identifies duplicates objects in the catalogue
        with a double check: by using their ID (TODO and after also their
        coordinates).
        """
        newid_obj, i = np.unique(self.id_obj, return_index=True)
        if len(newid_obj) < len(self.id_obj):
            print()
            print("I'm cleaning the catalogue for the presence of duplicate objects (by their ID)")
            delta = len(self.id_obj) - len(newid_obj)
            print(len(self.id_obj),"->",len(newid_obj))
            print(str(delta),"objects have been deleted ")
            print()
            self.id_obj, self.ra, self.dec, self.l, self.b, self.Jmag, self.Hmag, self.Kmag, self.E_jk = self.id_obj[i], self.ra[i], self.dec[i], self.l[i], self.b[i], self.Jmag[i], self.Hmag[i], self.Kmag[i], self.E_jk[i]
        else:
            print()
            print("There are not duplicate objects in this catalogue (by their ID)")
            print()
        


    def read_input_namelist(self, pathfile="./",from_line="no", print_info="no"):
        """
        This function read and import the namelist parameters.
        --------------------------
        pathfile : the entire path of the catalogue
        from_line : "yes" if you put variables by line, "no" if you want to read
        catalogues by the namelist file
        print_info : if "yes", function will print some information about the
        catalogue
        """
        do = McL.write_input(self, pathfile)
        if do == 'yes':
            keyw = np.genfromtxt(pathfile, usecols=(0), dtype='str')
            val = np.genfromtxt(pathfile, usecols=(1), dtype='str')
            if from_line=="no":
                path_cat = val[0]
                print("Catalogue path :",path_cat)
                McL.import_catalogues(self, path_cat)
            elif from_line=="yes":
                print("Data catalogue imported from line command")
                #McL.import_catalogue_fromline(self, id_obj, ra, dec, l, b, Jmag, Hmag, Kmag, E_jk)
                pass


            f = val[1]
            comma = 0
            for n in range(len(f)):
                if f[n]==",":
                    comma = n
            x1           =       float(f[1:comma])
            y1           =       float(f[comma+1:-1])

            val[2] = float(val[2])
            val[4] = float(val[4])


            f = val[6]
            comma = 0
            if f!='None':
              for n in range(len(f)):
                if f[n]==",":
                    comma = n
              x6           =       float(f[1:comma])
              y6           =       float(f[comma+1:-1])
            else:
              x6           =       None
              y6           =       None

            f = val[7]
            comma = 0
            if f!='None':
              for n in range(len(f)):
                if f[n]==",":
                    comma = n
              x7           =       float(f[1:comma])
              y7           =       float(f[comma+1:-1])
            else:
              x7           =       None
              y7           =       None


            self.PATH                =       str(val[0])
            self.CENTRE_CUT          =       np.array([x1,y1])
            self.RADIUS_CUT          =       float(val[2])/60.
            self.MAG_BAND_CUT        =       str(val[3])
            self.MAG_VAL_CUT         =       float(val[4])
            self.Type                =       str(val[5])
            self.COL_VAL             =       np.array([x6,y6])
            self.MAG_VAL             =       np.array([x7,y7])
            self.N_TELL              =       np.int(val[8])
            self.COL_TELL            =       float(val[9])
            self.MAG_TELL            =       float(val[10])
            self.N_SKY              =       int(val[11])
            self.N_CELL              =       int(val[12])
            self.N_SKY_XCELL         =       float(val[13])
            self.SKYSIZE             =       float(val[14])
            self.RADIUS_DEBL         =       str(val[15])
            self.EXT_CORR            =       str(val[16])
            self.hasGuideStar        =       str(val[17])
            self.guideStarRA         =       float(val[18])
            self.guideStarDEC        =       float(val[19])
            self.orientationGP       =       str(val[20])
            self.date                =       str(val[21])
            self.time                =       str(val[22])
            self.airmass_limit       =       float(val[23])
            self.observingMode       =       str(val[24])
            self.noddingSize         =       float(val[25])
            self.noddingDirection    =       float(val[26])
            self.n_fibres_on_sky     =       float(val[27])
            self.bi_policy           =       float(val[28])
            self.bi_iterations       =       float(val[29])
            self.bi_max_priority     =       float(val[30])
            self.T_initial           =       float(val[31])
            self.T_final             =       float(val[32])
            self.T_delta             =       float(val[33])
            self.iterations          =       float(val[34])
            self.random_seed         =       float(val[35])
        elif do == 'no':
            pass


        if print_info!='no':
            [print(keyw[n], val[n]) for n in range(len(val))]

        #McL.get_info(self)


    def get_data(self, index=None):
        """
        This function extrat data.
        --------------------------
        index : if index are a vector, data are extract according to thei index
        if it is None, all data are extract
        """
        if index==None:
            id_obj    =     self.id_obj
            ra        =     self.ra
            dec       =     self.dec
            l         =     self.l
            b         =     self.b
            Jmag      =     self.Jmag
            Hmag      =     self.Hmag
            Kmag      =     self.Kmag
            E_jk      =     self.E_jk

            return id_obj, ra, dec, l, b, Jmag, Hmag, Kmag, E_jk
        else:
            id_obj    =     self.id_obj[index]
            ra        =     self.ra[index]
            dec       =     self.dec[index]
            l         =     self.l[index]
            b         =     self.b[index]
            Jmag      =     self.Jmag[index]
            Hmag      =     self.Hmag[index]
            Kmag      =     self.Kmag[index]
            E_jk      =     self.E_jk[index]

            return id_obj, ra, dec, l, b, Jmag, Hmag, Kmag, E_jk


    def spherical_cut(self, CENTRE_CUT=None, RADIUS_CUT=None, option_self="yes"):
        """
        This function cut data into a circle
        --------------------------
        CENTRE_CUT : coordinates of the centre
        RADIUS_CUT : radius of the circle
        option_self : if "yes" it update the self.data, otherwise it extraxts
        cutted data
        """
        if option_self=="no":
            if (CENTRE_CUT==None) & (RADIUS_CUT==None) :
                l_centre    =    self.CENTRE_CUT[0]
                b_centre    =    self.CENTRE_CUT[1]

                radius = np.sqrt( (self.l-l_centre)**2+(self.b-b_centre)**2 )

                i = np.where(radius <= self.RADIUS_CUT)
            else:
                l_centre    =    CENTRE_CUT[0]
                b_centre    =    CENTRE_CUT[1]

                radius = np.sqrt( (self.l-l_centre)**2+(self.b-b_centre)**2 )

                i = np.where(radius <= RADIUS_CUT)


            id_obj    =     self.id_obj[i]
            ra        =     self.ra[i]
            dec       =     self.dec[i]
            l         =     self.l[i]
            b         =     self.b[i]
            Jmag      =     self.Jmag[i]
            Hmag      =     self.Hmag[i]
            Kmag      =     self.Kmag[i]
            E_jk      =     self.E_jk[i]

            return id_obj, ra, dec, l, b, Jmag, Hmag, Kmag, E_jk

        elif option_self=="yes":
            if (CENTRE_CUT==None) & (RADIUS_CUT==None) :
                l_centre    =    self.CENTRE_CUT[0]
                b_centre    =    self.CENTRE_CUT[1]

                radius = np.sqrt( (self.l-l_centre)**2+(self.b-b_centre)**2 )

                i = np.where(radius <= self.RADIUS_CUT)
            else:
                l_centre    =    CENTRE_CUT[0]
                b_centre    =    CENTRE_CUT[1]

                radius = np.sqrt( (self.l-l_centre)**2+(self.b-b_centre)**2 )

                i = np.where(radius <= RADIUS_CUT)


            self.id_obj    =     self.id_obj[i]
            self.ra        =     self.ra[i]
            self.dec       =     self.dec[i]
            self.l         =     self.l[i]
            self.b         =     self.b[i]
            self.Jmag      =     self.Jmag[i]
            self.Hmag      =     self.Hmag[i]
            self.Kmag      =     self.Kmag[i]
            self.E_jk      =     self.E_jk[i]


    def spherical_cut_radec(self, CENTRE_CUT=None, RADIUS_CUT=None, option_self="yes"):
        """
        This function cut data into a circle
        --------------------------
        CENTRE_CUT : coordinates of the centre
        RADIUS_CUT : radius of the circle
        option_self : if "yes" it update the self.data, otherwise it extraxts
        cutted data
        """
        if option_self=="no":
            if (CENTRE_CUT==None) & (RADIUS_CUT==None) :
                ra_centre    =    self.CENTRE_CUT[0]
                dec_centre    =    self.CENTRE_CUT[1]

                radius = np.sqrt( (self.ra-ra_centre)**2+(self.dec-dec_centre)**2 )

                i = np.where(radius <= self.RADIUS_CUT)
            else:
                ra_centre    =    CENTRE_CUT[0]
                dec_centre    =    CENTRE_CUT[1]

                radius = np.sqrt( (self.ra-ra_centre)**2+(self.dec-dec_centre)**2 )

                i = np.where(radius <= RADIUS_CUT)


            id_obj    =     self.id_obj[i]
            ra        =     self.ra[i]
            dec       =     self.dec[i]
            l         =     self.l[i]
            b         =     self.b[i]
            Jmag      =     self.Jmag[i]
            Hmag      =     self.Hmag[i]
            Kmag      =     self.Kmag[i]
            E_jk      =     self.E_jk[i]

            return id_obj, ra, dec, l, b, Jmag, Hmag, Kmag, E_jk

        elif option_self=="yes":
            if (CENTRE_CUT==None) & (RADIUS_CUT==None) :
                ra_centre    =    self.CENTRE_CUT[0]
                dec_centre   =    self.CENTRE_CUT[1]

                radius = np.sqrt( (self.ra-ra_centre)**2+(self.dec-dec_centre)**2 )

                i = np.where(radius <= self.RADIUS_CUT)
            else:
                ra_centre    =    CENTRE_CUT[0]
                dec_centre   =    CENTRE_CUT[1]

                radius = np.sqrt( (self.ra-ra_centre)**2+(self.dec-dec_centre)**2 )

                i = np.where(radius <= RADIUS_CUT)

            self.id_obj    =     self.id_obj[i]
            self.ra        =     self.ra[i]
            self.dec       =     self.dec[i]
            self.l         =     self.l[i]
            self.b         =     self.b[i]
            self.Jmag      =     self.Jmag[i]
            self.Hmag      =     self.Hmag[i]
            self.Kmag      =     self.Kmag[i]
            self.E_jk      =     self.E_jk[i]




    def select_tellurics(self):
        i = np.where((self.Hmag<self.MAG_TELL) & ((self.Jmag-self.Hmag)<self.COL_TELL))
        
        return list(i[0][:int(self.N_TELL)]), i





    def mag_cut(self, MAG_BAND_CUT=None, MAG_VAL_CUT=None, option_self="yes"):
        """
        This function cut data under a threshold magnitude
        --------------------------
        MAG_BAND_CUT : magnitude band (ex. "J", "H" or "K")
        MAG_VAL_CUT : values of threshold magnitude
        option_self : if "yes" it update the self.data, otherwise it extraxts
        cutted data
        """
        if option_self=="no":
            if (MAG_BAND_CUT==None) & (MAG_VAL_CUT==None):
                if self.MAG_BAND_CUT=='J':
                    mag = self.Jmag
                elif self.MAG_BAND_CUT=='H':
                    mag = self.Hmag
                elif self.MAG_BAND_CUT=='K':
                    mag = self.Kmag

                i = np.where(mag <= self.MAG_VAL_CUT)
            else:
                if MAG_BAND_CUT=='J':
                    mag = self.Jmag
                elif MAG_BAND_CUT=='H':
                    mag = self.Hmag
                elif (MAG_BAND_CUT=='K') or (MAG_BAND_CUT=='Ks'):
                    mag = self.Kmag

                i = np.where(mag <= self.MAG_VAL_CUT)


            id_obj    =     self.id_obj[i]
            ra        =     self.ra[i]
            dec       =     self.dec[i]
            l         =     self.l[i]
            b         =     self.b[i]
            Jmag      =     self.Jmag[i]
            Hmag      =     self.Hmag[i]
            Kmag      =     self.Kmag[i]
            E_jk      =     self.E_jk[i]

            return id_obj, ra, dec, l, b, Jmag, Hmag, Kmag, E_jk

        elif option_self=="yes":

            if (MAG_BAND_CUT==None) & (MAG_VAL_CUT==None):
                if self.MAG_BAND_CUT=='J':
                    mag = self.Jmag
                elif self.MAG_BAND_CUT=='H':
                    mag = self.Hmag
                elif self.MAG_BAND_CUT=='K':
                    mag = self.Kmag

                i = np.where(mag <= self.MAG_VAL_CUT)
            else:
                if MAG_BAND_CUT=='J':
                    mag = self.Jmag
                elif MAG_BAND_CUT=='H':
                    mag = self.Hmag
                elif MAG_BAND_CUT=='K':
                    mag = self.Kmag

                i = np.where(mag <= self.MAG_VAL_CUT)



            self.id_obj    =     self.id_obj[i]
            self.ra        =     self.ra[i]
            self.dec       =     self.dec[i]
            self.l         =     self.l[i]
            self.b         =     self.b[i]
            self.Jmag      =     self.Jmag[i]
            self.Hmag      =     self.Hmag[i]
            self.Kmag      =     self.Kmag[i]
            self.E_jk      =     self.E_jk[i]





    def get_info(self, return_dict=False):
        """
        This function print info of the imported catalogue(s).
        --------------------------
        return_dict : if "True" it extracts some information via
        python dictionary
        """
        raFIELD   =   np.mean(self.ra)
        decFIELD  =   np.mean(self.dec)
        raCENTER  =   np.min(self.ra) + (np.max(self.ra) - np.min(self.ra))/2.
        decCENTER =   np.min(self.dec) + (np.max(self.dec) - np.min(self.dec))/2.
        print()
        print("Number of objects", len(self.ra))
        print("FoV center RA", raFIELD)
        print("FoV center DEC", decFIELD)
        print("FIELD Ra",raCENTER)
        print("FIELD Dec",decCENTER)
        print()

        self.raFIELD        =       raFIELD
        self.decFIELD       =       decFIELD
        self.raCENTER        =       raCENTER
        self.decCENTER       =       decCENTER


        dict_info  =    {"FoV Center"   :  (raFIELD, decFIELD),"Field center" :  (raCENTER, decCENTER)}
        if return_dict==True:
            return dict_info




    def extinction_correction(self):
        """
        This function computes the de-reddened magnitudes.
        """
        # Extinction coefficient of Nishiyama et al. :

        self.Aj    =  1.526 * self.E_jk
        self.Ah    =  0.855 * self.E_jk
        self.Aks   =  0.528 * self.E_jk

        J0    =  self.Jmag - self.Aj
        H0    =  self.Hmag - self.Ah
        K0    =  self.Kmag - self.Aks

        self.J0mag    =   J0
        self.H0mag    =   H0
        self.K0mag    =   K0

        return self.J0mag, self.H0mag, self.K0mag






    def select_stars_fromCMD(self, COL="J-H", MAG="H", COL_VAL=None, MAG_VAL=None, return_index=False, print_columns="no"):
        """
        This function select stars from CMD. If more than one selection, use inline
        mode (put None in namelist)
        --------------------------
        COL             : color used
        MAG             : mag used
        COL_VAL         : color values [col_min, col_max]
        MAG_VAL         : mag values [mag_min, mag_max]
        return_index    : "yes" if you want only the index of data
        print_columns   : "yes" if you want informations about columns extracted
        """
        
        if ( (COL=="J0-K0") & (MAG=="K0") & (COL_VAL!=None) & (MAG_VAL!=None)):
            i = np.where( (self.J0mag-self.K0mag<=np.max(COL_VAL)) &
            (self.J0mag-self.K0mag>=np.min(COL_VAL)) &
            (self.K0mag<=np.max(MAG_VAL)) & (self.K0mag>=np.min(MAG_VAL)) )
        elif ( (COL=="J0-K0") & (MAG=="K0") & (COL_VAL==None) & (MAG_VAL==None)):
            i = np.where( (self.J0mag-self.K0mag<=np.max(self.COL_VAL)) &
            (self.J0mag-self.K0mag>=np.min(self.COL_VAL)) &
            (self.K0mag<=np.max(self.MAG_VAL)) & (self.K0mag>=np.min(self.MAG_VAL)) )

        elif ( (COL=="J-K") & (MAG=="K") & (COL_VAL!=None) & (MAG_VAL!=None)):
            i = np.where( (self.Jmag-self.Kmag<=np.max(COL_VAL)) &
            (self.Jmag-self.Kmag>=np.min(COL_VAL)) &
            (self.Kmag<=np.max(MAG_VAL)) & (self.Kmag>=np.min(MAG_VAL)) )
            j = np.where( ((self.Jmag-self.Kmag>np.max(COL_VAL)) |
            (self.Jmag-self.Kmag<np.min(COL_VAL)) |
            (self.Kmag>np.max(MAG_VAL)) | (self.Kmag<np.min(MAG_VAL))))
        elif ( (COL=="J-K") & (MAG=="K") & (COL_VAL==None) & (MAG_VAL==None)):
            i = np.where( (self.Jmag-self.Kmag<=np.max(self.COL_VAL)) &
            (self.Jmag-self.Kmag>=np.min(self.COL_VAL)) &
            (self.Kmag<=np.max(self.MAG_VAL)) & (self.Kmag>=np.min(self.MAG_VAL)) )

        elif ( (COL=="J0-H0") & (MAG=="H0") & (COL_VAL!=None) & (MAG_VAL!=None)):
            i = np.where( (self.J0mag-self.H0mag<=np.max(COL_VAL)) &
            (self.J0mag-self.H0mag>=np.min(COL_VAL)) &
            (self.H0mag<=np.max(MAG_VAL)) & (self.H0mag>=np.min(MAG_VAL)) )
        elif ( (COL=="J0-H0") & (MAG=="H0") & (COL_VAL==None) & (MAG_VAL==None)):
            i = np.where( (self.J0mag-self.H0mag<=np.max(self.COL_VAL)) &
            (self.J0mag-self.H0mag>=np.min(self.COL_VAL)) &
            (self.H0mag<=np.max(self.MAG_VAL)) & (self.H0mag>=np.min(self.MAG_VAL)) )

        elif ( (COL=="J-H") & (MAG=="H") & (COL_VAL!=None) & (MAG_VAL!=None)):
            i = np.where( (self.Jmag-self.Hmag<=np.max(COL_VAL)) &
            (self.Jmag-self.Hmag>=np.min(COL_VAL)) &
            (self.Hmag<=np.max(MAG_VAL)) & (self.Hmag>=np.min(MAG_VAL)) )
        elif ( (COL=="J-H") & (MAG=="H") & (COL_VAL==None) & (MAG_VAL==None)):
            i = np.where( (self.Jmag-self.Hmag<=np.max(self.COL_VAL)) &
            (self.Jmag-self.Hmag>=np.min(self.COL_VAL)) &
            (self.Hmag<=np.max(self.MAG_VAL)) & (self.Hmag>=np.min(self.MAG_VAL)) )


        if return_index==False:
            id_obj_cut    =     self.id_obj[i]
            ra_cut        =     self.ra[i]
            dec_cut       =     self.dec[i]
            l_cut         =     self.l[i]
            b_cut         =     self.b[i]
            Jmag_cut      =     self.Jmag[i]
            Hmag_cut      =     self.Hmag[i]
            Kmag_cut      =     self.Kmag[i]
            J0mag_cut     =     self.J0mag[i]
            H0mag_cut     =     self.H0mag[i]
            K0mag_cut     =     self.K0mag[i]
            E_jk_cut      =     self.E_jk[i]

            if ( (print_columns!="no")):
                print("Columns are: id_obj, ra, dec, l, b, Jmag, Hmag, Kmag, J0mag, H0mag, K0mag, E_jk")

            return id_obj_cut, ra_cut, dec_cut, l_cut, b_cut, Jmag_cut, Hmag_cut, Kmag_cut, J0mag_cut, H0mag_cut, K0mag_cut, E_jk_cut
        elif return_index==True:

            if ( (print_columns!="no")):
                print("Columns are: i (index of Target) and j (index acquisition)")

            return i




    def write_MOONLIGHT_param(self, savepath="./", cat_name=None, param_name=None):
        """
        This function select stars from CMD. If more than one selection, use inline
        mode (put None in namelist)
        --------------------------
        savepath       : path in which you want the innput parameter file of MOONLIGHT
        cat_name        : name of the catalogue which you will use for MOONLIGHT
        param_name      : name of the file without extention
        """
        #if ((self.guideStarRA==-999) & (self.guideStarDEC==-999)):
        #    self.guideStarRA, self.guideStarDEC = McL.find_guidestar(self)



        file=open(savepath+param_name,"w")
        file.write("#******************************************************************************    \n".format())
        file.write("#E.S.O. - VLT project                                                             \n".format())
        file.write("#                                                                                 \n".format())
        file.write("#An example Observation Preparation Software Parameter File.                      \n".format())
        file.write("#                                                                                 \n".format())
        file.write("#******************************************************************************    \n".format())
        file.write("                                                                                   \n".format())
        file.write("#Catalogue information                                                            \n".format())
        file.write("#                                                                                 \n".format())
        file.write("catalogueFilename  {}      #input source catalogue                          \n".format(str(cat_name)))
        file.write("                                                                                   \n".format())
        file.write("outputDirectory {}    #directory where OPS output will be written            \n".format(str(savepath)))
        file.write("                                                                                   \n".format())
        file.write("#THESE WILL BE PASSED ON FROM OB                                                   \n".format())
        file.write("                                                                                   \n".format())
        file.write("#************* Acquisition template parameters                                    \n".format())
        file.write("#                                                                                 \n".format())
        file.write("#Pointing  information  \n".format())
        file.write("#                                                                                 \n".format())
        file.write("FoVcenterRA         {:=6.2f}  #RA  pointing center in decimal degrees     \n".format(float(self.raFIELD)))
        file.write("FoVcenterDEC        {:=6.2f}  #DEC pointing center in decimal degrees     \n".format(float(self.decFIELD)))
        file.write("#                                                                                 \n ".format())
        file.write("#VLT guide star                                                                   \n".format())
        file.write("#                                                                                 \n ".format())
        file.write("hasGuideStar   {}      #if false the guide probe will not be considered (true or false)   \n".format(self.hasGuideStar))
        file.write("guideStarRA    {:=6.2f}     #RA of the VLT guide star in decimal degrees  \n".format(self.guideStarRA))
        file.write("guideStarDEC   {:=6.2f}     #DEC of the VLT guide star in decimal degrees \n".format(self.guideStarDEC))
        file.write("orientationGP  {}        #guide probe orientation (POS or NEG)            \n".format(self.orientationGP))
        file.write("                                                                                   \n".format())
        file.write("#******** Constraint set   \n".format())
        file.write("#                                                                                \n ".format())
        file.write("#Precession parameters ***** These will be substituted by hour angle range        \n".format())
        file.write("#                                                                                 \n".format())
        file.write("date          {}                      #Observing date, format mm:dd:yyyy   \n".format(str(self.date)))
        file.write("time          {}      #Observing time, format hh:mm:ss                     \n".format(str(self.time)))
        file.write("airmass_limit          {:=4.2f}      #Maximum airmass limit                \n".format(self.airmass_limit))
        file.write("                                                                                   \n".format())
        file.write("                                                                                   \n".format())
        file.write("#********  Observing template parameters                                         \n".format())
        file.write("observingMode     {}   #observing strategy: STARE, XSWITCH, STAREandNOD   \n".format(str(self.observingMode)))
        file.write("noddingSize       {:=4.1f}      #nodding size from 10 to 30 arcsec         \n".format(self.noddingSize))
        file.write("noddingDirection  {:=4.1f}      #nodding direction from 0 to 360 degrees   \n".format(self.noddingDirection))
        file.write("                                                                                   \n".format())
        file.write("                                                                                   \n".format())
        file.write("                                                                                   \n".format())
        file.write("                                                                                   \n".format())
        file.write("                                                                                   \n".format())
        file.write("# ********  MoonLight SPECIFIC   \n".format())
        file.write("n_fibres_on_sky  {:=3.0f}  #number of fibres dedicated to sky (STARE mode only) \n".format(self.n_fibres_on_sky))
        file.write("                                                                                   \n".format())
        file.write("#                                                                                 \n".format())
        file.write("#Fibre allocation strategy   \n".format())
        file.write("#first a burn-in (BI) with top priority targets then simulated annealing (SA)   \n".format())
        file.write("#during BI targets are allocated in order of increasing priority   \n".format())
        file.write("#if frozen, fibres allocated during BI are untouched in SA   \n".format())
        file.write("bi_policy         {:=2.0f}         #0 no burn-in; 1 burn-in; 2 burn-in and freeze   \n".format(self.bi_policy))
        file.write("bi_iterations     {:=2.0f}         #Number of iterations during burn-in   \n".format(self.bi_iterations))
        file.write("bi_max_priority   {:=2.0f}         #Max (lowest) priority for burn-in   \n".format(self.bi_max_priority))
        file.write("T_initial         {:=4.2f}       #Initial temperature of the system   \n".format(self.T_initial))
        file.write("T_final           {:=4.2f}       #Final temperature of the system   \n".format(self.T_final))
        file.write("T_delta           {:=4.2f}       #Variation in temperature   \n".format(self.T_delta))
        file.write("iterations        {:=2.0f}         #Iterations to reach thermal equilibrium   \n".format(self.iterations))
        file.write("random_seed       {:=2.0f}         #Random seed for different fibre configurations   \n".format(self.random_seed))
        file.write("   \n".format())
        file.write("   \n".format())
        file.write("   \n".format())
        file.close()

        print("Printed ->", param_name)




    def make_priority_and_flag0(self, i, prio_i, flag):
        """
        This function gives priority and flag
        --------------------------
        i, prio_i, flag       : index, priority and flag value for stars with this index
        """
        IDw, Raw, Decw = [], [], []
        flags, priority, mag_input = [], [], []
        for n in i[0]:
            IDw.append(self.id_obj[n])
            Raw.append(self.ra[n])
            Decw.append(self.dec[n])
            flags.append(flag)
            priority.append(prio_i)
            mag_input.append(self.Jmag[n])

        IDw, Raw, Decw, flag, priority, mag_input = np.asarray(IDw), np.asarray(Raw), np.asarray(Decw), np.asarray(flag), np.asarray(priority), np.asarray(mag_input)
        o = np.argsort(IDw)

        return IDw[o], Raw[o], Decw[o], flag[o], priority[o], mag_input[o]


    def make_priority_and_flag(self, id_obj, ra, dec, mag_input, i, prio_i, flag):
        """
        This function gives priority and flag
        --------------------------
        i, prio_i, flag       : index, priority and flag value for stars with this index
        """
        IDw, Raw, Decw = [], [], []
        flags, priority, mag = [], [], []

        for n in i[0]:
            IDw.append(id_obj[n])
            Raw.append(ra[n])
            Decw.append(dec[n])
            flags.append(flag)
            priority.append(prio_i)
            mag.append(mag_input[n])

        IDw, Raw, Decw, flag, priority, mag = np.asarray(IDw), np.asarray(Raw), np.asarray(Decw), np.asarray(flag), np.asarray(priority), np.asarray(mag)
        #o = np.argsort(IDw)

        return IDw, Raw, Decw, flag, priority, mag




    def merge_catalogues(self, ID1, Ra1, Dec1, flag1, priority1, mag_input1, ID2, Ra2, Dec2, flag2, priority2, mag_input2):
        """
        This function merge 2 MOONLIGHT catalogues from line
        """
        IDw = np.concatenate((ID1,ID2))
        Raw = np.concatenate((Ra1,Ra2))
        Decw = np.concatenate((Dec1,Dec2))
        flagw = np.concatenate((flag1,flag2))
        priorityw = np.concatenate((priority1,priority2))
        mag_inputw = np.concatenate((mag_input1,mag_input2))


        _, i = np.unique(IDw, return_index=True)
        IDw, Raw, Decw, flagw, priorityw, mag_inputw = IDw[i], Raw[i], Decw[i], flagw[i], priorityw[i], mag_inputw[i]

        return IDw, Raw, Decw, flagw, priorityw, mag_inputw



    def concatenate_catalogues2(vec1,vec2):
        """
        TODO
        """

        if len(vec1)==len(vec2):
            vec_conc = [np.concatenate((vec1[n],vec2[n])) for n in range(len(vec1))]

            return vec_conc
        else:
            print("arrays must be equal lenght")






    def write_MOONLIGHT_input(self, ID, Ra, Dec, flag, priority, mag, savepath="./", cat_name=None):
        """
        This function write the MOONLIGHT input catalogue from line
        --------------------------
        ID          :   id objects
        Ra          :
        Dec         :
        flag        :   target, acquiisition, sky
        priority    :
        mag         :   J or r magnitude
        savepath   :
        cat_name    :
        """
        _, i = np.unique(ID, return_index=True)
        file=open(savepath+cat_name,"w")
        file.write("# ID          RA        DEC        flag   priority   mag  \n".format())
        for val in zip(ID[i], Ra[i], Dec[i], flag[i], priority[i], mag[i]):
                file.write('  {:10}   {:=9.5f}   {:=9.5f}   {:1}   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4],val[5]))
        file.close()

        print("Printed ->", cat_name)


    def write_input(self, pathfile):
        """
        This function write an example of the McL namelist input parameters
        if it doesn't exist
        """
        try:
            fi = open(pathfile)
            fi.close()
            print('Namelist found...')
            do = 'yes'
        except:
            print('Namelist do not Exist. I write here an example for you...')
            do = 'no'



            file=open("./listname_example","w")
            file.write("#######################################################   \n".format())
            file.write("#             Select Stars for MOONS - GC             #   \n".format())
            file.write("#                 Cristiano Fanelli                   #   \n".format())
            file.write("#                   0.3.0  Version                    #   \n".format())
            file.write("#                   29 - 04 - 2022                    #   \n".format())
            file.write("#######################################################   \n".format())
            file.write("   \n".format())
            file.write("   \n".format())
            file.write("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   \n".format())
            file.write("# Comments:   \n".format())
            file.write("# For V-0.0.1 the files (both ascii and fits) must be arranged   \n".format())
            file.write("# in columns as: *"'ID, ra, dec, l, b, J, H, Ks, E(J-K)'"   \n".format())
            file.write("# TODO: Merging of more than 1 catalogue at the starting point   \n".format())
            file.write("# TODO: automatic recognition of columns file (maybe? Is it useful?)   \n".format())
            file.write("# of columns.   \n".format())
            file.write("# 25 - 10 - 2021        	Cristiano Fanelli   \n".format())
            file.write("# You can use more than 1 catalogue (Max 3 for now)   \n".format())
            file.write("# TODO: generalize for N catalogue, with N>3   \n".format())
            file.write("# 01 - 11 - 2021        	Cristiano Fanelli   \n".format())
            file.write("#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   \n".format())
            file.write("   \n".format())
            file.write("   \n".format())
            file.write("#The full path to the file name, with the extension (.fits, .dat, .txt, ...)   \n".format())
            file.write("#if more than one catalogue, adding it as a list path1,path2,...,pathN   \n".format())
            file.write("PATH		/Users/cristiano/phd/NICE/catalogues/333_d.fits,/Users/cristiano/phd/NICE/catalogues/334_d.fits   \n".format())
            file.write("   \n".format())
            file.write("#For spherical cut   \n".format())
            file.write("CENTRE_CUT	 (0,0)		  # FLOAT (X0,Y0)   \n".format())
            file.write("RADIUS_CUT	 0.10		  # FLOAT   \n".format())
            file.write("   \n".format())
            file.write("#For magnitude cut in J, H or K band   \n".format())
            file.write("MAG_BAND_CUT	H		  # STRING   \n".format())
            file.write("MAG_VAL_CUT	    17		  # FLOAT   \n".format())
            file.write("   \n".format())
            file.write("#For selecting a type of stars into a square   \n".format())
            file.write("#it gives values as in *Comments   \n".format())
            file.write("Type		RC		  # STRING   \n".format())
            file.write("J-K		   None		  # [(J-K)0, (J-K)1]   \n".format())
            file.write("K		   None		  # [K0, K1]   \n".format())
            file.write("   \n".format())
            file.write("   \n".format())
            file.write("   \n".format())
            file.write("#For the parameter file   \n".format())
            file.write("   \n".format())
            file.write("hasGuideStar          false       # if false the guide probe will not be considered (true or false)   \n".format())
            file.write("guideStarRA           266.41      # RA of the VLT guide star in decimal degrees   \n".format())
            file.write("guideStarDEC         -28.94       # DEC of the VLT guide star in decimal degrees   \n".format())
            file.write("orientationGP         POS         # guide probe orientation (POS or NEG)   \n".format())
            file.write("   \n".format())
            file.write("date                  06:11:2020  # Observing date, format mm:dd:yy   \n".format())
            file.write("time                  22:00:00    # Observing time, format hh:mm:ss   \n".format())
            file.write("airmass_limit         6           # Maximum airmass limit   \n".format())
            file.write("   \n".format())
            file.write("observingMode         XSWITCH     # observing strategy: STARE, XSWITCH, STAREandNOD   \n".format())
            file.write("noddingSize           30.0        # nodding size from 10 to 30 arcsec   \n".format())
            file.write("noddingDirection      90.0        # nodding direction from 0 to 360 degrees   \n".format())
            file.write("   \n".format())
            file.write("n_fibres_on_sky       100         # number of fibres dedicated to sky (STARE mode only)   \n".format())
            file.write("   \n".format())
            file.write("bi_policy             2           # 0 no burn-in; 1 burn-in; 2 burn-in and freeze   \n".format())
            file.write("bi_iterations         5           # Number of iterations during burn-in   \n".format())
            file.write("bi_max_priority       2           # Max (lowest) priority for burn-in   \n".format())
            file.write("T_initial             0.10        # Initial temperature of the system   \n".format())
            file.write("T_final               0.00        # Final temperature of the system   \n".format())
            file.write("T_delta               0.02        # Variation in temperature   \n".format())
            file.write("iterations            10          # Iterations to reach thermal equilibrium   \n".format())
            file.write("random_seed           0           # Random seed for different fibre configurations   \n".format())

            file.close()

        return do



    def plot_extinction(self, Ai = "Ah", namefig = None, savepath= './'):
        """
        This function shows an extintion map
        --------------------------
        Ai          : Ah, Aj or Ak
        namefig     : if a name is provided, then save the
        figure (inpdf)
        """
        AA='s'
        if Ai == "Aj":
            AA=self.Aj
        if Ai == "Ah":
            AA=self.Ah
        if (Ai == "Ak") or (Ai == "Aks"):
            AA=self.Aks

        fig=plt.figure(figsize=(9,7), dpi=300)
        #fig.subplots_adjust(left=0.1,bottom=0.12,right=0.95,top=0.98,hspace=0.24,wspace=0.28)
        ax =  fig.add_subplot(111)
        ax.set_xlabel("l", fontsize=20)
        ax.set_ylabel("b", fontsize=20)
        cm = plt.cm.get_cmap('jet')
        sc = ax.scatter(self.l,self.b,edgecolor="none", s=3, c=AA,cmap=cm, label=" "%np.sum(AA))
        cb=plt.colorbar(sc)
        cb.set_label(str(Ai[0])+"$_"+str(Ai[1])+"$", fontsize=20)
        fig.suptitle("    ", fontsize=20)
        fig.tight_layout()
        fig.savefig(savepath+namefig,dpi=120)
        plt.show()

        plt.close(fig)


    def plot_extinction_fibers(self, l, b, Ai = "Ah", namefig = None):
        """
        This function shows an extintion map
        --------------------------
        Ai          : Ah, Aj or Ak
        namefig     : if a name is provided, then save the
        figure (inpdf)
        """
        AA='s'
        if Ai == "Aj":
            AA=self.Aj
        if Ai == "Ah":
            AA=self.Ah
        if (Ai == "Ak") or (Ai == "Aks"):
            AA=self.Aks

        fig=plt.figure(figsize=(9,7), dpi=250)
        #fig.subplots_adjust(left=0.1,bottom=0.12,right=0.95,top=0.98,hspace=0.24,wspace=0.28)
        ax =  fig.add_subplot(111)
        ax.set_xlabel("l", fontsize=20)
        ax.set_ylabel("b", fontsize=20)
        cm = plt.cm.get_cmap('jet')
        sc = ax.scatter(self.l,self.b,edgecolor="none", s=3, c=AA,cmap=cm, label=" "%np.sum(AA))
        cb=plt.colorbar(sc)
        cb.set_label(str(Ai[0])+"$_"+str(Ai[1])+"$", fontsize=20)
        ax.scatter(l, b, color='k',zorder=5)
        fig.suptitle("    ", fontsize=20)
        fig.tight_layout()
        if namefig!=None:
            fig.savefig(namefig,dpi=120)
        else:
            plt.show()

        plt.close(fig)





    def get_data_dered(self, index=None):
        """
        This function extrat data (with alsothe de-renneded magnitudes)
        --------------------------
        index : if index are a vector, data are extract according to thei index
        if it is None, all data are extract
        """
        if (index==None) & (self.EXT_CORR=='yes'):
            id_obj    =     self.id_obj
            ra        =     self.ra
            dec       =     self.dec
            l         =     self.l
            b         =     self.b
            Jmag      =     self.Jmag
            Hmag      =     self.Hmag
            Kmag      =     self.Kmag
            J0mag     =     self.J0mag
            H0mag     =     self.H0mag
            K0mag     =     self.K0mag
            E_jk      =     self.E_jk

            return id_obj, ra, dec, l, b, Jmag, Hmag, Kmag, J0mag, H0mag, K0mag, E_jk

        elif (index==None) & (self.EXT_CORR=='yes'):
            id_obj    =     self.id_obj[index]
            ra        =     self.ra[index]
            dec       =     self.dec[index]
            l         =     self.l[index]
            b         =     self.b[index]
            Jmag      =     self.Jmag[index]
            Hmag      =     self.Hmag[index]
            Kmag      =     self.Kmag[index]
            J0mag     =     self.J0mag[index]
            H0mag     =     self.H0mag[index]
            K0mag     =     self.K0mag[index]
            E_jk      =     self.E_jk[index]
            return id_obj, ra, dec, l, b, Jmag, Hmag, Kmag, J0mag, H0mag, K0mag, E_jk


        if (index==None) & (self.EXT_CORR!='yes'):
            id_obj    =     self.id_obj
            ra        =     self.ra
            dec       =     self.dec
            l         =     self.l
            b         =     self.b
            Jmag      =     self.Jmag
            Hmag      =     self.Hmag
            Kmag      =     self.Kmag
            J0mag     =     -99*np.ones((len(id_obj)))
            H0mag     =     -99*np.ones((len(id_obj)))
            K0mag     =     -99*np.ones((len(id_obj)))
            E_jk      =     self.E_jk

            return id_obj, ra, dec, l, b, Jmag, Hmag, Kmag, J0mag, H0mag, K0mag, E_jk

        elif (index==None) & (self.EXT_CORR!='yes'):
            id_obj    =     self.id_obj[index]
            ra        =     self.ra[index]
            dec       =     self.dec[index]
            l         =     self.l[index]
            b         =     self.b[index]
            Jmag      =     self.Jmag[index]
            Hmag      =     self.Hmag[index]
            Kmag      =     self.Kmag[index]
            J0mag     =     -99*np.ones((len(self.id_obj[index])))
            H0mag     =     -99*np.ones((len(self.id_obj[index])))
            K0mag     =     -99*np.ones((len(self.id_obj[index])))
            E_jk      =     self.E_jk[index]
            return id_obj, ra, dec, l, b, Jmag, Hmag, Kmag, J0mag, H0mag, K0mag, E_jk




    def crossmatch(self, X1, X2, max_distance=np.inf):
            """Cross-match the values between X1 and X2

            By default, this uses a KD Tree for speed.

            Parameters
            ----------
            X1 : array_like
            first dataset, shape(N1, D)
            X2 : array_like
            second dataset, shape(N2, D)
            max_distance : float (optional)
            maximum radius of search.  If no point is within the given radius,
            then inf will be returned.

            Returns
            -------
            dist, ind: ndarrays
            The distance and index of the closest point in X2 to each point in X1
            Both arrays are length N1.
            Locations with no match are indicated by
            dist[i] = inf, ind[i] = N2
            """
            X1 = np.asarray(X1, dtype=float)
            X2 = np.asarray(X2, dtype=float)

            N1, D = X1.shape
            N2, D2 = X2.shape

            if D != D2:
                raise ValueError('Arrays must have the same second dimension')

            kdt = cKDTree(X2)

            dist, ind = kdt.query(X1, k=1, distance_upper_bound=max_distance, workers=-1)


            return dist, ind


    def crossmatch_angular(self, X1, X2, max_distance=np.inf):
            """Cross-match angular values between X1 and X2

            by default, this uses a KD Tree for speed.  Because the
            KD Tree only handles cartesian distances, the angles
            are projected onto a 3D sphere.

            Parameters
            ----------
            X1 : array_like
            first dataset, shape(N1, 2). X1[:, 0] is the RA, X1[:, 1] is the DEC,
            both measured in degrees
            X2 : array_like
            second dataset, shape(N2, 2). X2[:, 0] is the RA, X2[:, 1] is the DEC,
            both measured in degrees
            max_distance : float (optional)
            maximum radius of search, measured in degrees.
            If no point is within the given radius, then inf will be returned.

            Returns
            -------
            dist, ind: ndarrays
            The angular distance and index of the closest point in X2 to
            each point in X1.  Both arrays are length N1.
            Locations with no match are indicated by
            dist[i] = inf, ind[i] = N2
            """
            X1 = X1 * (np.pi / 180.)
            X2 = X2 * (np.pi / 180.)
            max_distance = max_distance * (np.pi / 180.)

            # Convert 2D RA/DEC to 3D cartesian coordinates
            Y1 = np.transpose(np.vstack([np.cos(X1[0]) * np.cos(X1[1]),
                                 np.sin(X1[0]) * np.cos(X1[1]),
                                 np.sin(X1[1])]))
            Y2 = np.transpose(np.vstack([np.cos(X2[:, 0]) * np.cos(X2[:, 1]),
                                 np.sin(X2[:, 0]) * np.cos(X2[:, 1]),
                                 np.sin(X2[:, 1])]))

            # law of cosines to compute 3D distance
            max_y = np.sqrt(2 - 2 * np.cos(max_distance))
            dist, ind = McL.crossmatch(self, Y1, Y2, max_y)

            # convert distances back to angles using the law of tangents
            not_inf = ~np.isinf(dist)
            x = 0.5 * dist[not_inf]
            dist[not_inf] = (180. / np.pi * 2 * np.arctan2(x,
                                 np.sqrt(np.maximum(0, 1 - x ** 2))))

            j = np.where(dist!=0.)

            return dist[j], ind[j]





    def single_malt_star(self, ID_all, ra_all, dec_all, mag_all,  ID, ra, dec, mag, rr=None, keyw=['a0','a1','a2','a3']):
        plt.hist(mag_all[mag_all<90])
        plt.show()
        plt.close()

        print("Searching in a subsample of "+str(len(ID))+" stars over the full sample of "+str(len(ID_all))+"...")
        if rr==None:
            rr = float(self.RADIUS_DEBL)
        rr = rr/3600.
        pts_all = np.array([ra_all, dec_all]).T
        pts = np.array([ra, dec]).T
        max_distance = rr

        Tree_cat = KDTree(pts_all)
        Tree_target = KDTree(pts)
        idx = Tree_target.query_ball_tree(Tree_cat, r=max_distance)


        bl = np.zeros((len(idx)))
        for n in range(len(idx)):
            n_blend = len(ID_all[idx[n]])
            if n_blend>1:
                for m in range(n_blend):
                    if (ID_all[idx[n]][m]!=ID[n]) & (mag_all[idx[n]][m]<=mag[n]):
                        bl[n] = 3
                        #print('3', mag_all[idx[n]][m],'<',mag[n])
                    elif (ID_all[idx[n]][m]!=ID[n]) & ((mag_all[idx[n]][m]>mag[n]) & (mag_all[idx[n]][m]<(mag[n]+2.5))):
                        bl[n] = 2
                        #print('2', mag_all[idx[n]][m],'...',mag[n]+2.5)
                    elif (ID_all[idx[n]][m]!=ID[n]) & (mag_all[idx[n]][m]>=(mag[n]+2.5)):
                        bl[n] = 1
                        #print('1', mag_all[idx[n]][m],'>',mag[n])
            elif (n_blend==1) & (mag[n]<17):
                bl[n] = 0
                print(mag[n])

        blend_type = []
        for n in range(len(bl)):
            if bl[n]==0:
                blend_type.append('a0')
            elif bl[n]==1:
                blend_type.append('a1')
            elif bl[n]==2:
                blend_type.append('a2')
            elif bl[n]==3:
                blend_type.append('a3')

        print('Priorities -> ',Counter(blend_type))


        return np.asarray(blend_type)




    def get_sky(self, ra_all, dec_all, rr = None, n_fiber_sky=10, show=True):
        if rr==None:
            rr = float(self.SKYSIZE)
        rr = rr/3600.
        sky = []
        k=-1

        while len(sky)<n_fiber_sky:
            k += 1

            rho = np.sqrt(np.random.uniform(0, self.RADIUS_CUT/4))
            phi = np.random.uniform(0, 2*np.pi)

            ra_random = rho * np.cos(phi)+self.raCENTER
            dec_random = rho * np.sin(phi)+self.decCENTER



            pts_all = np.array([ra_all, dec_all]).T
            pts = np.array([ra_random, dec_random]).T
            max_distance = (rr)*1.01

            Tree = KDTree(pts_all)
            idx = Tree.query_ball_point(pts, r=max_distance)
            if idx == []:
                sky.append(pts)
            if k%200 == 0:
                print('it is very crowded...')
            if k == 1000:
                print("infinite loop in finding sky. I found "+str(len(sky))+" for now..." )
                print("try to reduce the radius")
                break

        self.sky = np.array(sky)

        if show == True:
            self.sky = McL.show_skypoint(self, ra_all, dec_all, rr)

        print()
        print("Saved "+str(len(self.sky))+" sky positions")
        print()

        idd_sky = ['sky_'+str(n) for n in range(len(sky))]
        flag_sky = ['T']*len(sky)
        priority_sky = [1]*len(sky)
        mag_sky = [-99]*len(sky)

        return idd_sky, self.sky[:,0], self.sky[:,1], flag_sky, priority_sky, mag_sky
    
    
    
    def plot_hist_sky(self, ra_sky, dec_sky, gridx, gridy):
        fig, ax = plt.subplots(figsize=(5,5), dpi=500)
        #plt.pcolormesh(ra_sky, dec_sky, grid)
        plt.hist2d(ra_sky, dec_sky, bins =[gridx, gridy], vmin=0, vmax=self.N_SKY_XCELL)
        plt.plot(ra_sky, dec_sky, 'k.')
        plt.colorbar()
        plt.xlim(self.raCENTER-self.RADIUS_CUT*2, self.raCENTER+self.RADIUS_CUT*2)
        plt.ylim(self.decCENTER-self.RADIUS_CUT*2, self.decCENTER+self.RADIUS_CUT*2)
        ax.set_aspect("equal")
        plt.show()
        plt.close()

    
    def skygrid(self, sky_points, k):
        ra_sky, dec_sky = np.array(sky_points)[:,0], np.array(sky_points)[:,1]
        c = (self.raCENTER,self.decCENTER)
        rr = self.RADIUS_CUT
        r = rr#np.sqrt(2)*rr/2
        nbox = self.N_CELL 
        nbox += 1
        
        
        gridx = np.linspace(c[0]-r*1.2,c[0]+r*1.2, nbox)
        gridy = np.linspace(c[1]-r*1.2,c[1]+r*1.2, nbox)
        
        grid, _, _ = np.histogram2d(ra_sky, dec_sky, bins=[gridx, gridy])
        if (k%200 == 0) & (k!=0):
            #print(grid)
            McL.plot_hist_sky(self, ra_sky, dec_sky, gridx, gridy)
        #if len(ra_sky)>5:    
            #McL.plot_hist_sky(self, ra_sky, dec_sky, grid)

        grid_red = grid[1:-1,1:-1]
        conds = (np.min(grid_red)>=self.N_SKY_XCELL)
        
        if conds == True:
            McL.plot_hist_sky(self, ra_sky, dec_sky, gridx, gridy)
            return 1, grid
        elif conds == False:
            return 0, grid
        



    def get_sky_grid(self, ra_all, dec_all, rr = None, n_fiber_sky=None, show=False, allsky=True):
        if rr == None:
            rr = self.SKYSIZE
        rr = rr/3600.
        
        n_fiber_sky = self.N_SKY_XCELL*self.N_CELL**2
        sky = []
        k=-1
        
        cond = 0
        while cond != 1:
            k += 1

            rho = np.sqrt(np.random.uniform(0, self.RADIUS_CUT/4))
            phi = np.random.uniform(0, 2*np.pi)

            ra_random = rho * np.cos(phi)+self.raCENTER
            dec_random = rho * np.sin(phi)+self.decCENTER


            pts_all = np.array([ra_all, dec_all]).T
            pts = np.array([ra_random, dec_random]).T
            max_distance = rr*1.01

            Tree = KDTree(pts_all)
            idx = Tree.query_ball_point(pts, r=max_distance)
            if idx == []:
                sky.append(pts)
            if k%200 == 0:
                print('it is very crowded...')
            if k == 2000:
                print("infinite loop in finding sky. I found "+str(len(sky))+" for now..." )
                print("try to reduce the 'empty sky' radius (now is "+str(int(3600*rr))+" arcsec)")
                break
            
            
            #if len(np.array(sky)[:,0])>3:
            cond, grid = McL.skygrid(self, np.array(sky), k)



        self.sky = np.array(sky)

        if show == True:
            self.sky = McL.show_skypoint(self, ra_all, dec_all, rr)

        print()
        print("Saved "+str(len(self.sky))+" sky positions")
        print()        

        rra_sky, ddec_sky = McL.equisky(self, self.sky[:,0], self.sky[:,1], n_sky=None , show=False)

        idd_sky = ['sky_'+str(n) for n in range(len(rra_sky))]
        flag_sky = ['T']*len(rra_sky)
        priority_sky = [1]*len(rra_sky)
        mag_sky = [-99]*len(rra_sky)

        if allsky==True:
            idd_sky_all = ['sky_'+str(n) for n in range(len(self.sky[:,0]))]
            flag_sky_all = ['T']*len(self.sky[:,0])
            priority_sky_all = [1]*len(self.sky[:,0])
            mag_sky_all = [-99]*len(self.sky[:,0])

        print()
        print("Printed "+str(self.N_SKY)+" sky positions")
        print()        

        if allsky==True:
            return idd_sky, rra_sky, ddec_sky, flag_sky, priority_sky, mag_sky, idd_sky_all, self.sky[:,0], self.sky[:,1], flag_sky_all, priority_sky_all, mag_sky_all
        else:
            return idd_sky, rra_sky, ddec_sky, flag_sky, priority_sky, mag_sky



    def equisky(self, ra_sky, dec_sky, n_sky=None , show=True):
        if n_sky == None:
            n_sky=self.N_SKY
        k = n_sky
        points = np.array([ra_sky, dec_sky]).T
        pts2D = points

        kmeans = KMeans(n_clusters=k, random_state=0).fit(pts2D)
        labels = kmeans.predict(pts2D)
        cntr = kmeans.cluster_centers_




        # indices of nearest points to centres
        approx = []
    
        for j, c in enumerate(cntr):
            lab = np.where(labels == j)[0]
            pts = pts2D[lab]
            d = distance_matrix(c[None, ...], pts)
            idx1 = np.argmin(d, axis=1) + 1
            idx2 = np.searchsorted(np.cumsum(labels == j), idx1)[0]
            approx.append(idx2)
        
        c = (self.raCENTER,self.decCENTER)
        rr = self.RADIUS_CUT
        r = rr#np.sqrt(2)*rr/2
        nbox = self.N_CELL 
        nbox += 1

        gridx = np.linspace(c[0]-r*1.2,c[0]+r*1.2, nbox)
        gridy = np.linspace(c[1]-r*1.2,c[1]+r*1.2, nbox)


        fig, ax = plt.subplots(figsize=(5, 5), dpi=250)
        ax.hist2d(pts2D[:, 0], pts2D[:, 1], bins =[gridx, gridy], vmin=0, vmax=self.N_SKY_XCELL)
        ax.plot(pts2D[:, 0], pts2D[:, 1], 'k.')
        ax.plot(cntr[:, 0], cntr[:, 1], 'x')
        ax.plot(pts2D[approx, 0], pts2D[approx, 1], 'ro')
        ax.set_aspect("equal")
        #plt.colorbar(ax=a)
        fig.legend(["points", "centres", "selected"], loc=1)
        ax.set_xlim(self.raCENTER-self.RADIUS_CUT*2, self.raCENTER+self.RADIUS_CUT*2)
        ax.set_ylim(self.decCENTER-self.RADIUS_CUT*2, self.decCENTER+self.RADIUS_CUT*2)
        plt.show()
        plt.close()
        
        return pts2D[approx, 0], pts2D[approx, 1]





    def show_skypoint(self, ra_all, dec_all, rr):
        ind = ['y']*len(self.sky)
        ind = np.asarray(ind)
        for n in range(len(self.sky)):
            fig, ax = plt.subplots(figsize=(7,7), dpi=100)
            dd = (rr)*7
            ralim = [self.sky[n][0]-dd,self.sky[n][0]+dd]
            declim = [self.sky[n][1]-dd,self.sky[n][1]+dd]

            i = np.where(((ra_all>=ralim[0]) & (ra_all<=ralim[1])) & ((dec_all>=declim[0]) & (dec_all<=declim[1])))

            ax.scatter(ra_all[i], dec_all[i], marker='.', color='k', s=20)


            circle1 = plt.Circle((self.sky[n][0], self.sky[n][1]), rr/3600, color='r', fill=False)
            circle2 = plt.Circle((self.sky[n][0], self.sky[n][1]), 1.2/3600, color='r', fill=False)
            ax.add_artist(circle1)
            ax.add_artist(circle2)

            plt.xlim(ralim[0], ralim[1])
            plt.ylim(declim[0], declim[1])
            plt.show(block=False)
            ind[n] = str(input('Is it sky? [y,n]  '))
            plt.close()


        i = np.array(np.where(ind=='y')).reshape(-1)

        return self.sky[i]




    def control_doubles(namefile):
        a,b,c,e,f   =   np.genfromtxt(namefile, usecols=(0,1,2,4,5), unpack=True)
        d           =   np.genfromtxt(namefile, usecols=(3), dtype='str', unpack=True)
        j = len(a)
        _, i = np.unique(a, return_index=True)
        jj = len(a[i])
        a,b,c,d,e,f = a[i],b[i],c[i],d[i],e[i],f[i]

        file=open(namefile, "w")
        file.write("# ID          RA        DEC        flag   priority   mag  \n".format())
        for val in zip(a,b,c,d,e,f):
            file.write('  {:10}   {:=9.5f}   {:=9.5f}   {:=2}   {:=1.0f}   {:=7.4f}   \n'.format(val[0],val[1],val[2],val[3],val[4],val[5]))
        file.close()

        print("deleted "+str(j-jj)+" double elements")
