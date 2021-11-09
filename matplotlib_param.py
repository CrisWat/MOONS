#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pyramid as pyd
from matplotlib.pylab import *
import matplotlib.pyplot as plt
#matplotlib.use('TKAgg')


fnt = {'fontsize':20}; fnt_subp = {'fontsize':16}; fnt_leg = {'size':15}
def plot_layout(xt=16, yt=16):
    matplotlib.rc('xtick', labelsize=xt)
    matplotlib.rc('ytick', labelsize=yt)
    plt.rc('text', usetex=False); plt.rc('font', family='serif')
    plt.rcParams['mathtext.default']='regular'
    for string in ['xtick.direction', 'ytick.direction']:
        plt.rcParams[string]= 'in'
    for string in ['ytick.right', 'xtick.top']:
        plt.rcParams[string] = True
    for string in ['axes.linewidth', 'xtick.major.width', 'xtick.minor.width', 'ytick.major.width', 'ytick.minor.width']:
        plt.rcParams[string] = 1.5
    for string in ['xtick.major.size', 'ytick.major.size']:
        plt.rcParams[string]= 6
    for string in ['xtick.minor.size', 'ytick.minor.size']:
        plt.rcParams[string]= 3
