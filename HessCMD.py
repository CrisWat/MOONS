### Built with Python 2.6

### Import other packages
from pylab import *
from numpy import *
from random import *
from scipy import *
from matplotlib import *
rcdefaults()
matplotlib.rc('font',family='Bitstream Vera Serif')

#colormap=cm.gray_r

def plotHess(mag, color, magRC, colorRC, magRGB, colorRGB,levels='None', xlims=[-1,3], ylims=[21,5], colormap='rainbow', cbarr='None', cbarrtitle='',xlab='', ylab='', savepath='./', namefig='plot.pdf', ftitle=''):
    """
    plotHess(mag,color) is a function intended to create a contour plot of stellar density in the canonical Color-Magnitude Diagram (CMD). plotHess() plots all the stars in the CMD and then overplots a contour of their density.

    mag = magnitudes of the stars
    color = colors of the stars
    mag and color are required to run plotHess() and need to be the same length arrays.


    Optional keyword arguments:

    =========   =======================================================
    Keyword     Description
    =========   =======================================================
    levels:    levels to be used for density contour; if 'None', the defaults are defined by the contour() function, which may not be optimal for the users dataset and may be changed if necessary.
    ylims:      define the y-range of the plotted values; defaults will need to be changed to fit user's data.
    xlims:      define the x-range of the plotted values; defaults will need to be changed to fit user's data.
    colormap:   allows the user to choose a Python color map to use for the contour; the default is grayscale.
    cbarr:      'None' means the colorbar will not be plot as a default. To change this, set cbarr='Yes'; in a true Hess diagram, the cbarr values represent the stellar point density in the CMD plot.
    cbarrtitle: sets the title for the colorbar; should be string
    xlabel:     sets the xlabel for the plot; should be string
    ylabel:     sets the ylabel for the plot; should be string
    saveas:     pathway for saving the output plot. The default is to save in the same folder as "HessCMD.png"
    ftitle:     title

    """
    ### Set figure options
    fig = pylab.figure(1, figsize=(8, 7))

    ### Plot all stars in CMD as small points
    plot(color,mag,'ko',ms=0.2,zorder=1)
    #plot(colorRGB,magRGB,'ko',color='green',ms=0.1,zorder=2,alpha=0.75)
    #plot(colorRC,magRC,'ko',color='red',ms=0.1,zorder=3,alpha=0.75)
    ylabel(ylab, fontsize=18)
    xlabel(xlab, fontsize=18)
    title(ftitle, fontsize=18)
    ### Use color and magnitude data to split and average into a binned 2D histogram
    Z, xedges, yedges = np.histogram2d(color,mag,bins=(200,200))

    ### Find the cell centers from the cell edges
    x = 0.5*(xedges[:-1] + xedges[1:])
    y = 0.5*(yedges[:-1] + yedges[1:])

    ### Promote to 2D arrays
    Y, X = np.meshgrid(y, x)

    ### Mask out values without data
    Z = np.ma.array(Z, mask=(Z==0))

    ### Set levels options as default or user-defined option
    if levels=='None':
        cntr=plt.contourf(X,Y,Z,cmap=cm.gray_r,zorder=2)
    else:
        cntr=plt.contourf(X,Y,Z,levels=levels,cmap=colormap,zorder=2)
    cntr.cmap.set_under('white')



    ### Set plot limits as default or user-defined option
    xlim(xlims[0],xlims[1]);ylim(ylims[0],ylims[1])

    ### Determine whether to plot a color bar
    if cbarr!='None':
        cbar_ax = fig.add_axes([0.91, 0.14, 0.02, 0.7])
        d=fig.colorbar(cntr, cax=cbar_ax)#, title=cbarrtitle)
        d.set_label(label=cbarrtitle,size=18)

    #plot(colorRC,magRC,'ko',color='green',ms=0.2,zorder=1)
    #plot(colorRC,magRC,'ko',color='red',ms=0.2,zorder=1)
    ### Set final options and save figure to default (current folder) or to a user-defined location on disk
    #fig.suptitle(label=figtitle, fontsize=16)
    #fig.suptitle(label=ftitle,fontsize=14, fontweight='bold')
    #plt.text(0,0,"this is a test")
    savefig(savepath+namefig,dpi=300)
    plt.show()


def plotHess2(mag, color, mag_target, color_target, levels='None', xlims=[-1,3], ylims=[21,5], colormap='rainbow', cbarr='None', cbarrtitle='',xlab='', ylab='', saveas=' ', ftitle=''):
    """
    plotHess(mag,color) is a function intended to create a contour plot of stellar density in the canonical Color-Magnitude Diagram (CMD). plotHess() plots all the stars in the CMD and then overplots a contour of their density.

    mag = magnitudes of the stars
    color = colors of the stars
    mag and color are required to run plotHess() and need to be the same length arrays.


    Optional keyword arguments:

    =========   =======================================================
    Keyword     Description
    =========   =======================================================
    levels:    levels to be used for density contour; if 'None', the defaults are defined by the contour() function, which may not be optimal for the users dataset and may be changed if necessary.
    ylims:      define the y-range of the plotted values; defaults will need to be changed to fit user's data.
    xlims:      define the x-range of the plotted values; defaults will need to be changed to fit user's data.
    colormap:   allows the user to choose a Python color map to use for the contour; the default is grayscale.
    cbarr:      'None' means the colorbar will not be plot as a default. To change this, set cbarr='Yes'; in a true Hess diagram, the cbarr values represent the stellar point density in the CMD plot.
    cbarrtitle: sets the title for the colorbar; should be string
    xlabel:     sets the xlabel for the plot; should be string
    ylabel:     sets the ylabel for the plot; should be string
    saveas:     pathway for saving the output plot. The default is to save in the same folder as "HessCMD.png"
    ftitle:     title

    """
    ### Set figure options
    fig = pylab.figure(1, figsize=(8, 7))

    ### Plot all stars in CMD as small points
    plot(color,mag,'ko',ms=0.2,zorder=1)
    plot(color_target,mag_target,'ko',color='white',ms=3,zorder=3,alpha=0.5)
    ylabel(ylab, fontsize=18)
    xlabel(xlab, fontsize=18)
    title(ftitle, fontsize=18)
    ### Use color and magnitude data to split and average into a binned 2D histogram
    Z, xedges, yedges = np.histogram2d(color,mag,bins=(300,300))

    ### Find the cell centers from the cell edges
    x = 0.5*(xedges[:-1] + xedges[1:])
    y = 0.5*(yedges[:-1] + yedges[1:])

    ### Promote to 2D arrays
    Y, X = np.meshgrid(y, x)

    ### Mask out values without data
    Z = np.ma.array(Z, mask=(Z==0))

    ### Set levels options as default or user-defined option
    if levels=='None':
        cntr=plt.contourf(X,Y,Z,cmap=cm.gray_r,zorder=2)
    else:
        cntr=plt.contourf(X,Y,Z,levels=levels,cmap=colormap,zorder=2)
    cntr.cmap.set_under('white')



    ### Set plot limits as default or user-defined option
    xlim(xlims[0],xlims[1]);ylim(ylims[0],ylims[1])

    ### Determine whether to plot a color bar
    if cbarr!='None':
        cbar_ax = fig.add_axes([0.91, 0.14, 0.02, 0.7])
        d=fig.colorbar(cntr, cax=cbar_ax)#, title=cbarrtitle)
        d.set_label(label=cbarrtitle,size=18)

    #plot(colorRC,magRC,'ko',color='green',ms=0.2,zorder=1)
    #plot(colorRC,magRC,'ko',color='red',ms=0.2,zorder=1)
    ### Set final options and save figure to default (current folder) or to a user-defined location on disk
    #fig.suptitle(label=figtitle, fontsize=16)
    #fig.suptitle(label=ftitle,fontsize=14, fontweight='bold')
    #plt.text(0,0,"this is a test")
    if saveas!=None:
        savefig(saveas,dpi=300)
    plt.show()


def plotHess3(mag, color, levels='None', xlims=[-1,3], ylims=[21,5], colormap='rainbow', cbarr='None', cbarrtitle='',xlab='', ylab='', saveas=' ', ftitle=''):
    """
    plotHess(mag,color) is a function intended to create a contour plot of stellar density in the canonical Color-Magnitude Diagram (CMD). plotHess() plots all the stars in the CMD and then overplots a contour of their density.

    mag = magnitudes of the stars
    color = colors of the stars
    mag and color are required to run plotHess() and need to be the same length arrays.


    Optional keyword arguments:

    =========   =======================================================
    Keyword     Description
    =========   =======================================================
    levels:    levels to be used for density contour; if 'None', the defaults are defined by the contour() function, which may not be optimal for the users dataset and may be changed if necessary.
    ylims:      define the y-range of the plotted values; defaults will need to be changed to fit user's data.
    xlims:      define the x-range of the plotted values; defaults will need to be changed to fit user's data.
    colormap:   allows the user to choose a Python color map to use for the contour; the default is grayscale.
    cbarr:      'None' means the colorbar will not be plot as a default. To change this, set cbarr='Yes'; in a true Hess diagram, the cbarr values represent the stellar point density in the CMD plot.
    cbarrtitle: sets the title for the colorbar; should be string
    xlabel:     sets the xlabel for the plot; should be string
    ylabel:     sets the ylabel for the plot; should be string
    saveas:     pathway for saving the output plot. The default is to save in the same folder as "HessCMD.png"
    ftitle:     title

    """
    ### Set figure options
    fig = pylab.figure(1, figsize=(8, 7))

    ### Plot all stars in CMD as small points
    plot(color,mag,'ko',ms=0.2,zorder=1)
    ylabel(ylab, fontsize=18)
    xlabel(xlab, fontsize=18)
    #title(ftitle, fontsize=18)
    ### Use color and magnitude data to split and average into a binned 2D histogram
    Z, xedges, yedges = np.histogram2d(color,mag,bins=(150,150))

    ### Find the cell centers from the cell edges
    x = 0.5*(xedges[:-1] + xedges[1:])
    y = 0.5*(yedges[:-1] + yedges[1:])

    ### Promote to 2D arrays
    Y, X = np.meshgrid(y, x)

    ### Mask out values without data
    Z = np.ma.array(Z, mask=(Z==0))

    ### Set levels options as default or user-defined option
    if levels=='None':
        cntr=plt.contourf(X,Y,Z,cmap=cm.gray_r,zorder=2)
    else:
        cntr=plt.contourf(X,Y,Z,levels=levels,cmap=colormap,zorder=2)
    cntr.cmap.set_under('white')



    ### Set plot limits as default or user-defined option
    xlim(xlims[0],xlims[1]);ylim(ylims[0],ylims[1])

    ### Determine whether to plot a color bar
    if cbarr!='None':
        cbar_ax = fig.add_axes([0.91, 0.14, 0.02, 0.7])
        d=fig.colorbar(cntr, cax=cbar_ax)#, title=cbarrtitle)
        d.set_label(label=cbarrtitle,size=18)

    #plot(colorRC,magRC,'ko',color='green',ms=0.2,zorder=1)
    #plot(colorRC,magRC,'ko',color='red',ms=0.2,zorder=1)
    ### Set final options and save figure to default (current folder) or to a user-defined location on disk
    #fig.suptitle(label=figtitle, fontsize=16)
    #fig.suptitle(label=ftitle,fontsize=14, fontweight='bold')
    #plt.text(0,0,"this is a test")
    if saveas!=None:
        savefig(saveas,dpi=300)
    plt.show()


