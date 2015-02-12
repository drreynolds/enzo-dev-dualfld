#!/usr/bin/env python
# yt-based plotting script for Iliev et al. test #6
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
from yt.mods import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid


# set the total number of snapshots
te = 60

# set the graphics output type
#pictype = 'pdf'
pictype = 'png'

# set some constants
mp = 1.67262171e-24    # proton mass [g]
kpc = 3.0857e21        # 1 kpc in cm
megayear = 3.15576e13  # 1 Myr in sec


# Define some derived fields
#   Neutral Hydrogen fraction (log plot)
def _xHI(field, data):
    return (data["HI_Density"]/data["Density"])
add_field("xHI", take_log=True, function=_xHI, 
          display_name="Neutral\; Fraction")

#   Ionized Hydrogen fraction (log plot)
def _xHII(field, data):
    return (data["HII_Density"]/data["Density"])
add_field("xHII", take_log=True, function=_xHII, 
          display_name="Ionized fraction")

#   Radiation energy density (log plot)
def _logE(field, data):
    return (data["UV_Radiation"])
def _convertlogE(data):
    return ( data.convert("MassUnits")/data.convert("TimeUnits") / 
             data.convert("TimeUnits")/data.convert("cm"))
add_field("logE", take_log=True, function=_logE, 
          convert_function=_convertlogE, 
          display_name="Radiation\; Energy\; Density", 
          units=r"\rm{erg}/\rm{cm}^3")

#   Radius from domain center
def _radius(field, data):
    return (np.sqrt(data["x"]*data["x"] + data["y"]*data["y"] +
                    data["z"]*data["z"]))
add_field("radius", take_log=False, function=_radius, 
          display_name="radius", units=r"{r/L_{box}}")

#   Density (log plot)
def _logDens(field, data):
    return (data["Density"])
add_field("logDens", take_log=True, function=_logDens, 
          display_name="Density", units=r"\rm{g}/\rm{cm}^3")

#   Number Density (log plot)
def _logNDens(field, data):
    return (data["Density"]/mp)
add_field("logNDens", take_log=True, function=_logNDens, 
          display_name="number density [cm$^{-3}$]")

#   Velocity magnitude (log)
def _VelMag(field, data):
    return (np.sqrt(data["x-velocity"]*data["x-velocity"] + 
                    data["y-velocity"]*data["y-velocity"] +
                    data["z-velocity"]*data["z-velocity"] ))
add_field("VelMag", take_log=True, function=_VelMag, 
          display_name="|Velocity|", units=r"\rm{cm}/\rm{s}")

#   Sound speed (log)
def _VSound(field, data):
    #return (np.sqrt(np.divide(data["Pressure"],data["Density"])))
    return (np.sqrt(np.divide(data["logP"],data["Density"])))
add_field("Vsound", take_log=True, function=_VSound, 
          display_name="Sound speed", units=r"\rm{cm}/\rm{s}")

#   Mach number (log plot)
def _logMach(field, data):
    return (np.divide(data["VelMag"],data["Vsound"]))
add_field("logMach", take_log=True, function=_logMach, 
          display_name="Mach number")

#   Pressure (log plot)
def _logP(field, data):
    #return (data["Pressure"])
    return (2.0/3.0*data["Density"]*data["Internal_Energy"]*pf["LengthUnits"]**2/pf["TimeUnits"]**2)
add_field("logP", take_log=True, function=_logP, 
          display_name="Pressure", units=r"\rm{g}/\rm{cm}/\rm{s}^2")

#   Temperature (log plot)
def _logT(field, data):
    return (data["Temperature"])
add_field("logT", take_log=True, function=_logT, 
          display_name="Temperature [K]")




# initialize time-history outputs
#    row 1: time (t)
#    row 2: computed i-front radius
rdata = zeros( (2, te+1), dtype=float);


# loop over snapshots, loading values and times
for tstep in range(0,te+1):
    
    # load relevant information
    sdump = repr(tstep).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    pf = load(pfile)
    t = pf.current_time * pf["TimeUnits"]

    # determine if simulation was run with source in center or corner
    spherical = (2.0**(pf.domain_left_edge[0]/pf.domain_right_edge[0]+1.0) 
               * 2.0**(pf.domain_left_edge[1]/pf.domain_right_edge[1]+1.0) 
               * 2.0**(pf.domain_left_edge[2]/pf.domain_right_edge[2]+1.0))

    # compute I-front radius (assuming spherical)
    sp = pf.h.sphere([0.0, 0.0, 0.0], 1.0)
    HIIvolume = (sp["xHII"]*sp["CellVolumeCode"]*pf["cm"]**3).sum()*spherical
    radius = (3.0/4.0*HIIvolume/pi)**(1.0/3.0)
    
    # store data
    rdata[0][tstep] = t/megayear
    rdata[1][tstep] = radius
    
    # generate 2D plots at certain times
    if (tstep == 6) or (tstep == 20) or (tstep == 50):
        
        # set time label
        if (tstep == 6):
            Myr  = '3'
            Myr2 = '03'
            style = 'r:'
            lab = '3 Myr'
        elif (tstep == 20):
            Myr  = '10'
            Myr2 = '10'
            style = 'b--'
            lab = '10 Myr'
        else:
            Myr  = '25'
            Myr2 = '25'
            style = 'k-'
            lab = '25 Myr'
        
        # determine domain "center" for plotting
        xC = 0.5*(pf["DomainLeftEdge"][0] + pf["DomainRightEdge"][0])
        yC = 0.5*(pf["DomainLeftEdge"][1] + pf["DomainRightEdge"][1])
        zC = 0.5*(pf["DomainLeftEdge"][2] + pf["DomainRightEdge"][2])


        # increase the standard font sizes for the following plot
        plt.rcParams.update({'font.size': 20})


        fig = plt.figure(tstep)

        # See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
        grid = AxesGrid(fig, (0.08,0.075,0.85,0.85),
        #grid = AxesGrid(fig, (-0.02,-0.05,0.9,1.0),
                        nrows_ncols = (2, 2),
                        axes_pad = 1.0,
                        label_mode = "L",
                        share_all = True,
                        cbar_location="right",
                        cbar_mode="each",
                        cbar_pad="3%")

        fields = ['xHII', 'logMach', 'logNDens', 'logT']

        # Create the plot.  Since SlicePlot accepts a list of fields, we need only
        # do this once.
        p = SlicePlot(pf, 'z', fields, center=[xC,yC,zC], axes_unit='kpc', fontsize=20)

        # For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
        # axes.
        for i, field in enumerate(fields):
            slplot = p.plots[field]
            slplot.figure = fig
            slplot.axes = grid[i].axes
            slplot.axes.xaxis.set_ticklabels([])
            slplot.axes.yaxis.set_ticklabels([])
            slplot.cax = grid.cbar_axes[i]
            txt = slplot.cax.yaxis.label
            font = matplotlib.font_manager.FontProperties(size=20)
            txt.set_font_properties(font)

        # Finally, redraw the plot on the AxesGrid axes.
        p._setup_plots()
        plt.savefig('slices_' + Myr + 'Myr.' + pictype % pf)



        # generate profile plots by averaging results from multiple rays
        rays = np.array( ((1.0,0.0,0.0), (0.0,1.0,0.0), (0.0,0.0,1.0), 
                          (1.0,1.0,0.0), (1.0,0.0,1.0), (0.0,1.0,1.0), 
                          (1.0,1.0,1.0)) )
        nrays = 7
        nradii = 200
        #rvals = np.linspace(0,1*pf["cm"]/rs0,nradii)
        rvals = np.linspace(0,sqrt(3.0),nradii)
        HIprofile  = np.zeros(rvals.shape, dtype=float)
        HIIprofile = np.zeros(rvals.shape, dtype=float)
        Tprofile   = np.zeros(rvals.shape, dtype=float)
        nProfile   = np.zeros(rvals.shape, dtype=float)
        pProfile   = np.zeros(rvals.shape, dtype=float)
        mProfile   = np.zeros(rvals.shape, dtype=float)

        # generate 1D profiles from ray emanating out from box center to corner
        for iray in range(0,nrays):
            rvec = rays[iray,:]
            rvec = rvec / sqrt(rvec[0]**2 + rvec[1]**2 + rvec[2]**2)
            r = pf.h.ray([0.0, 0.0, 0.0], rvec)
            HIprof  = r["xHI"]
            HIIprof = r["xHII"]
            Tprof   = r["logT"]
            nProf   = r["Density"]/mp
            pProf   = r["logP"]
            mProf   = r["logMach"]
            Hradii  = r["radius"]
        
            # sort results by radius (since that isn't quite working correctly from yt)
            ptype = [('r', float), 
                     ('xHI', float), 
                     ('xHII', float), 
                     ('T', float),
                     ('n', float),
                     ('p', float),
                     ('m', float)]
            pdata = np.zeros(Hradii.shape, dtype=ptype);
            nrad = (Hradii.shape)[0]
            for irad in range(0,nrad):
                pdata[irad] = (Hradii[irad], 
                               HIprof[irad], 
                               HIIprof[irad], 
                               Tprof[irad],
                               nProf[irad],
                               pProf[irad],
                               mProf[irad])
            pdata = np.sort(pdata, order='r')
            
            # interpolate results into output arrays
            tmp = np.interp(rvals, pdata['r'], pdata['xHI'])
            HIprofile += tmp
            tmp = np.interp(rvals, pdata['r'], pdata['xHII'])
            HIIprofile += tmp
            tmp = np.interp(rvals, pdata['r'], pdata['T'])
            Tprofile += tmp
            tmp = np.interp(rvals, pdata['r'], pdata['n'])
            nProfile += tmp
            tmp = np.interp(rvals, pdata['r'], pdata['p'])
            pProfile += tmp
            tmp = np.interp(rvals, pdata['r'], pdata['m'])
            mProfile += tmp
        HIprofile  /= nrays
        HIIprofile /= nrays
        Tprofile   /= nrays
        nProfile   /= nrays
        pProfile   /= nrays
        mProfile   /= nrays

        # change the standard font sizes for the following plot
        plt.rcParams.update({'font.size': 20})

        # generate stored profile plots
        fig = plt.figure(1)
        plt.subplot(221)
        plt.semilogy(rvals,nProfile,style,hold=True)
        plt.ylabel('n [cm$^{-3}$]', fontsize=20)
        ax = gca()
        ax.yaxis.set_ticks_position("left")
        ax.yaxis.set_label_position("left")
        ax.xaxis.set_ticklabels([])
        plt.axis([ 0.0, 1.0, 2e-2, 1e1 ])
        plt.grid()
        plt.subplot(222)
        plt.semilogy(rvals,Tprofile,style,hold=True)
        plt.ylabel('T [K]', fontsize=20)
        ax = gca()
        ax.yaxis.set_ticks_position("right")
        ax.yaxis.set_label_position("right")
        ax.xaxis.set_ticklabels([])
        plt.axis(( 0.0, 1.0, 1e3, 1e5 ))
        plt.grid()
        plt.subplot(223)
        plt.semilogy(rvals,mProfile,style,hold=True,label=lab)
        plt.ylabel('Mach', fontsize=20)
        ax = gca()
        ax.yaxis.set_ticks_position("left")
        ax.yaxis.set_label_position("left")
        plt.axis(( 0.0, 1.0, 1e-4, 8e0 ))
        plt.xlabel('$r/L_{box}$', fontsize=20)
        plt.grid()
        plt.subplot(224)
        plt.semilogy(rvals,pProfile,style,hold=True)
        plt.ylabel('p [g cm$^{-1}$ s$^{-2}$]', fontsize=20)
        ax = gca()
        ax.yaxis.set_ticks_position("right")
        ax.yaxis.set_label_position("right")
        plt.axis(( 0.0, 1.0, 2e-16, 2e-11 ))
        plt.xlabel('$r/L_{box}$', fontsize=20)
        plt.grid()


# change the standard font sizes for the following plot
plt.rcParams.update({'font.size': 12})

# save accumulated profile plot
fig = plt.figure(1)
plt.subplot(223)
plt.legend( loc=3 )
plt.subplot(221).tick_params(axis='both', labelsize=20)
plt.subplot(222).tick_params(axis='both', labelsize=20)
plt.subplot(223).tick_params(axis='both', labelsize=20)
plt.subplot(224).tick_params(axis='both', labelsize=20)
plt.subplots_adjust(left=0.15)
plt.subplots_adjust(right=0.85)
plt.subplots_adjust(top=0.95)
plt.subplots_adjust(bottom=0.15)
plt.savefig('profiles.' + pictype)



# change the standard font sizes for the following plot
plt.rcParams.update({'font.size': 22})

# I-front radius/velocity
tdata = 0.5*(rdata[0,1:]+rdata[0,:-1])
vdata = (rdata[1,1:]-rdata[1,:-1])/(rdata[0,1:]-rdata[0,:-1])*(kpc/1e5/megayear)
fig = plt.figure(2)
plt.subplot(211)
plt.plot(rdata[0],rdata[1]/kpc,'b-')
plt.ylabel('$r_I$ [kpc]', fontsize=20)
ax = gca()
ax.xaxis.set_ticklabels([])
plt.axis([ 0.0, 25.0, 0.0, 0.59 ])
plt.grid()
plt.subplot(212)
plt.plot(tdata,vdata/kpc,'b-')
plt.xlabel('$t$ [Myr]', fontsize=20)
plt.ylabel('$v_I$ [km/s]', fontsize=20)
plt.axis([ 0.0, 25.0, 1e1, 3e1 ])
plt.grid()
plt.subplots_adjust(left=0.15)
plt.subplots_adjust(bottom=0.12)
plt.savefig('radius.' + pictype)


# end of script