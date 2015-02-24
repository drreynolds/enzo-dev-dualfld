#!/usr/bin/env python
# matplotlib-based plotting script for Iliev et al. test #2
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
from yt.mods import *
from os import *

# reduce log level of yt 
from yt.config import ytcfg
ytcfg["yt","loglevel"] = "50"

# set the total number of snapshots
te = 2

# set the graphics output type
#pictype = 'pdf'
pictype = 'png'

# set some constants
q0 = 0.5               # deceleration parameter
Nph = 5.0e48           # ionization source strength [photons/sec]
alpha2 = 2.52e-13      # recombination rate coefficient
mp = 1.67262171e-24    # proton mass [g]


##########
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
          display_name="Ionized\; Fraction")

#   Radiation energy densities (log plot)
def _logUV(field, data):
    return (data["UV_Radiation"])
add_field("logUV", take_log=True, function=_logUV, 
          display_name="UV\; Radiation\; Energy\; Density")
def _logXray(field, data):
    return (data["Xray_Radiation"])
add_field("logXray", take_log=True, function=_logXray, 
          display_name="X-ray\; Radiation\; Energy\; Density")

#   Radius from domain center
def _radius(field, data):
    return (np.sqrt(data["x"]*data["x"] + data["y"]*data["y"] +
                    data["z"]*data["z"]))
def _convertradius(data):
    return (data.convert("cm"))
add_field("radius", take_log=False, function=_radius, 
          convert_function=_convertradius, 
          display_name="radius", units=r"\rm{cm}")

##########



# loop over snapshots, loading values and times
for tstep in range(0,te+1):
    
    # load relevant information
    sdump = repr(tstep).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    pf = load(pfile)
    xR = pf["DomainRightEdge"][0]*pf["cm"]
    z  = pf["CosmologyCurrentRedshift"]

    # determine domain "center" for plotting
    xC = 0.5*(pf["DomainLeftEdge"][0] + pf["DomainRightEdge"][0])
    yC = 0.5*(pf["DomainLeftEdge"][1] + pf["DomainRightEdge"][1])
    zC = 0.5*(pf["DomainLeftEdge"][2] + pf["DomainRightEdge"][2])

    # set time label
    tout = repr(tstep).zfill(2)
        
    # begin plot collection (center at (xC,yC,0))
    pc = PlotCollection(pf, [xC,yC,0.0])
        
    # slices through z=0
    p = pc.add_slice("xHI",'z')
    p = pc.add_slice("xHII",'z')
    p = pc.add_slice("logUV",'z')
    p = pc.add_slice("logXray",'z')
    pc.save('tstep' + tout, format=pictype)

    # rename generated files
    f1 = 'tstep' + tout + '_Slice_z_logUV.' + pictype
    f2 = 'UVcontour_' + tout + '.' + pictype
    rename(f1,f2)
    f1 = 'tstep' + tout + '_Slice_z_logXray.' + pictype
    f2 = 'XrayContour_' + tout + '.' + pictype
    rename(f1,f2)
    f1 = 'tstep' + tout + '_Slice_z_xHI.' + pictype
    f2 = 'HIcontour_' + tout + '.' + pictype
    rename(f1,f2)
    f1 = 'tstep' + tout + '_Slice_z_xHII.' + pictype
    f2 = 'HIIcontour_' + tout + '.' + pictype
    rename(f1,f2)

