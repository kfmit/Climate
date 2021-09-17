#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Climate Variability and Diagnostics Utilities (cvd_utils)

Functions for use in the Python assignments for 12.860 CVD course, Fall 2021

"""

import numpy as np
import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import xarray as xr
import glob

#%%Tutorial 1 - Intro to Python.

def transect(x,d,v,returnplot=True):
    """
    Usage: [V,fig,ax,pcm,cbar]=transect(x,d,v)
    if returnplot is True, then V=transect(x,d,v)
    
    Transect produces a color-scaled transect diagram of oceanographic data from CTD profiles collected at different locations and times.
    This is a bare-bones, Pythonic translation of Chad A. Greene's "transect.m".
    
    By default, data is linearly interpolated between the maximum and minimum x and d values.
    This is performed first along depth (d) for each cast, then along each cast location (x). 
    Data outside the range is set to NaN.
    
    Inputs
    ------
    1) x : ARRAY[cast_number]
        Locations for each casts
    2) d : Nested ARRAY[cast_number][depth]
        Depths for each cast. This is a nested array where the outer index retrieves the cast, while the inner index retrieves the depth level
    3) v : Nested ARRAY[cast_number][depth]
        Corresponding variable value at each depth. Same structure as d.
    4) returnplot : BOOL
        Set to true to plot and return matplotlib objects (default). Set to false to just return the 2-D interpolated result.
    
    Output
    ------
    1) V : ARRAY[depth,x]
        Interpolated, 2-D array
    --- The following output applies if returnplot=True ---
    2) fig : matplotlib object
        Figure on which the plotting has occurred
    3) ax : matplotlib object
        Subplot axes where the plotting has occurred. (Use this to modify the axes and title labels!)
    4) pcm : matplotlib object
        The visualized, interpolated values
    5) cbar : matplotlib object
        Colorbar object
    
    Dependencies: 
        numpy as np
        matplotlib.pyplot as plt
    
    """
    # --------------------
    # Initial Error checks
    # --------------------
    assert np.all(len(x)==len(d)==len(v)),'Dimensions of x, d, and v must all agree.'
    
    # --------------
    # Ready the data
    # --------------
    # Get an array of d values so we know the depth ranges we're working with
    darray = []
    for n in range(len(v)):
        for z in range(len(d[n])):
            darray.append(d[n][z])
    darray = np.array(darray)
    
    # Common depths to interpolate to (using max and min depths)
    #if ~userdi:
    di = np.linspace(darray.min(),darray.max(),1000)
    
    # Final x locations to interpolate to
    #if ~userxi:
    xi = np.linspace(x.min(),x.max(),2000)
    
    # -----------
    # Interpolate
    # -----------
    # This section interpolates twice. First vertically, then horizontally.
    
    # Preallocate the first grid [depth, cast_number]
    V1 = np.ones((len(di),len(x)))*np.nan
    
    # Loop through each CTD location and interpolate in z
    for k in range(len(x)):
        
        z_within = ((di>d[k].min()) * (di<d[k].max())) # Only interpolate within data
        V1[z_within,k] = np.interp(di[z_within],d[k].squeeze(),v[k].squeeze()) # Squeeze if d and v are 2D
    
    # Preallocate the final grid [depth, distance]
    V = np.ones((len(di),len(xi)))*np.nan
    
    # Loop through each depth level
    for k in range(len(di)):
        V[k,:] = np.interp(xi,x,V1[k,:])
    
    # -------------
    # Make the plot
    # -------------
    # Initialize and plot the interpolated values
    fig,ax = plt.subplots(1,1)
    pcm = ax.pcolormesh(xi,di,V)
    
    # Plot a marker for each CTD cast
    for cast in range(len(d)):
        for z in range(len(d[cast])):
            ax.plot(x[cast],d[cast][z],marker=".",color="k")
    
    # Quick adjustments
    ax.invert_yaxis()
    cbar = fig.colorbar(pcm,ax=ax)
    
    return V,fig,ax,pcm,cbar

#%% Tutorial 2 - ENSO

def plot_anomaly(time,invar,
                 xlabfreq=None,ax=None,
                 barplot=False,lineplot=True):
    """
    Usage : ax =plot_anomaly(time,invar,
                 xlabfreq=None,ax=None,
                 barplot=False,lineplot=True)
    
    Plot anomaly of a timeseries onto a given axis
    where positive=red, negative=blue
    
    Parameters
    ----------
    time : 1-D Array or pd.Series [time,]
        The time-stamp for each observation
    invar : 1-D Array or pd.Series [time,]
        Variable to plot
    xlabfreq : INT, optional
        Frequency of x-ticks. The default is 10 ticks.
    ax : matplotlib object, optional
        Axis to plot on. The default grabs the current axis.
    barplot : BOOL, optional
        Set to True to make a barplot. The default is False.
    lineplot : BOOL, optional
        Set to False to exclude the original timeseries. The default is True.

    Returns
    -------
    ax : matplotlib object
        Returns the axis with the plotting result
        
    Dependencies:
        numpy as np
        matplotlib.pyplot as plt
    """
    
    # Get axis
    if ax is None:
        ax = ax.gca()
        
    # Set up time axis
    timeaxis = np.arange(0,len(time))  # Time axis for plotting
    if xlabfreq is None:
        xlabfreq = int(0.10*len(time)) # Defaults to ~10 ticks
    timeticks = timeaxis[::xlabfreq]
    timelabs = time[::xlabfreq]
    
    # Initialize and plot
    if barplot:
        id_pos = invar >= 0
        id_neg = invar <= 0
        bar_neg=ax.bar(timeaxis[id_neg],invar[id_neg],color='cornflowerblue',alpha=0.7,label="")
        bar_pos=ax.bar(timeaxis[id_pos],invar[id_pos],color='lightcoral',alpha=0.7,label="")
    else:
        ax.fill_between(timeaxis,0,invar,where=invar<=0,facecolor='cornflowerblue',interpolate=True,alpha=0.7)
        ax.fill_between(timeaxis,0,invar,where=invar>0,facecolor='lightcoral',interpolate=True,alpha=0.7)
    if lineplot: # include original timeseries in black
        ln = ax.plot(timeaxis,invar,color="black",lw=1) 
    
    ax.axhline(0,color='k',lw=0.75,ls='dashed') 
    
    # Time axis labels
    ax.set_xticks(timeticks)
    ax.set_xticklabels(timelabs)
    ax.set_xlim([0,timeaxis[-1]])
    ax.grid(True,ls="dotted")
    
    return ax

def numpy_to_da(invar,time,lat,lon,varname,savenetcdf=None):
    """
    Usage: da = numpy_to_da(invar,lon,lat,time,varname)
    
    Converts a NumPy array into an xr.DataArray with the same
    coordinates as the provided arrays.
    
    Parameters
    ----------
    invar : 3D ARRAY[time x lat x lon]
        Input variable
    lon :   1D ARRAY[lon]
        Longitude
    lat : 1D ARRAY[lat]
        Latitude
    time : 1D ARRAY[time]
        Time
    varname : STR
        Name of the variable
    savenetcdf : STR 
        If string argument is provided, saves as netcdf to the
        path indicated by the string. Default is None.

    Returns
    -------
    da : xr.DataArray
    """
    
    da = xr.DataArray(invar,
                dims={'time':time,'lat':lat,'lon':lon},
                coords={'time':time,'lat':lat,'lon':lon},
                name = varname
                )
    if savenetcdf is None:
        return da
    else:
        print("Saving netCDF to %s"%savenetcdf)
        da.to_netcdf(savenetcdf,
                 encoding={varname: {'zlib': True}})
        return da
    
def findnearestval(srchval,srcharray,bias=0):
    """
    Usage: i = findnearest(srchvalue,srcharray,bias)
    
    Find the nearest numerical value in an array to a search value
    All occurences are returned as array subscripts.
    Empty array is returned if nothing is found.
    
    Python translation of findnearest.m (2002) by
    Tom Benson of University College London 
    
    Inputs
    ------
    1) srchval : NUMERIC
        A numerical search value
    2) srcharray : ARRAY
        The array to be searched
    3) bias : -1, 0, or 1
         0 -- (default) no bias
        -1 -- bias output to lower search values
         1 -- bias output to higher search values
    
    Output
    ------
    1) i : INT
        Index of found value
        
    Dependencies:
        Numpy as np
    """
    # Find the Differences
    srcharray  = srcharray - srchval

    if bias == -1: # only choose values <= to the search value
        srcharray[srcharray>0] = np.nan
    elif bias == 1:# only choose values >= to the search value
        srcharray[srcharray<0] = np.nan

    if np.all(np.isnan(srcharray)):
        r = []   # None Found
    else:
        r = np.where(np.abs(srcharray) == np.nanmin(np.abs(srcharray)))
        r = r[0][0] # Just take the first value
    return r
    
#%% Tutorial 03 SSH

def plotmap(ax=None,bbox=None,proj=None,clon=0,figsize=(12,8),land_color=None,blabels=[1,0,0,1]):
    """
    Usage: fig,ax = plotmap(ax=None,bbox=None,proj=None,clon=0,figsize=(12,8),land_color=None,blabels=[1,0,0,1])
    if ax is provided: ax = plotmap(ax=ax,bbox=None,proj=None,clon=0,figsize=(12,8),land_color=None,blabels=[1,0,0,1])
    
    Initialize a figure or axes with coastlines and gridlines.

    Parameters (All Arguments are Optional!)
    ----------
    ax : Cartopy GeoAxes
        Axes to plot on. The default is to create both figure and axes within the function.
    bbox : LIST [west_bound,east_bound,south_bound,north_bound]
        Geographic bounding box/extent of the plot. First two elements are longitude bounds, last two
        are latitude bounds. The default is Global ([-180,180,-90,90]).
    proj : crs.projection
        The spatial projection. The default is ccrs.PlateCarree().
    clon : NUMERIC
        Central Longitude for the projection. The default is 0.
    figsize : LIST(Width,Height)
        Figure width and height in inches. The default is (12,8).
    land_color : STR,
        Color of the continents. The default is None.
    blabels : LIST of BOOL[Left, Right, Upper, Lower]
        Set to 1 to label the axis (for PlateCarree()). The default is [1,0,0,1].
    
    Returns
    -------
    if ax is None:
        returns fig,ax
    if ax is provided
        returns ax
        
    Dependencies
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    
    """
    # Set up projection
    if proj is None:
        proj = ccrs.PlateCarree(central_longitude=clon)

    # Initialize Figure
    init_fig=False
    if ax is None: # If no axis is provided, intialize figure
        init_fig=True
        fig = plt.figure(figsize=figsize)
        ax  = plt.axes(projection=proj)

    # Set plot extent
    if bbox is not None:
#         ax.set_extent([-180,180,-90,90],crs=ccrs.PlateCarree()) # Set the specified extent with bbox
#     else:
        ax.set_extent(bbox,crs=ccrs.PlateCarree()) # Set the specified extent with bbox

    # Add coastlines, continents
    ax.coastlines()
    if land_color:
        ax.add_feature(cfeature.LAND,facecolor=land_color)

    # Add gridlines
    gl = ax.gridlines(draw_labels=False,
                  linewidth=2, color='gray', alpha=1, linestyle="dotted",lw=0.75)

    # Remove the degree symbol
    gl.xformatter = LongitudeFormatter(direction_label=False,degree_symbol='')
    gl.yformatter = LatitudeFormatter(direction_label=False,degree_symbol='')
    
    # Turn off labels according to blabels
    gl.left_labels   = blabels[0]
    gl.right_labels  = blabels[1]
    gl.top_labels    = blabels[2]
    gl.bottom_labels = blabels[3]    
    
    if init_fig:
        return fig,ax
    else:
        return ax

def convert_datenum(matlab_datenum,datestr=False,fmt="%d-%b-%Y %H:%M:%S",autoreshape=True,return_datetimeobj=False,verbose=False):
    """
    Usage: python_datenum = convert_datenum(matlab_datenum,datestr=False,
                                fmt="%d-%b-%Y %H:%M%S",autoreshape=True,
                                return_datetimeobj=False,verbose=False)
    
    Converts an array of Matlab Datenumbers to either an array of np.datetime64 or strings (if datestr=True).
    This considers the offset for different origin date for matlab (Jan 1, year 0) vs. Python (Jan 1, 1970).
    If you prefer to work with datetime objects, set return_datetimeobj to True.
    
    Inputs
    ------
    1) matlab_datenum     : N-D ARRAY of floats
        Array containing matlab datenumbers to convert
    2) datestr            : BOOL
        Set to True to convert output to human-readable strings. Default returns np.datetime64[ns]
    3) fmt                : STR
        Format out output datestring. ex. of default format is "30-Jan-1990 23:59:00"
    4) autoreshape        : BOOL
        By default, the script flattens then reshapes the array to its original dimensions.
        Set to False to just return the flattened array.
    5) return_datetimeobj : BOOL
        By default, the script returns np.datetime64. To just return datetime, set this to true.
        NOTE: The dimensions will be flattened and autoreshape will not work.
    6) verbose            : BOOL
        Set to true to print messages
    
    Output
    ------
    1) python_datetime : N-D ARRAY of np.datetime64[64] or strings
        Array containing the result of the conversion.
    
    Dependencies
        import pandas as pd
        import numpy as np
        import datetime as dt
        
    """
    # Preprocess inputs (flattening n-d arrays)
    in_arr = np.array(matlab_datenum)
    dims   = in_arr.shape
    if len(dims)  > 1:
        if verbose:
            print("Warning! flattening %i-D array with dims %s"%(len(dims),str(dims)))
        in_arr = in_arr.flatten()

    # Calculate offset (In Python, reference date is Jan 1, year 1970 UTC, see "Unix Time")
    # Additional 366 days because reference date in matlab is January 0, year 0000 
    offset = dt.datetime(1970, 1, 1).toordinal() + 366

    # Convert to Python datetime, considering offset
    # Note that Matlab uses "days since", hence unit="D"
    python_datetime  = pd.to_datetime(in_arr-offset, unit='D')

    # Convert to datestring, and optionally convert to numpy array
    if datestr:
        if verbose:
            print("Converting to datestr")
        python_datetime = python_datetime.strftime(fmt)
        if return_datetimeobj:
            return python_datetime
        # Otherwise convert to string array
        python_datetime = np.array(python_datetime,dtype=str)
    else: # Convert to numpy array with datetime objects
        if return_datetimeobj:
            return python_datetime
        # Otherwise convert to np.datetime64 array
        python_datetime = np.array(python_datetime)

    # Reshape array if necessary (current works only for numpy arrays)
    if len(dims) > 1 and autoreshape:
        if verbose:
            print("Reshaping array to original dimensions!")
        # Reshape array to original dimensions
        python_datetime=python_datetime.reshape(dims) 
    return python_datetime

def read_ssh(region,verbose=False):
    """
    Usage: ssh, sshlon, sshlat, sshtime = read_ssh(region,verbose=False)
    
    This function loads sea level anomalies for a given region.
    Checks the current folder using glob.
    
    Inputs
    ------
    1) region : STR
        Region to load ssh for.
        "NA" is North Atlantic
        "PA" is Tropical Eastern Pacific
        
    Outputs
    -------
    1) ssh    : ARRAY[lat x lon x time]
        Sea surface height above geoid, in meters
    2) sshlon : ARRAY
        Longitude values
    3) sshlat : ARRAY
        Latitude values
    4) sshtime : ARRAY
        Timepoints
        
    Dependencies
        numpy as np
        xarray as xr
        glob
        
    """
    # Get list of regions to load
    if region == "NA":
        nclist = glob.glob("SSH_North_Atlantic*.nc")
    elif region == "PA":
        nclist = glob.glob("SSH_Tropical_Eastern_Pacific*.nc")
    else:
        print("Invalid region. This function accepts 'PA' or 'NA'")
    nclist.sort()
    if verbose:
        print("Found the following files!")
        print(nclist)

    # Load each region and concatenate after loading
    # Oddly, time is not monotonically increasing for portions of this dataset
    # so open_mfdataset is not working
    ssh = []
    sshtime =[]
    for i in range(len(nclist)):
        # Open dataset
        ds = xr.open_dataset(nclist[i])

        if i == 0: # Load lat/lon to numpy arrays
            sshlon = ds.longitude.values
            sshlat = ds.latitude.values

        # Append values and time, loading to np.arrays
        ssh.append(ds.adt.values)
        sshtime.append(ds.time.values)

    # Concatenate along time axis
    ssh = np.concatenate(ssh,axis=0) #[time x lon x lat]
    sshtime = np.hstack(sshtime)
    
    # Transpose to [lat x lon x time]
    ssh = ssh.transpose(2,1,0)
    return ssh, sshlon, sshlat, sshtime