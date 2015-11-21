### CHANGE: specify local python installation
#!home/shannon/usr/local/bin/python3.4

from netCDF4 import Dataset
from math import fsum, exp
from random import uniform
import numpy as np


# Example of offline cost function and adjoint forcing calculation 
# Shannon Capps | January 29, 2015

# To test the MDA8 calculation 
days = 7

# Iterate over this period of time
for dd in range(0, days):
    day = 1 + dd      # specific to benchmarking episode that starts on 10th
    # Ensure that the string for day will match the CMAQ file
    if (day < 10):
        # Add zero before 
        small = str('0')
        dStamp = small+str(day)
    else:
        dStamp = str(day)
    
    # For each day, open the hourly average concentration file
    ### CHANGE: specify location of checkpoint files
    chemchkname = '/data/shannon/cmaqadj/2007/output/ACONC.200705'+dStamp
    chemchkfile = Dataset(chemchkname, mode='r', open=True)
    
    # Read from output of model the ozone concentration
    cohr = np.squeeze(chemchkfile.variables['O3'][:][:][:])
    
    # Combine the 
    # Select only layer 1 (only one that exists in ACONC file)
    cohr = cohr[:,:,:]
    if dd == 0:
        cO3hr = cohr 
    else:
        cO3hr = np.concatenate((cO3hr,cohr),axis=0)
    
    concO3hr = cO3hr
    ltime, lrow, lcol = (concO3hr.shape)

# Read in file that contains hours to shift to local time
tzpath = '/data/cmaq/12us2_time.dat'
tzgrid = np.genfromtxt(tzpath, dtype=[('COL', '<i8'),('ROW', '<i8'),('CR', '<i8'),('TIMEZN', 'S8'),('ZONE','<f8')], skip_header=1)

# Initialize arrays for shifting concentrations in time
concO3hrLST = np.zeros([ltime,lrow,lcol])
tzshiftgrid = np.zeros([lrow,lcol])

# Write time zone file to check
### CHANGE: specify location of adjoint forcing file (possibly same as checkpoint files)
tzshiftpath = '/data/shannon/cmaqadj/2007/tzshift.nc'
tzshiftfile = Dataset(tzshiftpath, 'w', format='NETCDF3_64BIT')

setattr(tzshiftfile,"NVARS",1)
setattr(tzshiftfile,"NLAYS",1) 
setattr(tzshiftfile,"UPNAM","TZSHIFT_FILE")
setattr(tzshiftfile,"VAR-LIST","HR              ")
setattr(tzshiftfile,"FILEDESC","Timezone shift file")

tzshiftfile.createDimension("ROW", lrow)
tzshiftfile.createDimension("COL", lcol)

tzshiftfile.sync()

tzshiftvar = tzshiftfile.createVariable('HR','f4',('ROW','COL'))
tzshiftvar[:] = np.zeros([lrow,lcol])
tzshiftvar[:] = np.reshape(tzshiftgrid,[lrow,lcol])

tzshiftfile.close()


# Apply to time index
for colind in range(0,lcol):
        for rowind in range(0,lrow):
            ind = 3*(colind) + 9*lcol*(rowind) + 2
            #print(ind, 3*rowind, 3*colind)
            tzshift = int(tzgrid[ind][4])
            #print(tzshift)
            tzshiftgrid[rowind,colind] = tzshift
            for timeind in range(0,ltime-9):
                if timeind+tzshift > 0:
                    # tzshift is negative, not positive
                    concO3hrLST[timeind+tzshift,rowind,colind] = concO3hr[timeind,rowind,colind]

# Calculate 8-hour average O3 concentration exceeding an average of 60 ppb (0.060 ppm) across whole time period
ltime, lrow, lcol = (concO3hrLST.shape)
day8hravgO3 = np.zeros([ltime, lrow, lcol])
for hr in range(1, ltime):
    endhr = 7 + hr
    day8hravgO3[hr,:,:] = np.mean(concO3hrLST[hr:endhr,:,:],axis=0)

# Select the maximum of the 8-hr average while in LST
ltime, lrow, lcol = (day8hravgO3.shape)
ldays = int(ltime / 24)
max8hravgO3init = np.zeros([ldays, lrow, lcol])
max8hravgO3zeros = np.zeros([lrow, lcol])
max8hravgO3hrs = np.zeros([ltime, lrow, lcol])
for day in range(0, ldays):
    daystrthr = day*24 + 1
    dayendhr = day*24 + 24
    max8hravgO3init[day,:,:] = np.argmax(day8hravgO3[daystrthr:dayendhr,:,:],axis=0)

for dd in range(1, int(ldays)):
    firsthr = 24*(dd - 1)
    lasthr  = 24*( dd )
    day8hravgO3hrly = concO3hrLST[firsthr:lasthr,:,:]
    day8hravgO3max = np.zeros([24, lrow, lcol])
    bin8hravgO3max = np.zeros([24, lrow, lcol])
    print(day8hravgO3max.shape)
    O38hravgmaxind = np.argmax(day8hravgO3hrly, axis=0)  # Select indices of daily max values
    for rowind in range(0,lrow):
        for colind in range(0,lcol):
            for hrind in range(0,24):
                if hrind == O38hravgmaxind[rowind,colind]:
                    bin8hravgO3max[hrind,rowind,colind] = 1
    day8hravgO3max = bin8hravgO3max * day8hravgO3hrly
    if dd == 1:
        O3month8hravgmx = day8hravgO3max
        O3month8hravgbin = bin8hravgO3max
    else:
        print(O3mon8hravgmx.shape, day8hravgO3max.shape)
        print(O3mn8hravgbin.shape, bin8hravgO3max.shape)
        O3month8hravgmx = np.concatenate((O3mon8hravgmx,day8hravgO3max),axis=0)
        O3month8hravgbin = np.concatenate((O3mn8hravgbin, bin8hravgO3max), axis=0)
    O3mon8hravgmx = O3month8hravgmx
    O3mn8hravgbin = O3month8hravgbin

O3month8hrmxavg = np.sum(O3mon8hravgmx, axis=0)/ldays*1000

# Write MDA8 average to check
## CHANGE: specify location of adjoint forcing file (possibly same as checkpoint files)
mda8avgpath = '/data/shannon/cmaqadj/2007/WkMayMDA8avg.nc'
mda8avgfile = Dataset(mda8avgpath, 'w', format='NETCDF3_64BIT')

setattr(mda8avgfile,"NVARS",1)
setattr(mda8avgfile,"NLAYS",1) 
setattr(mda8avgfile,"UPNAM","MDA8avg_FILE")
setattr(mda8avgfile,"VAR-LIST","MDA8              ")
setattr(mda8avgfile,"FILEDESC","Sample average MDA8 file (ppb)")

mda8avgfile.createDimension("ROW", lrow)
mda8avgfile.createDimension("COL", lcol)

mda8avgfile.sync()

mda8avgvar = mda8avgfile.createVariable('MDA8','f4',('ROW','COL'))
mda8avgvar[:] = np.zeros([lrow,lcol])
mda8avgvar[:] = np.reshape(O3month8hrmxavg,[lrow,lcol])

mda8avgfile.close()

max8hravgO3frcLST = np.zeros([ltime, lrow, lcol])
for day in range(0, ldays):
    for rw in range(0,lrow):
        for cl in range(0,lcol):
            for hrdy in range(0,24):
                hr = day*24 + hrdy
                if hrdy == int(max8hravgO3init[day,rw,cl]):
                    int(max8hravgO3init[day,rw,cl])
                    endhr = hr + 8
                    # The scaling of 2 is to increase the steepness of the curve about 60 ppb so that the
                    #     the influence is close to zero at 57.5 ppb and one at 62.5 ppb.
                    max8hravgO3frcLST[hr:endhr,rw,cl] = 1/8

np.sum(max8hravgO3frcLST)

max8hravgO3frc = np.zeros([ltime,lrow,lcol])
# Convert adjoint forcing back to CMAQ time (from LST)
for colind in range(0,lcol):
        for rowind in range(0,lrow):
            ind = 3*(colind) + 9*lcol*(rowind) + 2
            #print(ind, 3*rowind, 3*colind)
            tzshift = int(tzgrid[ind][4])
            #print(tzshift)
            tzshiftgrid[rowind,colind] = tzshift
            for timeind in range(0,ltime-9):
                if timeind+tzshift > 0:
                    # tzshift is negative, not positive
                    max8hravgO3frc[timeind,rowind,colind] = max8hravgO3frcLST[timeind+tzshift,rowind,colind]

# Write adjoint forcing to one file per day
# Specifically, the forcing corresponds to the maximum 8 hr average ozone.

# Iterate over this period of time
dayincr = 0
days = 7
for dd in range(1, days):
    # Ensure that the string for day will match the CMAQ file
    if (dd < 10):
        # Add zero before 
        small = str('0')
        dStamp = small+str(dd)
    else:
        dStamp = str(dd)
    
    # For each day, open the hourly average concentration file
    chemchkname = '/data/shannon/cmaqadj/2007/output/ACONC.200705'+dStamp
    chemchkfile = Dataset(chemchkname, mode='r', open=True)
    
    # For each day, open the hourly average concentration file
    forcO38hmxpath = '/data/shannon/cmaqadj/2007/output/ADJ_FORCE.200705'+dStamp
    forcO38hmxfile = Dataset(forcO38hmxpath, mode='w', format='NETCDF3_64BIT')
    
    # I/O API attributes to be included in forcing file (originally coded by M.Russell / Carleton.U)
    attrs=["IOAPI_VERSION", "EXEC_ID", "FTYPE", "CDATE", "CTIME", "WDATE", "WTIME", "SDATE", "STIME", "TSTEP", "NTHIK", "NCOLS", "NROWS", "GDTYP", "P_ALP", "P_BET", "P_GAM", "XCENT", "YCENT", "XORIG", "YORIG", "XCELL", "YCELL", "VGTYP", "VGTOP", "VGLVLS", "GDNAM", "HISTORY"]
    
    for attr in attrs:
    # Read the attribute
        if hasattr(chemchkfile, attr): 
            attrVal = getattr(chemchkfile, attr);
            # Write it to the forcing file
            setattr(forcO38hmxfile, attr, attrVal)
    
    # I/O API attributes to be included in forcing file that do not match the output file (relevant to adjoint forcing specifically)
    setattr(forcO38hmxfile,"NVARS",1)
    setattr(forcO38hmxfile,"NLAYS",1) 
    setattr(forcO38hmxfile,"UPNAM","RD_FORCE_FILE")
    setattr(forcO38hmxfile,"VAR-LIST","O3              ")
    setattr(forcO38hmxfile,"FILEDESC","Adjoint forcing file. Cost function: daily maximum 8 hour average ozone (ppm)")
    
    ltime, lrow, lcol = (max8hravgO3frc.shape)
    llay = 1
    ldattim = 2
    ltime = 24
    strthr = dayincr*24
    endhr = (dayincr + 1)*24
    
    forcO38hmxfile.createDimension("TSTEP", None)
    forcO38hmxfile.createDimension("DATE-TIME", ldattim)
    forcO38hmxfile.createDimension("LAY", llay)
    forcO38hmxfile.createDimension("VAR", 1)
    forcO38hmxfile.createDimension("ROW", lrow)
    forcO38hmxfile.createDimension("COL", lcol)
    
    forcO38hmxfile.sync()
    
    dattim = np.squeeze(chemchkfile.variables['TFLAG'][:][:])
    dattim = dattim[0:ltime,1,:]
    
    # Create variables
    forc_O38hmx_tflag = forcO38hmxfile.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
    forc_O38hmx_tflag[:] = np.zeros([ltime,llay,ldattim])
    forc_O38hmx_tflag[:] = np.reshape(dattim,[ltime,llay,ldattim])
    
    varattrs=["long_name","units","var_desc"]
    for varattr in varattrs:
        # Read the attribute
        if hasattr(chemchkfile.variables['TFLAG'], varattr): 
            varattrVal = getattr(chemchkfile.variables['TFLAG'], varattr);
            # Write it to the forcing file
            setattr(forc_O38hmx_tflag, varattr, varattrVal)
    
    forc_O38hmx = forcO38hmxfile.createVariable('O3','f4',('TSTEP','LAY','ROW','COL'))
    forc_O38hmx[:] = np.zeros([ltime,llay,lrow,lcol])
    forc_O38hmx[:] = np.reshape(max8hravgO3frc[strthr:endhr,:,:],[ltime,llay,lrow,lcol])
    
    varattrs=["long_name","units","var_desc"]
    for varattr in varattrs:
        # Read the attribute
        if hasattr(chemchkfile.variables['O3'], varattr): 
            varattrVal = getattr(chemchkfile.variables['O3'], varattr);
            # Write it to the forcing file
            setattr(forc_O38hmx, varattr, varattrVal)
    
    forcO38hmxfile.close()
    dayincr = 1 + dayincr



