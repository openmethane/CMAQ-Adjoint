#! /apps/python/2.7.6/bin/python2.7

from netCDF4 import Dataset
import numpy as np
import os.path

import sys

dates = [ str(20070600 + i) for i in range(10,14) ]

for dstamp in dates:

    #?using CONC or ACONC?
    conc_name = '/home/563/spt563/cmaq_adj/CMAQadjBnmkData/GHG_output/CONC.' + dstamp
    if not os.path.isfile(conc_name):
        print 'cannot find', conc_name
        continue
    conc_file = Dataset(conc_name, mode='r', open=True)

    COhr = np.squeeze(conc_file.variables['CO2'][:][:][:])

    dayCOhrly = COhr[:, 0:2, :, :]
    #binCOmax = np.where(dayCOhrly[:,:,:] == np.max(dayCOhrly[:,:,:], axis=0), 1.0, 0.0)
    binCOmax = np.zeros(dayCOhrly.shape)
    if dstamp == "20070610":
#        binCOmax[12,12,12] = 1.0
#        binCOmax[1,12,12] = 1.0
#        s_hat = np.array([ [.25,.25,.25], [.25,.5,.25], [.25,.25,.25] ])
#        binCOmax[6, 1, 11:14, 11:14] = s_hat
#        binCOmax[5, 0, 11:14, 11:14] = s_hat
#        binCOmax[6, 0, 11:14, 11:14] = 2 * s_hat
#        binCOmax[7, 0, 11:14, 11:14] = s_hat
        binCOmax[6, 0, 11:14, 11:14] = 1.0
#        binCOmax[21, 0, :, :] = 1.0
    
    forcCOmxhrpath = '/home/563/spt563/cmaq_adj/CMAQadjBnmkData/GHG_output/ADJ_FORCE.' + dstamp
    forcCOmxhrfile = Dataset(forcCOmxhrpath, 'w', format='NETCDF3_64BIT')

    attrs = ["IOAPI_VERSION", "EXEC_ID", "FTYPE",
             "CDATE", "CTIME", "WDATE", "WTIME",
             "SDATE", "STIME", "TSTEP", "NTHIK",
             "NCOLS", "NROWS", "GDTYP", "P_ALP",
             "P_BET", "P_GAM", "XCENT", "YCENT",
             "XORIG", "YORIG", "XCELL", "YCELL",
             "VGTYP", "VGTOP", "VGLVLS", "GDNAM", "HISTORY"]

    for attr in attrs:
        if hasattr(conc_file, attr):
            attr_val = getattr(conc_file, attr)
            setattr(forcCOmxhrfile, attr, attr_val)

    setattr(forcCOmxhrfile, "NVARS", 1)
#    setattr(forcCOmxhrfile, "NLAYS", 1)
    setattr(forcCOmxhrfile, "NLAYS", 2)
    setattr(forcCOmxhrfile, "UPNAM", "RD_FORCE_FILE")
    setattr(forcCOmxhrfile, "VAR-LIST", "CO2            ")
    setattr(forcCOmxhrfile, "FILEDESC", "Adjoint forcing file. Cost function: sum of max CO2 in every cell, each day (ppm)")

    ltime, llay, lrow, lcol = (binCOmax.shape)
#    llay = 1
    ldattim = 2
    #ltime = 25

    forcCOmxhrfile.createDimension("TSTEP", None)
    forcCOmxhrfile.createDimension("DATE-TIME", ldattim)
    forcCOmxhrfile.createDimension("LAY", llay)
    forcCOmxhrfile.createDimension("VAR", 1)
    forcCOmxhrfile.createDimension("ROW", lrow)
    forcCOmxhrfile.createDimension("COL", lcol)

    forcCOmxhrfile.sync()

    dattim = np.squeeze(conc_file.variables['TFLAG'][:][:])
    #dattim = dattim[0:ltime, 1, :]
    dattime = dattim[0:ltime, :]

    forc_tflag = forcCOmxhrfile.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
    forc_tflag[:] = dattim.reshape(ltime, 1, ldattim)

    forc_COmxhr = forcCOmxhrfile.createVariable('CO2', 'f4', ('TSTEP', 'LAY', 'ROW', 'COL'))
    forc_COmxhr[:] = binCOmax.reshape(ltime, llay, lrow, lcol)

    varattrs = ["long_name", "units", "var_desc"]
    for varattr in varattrs:
        if hasattr(conc_file.variables['TFLAG'], varattr):
            varattr_val = getattr(conc_file.variables['TFLAG'], varattr)
            setattr(forc_tflag, varattr, varattr_val)
        if hasattr(conc_file.variables['CO2'], varattr):
            varattr_val = getattr(conc_file.variables['CO2'], varattr)
            setattr(forc_COmxhr, varattr, varattr_val)
        else:
            print 'CO2 has no attr', varattr

    forcCOmxhrfile.close()
    conc_file.close()
