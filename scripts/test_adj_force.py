#! /apps/python/2.7.6/bin/python2.7

from netCDF4 import Dataset
import numpy as np

import sys

dates = [ str(20070600 + i) for i in range(10,14) ]

for dstamp in dates:

    #?using CONC or ACONC?
    conc_name = '/home/563/spt563/cmaq_adj/CMAQadjBnmkData/spt_output/CONC.' + dstamp
    conc_file = Dataset(conc_name, mode='r', open=True)

    o3hr = np.squeeze(conc_file.variables['O3'][:][:][:])

    dayo3hrly = o3hr[0:25, 1, :, :]
    binO3max = np.where(dayo3hrly[:,:,:] == np.max(dayo3hrly[:,:,:], axis=0), 1.0, 0.0)
    
    forcO3mxhrpath = '/home/563/spt563/cmaq_adj/CMAQadjBnmkData/spt_output/ADJ_FORCE.' + dstamp
    forcO3mxhrfile = Dataset(forcO3mxhrpath, 'w', format='NETCDF3_64BIT')

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
            setattr(forcO3mxhrfile, attr, attr_val)

    setattr(forcO3mxhrfile, "NVARS", 1)
    setattr(forcO3mxhrfile, "NLAYS", 1)
    setattr(forcO3mxhrfile, "UPNAM", "RD_FORCE_FILE")
    setattr(forcO3mxhrfile, "VAR-LIST", "O3             ")
    setattr(forcO3mxhrfile, "FILEDESC", "Adjoint forcing file. Cost function: sum of max O# in every cell, each day (ppm)")

    ltime, lrow, lcol = (binO3max.shape)
    llay = 1
    ldattim = 2
    ltime = 25

    forcO3mxhrfile.createDimension("TSTEP", None)
    forcO3mxhrfile.createDimension("DATE-TIME", ldattim)
    forcO3mxhrfile.createDimension("LAY", llay)
    forcO3mxhrfile.createDimension("VAR", 1)
    forcO3mxhrfile.createDimension("ROW", lrow)
    forcO3mxhrfile.createDimension("COL", lcol)

    forcO3mxhrfile.sync()

    dattim = np.squeeze(conc_file.variables['TFLAG'][:][:])
    dattim = dattim[0:ltime, 1, :]

    forc_tflag = forcO3mxhrfile.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
    forc_tflag[:] = dattim.reshape(ltime, llay, ldattim)

    forc_O3mxhr = forcO3mxhrfile.createVariable('O3', 'f4', ('TSTEP', 'LAY', 'ROW', 'COL'))
    forc_O3mxhr[:] = binO3max.reshape(ltime, llay, lrow, lcol)

    varattrs = ["long_name", "units", "var_desc"]
    for varattr in varattrs:
        if hasattr(conc_file.variables['TFLAG'], varattr):
            varattr_val = getattr(conc_file.variables['TFLAG'], varattr)
            setattr(forc_tflag, varattr, varattr_val)
        if hasattr(conc_file.variables['O3'], varattr):
            varattr_val = getattr(conc_file.variables['O3'], varattr)
            setattr(forc_O3mxhr, varattr, varattr_val)
        else:
            print 'O3 has no attr', varattr

    forcO3mxhrfile.close()
    conc_file.close()
