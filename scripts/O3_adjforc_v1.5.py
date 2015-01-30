#!/usr/local/apps/Python/2.7.3/bin/python
###/Library/Frameworks/Python.framework/Versions/3.3/bin/python3

# Adjoint forcing calculator for CMAQ adjoint
#    Calculates derivative of ecosystem cost functions w.r.t. hourly ozone
#
# S. Capps - Initial commit July 2014 
# v.1.1 Include crop production and timber biomass in cost functions dependent on W126
# v.1.3 Changed the crop and timber forcing to depend on cumulative W126 rather than hourly O3
# v.1.4 Converted the timber biomass forcing to tons per grid cell
# v.1.5 Converted the crop biomass forcing calculation to match the timber method, eliminating NaNs


from netCDF4 import Dataset
from math import fsum, exp
from random import uniform
import numpy as np

class health:
    def O3window:
        # To use Jerrett et al. (2009), the O3 concentration of relevance is the
        #   maximum 6-month mean of 1-h daily max O3 in year
        #   Here, calculate which 6 month period in 2007 qualifies.
        
        seas = [5,6,7,8,9,10]
        # Calculate the hourly max O3 concentration
        for mn in seas:
            o3LST = mn_o3(mn)
            ltime, lrow, lcol = (o3LST.shape)
            mndays = ltime / 24 + 1
            for dd in range(1, mndays):
                firsthr = 24*(dd - 1)
                lasthr  = 24*( dd )
                dayO3hrly = o3LST[firsthr:lasthr,:,:]
                dayO3mx = np.zeros([24, lrow, lcol])
                binO3max = np.zeros([24,lrow, lcol])
                print dayO3hrly.shape
                O3maxind = np.argmax(dayO3hrly, axis=0)  # Select indices of daily max values
                for rowind in range(0,246):
                    for colind in range(0,396):
                        for hrind in range(0,24):
                            if hrind == O3maxind[rowind,colind]:
                                binO3max[hrind,rowind,colind] = 1
        
                dayO3mx = binO3max * dayO3hrly
                if dd == 1:
                    O3monthmx = dayO3mx
                else:
                    O3monthmx = np.concatenate((O3monmx,dayO3mx),axis=0)
                O3monmx = O3monthmx

            if mn == seas[0]:
                O3seasonhrmx = O3monmx
            else:
                O3seasonhrmx = np.concatenate((O3seashrmax,O3monmx),axis=0)
            O3seashrmax = O3seasonhrmx


        # Produce the tflag variable
        for mn in range(5,11):
            tflagseas = mn_tflag(mn)
            if mn == seas[0]:
                tflagall = tflagseas
            else:
                tflagall = np.concatenate((tflagseasall,tflagseas), axis = 0)
            tflagseasall = tflagall

        # Produce 6-month average
        # ppm - average accounting only for the filled values (1 / day)
        #sixmnmean = np.sum(O3seashrmax) / ( np.prod(O3seashrmax.shape) / 24 )
        lhrs, lrow, lcol = O3seashrmax.shape
        ldays = lhrs/24
        sixmnHrMaxMeanppb = np.zeros([lrow,lcol])
        # Need average of six month daily 1-hr max in each grid cell
        sixmnHrMaxMeanppb = 1000*np.sum(O3seashrmax,axis=0) / ldays

        # Also need to know which hours had maximum daily 1-hr concentration
        ltime, lrow, lcol = O3seashrmax.shape
        seasO3bin = np.zeros([ltime, lrow, lcol])
        seasones = np.ones([ltime, lrow, lcol])
        seaszeros = np.zeros([ltime, lrow, lcol])
        seasO3bin = np.where(O3seashrmax > 0.0, seasones, seaszeros)

        # Write to process separately
        seasO3binpath = '/work/CLIMSIM/slc/data/cost_fn/force_mort_07_max_bin.ncf'
        seasO3binfile = Dataset(seasO3binpath, 'w', format='NETCDF3_64BIT')

        # Make dimensions for file
        seasO3binfile.createDimension('rows',lrow)
        seasO3binfile.createDimension('cols',lcol)
        seasO3binfile.createDimension('allhrs',ltime)

        seasO3binfile.sync()

        # Create variables
        forc_mort_bin = seasO3binfile.createVariable('mort_bin','f8',('allhrs','rows','cols'))
        forc_mort_bin[:] = seasO3bin

        sixmeanmaxhr = seasO3binfile.createVariable('MeanMaxHr','f8',('rows','cols'))
        sixmeanmaxhr[:] = sixmnHrMaxMeanppb

        seasO3binfile.close()
  
        ### # Find index of maximum 6-month average; print time window associated.      
        # maxsixmnmean = np.argmax(sixmnmean)              
        ### Based on OAQPS hourly max values in LST, the six month mean of 1 hr max values is greatest May 7 - Oct 7 (days 127-307)
        ### This encompasses the JJA window for W126. 

    def mort_adj:
        # Calculate forcing based on ozone respiratory illness mortality from Jerrett et al., 2009:
        #    dJ = ( M_0 * PopOver30 ) * (B*exp(-B*C))

        # Load the array representing hours at which max occurred (0 - no max; 1 - max)
        # Shape of mort_bin array is the (hours.in.6months, rows, cols)
        O3Binfilepath = '/work/CLIMSIM/slc/data/cost_fn/force_mort_07_max_bin.ncf'
        O3Binfile = Dataset(O3Binfilepath, mode='r', open=True)
        O3Bin = O3Binfile.variables['mort_bin'][:]

        # Load the gridded baseline mortality * population data
        BMortPopfilepath = '/home/slc/cmaq_forcing/forcing/src/auxiliary_files/Mort_2010_12km.ncf'
        BMortPopfile = Dataset(BMortPopfilepath, mode='r', open=True)
        BMort = BMortPopfile.variables['MORT'][:]   # Rate of death / yr * population
        # Pop = BMortPopfile.variables['POP'][:]

        # Adjust the length of time over which the cost function will be applied as writing the forcing file
        #    Period over which the adjoint forcing will be applied
        forc_time = [6,7,8]
        forcdays = dys_in_mns( forc_time )
        forchrs = forcdays * 24

        #    Period over which the adjoint forcing will be zero  while the adjoint model continues to run
        neglect_time = [5]
        neglectdays = dys_in_mns( neglect_time )
        neglecthrs = neglectdays * 24

        tothrs = neglecthrs + forchrs

        # Beta determined from Jerrett et al., 2009 
        # (fractional increase in mortality due to respiratory illness / ppb O3 metric)
        beta = 0.04 / 10

        # Calculate the hourly adjoint forcing for mortality 
        #     Only force in hour that max value occurred (seasO3bin contains 0 for non-max hour, 1 for max hour)
        max_bin = readjust_time(O3Bin)
        ltime, lrow, lcol = max_bin.shape
        
        #     Determine divisor for specific 6-month season
        #     Note: this value does not change even if the period of adjoint forcing is shortened
        season = [5,6,7,8,9,10]
        seasdays = dys_in_mns( season )

        # Correct if MORT only were the rate of deaths and not multiplied by population
        # mort_adj_seas = ( BMort * Pop * np.exp(-beta * sixmnHrMaxMeanppb ) ) * beta / ( seasdays ) * max_bin
        # Use this since MORT includes the rate of death * population
        mort_adj_seas = ( BMort * np.exp(-beta * sixmnHrMaxMeanppb ) ) * beta / ( seasdays ) * max_bin

        # Reading I/O API attributes from this file
        aconcpath = '/work/ROMO/2007platform/CMAQv4.7.1/2007ee_ORDBC_v5_07c_v4.7.1_N5ao/12US2/extr/O3/2007ee_ORDBC_v5_07c_v4.7.1_N5ao.12US2_24.combine.O3.05'
        aconcfile = Dataset(aconcpath, mode='r', open=True)

        # Write to process separately (if desired)
        mortO3binpath = '/work/CLIMSIM/slc/data/cost_fn/force_mort_07.ncf'
        mortO3file = Dataset(mortO3binpath, 'w', format='NETCDF3_64BIT')

        # I/O Api attributes to be included in forcing file (originally coded by M.Russell / Carleton.U)
        attrs=["IOAPI_VERSION", "EXEC_ID", "FTYPE", "CDATE", "CTIME", "WDATE", "WTIME", "SDATE", "STIME", "TSTEP", "NTHIK", "NCOLS", "NROWS", "NLAYS", "NVARS", "GDTYP", "P_ALP", "P_BET", "P_GAM", "XCENT", "YCENT", "XORIG", "YORIG", "XCELL", "YCELL", "VGTYP", "VGTOP", "VGLVLS", "GDNAM", "UPNAM", "VAR-LIST", "FILEDESC", "HISTORY"]

        for attr in attrs:
            # Read the attribute
            if hasattr(aconcfile, attr): 
                attrVal = getattr(aconcfile, attr);
                # Write it to the forcing file
                setattr(mortO3file, attr, attrVal)
        
        # Make dimensions for IOAPI file creation
        llay = 1
        ldattim = 2
        ltime = tothrs-neglecthrs
        ltottim = tothrs

        # I/O API dimnensions
        mortO3file.createDimension("TSTEP", None)
        mortO3file.createDimension("DATE-TIME", ldattim)
        mortO3file.createDimension("LAY", 1)
        mortO3file.createDimension("VAR", 1)
        mortO3file.createDimension("ROW", lrow)
        mortO3file.createDimension("COL", lcol)

        mortO3file.sync()

        # Provide arrays for filling IOAPI variable 
        mort_adj_episode = np.zeros([1, tothrs, lrow, lcol])
        mort_adj_episode[0][neglecthrs:tothrs] = mort_adj_seas[:,neglecthrs:tothrs,:,:]
        tflag_episode = tflagseasall[0:tothrs,:]

        # Create TFLAG variable
        forc_mort_tflag = mortO3file.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))

        varattrs=["long_name","units","var_desc"]
        for varattr in varattrs:
            # Read the attribute
            if hasattr(aconcfile.variables['TFLAG'], varattr): 
                varattrVal = getattr(aconcfile.variables['TFLAG'], varattr);
                # Write it to the forcing file
                setattr(forc_mort_tflag, varattr, varattrVal)

        forc_mort_tflag[:] = np.reshape(tflag_episode,[tothrs,llay,ldattim])
        forc_mort_tflag.sync()
        # Create forcing variable
        forc_mort_bin = mortO3file.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        forc_mort_bin[:] = np.reshape(mort_adj_episode,[2952,1,246,396])

        varattrs=["long_name","units","var_desc"]
        for varattr in varattrs:
            # Read the attribute
            if hasattr(aconcfile.variables['O3'], varattr): 
                varattrVal = getattr(aconcfile.variables['O3'], varattr);
                # Write it to the forcing file
                setattr(forc_mort_bin, varattr, varattrVal)


        mortO3file.sync()
        mortO3file.close()

        mort_adj_episode_avg = (np.sum(mort_adj_episode,axis=0))/92
        mortadjepipath = '/work/CLIMSIM/slc/data/cost_fn/force_mort_episode_avg.ncf'
        mortadjepifile = Dataset(mortadjepipath, 'w', format='NETCDF3_64BIT')

        # For plotting purposes only
        for attr in attrs:
            # Read the attribute
            if hasattr(aconcfile, attr): 
                attrVal = getattr(aconcfile, attr);
                # Write it to the forcing file
                setattr(mortadjepifile, attr, attrVal)
         
        # I/O API dimnensions
        #dims = aconcfile.dimensions.keys()
        #for d in range(1,len(dims)):
        #    # Copy the dimension straight from the source
        #    v = aconcfile.dimensions[dims[d]]
        #    try:
        #        mortadjepifile.createDimension(dims[d], len(v))
        #    except IOError as ex:
        #        print "Cannot create dimension %s"%dims[d], ex

        mortadjepifile.createDimension(dims[0], 1)

        # Make dimensions for file
        mortadjepifile.createDimension('rows',lrow)
        mortadjepifile.createDimension('cols',lcol)

        # Create variables
        forc_mort_bin_avg = mortadjepifile.createVariable('mort_adj_epi_avg','f8',('rows','cols'))
        forc_mort_bin_avg[:] = mort_adj_episode_avg
        mortadjepifile.close()

class ecosystem:
    def crop_ryl_adj (seas):

        ### Adjoint forcing
    
        # Load file that contains crop yield data
        #    File produced from BELD land use data & NASS crop production data
        #       by Rob Pinder (US EPA) (August 2013) 
        yield_filepath = '/home/slc/cmaq_forcing/auxiliary/yield_BELD_NASS_12US2_396X246.ncf'
        yieldcropfile = Dataset(yield_filepath,'r')
    
        # Read in crop yield data
        yield_crn = np.squeeze(yieldcropfile.variables['corn'])
        yield_ctn = np.squeeze(yieldcropfile.variables['cotton'])
        yield_pto = np.squeeze(yieldcropfile.variables['potato'])
        yield_soy = np.squeeze(yieldcropfile.variables['soybean'])
        yield_wht = np.squeeze(yieldcropfile.variables['wheat'])
        nrows, ncols = yield_crn.shape 
        wholegrid = np.ones([nrows, ncols])
    
    	# Calculate hours in ozone season (highest W126 values)
    	#    W126 metric is based on 90-day accumulation of ozone concentrations 

        # Read O3 values and concatenate into one array
        seas = [5,6,7,8]
        grid_zeros = np.zeros([nrows,ncols])
        O3seas = seasO3(seas)

        # Shift O3 to proper CMAQ time and
        #    zero out hours other than 8 am - 8 pm 
        O3seasshift = w126prep ( O3seas )
    
        # Assign O3 values 
        ltime,lrow,lcol = (O3seasshift.shape)

        # Adjust the length of time over which the cost function will be applied as writing the forcing file
        forc_time = [6,7,8]
        forcdays = dys_in_mns( forc_time )
        forchrs = forcdays * 24

        neglect_time = [5]
        neglectdays = dys_in_mns( neglect_time )
        neglecthrs = neglectdays * 24

        tothrs = neglecthrs + forchrs
    
        # Produce the tflag variable
        for mn in range(5,11):
            tflagseas = mn_tflag(mn)
            if mn == seas[0]:
                tflagall = tflagseas
            else:
                tflagall = np.concatenate((tflagseasall,tflagseas), axis = 0)
            tflagseasall = tflagall

        tflagseasall = tflagseasall[0:tothrs,:]

        # Set parameters for making dimensions
        llay = 1
        ldattim = 2
        
        # Create file
        forcfilepath_crn = '/work/CLIMSIM/slc/data/cost_fn/force_crop_07_12US2_prod_wholeW126_crn.ncf'
        forccropdist_crn = Dataset(forcfilepath_crn, 'w', format='NETCDF3_64BIT')
        forcfilepath_ctn = '/work/CLIMSIM/slc/data/cost_fn/force_crop_07_12US2_prod_wholeW126_ctn.ncf'
        forccropdist_ctn = Dataset(forcfilepath_ctn, 'w', format='NETCDF3_64BIT')
        forcfilepath_pto = '/work/CLIMSIM/slc/data/cost_fn/force_crop_07_12US2_prod_wholeW126_pto.ncf'
        forccropdist_pto = Dataset(forcfilepath_pto, 'w', format='NETCDF3_64BIT')
        forcfilepath_soy = '/work/CLIMSIM/slc/data/cost_fn/force_crop_07_12US2_prod_wholeW126_soy.ncf'
        forccropdist_soy = Dataset(forcfilepath_soy, 'w', format='NETCDF3_64BIT')
        forcfilepath_wht = '/work/CLIMSIM/slc/data/cost_fn/force_crop_07_12US2_prod_wholeW126_wht.ncf'
        forccropdist_wht = Dataset(forcfilepath_wht, 'w', format='NETCDF3_64BIT')

        # Reading I/O API attributes from this file
        aconcpath = '/work/ROMO/2007platform/CMAQv4.7.1/2007ee_ORDBC_v5_07c_v4.7.1_N5ao/12US2/extr/O3/2007ee_ORDBC_v5_07c_v4.7.1_N5ao.12US2_24.combine.O3.05'
        aconcfile = Dataset(aconcpath, mode='r', open=True)

        # I/O API attributes (originally coded by M.Russell / Carleton.U)
        attrs=["IOAPI_VERSION", "EXEC_ID", "FTYPE", "CDATE", "CTIME", "WDATE", "WTIME", "SDATE", "STIME", "TSTEP", "NTHIK", "NCOLS", "NROWS", "NLAYS", "NVARS", "GDTYP", "P_ALP", "P_BET", "P_GAM", "XCENT", "YCENT", "XORIG", "YORIG", "XCELL", "YCELL", "VGTYP", "VGTOP", "VGLVLS", "GDNAM", "UPNAM", "VAR-LIST", "FILEDESC", "HISTORY"]


        for attr in attrs:
            # Read the attribute
            if hasattr(aconcfile, attr): 
                attrVal = getattr(aconcfile, attr);
                # Write it to the forcing file
                setattr(forccropdist_crn, attr, attrVal)
                setattr(forccropdist_ctn, attr, attrVal)
                setattr(forccropdist_pto, attr, attrVal)
                setattr(forccropdist_soy, attr, attrVal)
                setattr(forccropdist_wht, attr, attrVal)
        
#        # I/O API dimnensions
#        dims = aconcfile.dimensions.keys()
#        for d in range(1,len(dims)):
#            # Copy the dimension straight from the source
#            v = aconcfile.dimensions[dims[d]]
#            try:
#                forccropdist.createDimension(dims[d], len(v))
#            except IOError as ex:
#                print "Cannot create dimension %s"%dims[d], ex
#
#        # Add because no O3 files are of the correct length of time
#        forccropdist.createDimension(dims[0], tothrs)
#        
#        forccropdist.sync()

        # I/O API dimnensions
        forccropdist_crn.createDimension("TSTEP", None)
        forccropdist_ctn.createDimension("TSTEP", None)
        forccropdist_pto.createDimension("TSTEP", None)
        forccropdist_soy.createDimension("TSTEP", None)
        forccropdist_wht.createDimension("TSTEP", None)

        forccropdist_crn.createDimension("DATE-TIME", ldattim)
        forccropdist_ctn.createDimension("DATE-TIME", ldattim)
        forccropdist_pto.createDimension("DATE-TIME", ldattim)
        forccropdist_soy.createDimension("DATE-TIME", ldattim)
        forccropdist_wht.createDimension("DATE-TIME", ldattim)

        forccropdist_crn.createDimension("LAY", 1)
        forccropdist_ctn.createDimension("LAY", 1)
        forccropdist_pto.createDimension("LAY", 1)
        forccropdist_soy.createDimension("LAY", 1)
        forccropdist_wht.createDimension("LAY", 1)

        forccropdist_crn.createDimension("VAR", 1)
        forccropdist_ctn.createDimension("VAR", 1)
        forccropdist_pto.createDimension("VAR", 1)
        forccropdist_soy.createDimension("VAR", 1)
        forccropdist_wht.createDimension("VAR", 1)

        forccropdist_crn.createDimension("ROW", lrow)
        forccropdist_ctn.createDimension("ROW", lrow)
        forccropdist_pto.createDimension("ROW", lrow)
        forccropdist_soy.createDimension("ROW", lrow)
        forccropdist_wht.createDimension("ROW", lrow)

        forccropdist_crn.createDimension("COL", lcol)
        forccropdist_ctn.createDimension("COL", lcol)
        forccropdist_pto.createDimension("COL", lcol)
        forccropdist_soy.createDimension("COL", lcol)
        forccropdist_wht.createDimension("COL", lcol)

        forccropdist_crn.sync()
        forccropdist_ctn.sync()
        forccropdist_pto.sync()
        forccropdist_soy.sync()
        forccropdist_wht.sync()

        # Create TFLAG variable
        forc_crop_crn_tflag = forccropdist_crn.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        forc_crop_ctn_tflag = forccropdist_ctn.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        forc_crop_pto_tflag = forccropdist_pto.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        forc_crop_soy_tflag = forccropdist_soy.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        forc_crop_wht_tflag = forccropdist_wht.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))

        varattrs=["long_name","units","var_desc"]
        for varattr in varattrs:
            # Read the attribute
            if hasattr(aconcfile.variables['TFLAG'], varattr): 
                varattrVal = getattr(aconcfile.variables['TFLAG'], varattr);
                # Write it to the forcing file
                setattr(forc_crop_crn_tflag, varattr, varattrVal)
                setattr(forc_crop_ctn_tflag, varattr, varattrVal)
                setattr(forc_crop_pto_tflag, varattr, varattrVal)
                setattr(forc_crop_soy_tflag, varattr, varattrVal)
                setattr(forc_crop_wht_tflag, varattr, varattrVal)

        forc_crop_crn_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])
        forc_crop_ctn_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])
        forc_crop_pto_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])
        forc_crop_soy_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])
        forc_crop_wht_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])

        forccropdist_crn.sync()
        forccropdist_ctn.sync()
        forccropdist_pto.sync()
        forccropdist_soy.sync()
        forccropdist_wht.sync()
 
 
        # Algorithm for relative yield forcing each grid cell at each hour between between 8:00-20:00 LST
        #
        #    d(YL)        d(W126)    d(RYL)    d(YL)     
        #  ---------- = ---------- x ------- x -------
        #  d(O3 conc)   d(O3 conc)   d(W126)   d(RYL)
        #
        # Enter equation below in WolframAlpha where r, q, A, & B are defined as in code / comments
        #      based on the assumption that x, r, q, A, & B are positive.
        #
        #   dJ/dC_i = d/dx Y*(1-exp(-((x/(1+r exp(-q x)))/A)^B)) =
        #
        #        = Y*-(B A^(-B) x^(B-1) (q r x+e^(q x)+r) (r e^(-q x)+1)^(-B) e^(-A^(-B) x^B (r e^(-q x)+1)^(-B)))/(e^(q x)+r)
        #   
        #		x -> O3 conc (ppm) between 8:00-20:00 LST
        #       r -> W126 weighting
        #       q -> "
        #		A -> crop specific for relative yield calculation (see table above)
        #		B -> "
        #       Y -> Crop yield specific to each grid cell and type
        
        r = 4403.0  # Constants for W126 weighting
        q = 126.0   # where ozone concentration is in ppm
        
        w126seas = np.zeros([ltime, lrow, lcol])
        w126seas = (O3seasshift / (1.0+r*np.exp(-q*O3seasshift)))
        w126seas_cum = np.sum(w126seas[neglecthrs:tothrs], axis = 0)

        # Parameterizations of crop yield loss from Lehrer et al., EPA 452/R-07-002 (2007)
        Acr = np.array([94.4, 98.3, 99.5, 110.0, 53.7])
        Bcr = np.array([1.572, 2.973, 1.242, 1.367, 2.391])
            
        # Set parameter for number of crops being treated. 
        cropspecies = 5   
        forc_crop = np.zeros([cropspecies,tothrs, lrow, lcol])

        O3seasshift_forc = O3seasshift[neglecthrs:tothrs,:,:]
    
        # Calculate the adjoint forcing for all of the months
        delW126_delO3 = ((np.exp(q*O3seasshift_forc) * (r*q*O3seasshift_forc + r + np.exp(q*O3seasshift_forc))) / (r + np.exp(q*O3seasshift_forc))**2 )
        delW126_delO3 = np.where(O3seasshift_forc > 0.0, delW126_delO3, O3seasshift_forc)
        np.isnan(np.min(delW126_delO3))

        # Calculate the adjoint forcing for all of the months
        for ttyp in range(cropspecies):
            forc_crop[ttyp][neglecthrs:tothrs] = (Bcr[ttyp] * Acr[ttyp]**(-Bcr[ttyp]) * w126seas_cum**(Bcr[ttyp]-1.0) * np.exp(-(np.divide(Acr[ttyp],w126seas_cum)**(-Bcr[ttyp])))) * delW126_delO3
            
        # Scale impact by the crop yield in each cell
        
        # Create variables
        forc_crn = forccropdist_crn.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        # Scale by yield
        adj_crn = forc_crop[0]*yield_crn
        # Reshape to read in correctly as adjoint forcing
        forc_crn[:] = np.reshape(adj_crn,[2952,1,246,396])

        forc_ctn = forccropdist_ctn.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        adj_ctn = forc_crop[1]*yield_ctn
        forc_ctn[:] = np.reshape(adj_ctn,[2952,1,246,396])

        forc_pto = forccropdist_pto.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        adj_pto =  forc_crop[2]*yield_pto
        forc_pto[:] = np.reshape(adj_pto,[2952,1,246,396])

        forc_soy = forccropdist_soy.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        adj_soy = forc_crop[3]*yield_soy
        forc_soy[:] = np.reshape(adj_soy,[2952,1,246,396])

        forc_wht = forccropdist_wht.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        adj_wht = forc_crop[4]*yield_wht
        forc_wht[:] = np.reshape(adj_wht,[2952,1,246,396])

        varattrs=["long_name","units","var_desc"]
        for varattr in varattrs:
            # Read the attribute
            if hasattr(aconcfile.variables['O3'], varattr): 
                varattrVal = getattr(aconcfile.variables['O3'], varattr);
                # Write it to the forcing file
                setattr(forc_crn, varattr, varattrVal)
                setattr(forc_ctn, varattr, varattrVal)
                setattr(forc_pto, varattr, varattrVal)
                setattr(forc_soy, varattr, varattrVal)
                setattr(forc_wht, varattr, varattrVal)
        
        forccropdist_crn.close()
        forccropdist_ctn.close()
        forccropdist_pto.close()
        forccropdist_soy.close()
        forccropdist_wht.close()
            		
    	#### Checking forcing file - debug only
        ##forcckfilepath = '/work/CLIMSIM/slc/data/cost_fn/force_crop_jun07_12US2.ncf'
        ##forccropcheck = Dataset(forcckfilepath, 'r')
        ##frcchk_crn = forccropcheck.variables['crn_adj'][:]
        ##frcchk_ctn = forccropcheck.variables['ctn_adj'][:]
        ##frcchk_pto = forccropcheck.variables['pto_adj'][:]
        ##frcchk_soy = forccropcheck.variables['soy_adj'][:]
        ##frcchk_wht = forccropcheck.variables['wht_adj'][:]
        ##nrows, ncols = frcchk_crn.shape             
        ##forccropcheck.close()
    
    
        # For plotting purposes only
        cropadjepipath = '/work/CLIMSIM/slc/data/cost_fn/force_crop_episode_avg_W126based.ncf'
        cropadjepifile = Dataset(cropadjepipath, 'w', format='NETCDF3_64BIT')

        # Make dimensions for file
        cropadjepifile.createDimension('rows',lrow)
        cropadjepifile.createDimension('cols',lcol)

        # Create variables and write mean of episode to file
        forc_crop_crn = cropadjepifile.createVariable('crop_adj_epi_avg_crn','f8',('rows','cols'))
        forc_crop_crn[:] = np.mean(adj_crn,axis=0)

        forc_crop_ctn = cropadjepifile.createVariable('crop_adj_epi_avg_ctn','f8',('rows','cols'))
        forc_crop_ctn[:] = np.mean(adj_ctn,axis=0)

        forc_crop_pto = cropadjepifile.createVariable('crop_adj_epi_avg_pto','f8',('rows','cols'))
        forc_crop_pto[:] = np.mean(adj_pto,axis=0)

        forc_crop_soy = cropadjepifile.createVariable('crop_adj_epi_avg_soy','f8',('rows','cols'))
        forc_crop_soy[:] = np.mean(adj_soy,axis=0)

        forc_crop_wht = cropadjepifile.createVariable('crop_adj_epi_avg_wht','f8',('rows','cols'))
        forc_crop_wht[:] = np.mean(adj_wht,axis=0)

        varattrs=["long_name","units","var_desc"]
        for varattr in varattrs:
            # Read the attribute
            if hasattr(aconcfile.variables['O3'], varattr): 
                varattrVal = getattr(aconcfile.variables['O3'], varattr);
                # Write it to the forcing file
                setattr(forc_crop_crn, varattr, varattrVal)
                setattr(forc_crop_ctn, varattr, varattrVal)
                setattr(forc_crop_pto, varattr, varattrVal)
                setattr(forc_crop_soy, varattr, varattrVal)
                setattr(forc_crop_wht, varattr, varattrVal)

        cropadjepifile.close()


    def timb_ryl_adj (seas):

        ### Adjoint forcing

        # Load file that contains crop yield data
        #    File produced from BELD land use data & NASS crop production data
        #       by Jesse Bash (US EPA) (May 2014) 
        ### > Specify path
        timb_perarea_filepath = '/work/CLIMSIM/slc/data/cost_fn/beld4_total_species_biomass_Blackard_et_al_2008_12US2.ncf'
        timb_perareafile = Dataset(timb_perarea_filepath,'r')
        timb_scalefac_filepath = '/work/CLIMSIM/slc/data/summer2007/met/GRIDCRO2D_070610.w'
        timb_scalefacfile = Dataset(timb_scalefac_filepath,'r')
    
        # Read in timb per area data
        timb_perarea_rmp = np.squeeze(timb_perareafile.variables['Acer_rubrum'])        # red maple
        timb_perarea_smp = np.squeeze(timb_perareafile.variables['Acer_saccharum'])     # sugar maple
        timb_perarea_ral = np.squeeze(timb_perareafile.variables['Alnus_rubra'])        # red alder
        timb_perarea_tlp = np.squeeze(timb_perareafile.variables['Lir_tulipifera'])     # tulip poplar
        timb_perarea_pnp = np.squeeze(timb_perareafile.variables['Pinus_ponderosa'])    # ponderosa pine
        timb_perarea_ewp = np.squeeze(timb_perareafile.variables['Pinus_strobus'])      # eastern white pine
        timb_perarea_vap = np.squeeze(timb_perareafile.variables['Pinus_virginiana'])   # virginia pine
        timb_perarea_dfr = np.squeeze(timb_perareafile.variables['Pse_menziesii'])      # douglas fir
        timb_perarea_ect = np.squeeze(timb_perareafile.variables['Pop_deltoides'])      # eastern cottonwood
        timb_perarea_qas = np.squeeze(timb_perareafile.variables['Pop_tremuloides'])    # quaking aspen
        timb_perarea_bch = np.squeeze(timb_perareafile.variables['Prunus_serotina'])    # black cherry
        nrows, ncols = timb_perarea_rmp.shape 
        wholegrid = np.ones([nrows, ncols])
        grid_zeros = np.zeros([nrows,ncols])
    
        # Read in the scaling factor for the map projection
        mapscaling = np.squeeze(timb_scalefacfile.variables['MSFX2'])          # (m/m)^2 map scaling factor from MCIP

        # Calculate the biomass in grid cell
        # Multiply timber biomass per cell (tons/ha) by 12*12 km^2 * 100 ha / km^2 * scaling factor (m^2 / m^2)
        #   scaling factor is variable MSFX2 in GRIDCRO2D file
        timb_biompercell_rmp = timb_perarea_rmp*144*(mapscaling*100)                  # red maple
        timb_biompercell_smp = timb_perarea_smp*144*(mapscaling*100)                  # sugar maple
        timb_biompercell_ral = timb_perarea_ral*144*(mapscaling*100)                  # red alder
        timb_biompercell_tlp = timb_perarea_tlp*144*(mapscaling*100)                  # tulip poplar   <--- email J.Bash to find out equivalent
        timb_biompercell_pnp = timb_perarea_pnp*144*(mapscaling*100)                  # ponderosa pine
        timb_biompercell_ewp = timb_perarea_ewp*144*(mapscaling*100)                  # eastern white pine
        timb_biompercell_vap = timb_perarea_vap*144*(mapscaling*100)                  # virginia pine
        timb_biompercell_dfr = timb_perarea_dfr*144*(mapscaling*100)                  # douglas fir
        timb_biompercell_ect = timb_perarea_ect*144*(mapscaling*100)                  # eastern cottonwood - Populous is genus for this species (J.Herrick)
        timb_biompercell_qas = timb_perarea_qas*144*(mapscaling*100)                  # quaking aspen - Populous is genus for this species (J.Herrick)
        timb_biompercell_bch = timb_perarea_bch*144*(mapscaling*100)                  # black cherry - Prunus is genus for this species (J.Herrick)

        np.isnan(np.min(timb_biompercell_bch))

        
        # Calculate hours in ozone season (highest W126 values)
        #    W126 metric is based on 90-day accumulation of ozone concentrations 

        # Read O3 values and concatenate into one array
        seas = [5,6,7,8]
        grid_zeros = np.zeros([nrows,ncols])
        O3seas = seasO3(seas)

        # Shift O3 to proper CMAQ time and
        #    zero out hours other than 8 am - 8 pm 
        O3seasshift = w126prep ( O3seas )
    
        # Assign O3 values showq
        ltime,lrow,lcol = (O3seasshift.shape)

        # Adjust the length of time over which the cost function will be applied as writing the forcing file
        forc_time = [6,7,8]
        forcdays = dys_in_mns( forc_time )
        forchrs = forcdays * 24

        neglect_time = [5]
        neglectdays = dys_in_mns( neglect_time )
        neglecthrs = neglectdays * 24

        tothrs = neglecthrs + forchrs

        # Set parameters for making dimensions
        llay = 1
        ldattim = 2                

        # Reading I/O API attributes from this file
        aconcpath = '/work/ROMO/2007platform/CMAQv4.7.1/2007ee_ORDBC_v5_07c_v4.7.1_N5ao/12US2/extr/O3/2007ee_ORDBC_v5_07c_v4.7.1_N5ao.12US2_24.combine.O3.05'
        aconcfile = Dataset(aconcpath, mode='r', open=True)

        # I/O API attributes (originally coded by M.Russell / Carleton.U)
        attrs=["IOAPI_VERSION", "EXEC_ID", "FTYPE", "CDATE", "CTIME", "WDATE", "WTIME", "SDATE", "STIME", "TSTEP", "NTHIK", "NCOLS", "NROWS", "NLAYS", "NVARS", "GDTYP", "P_ALP", "P_BET", "P_GAM", "XCENT", "YCENT", "XORIG", "YORIG", "XCELL", "YCELL", "VGTYP", "VGTOP", "VGLVLS", "GDNAM", "UPNAM", "VAR-LIST", "FILEDESC", "HISTORY"]

        tflagseasall = np.squeeze(aconcfile.variables['TFLAG'][:][:])

        # Produce the tflag variable
        for mn in range(5,11):
            tflagseas = mn_tflag(mn)
            if mn == seas[0]:
                tflagall = tflagseas
            else:
                tflagall = np.concatenate((tflagseasall,tflagseas), axis = 0)
            tflagseasall = tflagall

        tflagseasall = tflagseasall[0:tothrs,:]

        # Algorithm for relative yield forcing each grid cell at each hour between between 8:00-20:00 LST
        #
        #    d(YL)        d(W126)    d(RYL)   
        #  ---------- = ---------- x ------- 
        #  d(O3 conc)   d(O3 conc)   d(W126) 
        #
        # Enter equation below in WolframAlpha where r, q, A, & B are defined as in code / comments
        #      based on the assumption that x, r, q, A, & B are positive.
        #
        #   dJ/dC_i = d/dx (1-exp(-((x/(1+r exp(-q x)))/A)^B)) =
        #
        #        = Y*-(B A^(-B) x^(B-1) (q r x+e^(q x)+r) (r e^(-q x)+1)^(-B) e^(-A^(-B) x^B (r e^(-q x)+1)^(-B)))/(e^(q x)+r)
        #   
        #       x -> O3 conc (ppm) between 8:00-20:00 LST
        #       r -> W126 weighting
        #       q -> "
        #       A -> Timber specific for relative yield calculation (see table above)
        #       B -> "
        #       Y -> Timber yield specific to each grid cell and type
        
        r = 4403.0  # Constants for W126 weighting
        q = 126.0   # where ozone concentration is in ppm
        
        # Parameterizations of timber relative yield loss from Lehrer et al., EPA 452/R-07-002 (2007)
        # index:
        #    0 - Acer_rurbrum - Red maple - rmp
        #    1 - Acer_saccharum - Sugar maple - smp
        #    2 - Alnus_rubra - Red alder - ral 
        #    3 - Lir_tulipifera - Tulip poplar - tlp
        #    4 - Pinus_ponderosa - Ponderosa pine - pnp
        #    5 - Pinus_strobus - Eastern white pine - ewp
        #    6 - Pinus_virginiana - Virginia pine - vap 
        #    7 - Pse_menziesii - Douglas fir - dfr
        #    8 - Pop_deltoides - Eastern cottonwood - ect
        #    9 - Pop_tremuloides - Quaking aspen - qas
        #    10 -Prunus_serotina - Black cherry - bch
        Atm = np.array([318.12, 36.35, 179.06, 51.38, 159.63, 63.23, 1714.64, 106.83, 10.10, 109.81, 38.92])
        Btm = np.array([1.3756, 5.7785, 1.2377, 2.0889, 1.190, 1.6582, 1.00, 5.9631, 1.7793, 1.2198, 0.9921])
            
        w126seas = np.zeros([ltime, lrow, lcol])
        w126seas = (O3seasshift / (1.0+r*np.exp(-q*O3seasshift)))
        w126seas_cum = np.sum(w126seas[neglecthrs:tothrs], axis = 0)

        # Set parameter for number of tree species being treated. 
        timbspecies = 11
        forc_timb = np.zeros([timbspecies, tothrs, lrow, lcol])

        O3seasshift_forc = O3seasshift[neglecthrs:tothrs,:,:]    
            
        # Calculate the adjoint forcing for all of the months
        delW126_delO3 = ((np.exp(q*O3seasshift_forc) * (r*q*O3seasshift_forc + r + np.exp(q*O3seasshift_forc))) / (r + np.exp(q*O3seasshift_forc))**2 )
        delW126_delO3 = np.where(O3seasshift_forc > 0.0, delW126_delO3, O3seasshift_forc)
        np.isnan(np.min(delW126_delO3))

        for ttyp in range(timbspecies):
            forc_timb[ttyp][neglecthrs:tothrs] = (Btm[ttyp] * Atm[ttyp]**(-Btm[ttyp]) * w126seas_cum**(Btm[ttyp]-1.0) * np.exp(-(np.divide(Atm[ttyp],w126seas_cum)**(-Btm[ttyp])))) * delW126_delO3
            
        np.isnan(np.min(forc_timb))


        # Create files
        forcfilepath_rmp = '/work/CLIMSIM/slc/data/cost_fn/force_timb_07_12US2_yield_wholeW126_rmp.ncf'
        forctimbdist_rmp = Dataset(forcfilepath_rmp, 'w', format='NETCDF3_64BIT')
        forcfilepath_smp = '/work/CLIMSIM/slc/data/cost_fn/force_timb_07_12US2_yield_wholeW126_smp.ncf'
        forctimbdist_smp = Dataset(forcfilepath_smp, 'w', format='NETCDF3_64BIT')
        forcfilepath_ral = '/work/CLIMSIM/slc/data/cost_fn/force_timb_07_12US2_yield_wholeW126_ral.ncf'
        forctimbdist_ral = Dataset(forcfilepath_ral, 'w', format='NETCDF3_64BIT')
        forcfilepath_tlp = '/work/CLIMSIM/slc/data/cost_fn/force_timb_07_12US2_yield_wholeW126_tlp.ncf'
        forctimbdist_tlp = Dataset(forcfilepath_tlp, 'w', format='NETCDF3_64BIT')
        forcfilepath_pnp = '/work/CLIMSIM/slc/data/cost_fn/force_timb_07_12US2_yield_wholeW126_pnp.ncf'
        forctimbdist_pnp = Dataset(forcfilepath_pnp, 'w', format='NETCDF3_64BIT')
        forcfilepath_ewp = '/work/CLIMSIM/slc/data/cost_fn/force_timb_07_12US2_yield_wholeW126_ewp.ncf'
        forctimbdist_ewp = Dataset(forcfilepath_ewp, 'w', format='NETCDF3_64BIT')
        forcfilepath_vap = '/work/CLIMSIM/slc/data/cost_fn/force_timb_07_12US2_yield_wholeW126_vap.ncf'
        forctimbdist_vap = Dataset(forcfilepath_vap, 'w', format='NETCDF3_64BIT')
        forcfilepath_dfr = '/work/CLIMSIM/slc/data/cost_fn/force_timb_07_12US2_yield_wholeW126_dfr.ncf'
        forctimbdist_dfr = Dataset(forcfilepath_dfr, 'w', format='NETCDF3_64BIT')
        forcfilepath_ect = '/work/CLIMSIM/slc/data/cost_fn/force_timb_07_12US2_yield_wholeW126_ect.ncf'
        forctimbdist_ect = Dataset(forcfilepath_ect, 'w', format='NETCDF3_64BIT')
        forcfilepath_qas = '/work/CLIMSIM/slc/data/cost_fn/force_timb_07_12US2_yield_wholeW126_qas.ncf'
        forctimbdist_qas = Dataset(forcfilepath_qas, 'w', format='NETCDF3_64BIT')
        forcfilepath_bch = '/work/CLIMSIM/slc/data/cost_fn/force_timb_07_12US2_yield_wholeW126_bch.ncf'
        forctimbdist_bch = Dataset(forcfilepath_bch, 'w', format='NETCDF3_64BIT')
        
        for attr in attrs:
            # Read the attribute
            if hasattr(aconcfile, attr): 
                attrVal = getattr(aconcfile, attr);
                # Write it to the forcing file
                setattr(forctimbdist_rmp, attr, attrVal)
                setattr(forctimbdist_smp, attr, attrVal)
                setattr(forctimbdist_ral, attr, attrVal)
                setattr(forctimbdist_tlp, attr, attrVal)
                setattr(forctimbdist_pnp, attr, attrVal)
                setattr(forctimbdist_ewp, attr, attrVal)
                setattr(forctimbdist_vap, attr, attrVal)
                setattr(forctimbdist_dfr, attr, attrVal)
                setattr(forctimbdist_ect, attr, attrVal)
                setattr(forctimbdist_qas, attr, attrVal)
                setattr(forctimbdist_bch, attr, attrVal)
        
        # I/O API dimnensions
        forctimbdist_rmp.createDimension("TSTEP", None)
        forctimbdist_smp.createDimension("TSTEP", None)
        forctimbdist_ral.createDimension("TSTEP", None)
        forctimbdist_tlp.createDimension("TSTEP", None)
        forctimbdist_pnp.createDimension("TSTEP", None)
        forctimbdist_ewp.createDimension("TSTEP", None)
        forctimbdist_vap.createDimension("TSTEP", None)
        forctimbdist_dfr.createDimension("TSTEP", None)
        forctimbdist_ect.createDimension("TSTEP", None)
        forctimbdist_qas.createDimension("TSTEP", None)
        forctimbdist_bch.createDimension("TSTEP", None)

        forctimbdist_rmp.createDimension("DATE-TIME", ldattim)
        forctimbdist_smp.createDimension("DATE-TIME", ldattim)
        forctimbdist_ral.createDimension("DATE-TIME", ldattim)
        forctimbdist_tlp.createDimension("DATE-TIME", ldattim)
        forctimbdist_pnp.createDimension("DATE-TIME", ldattim)
        forctimbdist_ewp.createDimension("DATE-TIME", ldattim)
        forctimbdist_vap.createDimension("DATE-TIME", ldattim)
        forctimbdist_dfr.createDimension("DATE-TIME", ldattim)
        forctimbdist_ect.createDimension("DATE-TIME", ldattim)
        forctimbdist_qas.createDimension("DATE-TIME", ldattim)
        forctimbdist_bch.createDimension("DATE-TIME", ldattim)

        forctimbdist_rmp.createDimension("LAY", 1)
        forctimbdist_smp.createDimension("LAY", 1)
        forctimbdist_ral.createDimension("LAY", 1)
        forctimbdist_tlp.createDimension("LAY", 1)
        forctimbdist_pnp.createDimension("LAY", 1)
        forctimbdist_ewp.createDimension("LAY", 1)
        forctimbdist_vap.createDimension("LAY", 1)
        forctimbdist_dfr.createDimension("LAY", 1)
        forctimbdist_ect.createDimension("LAY", 1)
        forctimbdist_qas.createDimension("LAY", 1)
        forctimbdist_bch.createDimension("LAY", 1)

        forctimbdist_rmp.createDimension("VAR", 1)
        forctimbdist_smp.createDimension("VAR", 1)
        forctimbdist_ral.createDimension("VAR", 1)
        forctimbdist_tlp.createDimension("VAR", 1)
        forctimbdist_pnp.createDimension("VAR", 1)
        forctimbdist_ewp.createDimension("VAR", 1)
        forctimbdist_vap.createDimension("VAR", 1)
        forctimbdist_dfr.createDimension("VAR", 1)
        forctimbdist_ect.createDimension("VAR", 1)
        forctimbdist_qas.createDimension("VAR", 1)
        forctimbdist_bch.createDimension("VAR", 1)

        forctimbdist_rmp.createDimension("ROW", lrow)
        forctimbdist_smp.createDimension("ROW", lrow)
        forctimbdist_ral.createDimension("ROW", lrow)
        forctimbdist_tlp.createDimension("ROW", lrow)
        forctimbdist_pnp.createDimension("ROW", lrow)
        forctimbdist_ewp.createDimension("ROW", lrow)
        forctimbdist_vap.createDimension("ROW", lrow)
        forctimbdist_dfr.createDimension("ROW", lrow)
        forctimbdist_ect.createDimension("ROW", lrow)
        forctimbdist_qas.createDimension("ROW", lrow)
        forctimbdist_bch.createDimension("ROW", lrow)

        forctimbdist_rmp.createDimension("COL", lcol)
        forctimbdist_smp.createDimension("COL", lcol)
        forctimbdist_ral.createDimension("COL", lcol)
        forctimbdist_tlp.createDimension("COL", lcol)
        forctimbdist_pnp.createDimension("COL", lcol)
        forctimbdist_ewp.createDimension("COL", lcol)
        forctimbdist_vap.createDimension("COL", lcol)
        forctimbdist_dfr.createDimension("COL", lcol)
        forctimbdist_ect.createDimension("COL", lcol)
        forctimbdist_qas.createDimension("COL", lcol)
        forctimbdist_bch.createDimension("COL", lcol)

        forctimbdist_rmp.sync()
        forctimbdist_smp.sync()
        forctimbdist_ral.sync()
        forctimbdist_tlp.sync()
        forctimbdist_pnp.sync()
        forctimbdist_ewp.sync()
        forctimbdist_vap.sync()
        forctimbdist_dfr.sync()
        forctimbdist_ect.sync()
        forctimbdist_qas.sync()
        forctimbdist_bch.sync()

        # Create TFLAG variable
        forc_timb_rmp_tflag = forctimbdist_rmp.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        forc_timb_smp_tflag = forctimbdist_smp.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        forc_timb_ral_tflag = forctimbdist_ral.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        forc_timb_tlp_tflag = forctimbdist_tlp.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        forc_timb_pnp_tflag = forctimbdist_pnp.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        forc_timb_ewp_tflag = forctimbdist_ewp.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        forc_timb_vap_tflag = forctimbdist_vap.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        forc_timb_dfr_tflag = forctimbdist_dfr.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        forc_timb_ect_tflag = forctimbdist_ect.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        forc_timb_qas_tflag = forctimbdist_qas.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        forc_timb_bch_tflag = forctimbdist_bch.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))

        varattrs=["long_name","units","var_desc"]
        for varattr in varattrs:
            # Read the attribute
            if hasattr(aconcfile.variables['TFLAG'], varattr): 
                varattrVal = getattr(aconcfile.variables['TFLAG'], varattr);
                # Write it to the forcing file
                setattr(forc_timb_rmp_tflag, varattr, varattrVal)
                setattr(forc_timb_smp_tflag, varattr, varattrVal)
                setattr(forc_timb_ral_tflag, varattr, varattrVal)
                setattr(forc_timb_tlp_tflag, varattr, varattrVal)
                setattr(forc_timb_pnp_tflag, varattr, varattrVal)
                setattr(forc_timb_ewp_tflag, varattr, varattrVal)
                setattr(forc_timb_vap_tflag, varattr, varattrVal)
                setattr(forc_timb_dfr_tflag, varattr, varattrVal)
                setattr(forc_timb_ect_tflag, varattr, varattrVal)
                setattr(forc_timb_qas_tflag, varattr, varattrVal)
                setattr(forc_timb_bch_tflag, varattr, varattrVal)

        #forc_timb_rmp_tflag[:] = np.zeros([tothrs,llay,ldattim])
        #forc_timb_smp_tflag[:] = np.zeros([tothrs,llay,ldattim])
        #forc_timb_ral_tflag[:] = np.zeros([tothrs,llay,ldattim])
        #forc_timb_tlp_tflag[:] = np.zeros([tothrs,llay,ldattim])
        #forc_timb_pnp_tflag[:] = np.zeros([tothrs,llay,ldattim])
        #forc_timb_ewp_tflag[:] = np.zeros([tothrs,llay,ldattim])
        #forc_timb_vap_tflag[:] = np.zeros([tothrs,llay,ldattim])
        #forc_timb_dfr_tflag[:] = np.zeros([tothrs,llay,ldattim])
        #forc_timb_ect_tflag[:] = np.zeros([tothrs,llay,ldattim])
        #forc_timb_qas_tflag[:] = np.zeros([tothrs,llay,ldattim])
        #forc_timb_bch_tflag[:] = np.zeros([tothrs,llay,ldattim])

        forc_timb_rmp_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])
        forc_timb_smp_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])
        forc_timb_ral_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])
        forc_timb_tlp_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])
        forc_timb_pnp_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])
        forc_timb_ewp_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])
        forc_timb_vap_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])
        forc_timb_dfr_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])
        forc_timb_ect_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])
        forc_timb_qas_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])
        forc_timb_bch_tflag[:] = np.reshape(tflagseasall,[tothrs,llay,ldattim])

        forctimbdist_rmp.sync()
        forctimbdist_smp.sync()
        forctimbdist_ral.sync()
        forctimbdist_tlp.sync()
        forctimbdist_pnp.sync()
        forctimbdist_ewp.sync()
        forctimbdist_vap.sync()
        forctimbdist_dfr.sync()
        forctimbdist_ect.sync()
        forctimbdist_qas.sync()
        forctimbdist_bch.sync()
 
        
        # Create variables
        forc_rmp = forctimbdist_rmp.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        # Scale by presence
        adj_rmp = forc_timb[0]*timb_biompercell_rmp
        # Reshape to read in correctly as adjoint forcing
        forc_rmp[:] = np.reshape(adj_rmp,[2952,1,246,396])

        forc_smp = forctimbdist_smp.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        adj_smp = forc_timb[1]*timb_biompercell_smp
        forc_smp[:] = np.reshape(adj_smp,[2952,1,246,396])

        forc_ral = forctimbdist_ral.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        adj_ral =  forc_timb[2]*timb_biompercell_ral
        forc_ral[:] = np.reshape(adj_ral,[2952,1,246,396])

        forc_tlp = forctimbdist_tlp.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        adj_tlp = forc_timb[3]*timb_biompercell_tlp
        forc_tlp[:] = np.reshape(adj_tlp,[2952,1,246,396])

        forc_pnp = forctimbdist_pnp.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        adj_pnp = forc_timb[4]*timb_biompercell_pnp
        forc_pnp[:] = np.reshape(adj_pnp,[2952,1,246,396])

        forc_ewp = forctimbdist_ewp.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        adj_ewp = forc_timb[5]*timb_biompercell_ewp
        forc_ewp[:] = np.reshape(adj_ewp,[2952,1,246,396])

        forc_vap = forctimbdist_vap.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        adj_vap =  forc_timb[6]*timb_biompercell_vap
        forc_vap[:] = np.reshape(adj_vap,[2952,1,246,396])

        forc_dfr = forctimbdist_dfr.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        adj_dfr = forc_timb[7]*timb_biompercell_dfr
        forc_dfr[:] = np.reshape(adj_dfr,[2952,1,246,396])

        forc_ect = forctimbdist_ect.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        adj_ect = forc_timb[8]*timb_biompercell_ect
        forc_ect[:] = np.reshape(adj_ect,[2952,1,246,396])

        forc_qas = forctimbdist_qas.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        adj_qas = forc_timb[9]*timb_biompercell_qas
        forc_qas[:] = np.reshape(adj_qas,[2952,1,246,396])

        forc_bch = forctimbdist_bch.createVariable('O3','f8',('TSTEP','LAY','ROW','COL'))
        adj_bch = forc_timb[10]*timb_biompercell_bch
        forc_bch[:] = np.reshape(adj_bch,[2952,1,246,396])

        varattrs=["long_name","units","var_desc"]
        for varattr in varattrs:
            # Read the attribute
            if hasattr(aconcfile.variables['O3'], varattr): 
                varattrVal = getattr(aconcfile.variables['O3'], varattr);
                # Write it to the forcing file
                setattr(forc_rmp, varattr, varattrVal)
                setattr(forc_smp, varattr, varattrVal)
                setattr(forc_ral, varattr, varattrVal)
                setattr(forc_tlp, varattr, varattrVal)
                setattr(forc_pnp, varattr, varattrVal)
                setattr(forc_ewp, varattr, varattrVal)
                setattr(forc_vap, varattr, varattrVal)
                setattr(forc_dfr, varattr, varattrVal)
                setattr(forc_ect, varattr, varattrVal)
                setattr(forc_qas, varattr, varattrVal)
                setattr(forc_bch, varattr, varattrVal)

        # For plotting purposes only
        
        # Read file
#        forcfilepath = '/work/CLIMSIM/slc/data/cost_fn/force_timb_07_12US2.ncf'
#        forctimbdist = Dataset(forcfilepath, 'r', format='NETCDF3_64BIT')
        
        forctimbdist_rmp.close()
        forctimbdist_smp.close()
        forctimbdist_ral.close()
        forctimbdist_tlp.close()
        forctimbdist_pnp.close()
        forctimbdist_ewp.close()
        forctimbdist_vap.close()
        forctimbdist_dfr.close()
        forctimbdist_ect.close()
        forctimbdist_qas.close()
        forctimbdist_bch.close()
        


        # Create files for visualizing and evaluating the data
        adj_rmp = forctimbdist_rmp.variables['O3'][:]
        adj_smp = forctimbdist_smp.variables['O3'][:]
        adj_ral = forctimbdist_ral.variables['O3'][:]
        adj_tlp = forctimbdist_tlp.variables['O3'][:]
        adj_pnp = forctimbdist_pnp.variables['O3'][:]
        adj_ewp = forctimbdist_ewp.variables['O3'][:]
        adj_vap = forctimbdist_vap.variables['O3'][:]
        adj_dfr = forctimbdist_dfr.variables['O3'][:]
        adj_ect = forctimbdist_ect.variables['O3'][:]
        adj_qas = forctimbdist_qas.variables['O3'][:]
        adj_bch = forctimbdist_bch.variables['O3'][:]
        
        timbadjepipath = '/work/CLIMSIM/slc/data/cost_fn/force_timb_episode_avg.ncf'
        timbadjepifile = Dataset(timbadjepipath, 'w', format='NETCDF3_64BIT')

        # Make dimensions for file
        lrow = 246
        lcol = 396
        timbadjepifile.createDimension('rows',lrow)
        timbadjepifile.createDimension('cols',lcol)

        # Create variables and write mean of episode to file
        forc_timb_rmp = timbadjepifile.createVariable('timb_adj_epi_avg_rmp','f8',('rows','cols'))
        forc_timb_rmp[:] = np.mean(adj_rmp,axis=0)

        forc_timb_smp = timbadjepifile.createVariable('timb_adj_epi_avg_smp','f8',('rows','cols'))
        forc_timb_smp[:] = np.mean(adj_smp,axis=0)

        forc_timb_ral = timbadjepifile.createVariable('timb_adj_epi_avg_ral','f8',('rows','cols'))
        forc_timb_ral[:] = np.mean(adj_ral,axis=0)

        forc_timb_tlp = timbadjepifile.createVariable('timb_adj_epi_avg_tlp','f8',('rows','cols'))
        forc_timb_tlp[:] = np.mean(adj_tlp,axis=0)

        forc_timb_pnp = timbadjepifile.createVariable('timb_adj_epi_avg_pnp','f8',('rows','cols'))
        forc_timb_pnp[:] = np.mean(adj_pnp,axis=0)

        forc_timb_ewp = timbadjepifile.createVariable('timb_adj_epi_avg_ewp','f8',('rows','cols'))
        forc_timb_ewp[:] = np.mean(adj_ewp,axis=0)

        forc_timb_vap = timbadjepifile.createVariable('timb_adj_epi_avg_vap','f8',('rows','cols'))
        forc_timb_vap[:] = np.mean(adj_vap,axis=0)

        forc_timb_dfr = timbadjepifile.createVariable('timb_adj_epi_avg_dfr','f8',('rows','cols'))
        forc_timb_dfr[:] = np.mean(adj_dfr,axis=0)

        forc_timb_ect = timbadjepifile.createVariable('timb_adj_epi_avg_ect','f8',('rows','cols'))
        forc_timb_ect[:] = np.mean(adj_ect,axis=0)

        forc_timb_qas = timbadjepifile.createVariable('timb_adj_epi_avg_qas','f8',('rows','cols'))
        forc_timb_qas[:] = np.mean(adj_qas,axis=0)

        forc_timb_bch = timbadjepifile.createVariable('timb_adj_epi_avg_bch','f8',('rows','cols'))
        forc_timb_bch[:] = np.mean(adj_bch,axis=0)

        timbadjepifile.close()

        # Create files
        ylfilepath = '/work/CLIMSIM/slc/data/cost_fn/yieldloss_timb_07_12US2.ncf'
        yltimbdist = Dataset(ylfilepath, 'w', format='NETCDF3_64BIT')

        # Make dimensions for file
        lrow = 246
        lcol = 396
        yltimbdist.createDimension('rows',lrow)
        yltimbdist.createDimension('cols',lcol)

        yl_timb_rmp=yltimbdist.createVariable('rmp_yl','f8',('rows','cols'))
        yl_timb_smp=yltimbdist.createVariable('smp_yl','f8',('rows','cols'))
        yl_timb_ral=yltimbdist.createVariable('ral_yl','f8',('rows','cols'))
        yl_timb_tlp=yltimbdist.createVariable('tlp_yl','f8',('rows','cols'))
        yl_timb_pnp=yltimbdist.createVariable('pnp_yl','f8',('rows','cols'))
        yl_timb_ewp=yltimbdist.createVariable('ewp_yl','f8',('rows','cols'))
        yl_timb_vap=yltimbdist.createVariable('vap_yl','f8',('rows','cols'))
        yl_timb_dfr=yltimbdist.createVariable('dfr_yl','f8',('rows','cols'))
        yl_timb_ect=yltimbdist.createVariable('ect_yl','f8',('rows','cols'))
        yl_timb_qas=yltimbdist.createVariable('qas_yl','f8',('rows','cols'))
        yl_timb_bch=yltimbdist.createVariable('bch_yl','f8',('rows','cols'))

        yl_timb_rmp[:] = timb_biompercell_rmp*(1.0-np.exp(-((w126seas_cum/Atm[0])**Btm[0])))
        yl_timb_smp[:] = timb_biompercell_smp*(1.0-np.exp(-((w126seas_cum/Atm[1])**Btm[1])))
        yl_timb_ral[:] = timb_biompercell_ral*(1.0-np.exp(-((w126seas_cum/Atm[2])**Btm[2])))
        yl_timb_tlp[:] = timb_biompercell_tlp*(1.0-np.exp(-((w126seas_cum/Atm[3])**Btm[3])))
        yl_timb_pnp[:] = timb_biompercell_pnp*(1.0-np.exp(-((w126seas_cum/Atm[4])**Btm[4])))
        yl_timb_ewp[:] = timb_biompercell_ewp*(1.0-np.exp(-((w126seas_cum/Atm[5])**Btm[5])))
        yl_timb_vap[:] = timb_biompercell_vap*(1.0-np.exp(-((w126seas_cum/Atm[6])**Btm[6])))
        yl_timb_dfr[:] = timb_biompercell_dfr*(1.0-np.exp(-((w126seas_cum/Atm[7])**Btm[7])))
        yl_timb_ect[:] = timb_biompercell_ect*(1.0-np.exp(-((w126seas_cum/Atm[8])**Btm[8])))
        yl_timb_qas[:] = timb_biompercell_qas*(1.0-np.exp(-((w126seas_cum/Atm[9])**Btm[9])))
        yl_timb_bch[:] = timb_biompercell_bch*(1.0-np.exp(-((w126seas_cum/Atm[10])**Btm[10])))

        mass_timb_rmp = yltimbdist.createVariable('rmp_mass','f8',('rows','cols'))
        mass_timb_smp = yltimbdist.createVariable('smp_mass','f8',('rows','cols'))
        mass_timb_ral = yltimbdist.createVariable('ral_mass','f8',('rows','cols'))
        mass_timb_tlp = yltimbdist.createVariable('tlp_mass','f8',('rows','cols'))
        mass_timb_pnp = yltimbdist.createVariable('pnp_mass','f8',('rows','cols'))
        mass_timb_ewp = yltimbdist.createVariable('ewp_mass','f8',('rows','cols'))
        mass_timb_vap = yltimbdist.createVariable('vap_mass','f8',('rows','cols'))
        mass_timb_dfr = yltimbdist.createVariable('dfr_mass','f8',('rows','cols'))
        mass_timb_ect = yltimbdist.createVariable('ect_mass','f8',('rows','cols'))
        mass_timb_qas = yltimbdist.createVariable('qas_mass','f8',('rows','cols'))
        mass_timb_bch = yltimbdist.createVariable('bch_mass','f8',('rows','cols'))

        mass_timb_rmp[:] = timb_biompercell_rmp
        mass_timb_smp[:] = timb_biompercell_smp
        mass_timb_ral[:] = timb_biompercell_ral
        mass_timb_tlp[:] = timb_biompercell_tlp
        mass_timb_pnp[:] = timb_biompercell_pnp
        mass_timb_ewp[:] = timb_biompercell_ewp
        mass_timb_vap[:] = timb_biompercell_vap
        mass_timb_dfr[:] = timb_biompercell_dfr
        mass_timb_ect[:] = timb_biompercell_ect
        mass_timb_qas[:] = timb_biompercell_qas
        mass_timb_bch[:] = timb_biompercell_bch

        yltimbdist.close()

class O3utility:
    
    # Concatenate months of O3 concentrations in local time together
    #    input: array of months (1-Jan, 2-Feb, etc.)
    #    output: O3 concentration in ppm in every hour at every grid cell in local time for all months in input array 
    def seasO3( seas ):
    	# Make the hour local time inside mn_O3 one month at a time
        #     then concatenate the months together.
        for mn in seas:
            print(mn)
            print(seas[0])
            O3mn = mn_o3(mn)
            print(O3mn.shape)
            if mn == seas[0]:
                O3season = O3mn
            else:
                O3season = np.concatenate((O3seas,O3mn),axis=0)
            O3seas = O3season
        
        print (O3seas.shape)
        return O3seas
    
    # Read ozone for month from file and return in array in local time
    #    input: variable representing month (1-Jan, 2-Feb, etc.)
    #    output: array of O3 at local time that is size of grid and length of time of the month in hours
    def mn_o3 ( mm ):
        # Returns the hourly ozone (ppm) in each grid cell in Local Standard Time
        yyyy = 2007
        dd30 = set([4,6,9,11])
        dd28 = set([2])
        if mm in dd28:
            dd = 28
        elif mm in dd30:
            dd = 30
        else:
            dd = 31
        if len(str(mm)) < 2:
            mStamp = '0' + str(mm)
        else:
            mStamp = str(mm)
        # For each day, open the hourly average concentration file
        aconcname = '/data/shannon/cmaqadj/2007/output/2007ee_ORDBC_v5_07c_v4.7.1_N5ao.12US2_24.combine.O3.'+mStamp
        aconcfile = Dataset(aconcname, mode='r', open=True)
        # Read in file that contains hours to shift to local time
        tzpath = '/Users/slc/Dropbox/Research/cmaqadjoint/PreProcessor/aux_files/12us2_time.dat'
        tzgrid = np.genfromtxt(tzpath, dtype=[('COL', '<i8'),('ROW', '<i8'),('CR', '<i8'),('TIMEZN', 'S8'),('ZONE','<f8')], skip_header=1)
        o3hr = np.squeeze(aconcfile.variables['O3'][:][:][:])
        ntime, nrows, ncols = o3hr.shape
        o3hr = o3hr/1000.0
        print (o3hr.shape)
        for rowind in range(0,nrows):
            for colind in range(0,ncols):
                ind = rowind*246 + colind
                tzshift = int(tzgrid[ind][4])
                for timeind in range(0,dd*24):
                    #print tzshift
                    if timeind + tzshift > 0:
                        o3hr[timeind+tzshift,rowind,colind] = o3hr[timeind,rowind,colind]
        return o3hr
    
    # Zero out hours that do not contribute to W126
    #    input: time series of hourly ozone concentrations in local standard time
    #    output: time series of hourly ozone concentrations zeroed out between 8 pm and 8 am LST and then shifted back to CMAQ time
    def w126prep ( O3series ): 
        ltime, lrow, lcol = O3series.shape
        
        # Make mask for filtering 8 pm - 8 am values
        nightfilter = np.zeros([ltime, lrow, lcol])
        for tind in range(0,ltime-1):
            for rowind in range(0,lrow):
                for colind in range(0,lcol):
                    hour = tind % 24
                    #print(hour)
                    if ((hour > 7) and (hour < 20)):
                        nightfilter[tind,rowind,colind] = 1
        
        print(sum(O3series))
        O3daytime = nightfilter * O3series
        print(sum(O3daytime))
        
        # Shift values back to CMAQ time
        
        # Read in file that contains hours to shift to local time
        tzpath = '/Users/slc/Dropbox/Research/cmaqadjoint/PreProcessor/aux_files/12us2_time.dat'
        tzgrid = np.genfromtxt(tzpath, dtype=[('COL', '<i8'),('ROW', '<i8'),('CR', '<i8'),('TIMEZN', 'S8'),('ZONE','<f8')], skip_header=1)
        
        O3w126ready = np.zeros([ltime,lrow,lcol])
        
        # Apply to time index
        for rowind in range(0,nrows):
            for colind in range(0,ncols):
                ind = rowind*246 + colind
                tzshift = int(tzgrid[ind][4])
                for timeind in range(0,ltime-9):
                    #print tzshift
                    if timeind-tzshift > 0:
                        O3w126ready[timeind-tzshift,rowind,colind] = O3daytime[timeind,rowind,colind]
        
        return O3w126ready
    
    # Determine days in month
    #    input: array of months (1-Jan, 2-Feb, etc.)
    #    output: sum of days in months input
    def dys_in_mns( mns ):
        days = 0
        for mm in mns:
            dd30 = set([4,6,9,11])
            dd28 = set([2])
            if mm in dd28:
                dd = 28
            elif mm in dd30:
                dd = 30
            else:
                dd = 31
            days += dd
        
        return days
        
    
    
    # Read binary array and shift to local time
    #    input: array in model time with 1's representing where the max 1-hr ozone occurred in the day and 0's elsewhere
    #    output: array in model time with 1's representing where the max 1-hr ozone occurred in the day and 0's elsewhere
    def readjust_time( binarray ):
        for mm in range(5,11):
            yyyy = 2007
            dd30 = set([4,6,9,11])
            dd28 = set([2])
            if mm in dd28:
                dd = 28
            elif mm in dd30:
                dd = 30
            else:
                dd = 31
            tzpath = '/Users/slc/Dropbox/Research/cmaqadjoint/PreProcessor/aux_files/12us2_time.dat'
            tzgrid = np.genfromtxt(tzpath, dtype=[('COL', '<i8'),('ROW', '<i8'),('CR', '<i8'),('TIMEZN', 'S8'),('ZONE','<f8')], skip_header=1)
            ntime, nrows, ncols = binarray.shape
            print (binarray.shape)
            for rowind in range(0,nrows):
                for colind in range(0,ncols):
                    ind = rowind*246 + colind
                    tzshift = int(tzgrid[ind][4])
                    for timeind in range(0,dd*24):
                        #print tzshift
                        if timeind - tzshift > 0:
                            binarray[timeind-tzshift,rowind,colind] = binarray[timeind,rowind,colind]
        return binarray
    
    # Read TFLAG variable from ACONC file for month indicated by input variable
    #    input: variable representing month (1-Jan, 2-Feb, etc.)
    #    output: the TFLAG variable for each month in dattim
    def mn_tflag ( mm ):
        ### Returns 
        # For each day, open the hourly average concentration file
        yyyy = 2007
        dd30 = set([4,6,9,11])
        dd28 = set([2])
        if mm in dd28:
            dd = 28
        elif mm in dd30:
            dd = 30
        else:
            dd = 31
        if len(str(mm)) < 2:
            mStamp = '0' + str(mm)
        else:
            mStamp = str(mm)
        aconcname = '/data/shannon/cmaqadj/2007/output/2007ee_ORDBC_v5_07c_v4.7.1_N5ao.12US2_24.combine.O3.'+mStamp
        aconcfile = Dataset(aconcname, mode='r', open=True)
        dattim = np.squeeze(aconcfile.variables['TFLAG'][:][:])
        ntime, ndattim = dattim.shape
        return dattim
    
