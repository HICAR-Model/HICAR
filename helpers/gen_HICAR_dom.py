import xarray as xr
import os
import sys
sys.path.append(os.getcwd())
import HICAR_Domain as hd

###########################################################################
#######################  User options, to be edited  ######################
###########################################################################

#The resolution of the domain
res = 250

# The target domain, including lat and lon variables named as "lat" and "lon", and
# a DEM labeled as "topo".
target_domain_fn = 'Target_domain.nc'
# Same as above, but for a domain with extent ~20km beyond the borders of the target
# domain.
large_domain_fn = 'Large_domain.nc'
# Name of output file
output_domain_fn = 'output_domain.nc'

# These are used in the calculation of ridelines, and can be tuned if the user
# is not satisfied with the deliniation of ridgelines in the output file
terr_filter = 10
TPI_thresh = 100

###########################################################################
############################  End user options  ###########################
###########################################################################

dom = xr.open_dataset(target_domain_fn)
dom_rad = xr.open_dataset(large_domain_fn)

dom_out = hd.wholeShebang(dom,dom_rad,res=res,terr_filter=terr_filter,TPI_thresh=TPI_thresh)

dom_out.to_netcdf(output_domain_fn)