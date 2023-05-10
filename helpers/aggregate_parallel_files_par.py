#!/usr/bin/env python

import glob

import multiprocessing as mp
from multiprocessing import shared_memory
import argparse
import numpy as np
import xarray as xr
import os
import time

pool = None
data_vars = None


def load_file(file_name,var_mems):
    '''Load a netcdf dataset into memory'''
    d = xr.open_dataset(file_name)
    subset(d,var_mems)
    return 'loaded '+file_name


def get_dims(dataset, section="d"):
    '''Get the global attributes defining the domain, memory, or tile space'''
    results = []
    for axis in ["i","j","k"]:
        for position in ["s","e"]:
            results.append(int(dataset.attrs[axis + section + position]))
    return results

def get_dim_offset(dims):
    '''Return x_offset, y_offset
    For the staggered dims, offset=1, otherwise offset=0'''
    x_off = 0
    if 'lon_u' in dims: x_off = 1

    y_off = 0
    if 'lat_v' in dims: y_off = 1

    return x_off, y_off

def set_up_data_vars(d,var_names):
    ids, ide, jds, jde, kds, kde = get_dims(d, section='d')
    nx = ide - ids + 1
    ny = jde - jds + 1
    nz = kde - kds + 1

    data_vars = dict()
    shms = dict()

    if var_names is None: var_names = d.variables

    for v in var_names:
        coords = [c for c in d[v].coords]
        dims   = d[v].dims
        name   = d[v].name
        attrs  = d[v].attrs

        x_off, y_off = get_dim_offset(dims)

        if len(dims) == 1:
            nt = d.dims[dims[0]]
            data = np.zeros((nt))
        if len(dims) == 2:
            data = np.zeros((ny + y_off, nx + x_off))
        if len(dims) == 3:
            data = np.zeros((d.dims[dims[0]], ny + y_off, nx + x_off))
        if len(dims) == 4:
            nt = d.dims[dims[0]]
            nz = d.dims[dims[1]]
            data = np.zeros((nt, nz, ny + y_off, nx + x_off))
        
        #make data into shared memory
        shm = mp.shared_memory.SharedMemory(create=True, size=data.nbytes)
        shms[v] = shm
        shm_data = np.ndarray(data.shape, dtype=data.dtype, buffer=shm.buf)
        # print(name, data.shape, dims, attrs)
        data_vars[v] = xr.DataArray(shm_data, dims=dims, name=name, attrs=attrs)#, coords=coords)
 
    return data_vars, shms

def set_up_dataset(d,data_vars):
    '''Create a dataset to cover the entire domain with the variables present in d

    d : an input dataset covering part of the domain
    d must have global attributes ids, ide, jds, jde, kds, kde that define the full domain

    A new dataset is created with all the variables+attributes in d covering the full domain
    '''
    ds = xr.Dataset(data_vars, attrs=d.attrs)
    ds.encoding = d.encoding
    ds["time"] = d["time"]
    cords = []
    
    for v in d.variables:
        if v in data_vars:
            for c in d[v].coords:
                if not(c in cords) and (c in data_vars):
                    cords.append(c)
    return ds.set_coords(cords)


def agg_file(first_file,var_names,verbose=True):
    '''Aggregated all files that come from the same time step as first_file

    first_file should have _001_ in the filename somewhere.  This will be replaced
    with * to search for all matching files from this date. Once files are found, a
    dataset containing the entire domain is created and the data from each file are
    added to the master dataset.

    Result: aggregated dataset is written to a netcdf file'''
    
    if verbose:print(first_file)
    date_search = first_file.replace("_000001_","*")
    outputfile = first_file.replace("000001_","_").replace("__","_")
    if os.path.isfile(outputfile):
        return
    this_date_files = glob.glob(date_search)
    this_date_files.sort()
    
    template = xr.open_dataset(this_date_files[0])
    data_vars, shms = set_up_data_vars(template,var_names)

    args = [(d,shms) for d in this_date_files]
    results = pool.starmap_async(load_file, args)

    #Just use get to wait for result
    message = results.get()

    data_set = set_up_dataset(template,data_vars)
    data_set.load().to_netcdf(outputfile)
    for key in shms:
        shms[key].close()
        shms[key].unlink()

def subset(d,var_mems):
        ids, ide, jds, jde, kds, kde = get_dims(d, section='d')
        ims, ime, jms, jme, kms, kme = get_dims(d, section='m')
        its, ite, jts, jte, kts, kte = get_dims(d, section='t')

        xts, xte = its - ims, ite - ims + 1
        yts, yte = jts - jms, jte - jms + 1
        zts, zte = kts - kms, kte - kms + 1

        xs, xe = its - ids, ite - ids + 1
        ys, ye = jts - jds, jte - jds + 1
        zs, ze = kts - kds, kte - kds + 1

        nx = ide - ids + 1
        ny = jde - jds + 1
        nz = kde - kds + 1

        if ims==ids:
            its = ids
        if ime==ide:
            ite = ide

        if jms==jds:
            jts = jds
        if jme==jde:
            jte = jde

        for v in var_mems:
            dims   = d[v].dims
            existing_mem = mp.shared_memory.SharedMemory(name=var_mems[v].name)
            x_off, y_off = get_dim_offset(dims)

            if len(dims) == 2:
                data = np.ndarray((ny + y_off, nx + x_off),buffer=existing_mem.buf)
                data[ys:ye, xs:xe] = d[v].values[yts:yte, xts:xte]
            if len(dims) == 3:
                data = np.ndarray((d.dims[dims[0]], ny + y_off, nx + x_off),buffer=existing_mem.buf)
                if dims[0] == "time":
                    data[:, ys:ye+y_off, xs:xe+x_off] = d[v].values[:, yts:yte+y_off, xts:xte+x_off]
                else:
                    data[zs:ze, ys:ye+y_off, xs:xe+x_off] = d[v].values[zts:zte, yts:yte+y_off, xts:xte+x_off]
            if len(dims) == 4:
                nt = d.dims[dims[0]]
                nz = d.dims[dims[1]]
                data = np.ndarray((nt, nz, ny + y_off, nx + x_off),buffer=existing_mem.buf)
                data[:,zs:ze, ys:ye+y_off, xs:xe+x_off] = d[v].values[:,zts:zte, yts:yte+y_off, xts:xte+x_off]
            existing_mem.close()
        return

def main(file_search,cpus,var_names):
    first_files = glob.glob(file_search.format(ens="000001"))
    first_files.sort()

    for f in first_files:
        agg_file(f,var_names,verbose=True)
    
def continuous(file_search,cpus,var_names):
    print("Running continuous aggregation, Ctrl-C to stop")
    while True:
        first_files = glob.glob(file_search.format(ens="000001"))
        first_files.sort()

        # skip the last file in the list as ICAR might still be running
        for f in first_files[:-1]:
            agg_file(f, var_names,verbose=False)
            
        time.sleep(10)

if __name__ == '__main__':

    # This should be an input, this is the search string that is assumed to match
    # the output files to be aggregated.
    file_search = "icar_out_{ens}_*"
 
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--n_cpus", type=int,
                    help="Number of cpus to use")
    parser.add_argument("--continuous", action="store_true",
                    help="Use continuous aggregation of files")
    parser.add_argument("-v", "--vars", type=str,
                    help="File containing var names to save (comma-delimited; .csv)")
    parser.add_argument("-s", "--search_string", type=str,
                    help="Format of search string")
    
    args = parser.parse_args()

    cpus = mp.cpu_count()    
    if ( args.n_cpus and args.n_cpus > 0 and args.n_cpus < cpus): cpus = args.n_cpus
    
    if ( args.vars): var_names = np.loadtxt(args.vars,dtype=str,delimiter=',')
    else: var_names = None
        
    if (args.search_string): file_search = args.search_string
       
    mp.set_start_method('spawn')
    pool = mp.Pool(cpus)
 
    if args.continuous:
        try:
            continuous(file_search,cpus,var_names)
        except KeyboardInterrupt:
            pass
    else:
        main(file_search,cpus,var_names)
    pool.close()
