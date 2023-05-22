#! /usr/bin/env python3
#
# SNAP: Servere Nuclear Accident Programme
# Copyright (C) 1992-2017   Norwegian Meteorological Institute
# 
# This file is part of SNAP. SNAP is free software: you can 
# redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the 
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import argparse
import netCDF4
import numpy
import re
import sys

def finishOutput(out_nc, isotope, times, conc, prob):
    probname = "{iso}_probability".format(iso=isotope)
    accprobname = "{iso}_acc_probability".format(iso=isotope)
    out_nc.variables[probname][:] = prob

    accprob = prob.copy()
    for idx, val in enumerate(times):
        if (idx != 0):
            accprob[idx,:,:] += accprob[idx-1,:,:]

    out_nc.variables[accprobname][:] = accprob

    out_nc.close()

def addFile(ncFile, isotope, times, conc, prob, isMeasurement=True):
    print(ncFile)
    with netCDF4.Dataset(ncFile, "r") as nc:
        ftimes = netCDF4.num2date(nc.variables['time'][:],units=nc.variables['time'].units).tolist()
        varname = "{iso}_avg_concentration_bl".format(iso=isotope)
        if not varname in nc.variables:
            raise Exception("variable {vn} not in {file}".format(vn=varname,file=ncFile))
        co = nc.variables[varname][:]
        pr = numpy.zeros(co.shape, dtype=int)
        nu = numpy.zeros(co.shape, dtype=int)
        if (isMeasurement):
            pr = numpy.select([co > 0], [numpy.ones(co.shape, dtype=numpy.int16)], numpy.zeros(1,dtype=numpy.int16))
        else:
            pr = numpy.select([co <= 0], [numpy.ones(co.shape, dtype=numpy.int16)], numpy.zeros(1,dtype=numpy.int16))
            nu += 1
        for i,t in enumerate(times):
            try:
                j = ftimes.index(t)
                conc[i,:] *= co[j,:]
                prob[i,:] *= pr[j,:]
            except ValueError as ve:
                print("{} not defined in file: {}".format(t,ncFile), file=sys.stderr)
                prob[i,:] *= nu[0,:]
    return (conc, prob)

def initializeNc(ncFile, outNcFile, isotope):
    print(ncFile)
    times = []
    conc = None
    prob = None
    out_nc = netCDF4.Dataset(outNcFile, "w", format="NETCDF4_CLASSIC")

    with netCDF4.Dataset(ncFile, "r") as nc:
        times = netCDF4.num2date(nc.variables['time'][:],units=nc.variables['time'].units)
        varname = "{iso}_concentration".format(iso=isotope)
        probname = "{iso}_probability".format(iso=isotope)
        accprobname = "{iso}_acc_probability".format(iso=isotope)
        if not varname in nc.variables:
            raise Exception("variable {vn} not in {file}".format(vn=varname,file=ncFile))
        conc = nc.variables[varname][:]
        prob = numpy.select([conc > 0], [numpy.ones(conc.shape, dtype=numpy.int16)], numpy.zeros(1,dtype=numpy.int16))

        # initialize the output file
        # create dimensions
        for dim in nc.dimensions.values():
            size = dim.size
            if (dim.isunlimited()):
                size=None
            out_nc.createDimension(dim.name, size)
        # global attributes
        for attr in nc.ncattrs():
            try:
                out_nc.setncattr(attr, nc.getncattr(attr))
            except:
                print("unable to set global attribute '{}', ignoring".format(attr), file=sys.stderr)
        for v in nc.variables.values():
            # avoid concentrations, depositions ...
            if re.search(r'(concentration|deposition|precipitation|height|air|arrival)', v.name):
                continue
            out_v= out_nc.createVariable(varname=v.name, datatype=v.dtype, dimensions=v.dimensions,
                                         zlib=True, shuffle=True)
            # variable-attributes
            for attr in v.ncattrs():
                try:
                    out_v.setncattr(attr, v.getncattr(attr))
                except:
                    print("unable to set attribute '{attr}' of variable '{var}'".format(attr=attr, var=v.name))

        # add extra variables for probabilities
        vn = nc.variables[varname]
        for pn in (probname, accprobname):
            out_v = out_nc.createVariable(varname=pn, datatype="i2", dimensions=vn.dimensions,
                                          zlib=True, shuffle=True)
            # variable-attributes
            for attr in vn.ncattrs():
                if attr != "units" and attr != "_FillValue":
                    try:
                        out_v.setncattr(attr, vn.getncattr(attr))
                    except:
                        print("unable to set attribute '{attr}' of variable '{var}'".format(attr=attr, var=pn))

        # data
        for out_v in out_nc.variables.values():
            if out_v.name in nc.variables:
                v = nc.variables[out_v.name]
                if (len(v.shape)):
                    try:
                        out_v[:] = v[:]
                    except:
                        print("unable to write data of variable '{var}'".format(var=out_v.name))


    return (times, out_nc, conc, prob)




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, action="append", help="one or more snap.nc file from an inverse run to a measurement")
    parser.add_argument("-n", action="append", help="zero or more snap.nc file from an inverse run to a stations without measurement")
    parser.add_argument("-o", required=True, help="output nc file")
    parser.add_argument("-I", default="I131", help="isotope (default: I131)")
    args = parser.parse_args()

    (times, out_nc, conc, prob) = initializeNc(args.i[0], args.o, args.I)
    # remove first input
    args.i.pop(0)
    for i in args.i:
        (conc, prob) = addFile(i, args.I, times, conc, prob)

    if args.n:
        for i in args.n:
            (conc, prob) = addFile(i, args.I, times, conc, prob, False)

    finishOutput(out_nc, args.I, times, conc, prob)

if __name__ == "__main__":
    main()
