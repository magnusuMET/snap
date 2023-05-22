#! /usr/bin/env python3
#
# SNAP: Servere Nuclear Accident Programme
# Copyright (C) 1992-2020   Norwegian Meteorological Institute
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
from collections import OrderedDict
import datetime
import os
import re
import subprocess
import sys

from Snappy.EcMeteorologyCalculator import EcMeteorologyCalculator
from Snappy.ICONMeteorologyCalculator import ICONMeteorologyCalculator
from Snappy.MeteorologyCalculator import MeteorologyCalculator
from Snappy.Resources import MetModel, Resources, snapNc_convert_to_grib
import xml.etree.ElementTree as ET


def parseIsoTimeDelta(period):
    """parse string like P2DT12H5M to timedelta-objects"""
    regex = re.compile(r'P((?P<days>\d+?)D)?T?((?P<hours>\d+?)H)?((?P<minutes>\d+?)M)?((?P<seconds>\d+?)S)?')
    match = regex.match(period)
    if not match:
        return
    time_params = {}
    for (name, param) in match.groupdict().items():
        if param:
            time_params[name] = int(param)
    return datetime.timedelta(**time_params)

def rimstermGetIsotopes(rimxml):
    """return ordered dictionary of isotopeNames (used by rimsterm) and isotopeIds"""
    isoItems = OrderedDict()
    for itemXML in rimxml.findall("ItemsTot/Item"):
        isoItems[itemXML.attrib['ItemName']] = int(itemXML.attrib['ItemNum'])
    return isoItems

def rimsterm2input(rimxml, argosxml, met_model):
    """create a full snap.input sourceterm from a rimsterm.xml file using the ec-nrpa meteo"""

    sourceTermTemplate = """
SET_RELEASE.POS= P=   {lat},   {lon}
TIME.START= {startTime}
TIME.RUN = {runTime:d}h
STEP.HOUR.OUTPUT.FIELDS= {outputTimestep:d}
"""
    res = Resources()
    lat = rimxml.find("Header/Place/Location/[@Latitude]").attrib['Latitude']
    lon = rimxml.find("Header/Place/Location/[@Longitude]").attrib['Longitude']

    startTime = rimxml.find("TimeOfInitialRelease").text
    runTime = 48
    outputTimestep = 3    
    if argosxml:
        runTime = int(argosxml.find("RunLength").text)
        outputTimestep = int(argosxml.find("OutputTimestep").text)
    
    releaseStartTime = datetime.datetime.strptime(startTime, '%Y-%m-%dT%H:%M:%SZ')
    startDTime = releaseStartTime.replace(minute=0, second=0, microsecond=0) # snap requires hourly start
    
    isoItems = rimstermGetIsotopes(rimxml)
    isos = res.getIsotopes()
    releases = []
    offsetTime = datetime.timedelta(hours=0)
    if releaseStartTime > startDTime:
        # add a 0 release from model start to release start
        # this handling does not cover non-hourly backward runs yet
        offsetTime = releaseStartTime-startDTime
        zeroIsoBqs = {}
        for isoName in isoItems.keys():
            zeroIsoBqs[isoName] = 0
        releaseDef =  {'endtime': offsetTime, 'lower': '0', 'upper': '1', 'isoBqs': zeroIsoBqs}
        releases.append(releaseDef)

    for relXML in rimxml.findall("Source/TimeDependent/ReleaseInterval"):
        # add real releases
        endtime = parseIsoTimeDelta(relXML.find("SourceTime/[@EndTime]").attrib['EndTime']) + offsetTime
        # offset endtime
        if runTime < 0: # backward run
            # StartTime is earliest measuremnt, Endtime end of measurement after starttime
            # snap needs to start at the end of the measurements, i.e. the last time 
            #     and run backupward in time
            startDTime = startDTime + endtime
        posXML = relXML.find("SourcePosition")
        lHeight = posXML.attrib['HeightAboveGround']
        uHeight = posXML.attrib['HeightAboveGroundMax']
        if uHeight < lHeight:
            uHeight = lHeight
        isoBqs = {}
        for isoName in isoItems.keys():
            xpath = 'SourceStrength[@ItemName="{}"]'.format(isoName)
            relValXML = relXML.find(xpath)
            if relValXML:
                isoBqs[isoName] = float(relValXML.find('BinStrength').attrib['Value'])
            else:
                isoBqs[isoName] = 0
        releaseDef =  {'endtime': endtime, 'lower': lHeight, 'upper': uHeight, 'isoBqs': isoBqs}
        if runTime >= 0:
            releases.append(releaseDef)
        else:
            releases.insert(0, releaseDef)

    timevals = ["0"]
    radiusvals = []
    lowervals = []
    uppervals = []
    for rel in releases:
        timevals.append("{:.2f}".format(rel['endtime'].total_seconds()/(60*60)))
        radiusvals.append("50")
        lowervals.append(rel['lower'])
        uppervals.append(rel['upper'])
    radiusvals.append("0")
    lowervals.append("0")
    uppervals.append("0")
    releaseTerm = "TIME.RELEASE.PROFILE.STEPS\n"
    releaseTerm += "MAX.PARTICLES.PER.RELEASE= {:d}\n".format(min(5000, 400*len(isoItems)))
    releaseTerm += "MAX.TOTALPARTICLES= {:d}\n".format(20000000) # rough estimate
    releaseTerm += "RELEASE.HOUR= {}\n".format(", ".join(timevals))
    releaseTerm += "RELEASE.RADIUS.M= {}\n".format(", ".join(radiusvals))
    releaseTerm += "RELEASE.LOWER.M= {}\n".format(", ".join(lowervals))
    releaseTerm += "RELEASE.UPPER.M= {}\n".format(', '.join(uppervals))
    for isoName, isoId in isoItems.items():
        vals = []
        for rel in releases:
            vals.append("{:10.3E}".format(rel['isoBqs'][isoName]))
        vals.append("0")
        vals.append("'{}'".format(isos[isoId]['isotope']))
        releaseTerm += "RELEASE.BQ/SEC.COMP= {}\n".format(','.join(vals))


    snapStartTime="{year} {month} {day} {hour}".format(year=startDTime.year,month=startDTime.month,day=startDTime.day, hour=startDTime.hour)
    sourceTerm = sourceTermTemplate.format(lat=lat,lon=lon,startTime=snapStartTime, runTime=runTime, outputTimestep=outputTimestep)
    sourceTerm += releaseTerm + "\n"

    sourceTerm += res.isotopes2snapinput(isoItems.values())

    metdef = res.getDefaultMetDefinitions(met_model)
    if (met_model == MetModel.NrpaEC0p1):
        files = res.getECMeteorologyFiles(dtime=startDTime, run_hours=int(runTime))
        if (len(files) == 0):
            raise Exception("no EC met-files found for {}, runtime {}".format(startDTime, runTime))
    elif met_model == MetModel.NrpaEC0p1Global:
        ecmet = EcMeteorologyCalculator(EcMeteorologyCalculator.getGlobalMeteoResources(), startDTime, float(lon), float(lat))
        ecmet.calc()
        if ecmet.must_calc():
            raise Exception("no EC met-files calculated for {}".format(startDTime))
        files = ecmet.get_meteorology_files()
        (metdef["startX"], metdef["startY"]) = ecmet.get_grid_startX_Y()
    elif met_model == MetModel.EC0p1Global:
        globalRes = EcMeteorologyCalculator.getGlobalMeteoResources()
        files = [x[1] for x in sorted(MeteorologyCalculator.findAllGlobalData(globalRes), key=lambda x: x[0])]
        lat0 = MeteorologyCalculator.getLat0(float(lat), globalRes.domainHeight)
        lon0 = MeteorologyCalculator.getLon0(float(lon), globalRes.domainWidth)
        sourceTerm += f"FIELD.TYPE=fimex\n"
        sourceTerm += f"FIMEX.FILE_TYPE=netcdf\n"
        sourceTerm += f"FIMEX.INTERPOLATION=nearest|+proj=latlon +R=6371000 +no_defs|{lon0},{lon0+0.2},...,{lon0+globalRes.domainWidth}|{lat0},{lat0+0.2},...,{lat0+globalRes.domainHeight}|degree\n"
    elif met_model == MetModel.Meps2p5:
        files = res.getMeteorologyFiles(met_model, startDTime, runTime, "best")
        if (len(files) == 0):
            raise Exception("no MEPS2_5 met-files found for {}, runtime {}".format(startDTime, runTime))
    elif met_model == MetModel.Icon0p25Global:
        metcalc = ICONMeteorologyCalculator(ICONMeteorologyCalculator.getGlobalMeteoResources(), startDTime, float(lon), float(lat))
        metcalc.calc()
        if metcalc.must_calc():
            raise Exception("no ICON met-files calculated for {}".format(startDTime))
        files = metcalc.get_meteorology_files()
        (metdef["startX"], metdef["startY"]) = metcalc.get_grid_startX_Y()
    elif met_model == MetModel.GfsGribFilter:
        files = res.getMeteorologyFiles(met_model, startDTime, runTime, "best")
        if (len(files) == 0):
            raise Exception("no GFS-grib-filter-fimex met-files found for {}, runtime {}".format(startDTime, runTime))
    else:
        raise Exception('unknown met_model: {}'.format(met_model))

    sourceTerm += res.getSnapInputMetDefinitions(met_model, files, **metdef)

    return sourceTerm

def snap4rimsterm(rimsterm, argosrequest, basedir, ident, met_model, bitmapCompress):
    '''run snap with rimsterm definition in basedir dir'''
    tree = ET.parse(rimsterm)
    root = tree.getroot()
    assert root.tag == 'CBRN_Sourceterm', "Not a rimsterm input file: {}".format(rimsterm)

    argosRoot = None
    if argosrequest:
        argosTree = ET.parse(argosrequest)
        argosRoot = argosTree.getroot()
        assert argosRoot.tag == 'Request', "Not a argos request xml file: {}". format(argosrequest)

    if (not os.path.isdir(basedir)):
        os.makedirs(basedir)
    snapNC = os.path.join(basedir, 'snap.nc')
    if os.path.exists(snapNC):
        os.remove(snapNC)

    snapInput = "title={}\n".format(ident) + rimsterm2input(root, argosRoot, met_model)
    with open(os.path.join(basedir, "snap.input"), "w") as fh:
        fh.write(snapInput)

    errlog = open(os.path.join(basedir, "snap.errlog"), "a")
    outlog = open(os.path.join(basedir, "snap.outlog"), "a")
    print("bsnap_naccident snap.input")
    proc = subprocess.Popen(['bsnap_naccident', 'snap.input'], cwd=basedir, stderr=errlog, stdout=outlog)
    if (proc.wait() != 0):
        if os.path.exists(snapNC):
            errlog.write("bsnap_naccident in {} did not finished completely. Continuing anyway.\n".format(basedir))
        else:
            raise Exception("bsnap_naccident in {} failed".format(basedir))

    print("snapAddToa snap.nc")
    proc2 = subprocess.Popen(['snapAddToa', 'snap.nc'], cwd=basedir, stderr=errlog, stdout=outlog)
    if (proc2.wait() != 0):
        errlog.write("snapAddToa snap.nc in {} failed\n".format(basedir))
    errlog.close()
    outlog.close()

    if argosRoot:
        for oformat in argosRoot.findall("OutputFormat"):
            print("output in format: {}".format(oformat.text))
            if oformat.text == 'GRIB':
                snapNc_convert_to_grib(snapNC, basedir, ident, rimstermGetIsotopes(root).values(), bitmapCompress=bitmapCompress)
            elif oformat.text == 'NetCDF':
                # create a filename which can be picked up by SMS
                os.symlink(snapNC, os.path.join(basedir,"{}_all.nc".format(ident)))
            else:
                raise Exception("unknown OutputFormat in request: {}".format(oformat.text))
    else:
        snapNc_convert_to_grib(snapNC, basedir, ident, rimstermGetIsotopes(root).values(), bitmapCompress=bitmapCompress)



def main():
    os.umask(0) # make sure files can be deleted later
    import argparse
    parser = argparse.ArgumentParser(description="run snap from a rimsterm.xml file and convert to grib-files named ident_conc, ident_depo ident_wetd ident_dose ident_tofa")
    parser.add_argument("--rimsterm", help="source-term in rimsterm format", required=False)
    parser.add_argument("--argosrequest", help="optional argos-request xml file in addition to --rimsterm file", required=False)
    parser.add_argument("--metmodel", help="select the NWP input model, nrpa_ec_0p1, nrpa_ec_0p1_global, meps_2_5km, ec_0p1_global, gfs_grib_filter_fimex, icon_0p25_global", required=True, default='nrpa_ec_0p1')
    parser.add_argument("--dir", help="output-directory", required=True)
    parser.add_argument("--ident", help="output-file identifier", required=True)
    parser.add_argument("--bitmapCompress", help="enable grib bitmap-compression", action='store_true')

    args = parser.parse_args()
    if args.rimsterm is not None:
        snap4rimsterm(args.rimsterm, args.argosrequest, args.dir, args.ident, args.metmodel, bitmapCompress=args.bitmapCompress)
    else:
        print("need rimsterm option", file=sys.stderr)


if __name__ == "__main__":
    main()
