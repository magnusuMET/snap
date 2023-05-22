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
'''
Created on Nov 23, 2016

The eemep model runner parses a input directory for model-setup files (volcano.xml/npp.xml)
and runs the model on a super-computer. eemepModelRunner should be started from an
environment which starts the model-runner regularly, e.g. crontab

@author: heikok
'''

import atexit
import datetime
import os
import pwd
import sys
import traceback
from pathlib import Path

from METNO.HPC import typed_property, HPC
from Snappy.EEMEP.ModelRunner import ModelRunner
from Snappy.Utils import dirIsWritable, delete_oldfiles


def _cleanupFileCallable(filename):
    '''closure for cleaning up a file atexit'''
    def cleanup():
        os.unlink(filename)
    return cleanup

class EemepModelRunner():
    hpc = typed_property("hpc", HPC)
    directory = typed_property("directory", str)
    directory2 = typed_property("directory2", str)
    statusfile = typed_property("statusfile", str)
    dryrun = typed_property("dryrun", bool)

    def _getLogger(self, dirname):
        if not hasattr(self, '_logger'):
            self._logger = ModelRunner.getLogger(path=dirname)
        return self._logger

    def _log(self, dirname, msg):
        self._getLogger(dirname).debug(msg)

    def _findVolcanoFiles(self):
        '''find volcano.xml/npp.xml files in directory and directory2'''
        volcanoes = []
        dirs = [self.directory]
        if self.directory2:
            dirs.append(self.directory2)
        for directory in dirs:
            if directory:
                for (dirpath, dirnames, filenames) in os.walk(directory):
                    for x in filenames:
                        if x in self.xml_filesnames:
                            file = os.path.join(dirpath, x)
                            volcanoes.append((file, dirpath))
        return volcanoes

    def _copyResults(self, dirPath):
        '''make a copy from dirPath to directory2 (or directory)'''
        if self.directory2:
            outdir = self.directory2
            if dirPath.startswith(self.directory2):
                outdir = self.directory
            indir = dirPath.rstrip('/') #remove trailing slash to copy directory-name, too
            os.system("rsync -a -O --no-perms {indir} {outdir}".format(indir=indir, outdir=outdir))


    def __init__(self, directory, hpc, directory2, dryrun=False):
        # self.xml_filename = "npp.xml" if npp else "volcano.xml"
        self.npp_xml_filename = "npp.xml"
        self.volcano_xml_filename = "volcano.xml"
        self.xml_filesnames = [self.volcano_xml_filename, self.npp_xml_filename]

        self.dryrun = dryrun
        if dirIsWritable(directory):
            self.directory = directory
            if dirIsWritable(directory2):
                self.directory2 = directory2
            else:
                if (self.dryrun):
                    print("directory2: '{}' not writable and disabled".format(directory2), file=sys.stderr)
                self.directory2 = ""
        elif dirIsWritable(directory2):
            if (self.dryrun):
                print("directory: '{}' not writable and disabled, using '{}' as default ".format(directory, directory2), file=sys.stderr)
            self.directory = directory2
            self.directory2 = ""
        else:
            raise Exception("{dir1} and {dir2} not writable".format(dir1=directory, dir2=directory2))

        self.statusfile = os.path.join(self.directory, "eemepModelRunner_working")
        # make sure only one instance is running, not failsafe (no flock on lustre)
        if (os.path.exists(self.statusfile)):
            file_modified = datetime.datetime.fromtimestamp(os.lstat(self.statusfile).st_mtime)
            if (self.dryrun):
                with open(self.statusfile, 'rt') as fh: msg = fh.read()
                print("status-file exists at '{}' with:".format(self.statusfile), file=sys.stderr)
                print(msg, file=sys.stderr)
            else:
                if datetime.datetime.now() - file_modified > datetime.timedelta(hours=3):
                    # return statusfile if hanging for more than 3 hours
                    print("cleaning up {} after 3 hours".format(self.statusfile), file=sys.stderr)
                    _cleanupFileCallable(self.statusfile)()
                return
        else:
            if not self.dryrun:
                with open(self.statusfile, 'wt') as fh:
                    atexit.register(_cleanupFileCallable(self.statusfile))
                    fh.write("working pid: {} on node: {}\n".format(os.getpid(), os.uname().nodename))
        self.hpc = HPC.by_name(hpc)


        volcanoes = self._findVolcanoFiles()
        if self.dryrun:
            print("volcano stack:", file=sys.stderr)
            for (ifile,dpath) in sorted(volcanoes, key=lambda x: os.lstat(x[0]).st_mtime, reverse=True):
                print(ifile, sys.stderr)
            return
        if (len(volcanoes) == 0):
            return

        # work only on newest volcano.xml/npp.xml file and exit, avoiding timeouts
        (infile, dirpath) = max(volcanoes, key=lambda x: os.lstat(x[0]).st_mtime)
        with open(self.statusfile, 'at') as fh:
            fh.write("current volcano: '{}'\n".format(infile))
            try:
                pwstruct = pwd.getpwuid(os.lstat(infile).st_uid)
                fh.write("  started by: {} ({})\n".format(pwstruct.pw_gecos, pwstruct.pw_name))
            except:
                fh. write("  started by: unknown\n")


        xml_filename = Path(infile).name

        try:
            self._log(dirpath, 'Starting run for: {}'.format(infile))
            run_npp = xml_filename == self.npp_xml_filename
            mr = ModelRunner(dirpath, hpc, npp=run_npp)
            os.rename(infile,
                      os.path.join(dirpath,"{}.{}".format(xml_filename, datetime.datetime.now().strftime("%Y%m%dT%H%M%SZ"))))
            mr.work()
            self._log(dirpath, 'Finished run')
        except Exception as ex:
            print('Run failed in {}: {}'.format(dirpath,ex.args),file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            self._log(dirpath, 'Run failed: {}'.format(ex.args))
            if (os.path.exists(infile)):
                os.rename(infile,
                          os.path.join(dirpath,"{}.{}".format(xml_filename, datetime.datetime.now().strftime("%Y%m%dT%H%M%SZ"))))
        try:
            self._copyResults(dirpath)
        except:
            pass
        return



def main():
    os.umask(0)
    import argparse
    parser = argparse.ArgumentParser(description="Search for volcano.xml/npp.xml files in the top-level directory, take the newest one and run the model on a HPC machine")
    parser.add_argument("--dir", help="top-level dir to search for volcano.xml/npp.xml files", required=True)
    parser.add_argument("--dir2", help="second top-level dir, main-dir if --dir is not writable")
    parser.add_argument("--hpc", help="HPC-machine to run job on", required=True)
    parser.add_argument("--cleanup", type=int, default=0, help="remove files in dir older than cleanup days")
    parser.add_argument("--dryrun", action='store_true', default=False, help="just test what would be done, don't do anything but write debug messages")
    args = parser.parse_args()
    if "dir2" in args:
        dir2 = args.dir2
    else:
        dir2 = ""
    EemepModelRunner(directory=args.dir, hpc=args.hpc, directory2=dir2, dryrun=args.dryrun)

    if (args.cleanup > 0) and not args.dryrun:
        if dirIsWritable(args.dir):
            delete_oldfiles(args.dir, args.cleanup)
        if dir2 and dirIsWritable(dir2):
            delete_oldfiles(dir2, args.cleanup)


if __name__ == '__main__':
    main()
