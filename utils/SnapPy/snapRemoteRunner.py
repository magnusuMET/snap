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
Created on Mar 2, 2018

The snap remote runner reads input files from the ARGOS DSS, uploaded to a input directory on
a remote machine and runs the model on a HPC-like machine with queue system. snap remote runner 
should be started from an environment which starts the runner regularly, e.g. crontab

Structure below directory is:   
     upload/ : files uploaded to remote machine upload directory
     results/: files delivered to remote machine root directory
     runs/:    job-directory, with each job in an own directory
     rejected/ : files rejected 

@author: heikok
'''

from METNO.HPC import typed_property, HPC
from Snappy.Utils import delete_oldfiles, dirIsWritable
import atexit
import datetime
import os
import re
import sys
import traceback
import zipfile

from Snappy.SnapJobEC import SnapJobEC


def _cleanupFileCallable(filename):
    '''closure for cleaning up a file atexit'''
    def cleanup():
        os.unlink(filename)
    return cleanup

class SnapTask():
    topdir = typed_property('topdir', str)
    zipfile = typed_property('zipfile', str)
    model = typed_property('model', str)
    id = typed_property('id', str)
    scpdestination = typed_property('scpdestination', str)
    timestamp = typed_property('timestamp', datetime.datetime)
    rundir = typed_property('rundir', str)

    def __init__(self, topdir, zip_file, model, ident, scpdestination):
        self.topdir = topdir
        self.zipfile = zip_file
        self.model = model
        self.id = ident
        self.scpdestination = scpdestination
        self.timestamp = datetime.datetime.now()

    def status_filename(self):
        return "{ident}_{model}_status".format(ident=self.id, model=self.model)

    def is_complete(self, reldir):
        infile = os.path.join(self.topdir, reldir, self.zipfile)
        try:
            with zipfile.ZipFile(infile, 'r') as zf:
                if zf.testzip() == None:
                    return True
        except:
            pass
        return False
    
    def handle(self, hpc):
        ''' Handle the job on the hpc. Return True if job is submitted '''
        retval = False
        try:
            retval = self.handle(hpc)
        except:
            traceback.print_exc()
        return retval
            
    def _handle(self, hpc): 
        top_rundir = os.path.join(self.topdir, SnapRemoteRunner.RUN_DIR)
        if not os.path.isdir(top_rundir):
            os.mkdir(top_rundir)
        self.rundir = os.path.join(top_rundir, "{dt}_{ident}".format(dt=self.timestamp.strftime('%Y-%m-%dT%M%H%S'),
                                                                     ident=self.id))
        os.mkdir(self.rundir)
        infile = os.path.join(self.topdir, SnapRemoteRunner.UPLOAD_DIR, self.zipfile)
        workfile = os.path.join(self.rundir, self.zipfile)
        os.rename(infile, workfile)
        
        with zipfile.ZipFile(workfile, 'r') as zf:
            # files = zf.namelist()
            zf.extractall(path=self.rundir)
        
        # start a remote detached qsub job
        snapJob = SnapJobEC(self, hpc)
        jobscript = snapJob.job_script()
        jobfile = os.path.join(self.rundir, 'snap.job')
        with open(jobfile, 'w') as jh:
            jh.write(jobscript)
        # push the job into the queue, no feedback
        qjob = hpc.submit_job(jobfile)
        if (qjob == None):
            return False
        return True
        


class SnapRemoteRunner():
    UPLOAD_DIR = 'upload'
    RUN_DIR = 'run'
    REJECTED_DIR = 'rejected'
    WORK_DIR = 'work'
    
    
    hpc = typed_property("hpc", HPC)
    directory = typed_property("directory", str)
    directory2 = typed_property("directory2", str)
    dryrun = typed_property("dryrun", bool)
    remote = typed_property("remote", str)
    remote_dir = typed_property("remote_dir", str)
    remote_user = typed_property("remote_user", str)
    statusfile = typed_property("statusfile", str)

    def __init__(self, directory, hpc, directory2, remote, remoteUser, remoteDir, dryrun=False):
        self.dryrun = dryrun
        self.hpc = HPC.by_name(hpc)
        self.remote = remote
        self.remote_user = remoteUser
        self.remote_dir = remoteDir
        
        self.scpdestination = "{remote}:{remoteDir}".format(remote=self.remote, remoteDir=self.remote_dir)
        if self.remote_user:
            self.scpdestination = self.remote_user + '@' + self.scpdestination
         
        
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

        self.statusfile = os.path.join(self.directory, "snapRemoteRunner_working")
        # make sure only one instance is running, not failsafe (no flock on lustre, eventually in different directories, but good enough)
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


    def write_status(self, task, tag, msg=""):
        filename = task.status_filename()
        work_file = os.path.join(self.directory, self.WORK_DIR, filename)
        timestamp = datetime.datetime.now().strftime('%Y%m%d%H%m')
        with open(work_file, 'a') as fh:
            if (tag == 'downloading'):
                fh.write(":{ts} downloading run-description {file}\n".format(ts=timestamp, file=task.zipfile))
            elif (tag == 'success'):
                fh.write(":Finished extracting {model} data for ARGOS\n".format(model=task.model))
            elif (tag == 'error'):
                fh.write(":{model} output data do not exist\n".format(model=task.model))
            elif (tag == 'running'):
                fh.write("{ts} running\n".format(ts=timestamp))
            elif (tag == 'internal'):
                fh.write(":internal error, cannot start job in queue in dir '{rundir}'".format(rundir=task.rundir))
            else:
                fh.write("{tag}:{ts} {msg}\n".format(ts=timestamp, tag=tag, msg=msg))
        self.hpc.put_files([work_file], self.remote_dir, 30)

    def _check_and_unpack_new_files(self):
        '''Download new files from the remote machine to the upload directory.
        - Move invalid files to rejected. (Wrong name, not containing *ARGOS2*.zip)
        - Unpack zip-files in project-folder / delete ignore incomplete files.
            - status for complete and incomplete files
        - Remove complete files from remote-upload and local upload 
        - create modelruns
        
        throws an exception when download / unpack failed unexpectedly
        '''
        remote_files = os.path.join(self.UPLOAD_DIR, '*')
        local_upload = os.path.join(self.directory, self.UPLOAD_DIR)
        if not os.path.isdir(local_upload):
            os.mkdir(local_upload)
        local_rejected = os.path.join(self.directory, self.REJECTED_DIR)
        if not os.path.isdir(local_rejected):
            os.mkdir(local_rejected)
        self.hpc.get_files(remote_files, local_upload, 30)
        
        delete_in_upload = []
        for f in os.listdir(local_upload):
            if os.path.isfile(os.path.join(local_upload, f)):
                m = re.match(r'([\w\-\.:]*)_ARGOS2(.*)\.zip', f)
                if m:
                    task = SnapTask(topdir=self.directory,
                                    zipfile=f,
                                    ident=m.group(1),
                                    model=m.group(2),
                                    scpdestination=self.scpdestination)
                    if task.is_complete(reldir=self.UPLOAD_DIR):
                        if task.handle(self.hpc):
                            self.write_status(task, tag='running')
                        else:
                            self.write_status(task, tag='internal')
                        delete_in_upload.append(f)
                    else:
                        self.write_status(task, tag='downloading')
                else:                    
                    os.rename(os.path.join(local_upload, f), os.path.join(local_rejected, f))
                    delete_in_upload.append(f)

        delete_upload_files = [os.path.join(self.UPLOAD_DIR, f) for f in delete_in_upload]
        self.hpc.syscall('rm', delete_upload_files, 30)


if __name__ == '__main__':
    os.umask(0)
    import argparse
    parser = argparse.ArgumentParser(description="Read snap-job description files from a remote machine, and run the model on a HPC-like machine")
    parser.add_argument("--dir", help="top-level working dir", required=True)
    parser.add_argument("--dir2", help="secondary top-level dir, main-dir if --dir is not writable")
    parser.add_argument("--remote", help="remote machine name", required=True)
    parser.add_argument("--remoteUser", help="remote username", required=True)
    parser.add_argument("--remoteDir", help="remote main directory", default="")
    parser.add_argument("--hpc", help="HPC-machine to run job on", required=True)
    parser.add_argument("--cleanup", type=int, default=0, help="remove files in dir older than cleanup days")
    parser.add_argument("--dryrun", action='store_true', default=False, help="just test what would be done, don't do anything but write debug messages")
    args = parser.parse_args()
    if "dir2" in args:
        dir2 = args.dir2
    else:
        dir2 = ""
    SnapRemoteRunner(directory=args.dir, hpc=args.hpc, directory2=dir2, 
                     remoteDir=args.remoteDir, remoteUser=args.remoteUser, remote=args.remote,
                     dryrun=args.dryrun)

    if (args.cleanup > 0) and not args.dryrun:
        if dirIsWritable(args.dir):
            delete_oldfiles(args.dir, args.cleanup)
        if dir2 and dirIsWritable(dir2):
            delete_oldfiles(dir2, args.cleanup)
