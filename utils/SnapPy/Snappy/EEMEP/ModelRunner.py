'''
Created on Sep 2, 2016

@author: heikok
'''
import datetime
import os
import re
import subprocess
from time import sleep
import unittest

from METNO.HPC import HPC, StatusFile, QJobStatus
from Snappy.EEMEP.Resources import Resources
from Snappy.EEMEP.VolcanoRun import VolcanoRun
from posix import R_OK


class ModelRunner():


    def __init__(self, path, hpcMachine):
        self.upload_files = set()
        self.timestamp = datetime.datetime.now()
        self.jobscript = "eemep_script.job"
        self.statusfile = "status"
        self.res = Resources()
        self.hpc = HPC.by_name(hpcMachine)
        self.hpcMachine = hpcMachine
        self.inpath = path

        self.rundir = self.res.HPC[self.hpcMachine]["RUNDIR"]
        volcano_path = os.path.join(path, "volcano.xml")
        if not os.path.exists(volcano_path):
            raise Exception("no such file or directory: {}".format(volcano_path))
        self.volcano = VolcanoRun(volcano_path)
        self.runtag = "eemep_{}".format(os.path.basename(self.volcano.outputDir))
        self.hpc_outdir = os.path.join(self.rundir, self.runtag)

        self.path = self.volcano.outputDir
        os.makedirs(name=self.path, exist_ok=True)
        self._volcano_to_column_source()
        self._get_meteo_files()
        self._get_restart_file()
        self._create_job_script()

    def _write_log(self, msg):
        with open(os.path.join(self.inpath, 'volcano.log'), 'a') as fh:
            print("{}: {}".format(datetime.datetime.now().strftime("%Y%m%dT%H%M%SZ"), msg), file=fh)

    def _volcano_to_column_source(self):
        '''write columnsource_location.csv and columnsource_emissions.csv from volcano.xml'''
        location = os.path.join(self.path, "columnsource_location.csv")
        with open(location, 'wt') as lh:
            lh.write(self.volcano.get_columnsource_location())
        self.upload_files.add(location)

        emission = os.path.join(self.path, "columnsource_emission.csv")
        with open(emission, 'wt') as eh:
            eh.write(self.volcano.get_columnsource_emission())
        self.upload_files.add(emission)

    def _generate_meteo_file(self, outfile, date_files):
        '''Generate a meteo file, eventually by concatenting several input-files to get a file with 8 timesteps

        Args:
           outfile: output filename
           date_files: list of pairs, each pair consisting of a input-file and the number of time-steps (3hourly) containted in this timestep
        '''
        if os.path.islink(outfile):
            os.unlink(outfile)
        if os.path.isfile(outfile):
            if os.access(outfile, os.W_OK):
                os.unlink(outfile)
            elif os.access(outfile, os.R_OK):
                # file only readable, forced to that file
                return
        if (date_files[0][1] == 8):
            # no file-concatenation needed, just create a link
            if not os.path.lexists(outfile):
                os.symlink(date_files[0][0], outfile)
        else:
            # find the number of steps needed for which file (latest date first)
            timesteps_in_file = 0
            totalsteps = 0
            use_steps = []
            for (file, tsteps) in date_files:
                newsteps = tsteps - timesteps_in_file
                if newsteps <= 0:
                    continue
                use_steps.append((file, newsteps, timesteps_in_file))
                timesteps_in_file = timesteps_in_file + newsteps
                assert(timesteps_in_file <= 8)
                if (timesteps_in_file == 8):
                    break
                step = 0
            # create a list of all files needed (first date first) and
            # find the timesteps to select from the joined files. (from first date to last date)
            steps = []
            files = []
            pos = 0
            for (file, use_steps, skip_steps) in reversed(use_steps):
                for i in range(0, use_steps):
                    steps.append(pos)
                    pos = pos + 1
                pos = pos + skip_steps
                files.append('<netcdf location="{file}" />'.format(file=file))

            # run fimex on files/steps
            joinfiles = '''<?xml version="1.0" encoding="UTF-8"?>
<netcdf xmlns="http://www.unidata.ucar.edu/namespaces/netcdf/ncml-2.2"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
<aggregation type="joinExisting">
    {files}
</aggregation>
</netcdf>
'''
            ncml_file = os.path.join(self.path, 'joinMeteo.ncml')
            with open(ncml_file, 'wt') as fh:
                fh.write(joinfiles.format(files="\n".join(files)))
            subprocess.call(args=['fimex', '--input.file', ncml_file,
                                  '--output.file', outfile,
                                  '--output.type=nc4',
                                  '--extract.pickDimension.name=time',
                                  '--extract.pickDimension.list={}'.format(','.join(str(x) for x in steps))])


    def _get_meteo_files(self):
        '''Create meteorology files in the output-directory of volcano.xml.
        This involves linking and copying of needed meteorology. and eventually
        addition of a few timesteps at the beginning of the run

        Args:
            overwrite_model_start_time: don't used the volcano meteo, but the provided, used only for testing

        Returns: list of meteorology files
        '''
        (ref_date, model_start_time) = self.volcano.get_meteo_dates()

        files = self.res.getECMeteorologyFiles(model_start_time, 72, ref_date)
        for i, date_files in enumerate(files):
            file_date = model_start_time + datetime.timedelta(days=i)
            outfile = os.path.join(self.path, "meteo{date}.nc".format(date=file_date.strftime("%Y%m%d")))
            self._generate_meteo_file(outfile, date_files)
            self.upload_files.add(outfile)
        # vlevel-definition
        vlevels = self.res.getVerticalLevelDefinition()
        vfile = os.path.join(self.path, "Vertical_levels.txt")
        with open(vfile, 'w') as vh:
            vh.write(vlevels)
            self.upload_files.add(vfile)

    def _get_restart_file(self):
        if (self.volcano.run_as_restart()):
            model_start_time = self.volcano.get_meteo_dates()[1]
            restart_file = os.path.join(self.path, "EMEP_IN_{date}.nc".format(date=model_start_time.strftime("%Y%m%d")))
            if (os.path.exists(restart_file)):
                self.upload_files.add(restart_file)

    def _create_job_script(self):
        job = self.res.get_job_script(self.hpcMachine)
        defs = {"rundir": self.rundir,
                "runtag": self.runtag} # year, month, day, hour, runhour
        defs["runhour"] = "{}".format(int(self.volcano.runTimeHours))
        start_time = self.volcano.get_meteo_dates()[1]
        defs["year"] = start_time.year
        defs["month"] = "{:02d}".format(start_time.month)
        defs["day"] = "{:02d}".format(start_time.day)
        defs["hour"] = "{:02d}".format(start_time.hour)

        filename = os.path.join(self.path, self.jobscript)
        with open(filename, 'wt') as jh:
            jh.write(job.format(**defs))
        self.upload_files.add(filename)

    def do_upload_files(self):
        self._write_log("uploading to {}:{}".format(self.hpcMachine, self.hpc_outdir))
        self.hpc.syscall("mkdir", ["-p", self.hpc_outdir])
        for f in self.upload_files:
            self._write_log("uploading '{}'".format(f))
            self.hpc.put_files([f], self.hpc_outdir, 600)

    def run_and_wait(self):
        '''Start the model and wait for it to finish

        Returns QJobStatus code
        '''
        self._write_log("removing old status on {}: {}".format(self.hpcMachine, self.hpc_outdir))
        self.hpc.syscall("rm", ["-f", os.path.join(self.hpc_outdir, self.statusfile)])
        self._write_log("starting run on hpc {}: {}".format(self.hpcMachine, self.hpc_outdir))
        remote_jobscript = os.path.join(self.hpc_outdir, self.jobscript)
        qjob = self.hpc.submit_job(remote_jobscript, [])
        qjob.status_file = StatusFile(os.path.join(self.hpc_outdir, self.statusfile), "finished")

        # wait for 60 minutes to finish, check every minute
        sleep_time = 60 # seconds
        count = 60 # * sleep_time
        status = self.hpc.get_status(qjob)
        while not (status == QJobStatus.finished or status == QJobStatus.failed):
            sleep(sleep_time)
            count -= 1
            if count == 0:
                break
            status = self.hpc.get_status(qjob)
            self._write_log("jobstatus on hpc {} jobid={}: {}".format(self.hpcMachine, qjob.jobid, status))

        return status

    def download_results(self):
        '''download the result-files, and rename them as appropriate'''
        start_time = self.volcano.get_meteo_dates()[1]
        tomorrow = (start_time + datetime.timedelta(days=1)).strftime("%Y%m%d")
        timestamp = self.timestamp.strftime("%Y%m%dT%H%M%S")
        simulationstart = self.timestamp.strftime("%Y-%m-%d_%H:%M:%S")

        file = 'eemep_hourInst.nc'
        self._write_log("downloading {}:{} to {}".format(file, self.hpcMachine, self.path))
        self.hpc.get_files([os.path.join(self.hpc_outdir, file)], self.path, 1200)
        subprocess.call(args=['ncatted', '-a', 'SIMULATION_START_DATE,global,o,c,{}'.format(simulationstart),
                              os.path.join(self.path, file)])
        os.rename(os.path.join(self.path, file),
                  os.path.join(self.path, 'eemep_hourInst_{}.nc'.format(timestamp)))

        file = 'eemep_hour.nc'
        self._write_log("downloading {}:{} to {}".format(file, self.hpcMachine, self.path))
        self.hpc.get_files([os.path.join(self.hpc_outdir, file)], self.path, 1200)
        subprocess.call(args=['ncatted', '-a', 'SIMULATION_START_DATE,global,o,c,{}'.format(simulationstart),
                              os.path.join(self.path, file)])
        os.rename(os.path.join(self.path, file),
                  os.path.join(self.path, 'eemep_hour_{}.nc'.format(timestamp)))


        file = 'EMEP_OUT_{}.nc'.format(tomorrow)
        self._write_log("downloading {}:{} to {}".format(file, self.hpcMachine, self.path))
        try :
             self.hpc.get_files([os.path.join(self.hpc_outdir, file)], self.path, 1200)
        except Exception as ex:
            # not dangerous if it fail, but remove file
            self._write_log("couldn't download '{}', ignoring: {}".format(file, ex.args))
            filename = os.path.join(self.hpc_outdir, file)
            if os.path.lexists(filename): os.unlink(filename)
        else:
            os.rename(os.path.join(self.path, file),
                      os.path.join(self.path, 'EMEP_IN_{}.nc'.format(tomorrow)))

    def work(self):
        '''do the complete work, e.g. upload, run, wait and download'''
        self.do_upload_files()
        self.run_and_wait()
        self.download_results()

class TestModelRunner(unittest.TestCase):
    hpcMachine = "frost"
    doRun = False

    def setUp(self):
        unittest.TestCase.setUp(self)
        self.indir = os.path.join(os.path.dirname(__file__),"test")
        volcanoFile = os.path.join(self.indir, "volcano.xml")
        volc = VolcanoRun(volcanoFile)
        self.dir = volc.outputDir
        os.makedirs(name=self.dir, exist_ok=True)
        yesterday = datetime.datetime.now() - datetime.timedelta(days=1)
        with open(os.path.join(self.dir, "volcano.xml"), "wt") as oh:
            with open(volcanoFile, "rt") as ih:
                for line in ih:
                    oh.write(re.sub('2016-11-03',yesterday.strftime('%Y-%m-%d'), line))
        self.files = ('columnsource_location.csv', 'columnsource_emission.csv', 'eemep_script.job')
        for file in self.files:
            if (os.path.lexists(file)):
                os.unlink(file)
        for f in os.scandir(self.dir):
            if f.is_symlink():
                os.unlink(f.path)


    def testModelRunner(self):
        mr = ModelRunner(self.dir, TestModelRunner.hpcMachine)
        self.assertTrue(len(self.files) == 3)
        #self.assertTrue(os.path.exists(self.files[0]))
        for x in self.files:
            self.assertTrue(os.path.exists(os.path.join(mr.path, x)), "file created: {}".format(x))
        meteo_count = 0
        for x in os.scandir(self.dir):
            if re.search(r'meteo\d{8}.nc', x.path):
                meteo_count += 1
        self.assertEqual(meteo_count, 4, "meteo files created")


    @unittest.skipIf(doRun, "testet in upload and run")
    def test_upload_files(self):
        mr = ModelRunner(self.dir, TestModelRunner.hpcMachine)
        mr.do_upload_files()

    @unittest.skipIf(doRun == False, "no run")
    def test_upload_files_and_run(self):
        mr = ModelRunner(self.dir, TestModelRunner.hpcMachine)
        mr.do_upload_files()
        status = mr.run_and_wait()
        print(status)

    @unittest.skipIf(doRun == False, "no run")
    def test_download(self):
        mr = ModelRunner(self.dir, TestModelRunner.hpcMachine)
        mr.download_results()


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
