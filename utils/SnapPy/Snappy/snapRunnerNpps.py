#! /usr/bin/env python3

import argparse
import datetime
import multiprocessing
import subprocess
from Snappy.Resources import Resources
import os
import sys
from time import gmtime, strftime


def get_isotope_release(res, reldict):
    errors = ""
    for tag in ("releaseTime", "radius", "lowerHeight", "upperHeight"):
        if tag not in reldict:
            errors += "Cannot interpret {}: {}".format(tag, reldict[tag])

    source_tmpl = """
MAX.PARTICLES.PER.RELEASE= 2000
TIME.RELEASE.PROFILE.STEPS
RELEASE.HOUR= 0, {releaseTime}
RELEASE.RADIUS.M= {radius}, {radius}
RELEASE.LOWER.M= {lowerHeight}, {lowerHeight}
RELEASE.UPPER.M= {upperHeight}, {upperHeight}
"""
    source_term = source_tmpl.format(
        releaseTime=reldict["releaseTime"],
        radius=reldict["radius"],
        lowerHeight=reldict["lowerHeight"],
        upperHeight=reldict["upperHeight"],
    )

    isotopes = {"relI131": "I131", "relXE133": "Xe133", "relCS137": "Cs137"}
    for rel, iso in isotopes.items():
        emis = 0.0
        try:
            emis = float(reldict[rel])
        except Exception:
            pass
        if emis > 0.0:
            source_term += "RELEASE.BQ/SEC.COMP= {rel}, {rel}, '{iso}'\n".format(
                rel=reldict[rel], iso=iso
            )

    # add Cs137, I131 and Xe133
    source_term += res.isotopes2snapinput([169, 158, 148])

    return (source_term, errors)


def runsnapsingle(npp, release_dt, tag_time, dirname, res):
    nPPs = res.readNPPs()
    reldict = {
        "releaseTime": 96,
        "radius": 50,
        "lowerHeight": 50,
        "upperHeight": 500,
        "relI131": 1.39e13,
        "relXE133": 1.0e13,
        "relCS137": 2.6e11,
    }
    (term, errors) = get_isotope_release(res, reldict)
    if len(errors) > 0:
        print(errors, file=sys.stderr)
        return False
    lat = float(nPPs[npp]["lat"])
    lon = float(nPPs[npp]["lon"])
    simStart = strftime("%Y-%m-%d_%H:%M:%S", tag_time)
    print(f"outputdir for {simStart}: {dirname}")
    os.mkdir(dirname)
    sourceTerm = f"""
TITLE={npp} {release_dt:%Y %m %d %H}
SIMULATION.START.DATE={simStart}
SET_RELEASE.POS= P=   {lat},   {lon}
TIME.START= {release_dt:%Y %m %d %H}
TIME.RUN = 52h
STEP.HOUR.OUTPUT.FIELDS= 3
"""
    sourceTerm += term
    with open(os.path.join(dirname, "snap.input"), "w") as fh:
        fh.write(sourceTerm)
        files = res.getECMeteorologyFiles(release_dt, 52, "best")
        fh.write(res.getSnapInputMetDefinitions("nrpa_ec_0p1", files))
    with open(f"{dirname}/snap.stdout.log", "wt") as stdout:
        with open(f"{dirname}/snap.stderr.log", "wt") as stderr:
            os.environ["OMP_NUM_THREADS"] = "2"  # moderate thread usage
            try:
                subprocess.run(
                    ["bsnap_naccident", "snap.input"],
                    cwd=dirname,
                    stdout=stdout,
                    stderr=stderr,
                    check=True,
                )
            except subprocess.CalledProcessError as e:
                print(f"{e}\nContinuing despite errors")
            # add time of arrival
            try:
                subprocess.run(
                    ["snapAddToa", "snap.nc"],
                    cwd=dirname,
                    stdout=stdout,
                    stderr=stderr,
                    check=true,
                )
            except subprocess.CalledProcessError as e:
                print(f"\nContinuing despite errors")
    # and plot
    prod_dir = os.path.join(dirname, "prod")
    os.mkdir(prod_dir)
    with open(os.path.join(prod_dir, "diana.setup"), "wt") as fh:
        di_setup = f"""
%include /etc/diana/setup/diana.setup-COMMON
<COLOURS>
seablueOSM=174,207,224
landOSM=242,239,233
</COLOURS>

<PALETTES>
dsa_toa=255:0:255,255:128:255,128:0:128,128:0:255,128:128:255,128:128:192,192:192:192
</PALETTES>

<OBSERVATION_FILES>
PROD=ascii:EuropeISO3
file={res.directory}/europe3.txt
</OBSERVATION_FILES>

<FIELD_FILES>
filegroup=SNAP
m=SNAP.current t=fimex format=netcdf f={dirname}/snap.nc
</FIELD_FILES>
"""
        fh.write(di_setup)

    hours = ",".join(["{x} | {y}".format(x=x, y=x + 3) for x in range(0, 48, 3)])
    component = "Cs137"
    dianaIn = os.path.join(prod_dir, "diana.in")
    dianaInTmpl = ""
    with open(res.getBSnapInputFile(), "rt") as fh:
        for line in fh:
            if line.startswith("AREA"):
                # zoom into area (lon_min:lon_max:lat_min:lat_max in radians)
                dianaInTmpl += 'AREA proj4string="+proj=latlon +R=6371000 +no_defs" rectangle=-0.38:1.0:0.55:1.35\n'
            else:
                dianaInTmpl += line
    with open(dianaIn, "wt") as fh:
        fh.write(dianaInTmpl.format(hours=hours, component=component))
    # env-vars for bdiana without display
    os.environ["QT_QPA_PLATFORMTHEME"] = ""
    os.environ["QT_QPA_FONTDIR"] = "/usr/share/fonts/truetype"
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    try:
        subprocess.run(["bdiana", "-i", "diana.in", "-s", "diana.setup"], cwd=prod_dir, check=True)
    except subprocess.CalledProcessError as e:
        print(f"{e}")
        return False
    return True


def runsnap(npps, outdir):
    release_dt = datetime.datetime.now() + datetime.timedelta(hours=1)

    nppdir = {}
    runparams = []
    res = Resources()
    res.getLustreDir()  # set lustredir within res before parallelization - mapp-services.sh not parallelizable yet
    nPPs = res.readNPPs()
    for i, tag in enumerate(npps):
        if tag not in nPPs:
            print(f"unknown NPP: {tag}", file=sys.stderr)
            sys.exit(1)
        tag_time = gmtime(
            release_dt.timestamp() + i
        )  # ensure a few secs difference for tags
        nppdir[tag] = os.path.join(
            res.getSnapOutputDir(),
            "{0}_{1}".format(tag, strftime("%Y-%m-%dT%H%M%S", tag_time)),
        )
        runparams.append([tag, release_dt, tag_time, nppdir[tag], res])

    with multiprocessing.Pool(4) as p:
        print(p.starmap(runsnapsingle, runparams))

    # montage
    opts = ["montage", "-mode", "concatenate", "-tile", f"4x{len(npps)}"]
    for tag in npps:
        opts += ["-label", tag]
        for t in [12, 24, 36, 48]:
            opts += [f"{nppdir[tag]}/prod/snap_{t}h.png"]
    opts.append(f"{outdir}/plots_{release_dt:%Y-%m-%dT%H}0000Z.png")
    print("running montage with:\n" + " ".join(opts))
    subprocess.run(opts)

    # montage toa
    opts = ["montage", "-mode", "concatenate", "-tile", f"{len(npps)}x1"]
    for tag in npps:
        opts += ["-label", tag]
        opts += [f"{nppdir[tag]}/prod/snap_toa.png"]
    opts.append(f"{outdir}/plots_{release_dt:%Y-%m-%dT%H}0000Z_toa.png")
    print("running montage for toa with:\n" + " ".join(opts))
    subprocess.run(opts)


def main():
    os.umask(0o002)

    parser = argparse.ArgumentParser(
        description="Run snap for several Nuclear Power Plants and receive figures of dispersion",
        usage=f"{sys.argv[0]} NPP [NPP NPP ...]",
    )
    parser.add_argument(
        "--dir", help="output directory", default="/home/VGLSHARE/SNAP4DSA/"
    )
    parser.add_argument(
        "--store",
        help="storeA or storeB, meteo and runtime-datastore, default used from MAPP-system",
    )
    parser.add_argument(
        "NPP",
        help="name of Nuclear Power Plant, e.g. CHERNOBYL (see list in SnapPy)",
        nargs="+",
    )
    args = parser.parse_args()

    if args.store:
        os.environ["STORE"] = args.store
    runsnap(args.NPP, args.dir)


if __name__ == "__main__":
    main()
