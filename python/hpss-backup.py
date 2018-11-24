# !/usr/bin/python
#SBATCH -q xfer
#SBATCH -J hpss-xfer
#SBATCH -o hpss-xfer.%j
#SBATCH -n 1
#SBATCH -t 24:00:00

"""
A script to backup / restore FastPM snapshots to / from HPSS.

Assumes the system has hsi and htar utilities.

The script can be directly submitted to the Slurm queue.

Currently each dataset (for particle type) is combined into a single tar file.
Each leaf directory is combined into a single tar file.
Other files are directly copied as files to the HPSS.

"""

from argparse import ArgumentParser

ap = ArgumentParser()
ap.add_argument("-r", "--restore", action='store_true', default=False,
            help="Reverse the operation; recover the files instead of backup.")
ap.add_argument("-g", "--debug", action='store_true', default=False, help="print debugging information.")

# use these options to recover a failed / partial operation.
ap.add_argument("-f", "--file-start", type=int, default=0, help="start the retrieval of files from this")
ap.add_argument("-F", "--file-end", type=int, default=None,
            help="end the retrieval of files before this; set to equal to -f to skip all files")
ap.add_argument("-t", "--tar-start", type=int, default=0, help="start the retrieval of tars from this")
ap.add_argument("-T", "--tar-end", type=int, default=None,
            help="end the retrieval of files before this; set to equal to -t to skip all tar balls")

ap.add_argument("src", help="location on local file system that contains multiple fastpm snapshots; ")
ap.add_argument("dest", help="location on hpss system that will receive the backups;")

import os

from fnmatch import fnmatch
from subprocess import check_call, check_output
import subprocess

def find_tarables(rootdir, relpath=True, leaflist=[]):
    """ discover directories and files that can be htar or hsi put.
        
        A directory will be htarred if it is a leaf directory, or if
        it matches any of the patterns in leaflist.

        Any files that is not included in a tar file is will be copied with hsi put

        relpath : if True, return relative to rootdir
        leaflist :
        returns

            tarables, extrafiles

        tarables is a list of directories ready to be tarred.

        extrafiles is a list of files ready to be copied.

    """

    tarables = []
    extrafiles = []

    for root, dirs, files in os.walk(rootdir, topdown=True):
        if len(dirs) == 0: # leaf
            tarables.append(root)
        else:
            # need to include the files
            extrafiles.extend(
                [os.path.join(root, file) for file in files])

            newdirs = [ ]

            for dir in dirs:
                if not any(fnmatch(dir, glob) for glob in leaflist):
                    newdirs.append(dir)
                else:
                    tarables.append(os.path.join(root, dir))

            while len(dirs):
                dirs.pop() 

            dirs.extend(newdirs)
    if relpath:
        tarables = [os.path.relpath(t, rootdir) for t in tarables]
        extrafiles = [os.path.relpath(t, rootdir) for t in extrafiles]

    return sorted(tarables), sorted(extrafiles)


def hput(workdir, targetdir, files, verbose=False):
    if verbose:
        run = check_call
    else:
        run = check_output

    return run(["hsi" , "-q", "mkdir -p %(hpss)s; cd %(hpss)s; put -P %(files)s" % dict(
        hpss=targetdir,
        files=' '.join(files))
        ], cwd=workdir, stderr=subprocess.STDOUT)

def hget(workdir, targetdir, files, verbose=False):
    if verbose:
        run = check_call
    else:
        run = check_output

    return run(["hsi" , "-q", "cd %(hpss)s; get %(files)s" % dict(
        hpss=targetdir,
        files=' '.join(files))
        ], cwd=workdir, stderr=subprocess.STDOUT)

def htar(workdir, target, src, verbose=False):
    if verbose:
        run = check_call
    else:
        run = check_output
    return run(["htar", "-q", "-P", "-cf" if not verbose else "-cvf", target, src], cwd=workdir, stderr=subprocess.STDOUT)

def huntar(workdir, target, src, verbose=False):
    if verbose:
        run = check_call
    else:
        run = check_output
    return run(["htar", "-q", "-xf" if not verbose else "-xvf", target, src], cwd=workdir, stderr=subprocess.STDOUT)

def backup(ns):
    print("# FastPM snapshots on local systems at:")
    print(ns.src)
    print("# HPSS target location:")
    print(ns.dest )

    if ns.restore:
        try:
            os.makedirs(ns.src)
        except OSError:
            pass

        hget(ns.src, ns.dest, ["backup.files", "backup.tars"], verbose=ns.debug)

        with open(os.path.join(ns.src, "backup.files"), 'r') as ff:
            extrafiles = [s.rstrip('\n') for s in ff.readlines()]

        with open(os.path.join(ns.src, "backup.tars"), 'r') as ff:
            tarables = [s.rstrip('\n') for s in ff.readlines()]

    else:
        tarables, extrafiles = find_tarables(ns.src,
                relpath=True,
                leaflist=['[0-9]', 'LL-*'])

        extrafiles.remove("backup.files")
        extrafiles.remove("backup.tars")

        print("# saving meta data for recovery:")

        with open(os.path.join(ns.src, "backup.files"), 'w') as ff:
            for file in extrafiles:
                ff.write("%s\n" % file)

        with open(os.path.join(ns.src, "backup.tars"), 'w') as ff:
            for file in tarables:
                ff.write("%s\n" % file)

        hput(ns.src, ns.dest, ["backup.files", "backup.tars"], verbose=ns.debug)

    print("Found %d extra files" % len(extrafiles))
    print("Found %d tarable directories" % len(tarables))

    if ns.file_end is None or ns.file_end > len(extrafiles):
        ns.file_end = len(extrafiles)
    if ns.tar_end is None or ns.tar_end > len(tarables):
        ns.tar_end = len(tarables)

    chunksize=32
    for i in range(ns.file_start, ns.file_end, chunksize):
        i1 = i + chunksize
        if i1 > ns.file_end:
            i1 = ns.file_end
        if ns.restore:
            print("Retriving extra files (%d / %d) ..." %( i, len(extrafiles)))
            hget(ns.src, ns.dest, extrafiles[i:i1], verbose=ns.debug)
        else:
            print("Storing extra files (%d / %d) ..." %( i, len(extrafiles)))
            hput(ns.src, ns.dest, extrafiles[i:i1], verbose=ns.debug)

    for i in range(ns.tar_start, ns.tar_end):
        tarable = tarables[i]
        if ns.restore:
            print("Retreving tar for dataset (%d / %d) %s : " % (i, len(tarables), tarable))

            huntar(ns.src, os.path.join(ns.dest, "%s.tar" % tarable), tarable, verbose=ns.debug)
        else:
            print("Creating tar for dataset (%d / %d) %s : " % (i, len(tarables), tarable))

            htar(ns.src, os.path.join(ns.dest, "%s.tar" % tarable), tarable, verbose=ns.debug)

if __name__ == "__main__":
    ns = ap.parse_args()

    backup(ns)

#print('\n'.join(tarables))
#print('\n'.join(extrafiles))

