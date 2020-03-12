#! /usr/bin/python -u
#SBATCH -q xfer
#SBATCH -o %x.%j
#SBATCH -t 48:00:00

from __future__ import print_function
"""
A script to backup / restore FastPM snapshots to / from HPSS.

Assumes the system has hsi and htar utilities.

The script can be directly submitted to the Slurm queue.

Currently each dataset (for particle type) is combined into a single tar file.
Each leaf directory is combined into a single tar file.
Other files are directly copied as files to the HPSS.

"""

from argparse import ArgumentParser

def print(*args, **kwargs):
    from __builtin__ import print
    return print('hpss-backup: ', *args, **kwargs)

ap = ArgumentParser()
ap.add_argument("-r", "--restore", action='store_true', default=False,
            help="Reverse the operation; recover the files instead of backup.")
ap.add_argument("-g", "--debug", action='store_true', default=False, help="print debugging information.")
ap.add_argument("--dry-run", action='store_true', default=False, help="Only print the commands.")

ap.add_argument("-p", "--pattern", type=str, action='append',
                default=[], help="only include files and tars that fnmatch patterns in this list. "
                                 "Relative to src. '/' has no special meaning.")
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
from subprocess import check_call, check_output, CalledProcessError
import subprocess

def print_call(args, cwd, stderr):
   print("( cd %(cwd)s; %(cmd)s )" % dict(cwd=cwd, cmd=" ".join(args)))

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


def hpathnames(files):
    if isinstance(files, (list, tuple)):
        return [hpathnames(f) for f in files]
    else:
        return '%s : %s' % (files, files)

def hput(workdir, targetdir, files, verbose=False, dry_run=False):
    if dry_run:
        run = print_call
    elif verbose:
        run = check_call
    else:
        run = check_output

    return run(["hsi" , "-q", "mkdir -p %(hpss)s; cd %(hpss)s; put -P %(files)s" % dict(
        hpss=targetdir,
        files=' '.join(hpathnames(files)))
        ], cwd=workdir, stderr=subprocess.STDOUT)

def hget(workdir, targetdir, files, verbose=False, dry_run=False):
    if dry_run:
        run = print_call
    elif verbose:
        run = check_call
    else:
        run = check_output

    if len(files) == 0:
        return

    return run(["hsi" , "-q", "cd %(hpss)s; get %(files)s" % dict(
        hpss=targetdir,
        files=' '.join(hpathnames(files)))
        ], cwd=workdir, stderr=subprocess.STDOUT)

def hexists(workdir, file, verbose=False):
    if verbose:
        run = check_call
    else:
        run = check_output

    try:
        run(["hsi" , "-q", "ls %(file)s" % dict(
            file=file)
            ], cwd=workdir, stderr=subprocess.STDOUT)
    except CalledProcessError as e:
        return False

    return True

def htar(mode, workdir, target, src, verbose=False, dry_run=False):
    """ mode can be c, t, x """
    if dry_run:
        run = print_call
    elif verbose:
        run = check_call
    else:
        run = check_output

    if verbose:
        verbose = "v"
    else:
        verbose = ""

    return run(["htar", "-q", "-P", "-%s%sf" % (mode, verbose), target, src], cwd=workdir, stderr=subprocess.STDOUT)

def match_filename(ns, filename):
    if len(ns.pattern) == 0:
        return True
    for pattern in ns.pattern:
        if fnmatch(filename, pattern):
            return True
    return False

def backup(ns):
    print("FastPM snapshots on local systems at:")
    print(ns.src)
    print("HPSS target location:")
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

        if "backup.files" in extrafiles:
            extrafiles.remove("backup.files")
        if "backup.tars" in extrafiles:
            extrafiles.remove("backup.tars")

        print("saving meta data for recovery ...")

        with open(os.path.join(ns.src, "backup.files"), 'w') as ff:
            for file in extrafiles:
                ff.write("%s\n" % file)

        with open(os.path.join(ns.src, "backup.tars"), 'w') as ff:
            for file in tarables:
                ff.write("%s\n" % file)

        hput(ns.src, ns.dest, ["backup.files", "backup.tars"], verbose=True)

    print("Found %d extra files" % len(extrafiles))
    print("Found %d tarable directories" % len(tarables))

    if ns.file_end is None or ns.file_end > len(extrafiles):
        ns.file_end = len(extrafiles)
    if ns.tar_end is None or ns.tar_end > len(tarables):
        ns.tar_end = len(tarables)

    chunksize = 32
    for i in range(ns.file_start, ns.file_end, chunksize):
        i1 = i + chunksize
        if i1 > ns.file_end:
            i1 = ns.file_end
        matched = [fn for fn in extrafiles[i:i1] if match_filename(ns, fn)]

        if ns.restore:
            print("Retriving extra file (%d / %d) ..." %(i, len(extrafiles)))
            print(" ".join(matched))
            hget(ns.src, ns.dest, matched, verbose=ns.debug, dry_run=ns.dry_run)
        else:
            print("Storing extra files (%d / %d) ..." %( i, len(extrafiles)))
            print(" ".join(matched))
            hput(ns.src, ns.dest, matched, verbose=ns.debug, dry_run=ns.dry_run)

    for i in range(ns.tar_start, ns.tar_end):
        tarable = tarables[i]
        if not match_filename(ns, tarable):
            print("Skipped tar for dataset (%d / %d) %s : " % (i, len(tarables), tarable))
            continue

        if ns.restore:
            print("Retreving tar for dataset (%d / %d) %s : " % (i, len(tarables), tarable))
            htar('x', ns.src, os.path.join(ns.dest, "%s.tar" % tarable), tarable,
                verbose=ns.debug, dry_run=ns.dry_run)
        else:
            print("Creating tar for dataset (%d / %d) %s : " % (i, len(tarables), tarable))

            if hexists(ns.src, os.path.join(ns.dest, "%s.tar.idx" % tarable)):
                print("File %s already exists with the following contenxt. Delete it if you want to overwrite it." % os.path.join(ns.dest, "%s.tar" % tarable))
                # use verbose to obtain list of files on stdout
                htar('t', ns.src, os.path.join(ns.dest, "%s.tar" % tarable), tarable,
                    verbose=True, dry_run=ns.dry_run)
            else:
                htar('c', ns.src, os.path.join(ns.dest, "%s.tar" % tarable), tarable,
                    verbose=ns.debug, dry_run=ns.dry_run)

if __name__ == "__main__":
    ns = ap.parse_args()

    backup(ns)

