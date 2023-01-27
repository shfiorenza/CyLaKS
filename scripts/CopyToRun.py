import subprocess
import os
import argparse

# (Shamelessly stolen from Wen Yan's AMSOS code)

# this script is supposed to be called from the base folder of CyLaKS:
#   python3 ./scripts/CopyTuRun.py /path/to/run/folder


def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


parser = argparse.ArgumentParser()
parser.add_argument('path', type=dir_path)

args = parser.parse_args()

exePath = args.path
print('copy executables to ', exePath)

os.system('rsync -avP cylaks_latest.sif '+exePath+'/')
os.system('rsync -avP params/kif4a.yaml '+exePath+'/')
# os.system('rsync -avP analysis_/get_motor_stats.m '+exePath+'/')
# os.system('rsync -avP gradOpt_lattice_coop.py '+exePath+'/')
os.system('rsync -avP SubmitJob.slurm '+exePath+'/')

# this script then creates a txt file to the executable folder:
# commit.txt
# which is a record of the git commit hash string where
# sim is compiled from

process = subprocess.Popen(
    ['git', 'rev-parse', 'HEAD'], shell=False, stdout=subprocess.PIPE)
git_current_hash = process.communicate()[0].strip()
print('current git commit: ', str(git_current_hash))
os.system('echo '+str(git_current_hash)+' > '+exePath+'/commit.txt')
