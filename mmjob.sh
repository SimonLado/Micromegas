#!/bin/sh -f
#
#          Job name (default is name of pbs script file)
#PBS -N micromegas
#
#          Resource limits: max. wall clock time during which job can be running
#PBS -l select=1:vnode=^wn-01-01-06.cluster.roma3
#PBS -l walltime=17:00:00
#
#          The directive below directs that the standard output and
#          error streams are to be merged, intermixed, as standard
#          output.
#PBS -j oe
#
##########################################

echo ------------------------------------------------------
echo 'Job is running on node(s): '
cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: executing queue is $PBS_QUEUE
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

#WORKDIR=/storage/DATA-01/ATLAS/LOCAL/Users/dinardo/RHUM_SIM/Amber_sim/build/jobs/
WORKDIR=/storage/scratch/ladogana/results

cd ${WORKDIR}

# COMMAND to EXECUTE:

#echo "$PBS_JOBID" > ${PBS_JOBID}.txt

ls  /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase

#export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
#alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
#setupATLAS

source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
lsetup "views LCG_104 x86_64-centos7-gcc11-opt"
source /storage/DATA-01/ATLAS/LOCAL/Users/dinardo/RHUM_SIM/garf_install/share/Garfield/setupGarfield.sh

#/storage/DATA-01/ATLAS/LOCAL/Users/dinardo/RHUM_SIM/Amber_sim/build/./MM ${var1} ${var2} ${var3}

/storage/DATA-01/ATLAS/LOCAL/Users/dinardo/RHUM_SIM/tesi_ladogana/build/./MM ${var1} ${var2} ${var3}

#cp mm.*${var1}*.root
#/storage/DATA-01/ATLAS/LOCAL/Users/dinardo/RHUM_SIM/Amber_sim/build/jobs/
#cp ${PBS_JOBID}.txt
#/storage/DATA-01/ATLAS/LOCAL/Users/dinardo/RHUM_SIM/Amber_sim/build/jobs/

exit
# End of job script