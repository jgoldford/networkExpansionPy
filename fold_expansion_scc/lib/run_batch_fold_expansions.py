import glob
import subprocess as sb

def make_qsub(job_name,pythonScript,out_id):
    out_string = f"""#!/bin/bash -l
# Set SCC project
#$ -P bioinfor

# Specify hard time limit for the job. 
#   The job will be aborted if it runs longer than this time.
#   The default time, also selected here, is 12 hours.  You can increase this up to 720:00:00 for single processor jobs but your job will take longer to start.
#$ -l h_rt=04:00:00

# Give job a name
#$ -N {job_name}

# Combine output and error files into a single file
#$ -j y

module load miniconda
conda activate network_expansion

"""

    python_call = 'python ' + pythonScript + ' -s ' + out_id ;
    out_string = out_string + python_call
    return out_string

python_call = '/projectnb2/bioinfor/SEGRE/goldford/network_expansion/networkExpansionPy/fold_expansion_scc/lib/run_fold_expansion.py'
qsub_dir = '/projectnb2/bioinfor/SEGRE/goldford/network_expansion/networkExpansionPy/fold_expansion_scc/qsub'

numJobs = 1000
for i in range(numJobs):
    sid = 's' + str(i)
    job_name = "fold_expansion_" + sid
    qsub_call = make_qsub(job_name,python_call,sid)
    sh_file = qsub_dir + '/foldExpansion.randOrder.' +  sid + '.qsub'
    with open(sh_file,'w') as f:
        f.write(qsub_call)
        
    # call qsub file to submit job
    sb.call('qsub ' + sh_file,shell=True)