
function filename=jobs_Lonestar5(core,TaskNum,RunTime,TaskPerNode)
if nargin<3
    RunTime=48;
end
if nargin<4
    TaskPerNode = 10;
end
filename = [core,'.LS5'];
fid = fopen(filename,'W');



fprintf(fid,'#!/bin/bash\n');
fprintf(fid,'#\n');
fprintf(fid,'# Simple SLURM script for submitting multiple serial\n');
fprintf(fid,'# jobs (e.g. parametric studies) using a script wrapper\n');
fprintf(fid,'# to launch the jobs.\n');
fprintf(fid,'#\n');
fprintf(fid,'# To use, build the launcher executable and your\n');
fprintf(fid,'# serial application(s) and place them in your WORKDIR\n');
fprintf(fid,'# directory.  Then, edit the CONTROL_FILE to specify\n');
fprintf(fid,'# each executable per process.\n');
fprintf(fid,'#-------------------------------------------------------\n');
fprintf(fid,'#-------------------------------------------------------\n');
fprintf(fid,'#\n'); 




fprintf(fid,'# \n');
fprintf(fid,'#------------------Scheduler Options--------------------\n');
fprintf(fid,'#SBATCH -J %s          # Job name\n',core);
%fprintf(fid,'#SBATCH -N %d                   # Total number of nodes (16 cores/node)\n',ceil(TaskNum/TaskPerNode));
fprintf(fid,'#SBATCH -n %d                  # Total number of tasks\n',TaskNum);
fprintf(fid,'#SBATCH -p normal          # Queue name\n');
fprintf(fid,'#SBATCH -o %s.o%%j       # Name of stdout output file (%%j expands to jobid) \n',core);
fprintf(fid,'#SBATCH -t %d:00:00            # Run time (hh:mm:ss)\n',RunTime);
fprintf(fid,'#SBATCH --mail-user=mzhou@utexas.edu\n');
fprintf(fid,'#SBATCH --mail-type=begin  # email me when the job starts\n');
fprintf(fid,'#SBATCH --mail-type=end    # email me when the job finishes\n');
fprintf(fid,'#      <------------ Account String ------------>\n');
fprintf(fid,'# <--- (Use this ONLY if you have MULTIPLE accounts) --->\n');
fprintf(fid,'#SBATCH -A ParNBP\n'); 

fprintf(fid,'module load launcher\n');
%fprintf(fid,'module load matlab/2014b\n');
fprintf(fid,'module load matlab\n');
fprintf(fid,'export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins\n');
fprintf(fid,'export LAUNCHER_RMI=SLURM\n');
%fprintf(fid,'export EXECUTABLE=init_launcher\n');
%fprintf(fid,'export WORKDIR=.\n');
fprintf(fid,'export LAUNCHER_JOB_FILE=%s.q\n',core);
fprintf(fid,'\n');
fprintf(fid,'$LAUNCHER_DIR/paramrun\n');
fclose(fid);