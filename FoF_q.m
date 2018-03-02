%function submit = NBP_Q_files(I0,I1,I2,I3,I4)
I1=2
I2=5
I3=5
core='FoF_extrapolate';
% mode is the number of different parameter settings.
% core is the core function to be called.
% Example: mode = 23; core = 'FA_Gibbs'; submit = Q_files(mode,core);
submit =[];
submit1 = [];
%for i0 = 1:1
    %for i1=[1,9]
    for i1=11
        for i2 = 1:5
            for i3 = 1:5
               % for i4 = 1:I4
                   
                        submit = [submit 'sbatch ' core '_' num2str(i1) '_' num2str(i2) '_'  num2str(i3) '.q; '];
                        fid = fopen([core '_' num2str(i1) '_' num2str(i2) '_'  num2str(i3) '.q'],'W');
                        fprintf(fid,'#!/bin/bash\n');
                        fprintf(fid,'#SBATCH -J %s      # job name\n', [core '_' num2str(i1) '_' num2str(i2) '_' num2str(i3) ]);
                        fprintf(fid,'#SBATCH -o %s\n', [core '_' num2str(i1) '_' num2str(i2) '_' num2str(i3) '.txt']);
                        %fprintf(fid,'#SBATCH -o NBP.o%j         # output and error file name (%j expands to jobID)\n');
                        fprintf(fid,'#SBATCH -N 1\n');
                        fprintf(fid,'#SBATCH -n 16              # total number of mpi tasks requested\n');
                        fprintf(fid,'#SBATCH -p normal     # queue (partition) -- normal, development, etc.\n');
                        fprintf(fid,'#SBATCH -t 48:00:00        # run time (hh:mm:ss) - 12 hours.\n');
                        fprintf(fid,'#SBATCH --mail-user=mzhou@utexas.edu\n');
                        fprintf(fid,'#SBATCH --mail-type=begin  # email me when the job starts\n');
                        fprintf(fid,'#SBATCH --mail-type=end    # email me when the job finishes\n');
                        fprintf(fid,'module load matlab\n');
                        fprintf(fid,'matlab  -nodisplay -nosplash -nodesktop -r "Call_%s(%d,%d,%d)"\n', core,i1,i2,i3);
                        fclose(fid);
                
              %  end
            end
        end
    end
%end
disp(submit);

    
    
