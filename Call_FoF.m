datasets = {'text_tomsawyer','text_HOLMES','text_tale','text_WarAndPeace', 'text_MOBY_DICK','text_pride','text_HuckleberryFinn','gene_core','gene_sultan','gene_yang','microdata'};
submit=[];
submit1=[];
%for dataset=[1,2,3,4,5,6,7,8,9,11]
for dataset=[1,8,9,11,2,3,4,5,6,7]
    
    corejob=['FoF',datasets{dataset}];
    fid = fopen([corejob,'.q'],'W');
    TaskNum=0;
    RunTime = 1;
    for randtry=1 %1:5
        for ttt=1 %1:5
            for model = 1:6
                TaskNum = TaskNum+1;
                fprintf(fid,'matlab  -nodisplay -nosplash -nodesktop -r "FoF_extrapolate_v2(%d,%d,%d,%d,false);" -logfile FoF_%d_%d_%d_%d.txt\n', dataset,randtry,ttt,model, dataset,randtry,ttt,model);
            end
            model = 8;
            TaskNum = TaskNum+1;
            fprintf(fid,'matlab  -nodisplay -nosplash -nodesktop -r "FoF_extrapolate_PitmanYor_v2(%d,%d,%d,%d,false);" -logfile FoF_%d_%d_%d_%d.txt\n', dataset,randtry,ttt,model, dataset,randtry,ttt,model);
        end
    end
    
    fclose(fid);
    TaskPerNode = 16;
    %filename=jobs_Stampede(corejob,TaskNum,RunTime,TaskPerNode);
    TaskPerNode = 24;
    filename1=jobs_Lonestar5(corejob,TaskNum,RunTime,TaskPerNode);
    %submit=[submit,'sbatch ',filename,'; '];
    submit1=[submit1,'sbatch ',filename1,'; '];
end
