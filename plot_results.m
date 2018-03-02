Chi2_all=zeros(5,2);
[z0,n_k,M,name_dataset]=read_data(11);
for i=1:5
    figure;
    load (['FoF_PY_microdata_78' num2str(i) '.mat'])
    Chi2_all(i,1) = Chi2;
    loglog(Mave/500,'b');hold on
    load (['FoF_microdata_71' num2str(i) '.mat'])
    Chi2_all(i,2) = Chi2;
    loglog(Mave/500,'r');hold on;loglog(M,'o')
end

Chi2_all=zeros(8,2);
[z0,n_k,M,name_dataset]=read_data(11);
for t=1:5
    figure;
    %load (['FoF1_PY_microdata_' num2str(t) '81.mat'])
    load (['FoF1_microdata_' num2str(t) '11.mat'])
    Chi2_all(t,1) = Chi2;
    loglog(Mave/500,'b');hold on
    load (['FoF_NGG_microdata_' num2str(t) '51.mat'])
   % load (['FoF1_microdata_' num2str(t) '11.mat'])
    Chi2_all(t,2) = Chi2;
    loglog(Mave/500,'r');hold on;loglog(M,'o')
end



Chi2_all=zeros(5,2);
[z0,n_k,M,name_dataset]=read_data(1);
for i=1:5
    figure;
    load (['FoF_PY_text_tomsawyer_88' num2str(i) '.mat'])
    Chi2_all(i,1) = Chi2;
    loglog(Mave/500,'b');hold on
    load (['FoF_text_tomsawyer_81' num2str(i) '.mat'])
    Chi2_all(i,2) = Chi2;
    loglog(Mave/500,'r');hold on;loglog(M,'o')
end


Chi2_all=zeros(8,2);
[z0,n_k,M,name_dataset]=read_data(1);
for t=1:5
    figure;
    load (['FoF1_PY_text_tomsawyer_' num2str(t) '81.mat'])
    Chi2_all(t,1) = Chi2;
    loglog(Mave/500,'b');hold on
    load (['FoF_text_tomsawyer_' num2str(t) '11.mat'])
    Chi2_all(t,2) = Chi2;
    loglog(Mave/500,'r');hold on;loglog(M,'o')
end


[z0,n_k,M,name_dataset]=read_data(9);

Chi2_all=zeros(5,2);
for i=1:5
    figure;
    load (['FoF_PY_gene_sultan_68' num2str(i) '.mat'])
    Chi2_all(i,1) = Chi2;
    loglog(Mave/500,'b');hold on
    load (['FoF_gene_sultan_61' num2str(i) '.mat'])
    Chi2_all(i,2) = Chi2;
    loglog(Mave/500,'r');hold on;loglog(M,'o')
end