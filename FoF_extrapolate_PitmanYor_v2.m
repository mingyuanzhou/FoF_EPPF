function [samples,Chi2,RMSE,Mave,Chi2_extrapolate,RMSE_extrapolate,Mave_extrapolate]=FoF_extrapolate_PitmanYor_v2(dataset,randtry,ttt,model,IsDisplay)
%%Demo
% dataset = 11; %9
% randtry = 1; %2,3,4,5
% ttt = 1; %2,3,4,5
% model =8; %2,3,4,5,6
%%
if nargin<5
    IsDisplay = true;
end
Burnin = 500;
Collections = 500;
Ratios=[1/32,1/16,1/8,1/4,1/2];
SampleRatio = Ratios(ttt);

[z0,n_k,M,name_dataset]=read_data(dataset);
dexM = find(M(1:100)>0);

figure %(model)
rng(randtry,'twister')
[z,idx] = datasample(z0,round(length(z0)*SampleRatio),'Replace',false);

z_extrapolate = z0;
z_extrapolate(idx)=[];
[~,~,z_extrapolate]=unique(z_extrapolate);
n_k_extrapolate=full(sparse(z_extrapolate,1,1));
M_extrapolate = nk_to_m(n_k_extrapolate);
dexM_extrapolate = find(M_extrapolate(1:100)>0);

n_k0=sparse(z0,1,1);
[~,~,z]=unique(z);
n_k=full(sparse(z,1,1));

e_0=1e-2; f_0=1e-2;
a=0.5;
p=0.9;

samples = zeros(7, Collections);
Mave = 0;
Mave_extrapolate = 0;
for iter=1:Burnin+Collections
    %iter
    
    if iter==1
        N = length(z);
        ell = length(n_k);
        ellz = ell;
        Nz = N;
        r=N;
        nz_k=n_k;
    end
    
  
    if model==7
        %Chinese restaurant process
        a=0;
        %Sample p
        p = betarnd(Nz-1,r+1);
        %Sample r
        r = gamrnd(e_0 + ellz-1, 1./(f_0-log1p(-p)));
    elseif model==8
        %Pitman-Yor process
        
%         %Sample p
%         p = betarnd(Nz,r);
%         %Sample r
%         sumy = sum(double(rand(1,ellz)<r./(r+a*(1:ellz))));
%         r = gamrnd(e_0+ sumy,1/(f_0 -log1p(-p)));
%         %Sample a
%         sumb = 0;
%         for k=1:ellz
%             sumb = sumb+sum(double(rand(1,nz_k(k))<(0:nz_k(k)-1)./((1:nz_k(k))-a)));
%         end
%         a = betarnd(1e-0+ellz-sumy,1e-0+Nz-sumb);
        
        
         %Sample p
        p = betarnd(Nz-1,r+1);
        %Sample r
        sumy = sum(double(rand(1,ellz-1)<r./(r+a*(1:ellz-1))));
        r = gamrnd(e_0+ sumy,1/(f_0 -log1p(-p)));
        %Sample a
        sumb = 0;
        for k=1:ellz
            if nz_k(k)>=2
                sumb = sumb+sum(1-double(rand(1,nz_k(k)-1)<(0:nz_k(k)-2)./((1:nz_k(k)-1)-a)));
            end
        end
        a = betarnd(1e-0+ellz-1-sumy,1e-0+sumb);
    elseif model==9
        a = 0.5;
        p = betarnd(Nz-1,r+1);
        %Sample r
        sumy = sum(double(rand(1,ellz-1)<r./(r+a*(1:ellz-1))));
        r = gamrnd(e_0+ sumy,1/(f_0 -log1p(-p)));
    end
    %[iter,r/1000,a,p,e_0,f_0]
    

    if iter==1
        zall = zeros(1,length(z0));
        zall(1:length(z))=z;
        n_k=full(sparse(z,1,1));
        SampleDex=find(zall==0);
    end

    if 1
        [zall,n_k] = PitmanYor_gibbs(zall,n_k,SampleDex(randperm(length(SampleDex))),r,a);
        [kk,kki,kkj] = unique(zall);
        zall=kkj;
        n_k=full(sparse(zall,1,1));
    else
        for iz=length(z)+1:length(z0)
            if zall(iz)>0
                n_k(zall(iz))=n_k(zall(iz))-1;
            end
            dex=find(n_k==0);
            if ~isempty(dex)
                zall(zall>dex)=zall(zall>dex)-1;
                n_k(dex)=[];
            end
            prob = [max(n_k-a,0);r+a*nnz(n_k)];
            zall(iz)= multrnd_unnormalized(prob);
            if zall(iz)==length(prob);
                n_k(end+1)=0;
                %eta(end+1,:)=gamrnd(0.05,1);
            end
            n_k(zall(iz))=n_k(zall(iz))+1;
        end
    end
    ell = length(n_k);

    
    
    if iter>Burnin %&& mod(iter,5)==0
        

        z_extrapolate = zall(SampleDex);
        [~,~,z_extrapolate]=unique(z_extrapolate);
        n_k_extrapolate=full(sparse(z_extrapolate,1,1));
        
        M1 = nk_to_m(n_k_extrapolate);
        len=max(length(M1), length(Mave_extrapolate));
        temp1 = zeros(1,len);
        temp2 = zeros(1,len);
        temp1(1:length(M1))=M1;
        temp2(1:length(Mave_extrapolate))=Mave_extrapolate;
        Mave_extrapolate = temp1+temp2;

        M1 = nk_to_m(n_k);
        len=max(length(M1), length(Mave));
        temp1 = zeros(1,len);
        temp2 = zeros(1,len);
        temp1(1:length(M1))=M1;
        temp2(1:length(Mave))=Mave;
        Mave = temp1+temp2;
        
        
        RMSE_extrapolate = log(M_extrapolate(dexM_extrapolate))-log(Mave_extrapolate(dexM_extrapolate)/(iter-Burnin));
        RMSE_extrapolate = sqrt(sum(RMSE_extrapolate.^2)/length(dexM_extrapolate));
        Chi2_extrapolate = chi2_rate(M_extrapolate, Mave_extrapolate/(iter-Burnin),50);
        
        RMSE = log(M(dexM))-log(Mave(dexM)/(iter-Burnin));
        RMSE = sqrt(sum(RMSE.^2)/length(dexM));
        Chi2 = chi2_rate(M, Mave/(iter-Burnin),50);
        
        
        samples(:,iter)=[r,a,p,Chi2,RMSE,Chi2_extrapolate,RMSE_extrapolate]';
        if IsDisplay
            subplot(3,2,1);plot(Burnin+1:iter, samples(5,Burnin+1:iter)');title(num2str(RMSE))
            subplot(3,2,2); loglog((1:length(M)),M,'o',(1:length(Mave)),max(Mave/(iter-Burnin),0.1));title(num2str(Chi2))
            subplot(3,2,3);plot(Burnin+1:iter, samples(7,Burnin+1:iter)');title(num2str(RMSE_extrapolate))
            subplot(3,2,4); loglog((1:length(M_extrapolate)),M_extrapolate,'o',(1:length(Mave_extrapolate)),max(Mave_extrapolate/(iter-Burnin),0.1));title(num2str(Chi2_extrapolate))
            subplot(3,2,5);plot(samples(2,1:iter));title(num2str(a))
            drawnow
        end
        %[Chi2,RMSE,Chi2_Pois,RMSE_Pois]
        %sum(log(poisspdf(M,PoissonRate/(iter-Burnin))))
    else
        samples(:,iter)=[r,a,p,0,0,0,0]';
    end
    %a
end
save(['FoFresults/FoFPYv2_',name_dataset,'_',num2str(ttt),num2str(model),num2str(randtry),'.mat'],'samples','Chi2','RMSE','Mave','Chi2_extrapolate','RMSE_extrapolate','Mave_extrapolate');