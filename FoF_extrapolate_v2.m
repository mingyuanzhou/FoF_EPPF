function [samples,Chi2,RMSE,Mave,Chi2_extrapolate,RMSE_extrapolate,Mave_extrapolate]=FoF_extrapolate_v2(dataset,randtry,ttt,model,IsDisplay)
% %%Demo
%  dataset = 9; %9
%  randtry = 1; %2,3,4,5
%  ttt = 1; %2,3,4,5
%  model =3; %2,3,4,5,6
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
a=0;
p=0.9;

samples = zeros(7, Collections);
Mave = 0;
Mave_extrapolate = 0;
for iter=1:Burnin+Collections
    %iter
    
    if iter==1
        N = length(z);
        ell = length(n_k);
    else
        N = length(z0);
    end
    
    if a==0
        %Sample r
        r = gamrnd(e_0 + ell, 1./(f_0-log(max(1-p,realmin))));
        %Sample p
        p = betarnd(N+1,r+1);
    else        
        %Sample r
        r = gamrnd(e_0+ ell,1/(f_0 + (1-(1-p)^a)/(a*p^a)));
        %Sample p
        p = 0.0001:0.0001:0.9999;
        %p = min(betarnd(0.1,0.1,1,10000),1-eps);
        logprob = -r*(1-(1-p).^a)./(a.*p.^a) + (N-a*ell)*log(p);
        prob = exp(logprob - max(logprob));
        prob(isnan(prob))=0;
        cdf =cumsum(prob);
        p = p(sum(rand(1)*cdf(end)>cdf)+1);
    end
    
    %Sample a
    if model==2
        a=0;
    elseif model==3
        a=0.5;
    elseif model==4
        a=-1;
    else
        if model==1
            %a<1
            at = 0.0001:0.0001:0.9999;
            a = 2-1./at;
        elseif model==5
            %0<=a<1
            at = 0:0.0001:0.9999;
            a = at;
        elseif model==6
            %a<0
            at = 0.0001:0.0001:0.9999;
            a = 1-1./at;
        end
        temp = (1-(1-p).^a)./(a.*p.^a);
        temp(a==0) = -log(1-p);
        logprob = -r*temp -a.*ell*log(p) + ...
            sum(gammaln(bsxfun(@minus,n_k,a)),1)- length(n_k)*gammaln(1-a);
        prob = exp(logprob - max(logprob));
        prob(isnan(prob))=0;
        cdf =cumsum(prob);
        a = a(sum(rand(1)*cdf(end)>cdf)+1);
    end
    %[iter,r/1000,a,p,e_0,f_0]

    if a==0
        IIII=5;
    else
        IIII=5;
    end

    if iter==1
        zall = zeros(1,length(z0));
        zall(1:length(z))=z;
        n_k=full(sparse(z,1,1));
        SampleDex=find(zall==0);
    end
    for iiii=1:IIII
        if 1
            [zall,n_k] = gCRSF_gibbs(zall,n_k,SampleDex(randperm(length(SampleDex))),r,a,p);
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
                prob = [max(n_k-a,0);r*p^(-a)];
                zall(iz)= multrnd_unnormalized(prob);
                if zall(iz)==length(prob);
                    n_k(end+1)=0;
                    %eta(end+1,:)=gamrnd(0.05,1);
                end
                n_k(zall(iz))=n_k(zall(iz))+1;
            end
        end
        ell = length(n_k);
    end
    
    
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
        %  [r/100,a,p,RMSE]
        %  plot(samples(5,Burnin+1:iter)'); drawnow
        %subplot(1,2,1);plot(1:iter, samples(4,1:iter)'-S0, 1:iter, samples(5,1:iter)'-S0);
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
end
save(['FoFresults/FoFv2_',name_dataset,'_',num2str(ttt),num2str(model),num2str(randtry),'.mat'],'samples','Chi2','RMSE','Mave','Chi2_extrapolate','RMSE_extrapolate','Mave_extrapolate');