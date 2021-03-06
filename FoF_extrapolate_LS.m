function [Chi2,RMSE,Mave,Chi2_1,RMSE1,Mave1]=FoF_extrapolate_LS(dataset,randtry,ttt)
%%Demo
% dataset = 11; %9
% randtry = 1; %2,3,4,5
% ttt = 1; %2,3,4,5
% model =8; %2,3,4,5,6
%%
Burnin = 500;
Collections = 500;
Ratios=[1/32,1/16,1/8,1/4,1/2];
SampleRatio = Ratios(ttt);

[z0,n_k,M,name_dataset]=read_data(dataset);
dexM = find(M(1:100)>0);

%figure %(model)
rng(randtry,'twister')
z = datasample(z0,round(length(z0)*SampleRatio),'Replace',false);

n_k0=sparse(z0,1,1);
[~,~,z]=unique(z);
n_k=full(sparse(z,1,1));
S0_sample= 1- sum(n_k.*(n_k-1)/(length(z)*(length(z)-1)));

e_0=1e-2; f_0=1e-2;
a=0.5;
p=0.9;

samples = zeros(5, Collections);
PoissonRate=0;
Mave = 0;

M_subsample = nk_to_m(n_k);
[alpha, xmin, L]=plfit(n_k);


%Kratio = length(n_k0)/length(n_k);
%plplot_FoF(n_k0,xmin,alpha,Kratio);
Slope = -alpha;
%Inter  = mean(log(M_subsample(M_subsample>=10)))- Slope*mean(log(1:length(M_subsample(M_subsample>=10))));
idxM=find(M>=3);
%idxM(idxM>100)=[];
Inter  = mean(log(M(idxM)))- Slope*mean(log(idxM));

Mave=zeros(1,length(z0));
Mave(xmin:end) = exp(Inter + Slope*log(xmin:length(Mave)));

% pmf = (xmin:length(z0)).^(-alpha)/(zeta(alpha)- sum((1:xmin-1).^-alpha));
% %Mave(xmin:end) = sum(M(xmin:end))*pmf;
% Mave(xmin:end) = sum(length(z0)/length(z)*M_subsample(xmin:end))*pmf;
% %Mave(1:xmin-1)=M_subsample(1:xmin-1)*length(z0)/length(z);

if xmin>2
    idx = find(M_subsample>=1);
    %mdl = fitlm(log(1:length(M_subsample(M_subsample<xmin))),log(M_subsample(M_subsample<xmin)));
    mdl = fitlm(log(idx(idx<xmin)),log(M_subsample(idx(idx<xmin))));
    %mdl = fitlm(log(1:length(M_subsample(M_subsample>=10))),log(M_subsample(M_subsample>=10)));
    Slope = mdl.Coefficients.Estimate(2);
    Inter = mdl.Coefficients.Estimate(1);
    
    idx = find(M>=1);
   Inter = mean(log(M(idx<xmin)))- Slope*mean(log(idx(idx<xmin)));
    %Inter = Inter + log(length(n_k0)/length(n_k));
   % Inter = Inter + log(length(z0)/length(z));
    Mave(1:xmin-1) = exp(Inter + Slope*log(1:xmin-1));
else
    Mave(1:xmin-1) = M(1:xmin-1);
end
%Mave(1:xmin-1)=M(1:xmin-1);
%Mave = exp(Inter + log(length(z0)/length(z)) + Slope*xx);

RMSE = log(M(dexM))-log(Mave(dexM));
RMSE = sqrt(sum(RMSE.^2)/length(dexM));

Chi2 = chi2_rate(M, Mave,50);

idx = find(M_subsample>=10);
mdl = fitlm(log(idx),log(M_subsample(idx)));
%mdl = fitlm(log(1:length(M_subsample(M_subsample>=10))),log(M_subsample(M_subsample>=10)));
Slope = mdl.Coefficients.Estimate(2);
Inter = mdl.Coefficients.Estimate(1);

xx=log(1:length(z0));
Mave1 = exp(Inter + log(length(n_k0)/length(n_k)) + Slope*xx);
%Mave1 = exp(mdl.Coefficients.Estimate(1) + log(length(unique(z0))/length(unique(z))) + mdl.Coefficients.Estimate(2)*xx);

RMSE1 = log(M(dexM))-log(Mave1(dexM));
RMSE1 = sqrt(sum(RMSE1.^2)/length(dexM));

Chi2_1 = chi2_rate(M, Mave1,50);

figure;
loglog(M,'o');
hold on;loglog(Mave,'r');
loglog(Mave1,'b');

save(['FoF_LS_',name_dataset,'_',num2str(ttt),num2str(randtry),'.mat'],'Chi2','RMSE','Mave','Chi2_1','RMSE1','Mave1');
