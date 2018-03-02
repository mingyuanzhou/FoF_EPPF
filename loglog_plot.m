r=10;
for a=[-2:0.5:0,0.1:0.1:0.9]
p=0.99;
PoissonRate_log = (log(r)+gammaln((1:length(M))-a)-gammaln(1-a)...
                    -gammaln(1+(1:length(M))) + ((1:length(M))*0-a)*log(p));
loglog(1:length(M),exp(PoissonRate_log))  ;hold on
loglog(1:length(M),(1:length(M)).^(-a-1)*r*p^(-a)/gamma(1-a))  ;hold on
title(num2str(a))
pause

end

alpha0=2;
n_k = round((0.8-1/2)*(1-rand(round(100000),1)).^(-1/(alpha0-1))+0.5);
M = nk_to_m(n_k);
[z,n_k] = m_i_to_z(M(M>0),find(M));


figure
alpha_mean=zeros(10,8);
x_mean=zeros(10,8);
count=0;
Simpson =zeros(10,8);
Slope=zeros(10,8);
for t=4.^[7:-1:0]
    count=count+1;
    for iter=1:10
        
        %dex=randperm(length(z));
        
        %n_k=full(sparse(1,z(dex(1:fix(end/t))),1));
        if t==1
            y = datasample(z,fix(length(z)/t));
        else
            y = datasample(z,fix(length(z)/t),'Replace',false);
        end
        n_k = full(sparse(1,y,1));
        
      % n_k = round((1-1/2)*(1-rand(round(50000/t),1)).^(-1/(alpha0-1))+0.5);
        
        n_k=n_k(n_k>0);
        M = nk_to_m(n_k);
        
        
        %figure(1)
        loglog(1:length(M),M,'-.');hold on
        drawnow
        %pause(0.1)
        
        %  if count<6
        %      [alpha, xmin, L]=plfit(n_k);
        %  else
        %      xmin=9-log2(t);
        %     [alpha, xmin, L]=plfit(n_k, 'xmin',xmin);
        %  end
       % xmin=3;
        [alpha, xmin, L]=plfit(n_k); %, 'xmin',xmin);
        plplot(n_k,xmin,alpha);alpha
        
        
        Simpson(iter,count)=1- sum(n_k.*(n_k-1)/(sum(n_k)*(sum(n_k)-1)));
        close
        alpha_mean(iter,count)=alpha;
        x_mean(iter,count)=xmin;
        
        mdl = fitlm(log(1:length(M(M>=10))),log(max(M(M>=10),0.5)));
        Slope(iter,count)=abs(mdl.Coefficients.Estimate(2));
        xx=log(1:length(M));
        loglog(1:length(M),exp(mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*xx) );hold on
    end
end
figure
alpha_mean
x_mean
subplot(2,2,1);boxplot(alpha_mean)
subplot(2,2,2);boxplot(x_mean);
subplot(2,2,3);boxplot(Simpson);
subplot(2,2,4);boxplot(Slope);
%[alpha, xmin, L]=plfit(n_k,'xmin',1);plplot(n_k,xmin,alpha);alpha


M1=M(M>0);
figure;loglog(1:length(M),M);hold on;loglog(1:length(M1),M1)

x=log(1:length(M(1:end)));
y=log(max(M(1:end),0.5));
xx=log(1:length(M));
figure;plot(xx,log(max(M,0.5)));hold on
mdl = fitlm(log(1:length(M(M>10))),log(max(M(M>10),0.5)))
plot(xx,mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*xx );hold on
mdl.Coefficients.Estimate(2)
%mdl.plot
a = -1-mdl.Coefficients.Estimate(2)


Rate=(log(r)+gammaln((1:length(M))-a)-gammaln(1-a)-gammaln(1+(1:length(M))) + ((1:length(M))-a)*log(p));
plot(xx,Rate)
n_k_sort=sort(n_k,'descend');
figure;loglog(n_k_sort)