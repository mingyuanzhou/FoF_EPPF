%Power_law
figure
z=read_data(9);

Simpson =zeros(20,5);
Slope=zeros(20,5);
UnitRatio = zeros(20,5);
%colors = {'r.','b.','k.','c.','m.'};
%colors1 = {'r-','b-','k-','c-','m-'};
rng(0,'twister')
for iter=1:20
    count=0;
    for t=4.^(0:1:4)
        count=count+1;
        
        
        if t==1
            if iter==1
                y=z;
            else
                y = datasample(z,fix(length(z)/t));
            end
        else
            y = datasample(z,fix(length(z)/t),'Replace',false);
        end
        n_k = full(sparse(1,y,1));
        n_k=n_k(n_k>0);
        M = nk_to_m(n_k);
        
        UnitRatio(iter,count)= nnz(n_k==1)/length(n_k);
        
        % Mmax=10^(ceil(log10(max(M)))+0);
       % Mmax=max(M);
        Mmax=10^(ceil(log10(max(max(M),length(M))))+0);
        Simpson(iter,count)=1- sum(n_k.*(n_k-1)/(sum(n_k)*(sum(n_k)-1)));
        %mdl = fitlm(log(1:length(M(M>=10))),log(M(M>=10)));
        %bug fix Feb 13,2016
        idx = find(M>=10);
        mdl = fitlm(log(idx),log(M(idx)));
       
        Slope(iter,count)=(mdl.Coefficients.Estimate(2));
        
        
        [alpha, xmin, L]=plfit(n_k);
        
        idx = find(M>=3);
        %idx = idx(idx>=min(xmin,10));
        idx=idx(idx>=xmin);
        Slope(iter,count) = -alpha;
        Intercept = mean( log(M(idx)))-Slope(iter,count)*mean(log(idx));
        
        subplot(1,5,count);
        %figure(1)
        %if iter==1
        loglog(1:length(M),M,'o');hold on
        
        xx=log(1:length(M));
        loglog(1:length(M),exp(Intercept + Slope(iter,count)*xx) );hold on
        ylim([0.5,Mmax])
        xlim([1,Mmax])
        drawnow
        %end
        
    end
end

texts={'(a) Sampling ratio = 1','(b) Sampling ratio = 1/4','(c) Sampling ratio = 1/16','(d) Sampling ratio = 1/64','(e) Sampling ratio = 1/256'}
for i=1:5
    subplot(1,5,i)
    xlabel('ln(i)')
    ylabel('ln(m_i)')
    title(texts{i})
end

figure
subplot(1,3,1);boxplot(Slope);
ax = gca;
ax.XTickLabel = {'1','1/4','1/16','1/64','1/256'};
xlabel('Sampling ratio')
ylabel('Slope')
title('(a)')

subplot(1,3,2);boxplot(UnitRatio);
ax = gca;
ax.XTickLabel = {'1','1/4','1/16','1/64','1/256'};
xlabel('Sampling ratio')
ylabel('Ratio of unit-size clusters')
title('(b)')

subplot(1,3,3);boxplot(Simpson);
ax = gca;
ax.XTickLabel = {'1','1/4','1/16','1/64','1/256'};
xlabel('Sampling ratio')
ylabel('Simpson''s index of diversity')
title('(c)')




