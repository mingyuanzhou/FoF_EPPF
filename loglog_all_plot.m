%Power_law
figure


datasets = {'text_WarAndPeace','text_MOBY_DICK', 'text_tale','text_pride','text_tomsawyer','text_HuckleberryFinn','text_HOLMES','gene_core','gene_sultan','gene_yang','microdata'};
count=0;
figure
texts={'(a) The adventures of Tom Sawyer','(b) Adventures of Sherlock Holmes', '(c) A tale of two cities','(d) War and peace', ...
    '(e) Gene Core', '(f) Gene Sultan', '(g) Gene Yang', '(h) Microdata' }

Simpson =zeros(1,8);
    Slope=zeros(1,8);
for dataset = [1,2,3,4,8,9,10,11]
    count=count+1;
    z=read_data(dataset);
    length(z)
    %Power_law
    subplot(2,4,count);
    
    %colors = {'r.','b.','k.','c.','m.'};
    %colors1 = {'r-','b-','k-','c-','m-'};
    rng(0,'twister')
    for iter=1:1

        for t=1
            
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
            
            [alpha, xmin, L]=plfit(n_k);
            
            M = nk_to_m(n_k);
            
            % Mmax=10^(ceil(log10(max(M)))+0);
            % Mmax=max(M);
            Mmax=10^(ceil(log10(length(M)))+0);
            Simpson(iter,count)=1- sum(n_k.*(n_k-1)/(sum(n_k)*(sum(n_k)-1)));
            %mdl = fitlm(log(1:length(M(M>=10))),log(M(M>=10)));
            %bug fix Feb 13,2016
            idx = find(M>=10);
            
            idx = find(M>=3);
            %idx = idx(idx>=min(xmin,10));
            idx=idx(idx>=xmin);
            mdl = fitlm(log(idx),log(M(idx)));
            Slope(iter,count)=(mdl.Coefficients.Estimate(2));
            Slope(iter,count)=-(alpha);
            Intercept = mean( log(M(idx)))-Slope(iter,count)*mean(log(idx));
            subplot(2,4,count);
            %figure(1)
            %if iter==1
            
            
            xx=log(1:length(M));
            
            %loglog(1:length(M),exp(mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*xx),'LineWidth',2 );hold on
            loglog(1:length(M),exp(Intercept + Slope(iter,count)*xx),'LineWidth',2 );hold on
            
            ylim([0.5,Mmax])
            xlim([1,Mmax])
            %loglog(xmin*ones(1,round(Mmax/2)),1:round(Mmax/2),'k');hold on;
            loglog(xmin,exp(Intercept + Slope(iter,count)*log(xmin)),'ks','MarkerSize',16);hold on;
            loglog(1:length(M),M,'o');hold on
            
            legend(['Slope = ',num2str(Slope(iter,count))], ['i_{min} = ',num2str(xmin)]);
            drawnow
            %end
            
        end
    end
    
    
    
    xlabel('ln(i)')
    ylabel('ln(m_i)')
    title(texts{count})
    
end

figure
%subplot(1,2,1);
%boxplot(Slope);
stem(Slope,'o')
ax = gca;
ax.XTickLabel = {'(a)','(b)','(c)','(d)','(e)', '(f)', '(g)', '(h)'};
xlabel('Sampling ratio')
ylabel('Slope')
%title('(a)')

% subplot(1,2,2);boxplot(Simpson);
% %ax = gca;
% %ax.XTickLabel = {'1','1/4','1/16','1/64','1/256'};
% xlabel('Sampling ratio')
% ylabel('Simpson''s index of diversity')
% title('(b)')


