
addpath results
%name_dataset='text_tomsawyer';
%[z0,n_k,M,name_dataset]=read_data(1);

%name_dataset='gene_sultan';
[z0,n_k,M,name_dataset]=read_data(9);

%name_dataset='microdata';
%[z0,n_k,M,name_dataset]=read_data(11);

SSS1 = zeros(5,5,9);
SSS2 = zeros(5,5,9);
for model=[1:6,7,8,9]
    for ttt=1:5
        for randtry=1:5
            if model<=6
                load(['FoF_',name_dataset,'_',num2str(ttt),num2str(model),num2str(randtry),'.mat'],'Chi2','RMSE');
            else
                load(['FoF_PY_',name_dataset,'_',num2str(ttt),num2str(model),num2str(randtry),'.mat'],'Chi2','RMSE');
            end
            SSS1(randtry,ttt,model)=RMSE;
            %Chi2=chi2_rate(M, Mave/500,51);
            SSS2(randtry,ttt,model)=Chi2;
        end
    end
end

figure;
Colors={'r-o','g-x','b-->','k--d','c:s','m:*','-.^','-.p','-.+'};
%Colors={'-o','-x','-->','--d',':s',':*','-.h','-.p','-.+'};

%Colors={'r-','g-','b-.','k-.','c:','m:'};
subplot(1,2,1);
count=0;
%for model=[4,2,3,6,5,1,9,8]
for model=[4,2,3,9,6,5,8,1]
    count=count+1;
    h=errorbar((1:5)+(count-3)*0.04, mean(SSS1(:,:,model),1), std(SSS1(:,:,model),1),Colors{count},'LineWidth',1,'MarkerSize',8); hold on;
end
%set(get(h,'Parent'), 'YScale', 'log')
set(get(h,'Parent'), 'XTick', 1:5, 'XTickLabel',{'1/32','1/16','1/8','1/4','1/2'})
legend('a = --1','a = 0','a = 0.5','PY a=0.5','a < 0','0<=a<1','PY','a < 1')
xlabel('Sampling ratio');
ylabel('RMSE')
title('(a)')
subplot(1,2,2);
count=0;
for model=[4,2,3,9,6,5,8,1]
    %model=[4,2,3,6,5,1,9,8]
    count=count+1;
    h=errorbar((1:5)+(count-3)*0.04, mean(SSS2(:,:,model),1), std(SSS2(:,:,model),1),Colors{count},'LineWidth',1,'MarkerSize',8); hold on;
end
set(get(h,'Parent'), 'YScale', 'log')
set(get(h,'Parent'), 'XTick', 1:5, 'XTickLabel',{'1/32','1/16','1/8','1/4','1/2'})
%legend('a = --1','a = 0','a = 0.5','a < 0','0<=a<1','a < 1','PY a=0.5','PY')
legend('a = --1','a = 0','a = 0.5','PY a=0.5','a < 0','0<=a<1','PY','a < 1')
xlabel('Sampling ratio');
ylabel('\chi^2')
title('(b)')





%%


%name_dataset='text_tomsawyer';
[z0,n_k,M,name_dataset]=read_data(1);

%name_dataset='gene_sultan';
[z0,n_k,M,name_dataset]=read_data(9);


%name_dataset='microdata';
%[z0,n_k,M,name_dataset]=read_data(9);
figure;
N = sum((1:length(M)).*M);
ii=1:length(M);
Entropy0=sum(-(ii/N.*log(ii./N)).*M)

Colors={'r-o','g-x','b-->','k--d','c:s','m:*','-.^','-.p','-.+'};
Colors1={'r-','g-','b--','k--','c:','m:','-.','-.','-.'};
texts={'(a) Sampling ratio = 1/32','(b) Sampling ratio = 1/16','(c) Sampling ratio = 1/8','(d) Sampling ratio = 1/4','(e) Sampling ratio = 1/2'}
Entropy = zeros(5,8);
for ttt=1:5
    ttt
    subplot(1,5,ttt)
    loglog(M,'o');hold on
    count=0;
    temp=0;
    for model=[4,2,3,9,6,5,8,1]
        randtry=1;
        count=count+1;
        if model<=6
            load(['FoF_',name_dataset,'_',num2str(ttt),num2str(model),num2str(randtry),'.mat'],'Mave','PoissonRate'); 
        else
            load(['FoF_PY_',name_dataset,'_',num2str(ttt),num2str(model),num2str(randtry),'.mat'],'Mave');
        end
        
        loglog(Mave/500,Colors1{count},'LineWidth',2);
        if max(Mave/500)>temp
            temp=10^(ceil(log10(max(Mave/500))));
            ylim([0.1,temp])
        end
        xlim([0,length(M)])
        %loglog(PoissonRate/500,'.');
        ii=1:length(Mave/500);
        N = sum((1:length(Mave)).*Mave/500);
        Entropy(ttt,count)=sum(-(ii/N.*log(ii./N)).*Mave/500);
    end
    % legend('FoF','a = --1','a = 0','a = 0.5','a < 0','0<=a<1','a < 1')
    legend('FoF','a = --1','a = 0','a = 0.5','PY a=0.5','a < 0','0<=a<1','PY','a < 1')
    
    xlabel('ln(i)');
    ylabel('log(m_i)');
    title(texts{ttt})
end

Colors1={'r-','g-','b--','k--','c:','m:'};
figure;
for count=1:8
    plot(Entropy(:,count)-Entropy0,Colors{count},'LineWidth',2);
    hold on
end
%legend('a = --1','a = 0','a = 0.5','a < 0','0<=a<1','a < 1')
legend('a = --1','a = 0','a = 0.5','PY a=0.5','a < 0','0<=a<1','PY','a < 1')
    


%
%
%
% Entropy = zeros(5,5,6);
% for model=1:6
%     for ttt=1:5
%         for randtry=1:5
%             load(['FoF_',name_dataset,'_',num2str(ttt),num2str(model),num2str(randtry),'.mat'],'Mave');
%             ii=1:length(Mave/500);
%             N = sum((1:length(Mave)).*Mave/500);
%             Entropy(randtry,ttt,model)=sum(-(ii/N.*log(ii./N)).*Mave/500)-Entropy0;
%         end
%     end
% end
% figure;
% Colors={'r-o','g-x','b-->','k--d','c:s','m:*'};
% %Colors={'r-','g-','b-.','k-.','c:','m:'};
% subplot(1,2,1);
% count=0;
% for model=[4,2,3,6,5,1]
%     count=count+1;
%     h=errorbar((1:5)+(count-3)*0.04, mean( Entropy(:,:,model),1), std( Entropy(:,:,model),1),Colors{count},'LineWidth',1,'MarkerSize',8); hold on;
% end
% %set(get(h,'Parent'), 'YScale', 'log')
% set(get(h,'Parent'), 'XTick', 1:5, 'XTickLabel',{'1/32','1/16','1/8','1/4','1/2'})
% legend('a = --1','a = 0','a = 0.5','a < 0','0<=a<1','a < 1')
% xlabel('Sampling ratio');
% ylabel('RMSE')
% title('(a)')
