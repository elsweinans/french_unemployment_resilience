% %% the data
% data=readmatrix('swiss_data2020_v3.txt');
% data=data(:,3:end);
% figure
% plot(data)
% xticks([1:12:537])
% xticklabels([1975:2019])

%% alternative data
%T = readtable('swiss_data2020_v3.txt','TreatAsEmpty',{'.','NA'});
T = readtable('french_data2020.txt','TreatAsEmpty',{'.','NA'});
data=T{:,3:end};

crisis1 = 1990;
crisis2 = 2000;
crisis3 = 2008;
crisis4 = 2012;

crisis1xa = [(crisis1-1982)*4+1 (crisis1-1982)*4+5 (crisis1-1982)*4+5 (crisis1-1982)*4+1];
crisis1xb = [(crisis1-1982)*4-39 (crisis1-1982)*4-35 (crisis1-1982)*4-35 (crisis1-1982)*4-39];
crisis2xa = [(crisis2-1982)*4+1 (crisis2-1982)*4+5 (crisis2-1982)*4+5 (crisis2-1982)*4+1];
crisis2xb = [(crisis2-1982)*4-39 (crisis2-1982)*4-35 (crisis2-1982)*4-35 (crisis2-1982)*4-39];
crisis3xa = [(crisis3-1982)*4+2 (crisis3-1982)*4+6 (crisis3-1982)*4+6 (crisis3-1982)*4+2];
crisis3xb = [(crisis3-1982)*4-38 (crisis3-1982)*4-34 (crisis3-1982)*4-34 (crisis3-1982)*4-38];
crisis4xa = [(crisis4-1982)*4+4 (crisis4-1982)*4+5 (crisis4-1982)*4+5 (crisis4-1982)*4+4];
crisis4xb = [(crisis4-1982)*4-36 (crisis4-1982)*4-35 (crisis4-1982)*4-35 (crisis4-1982)*4-36];

crisisypl1 = [0 0 16 16];
crisisypl2 = [-1 -1 0.3 0.3];
crisisypl3 = [60 60 100 100];
crisisypl4 = [0.05 0.05 0.35 0.35];


figure('position',[100,50,900,600])
subplot(4,1,1)
plot(data)
xticks([1:4:150]) %xticks([1:12:537])
xticklabels([]) %xticklabels([1982:2019])
xlim([1,150])
ylim([min(min(data)),max(max(data))])
ylabel('unemployment')
patch(crisis1xa,crisisypl1,'b','FaceAlpha',0.1)
patch(crisis2xa,crisisypl1,'b','FaceAlpha',0.1)
patch(crisis3xa,crisisypl1,'b','FaceAlpha',0.1)
patch(crisis4xa,crisisypl1,'b','FaceAlpha',0.1)

%% moving window PCA
windowyears = 10;
windowsize = windowyears * 4; %windowsize = windowyears * 12;
datasize = length(data);
nrsteps = datasize - windowsize + 1;
PCs = zeros(nrsteps,96);
explvars = zeros(nrsteps,1);

for i = 1:nrsteps
   
    data2 = data(i:i+windowsize - 1,:);
    [coeff,~,~,~,explained] = pca(data2);
    PCs(i,:) = coeff(:,1); % * explained(1);
    explvars(i) = explained(1);
    
end
   
subplot(4,1,3)
plot(PCs)
xticks([1:4:nrsteps])
xticklabels([])
xlim([-windowyears*4+1,nrsteps-1])
ylim([min(min(PCs)),max(max(PCs))])
ylabel('PC1 coefficients')
patch(crisis1xb,crisisypl2,'b','FaceAlpha',0.1)
patch(crisis2xb,crisisypl2,'b','FaceAlpha',0.1)
patch(crisis3xb,crisisypl2,'b','FaceAlpha',0.1)
patch(crisis4xb,crisisypl2,'b','FaceAlpha',0.1)


subplot(4,1,4)
plot(explvars,'-k','linewidth',1.5)
xticks([1:4:nrsteps])
xticklabels([1992:2019]) %xticklabels([1982 + windowyears:2019])
xlim([-windowyears*4+1,nrsteps-1])
ylim([min(explvars),max(explvars)])
ylabel('explained variance')
xlabel('year')
patch(crisis1xb,crisisypl3,'b','FaceAlpha',0.1)
patch(crisis2xb,crisisypl3,'b','FaceAlpha',0.1)
patch(crisis3xb,crisisypl3,'b','FaceAlpha',0.1)
patch(crisis4xb,crisisypl3,'b','FaceAlpha',0.1)

tsMI=importdata('tsMIS_France.txt')
subplot(4,1,2)
plot(tsMI.data,'-k','linewidth',1.5)
xticks([1:4:150]) %xticks([1:12:537])
xticklabels([1982:2019])
ylabel('Moran''s I')
xlim([1,150])
ylim([min(tsMI.data),max(tsMI.data)])
patch(crisis1xa,crisisypl4,'b','FaceAlpha',0.1)
patch(crisis2xa,crisisypl4,'b','FaceAlpha',0.1)
patch(crisis3xa,crisisypl4,'b','FaceAlpha',0.1)
patch(crisis4xa,crisisypl4,'b','FaceAlpha',0.1)

%exportgraphics(gcf,'TS_france.png','Resolution',600)

%writematrix(PCs, "PCs_frenchdata.txt")

%% significance for moran's I

moransmax=zeros(150-11,1);
kendalltaus=zeros(150-11,1);
kendallp=zeros(150-11,1);
alarm=zeros(150-11,1);

for i=1:150-11
    
    moranswindow=tsMI.data(i:i+11);
    moransmax(i)=mean(moranswindow)+std(moranswindow); 
    [rho,p]=corr([1:12]',moranswindow,'Type','Kendall');
    kendalltaus(i)=rho;
    kendallp(i)=p;
    if kendalltaus(i)>0.9 && kendallp(i)<0.05 && moransmax(i)<moranswindow(12)
        alarm(i)=1
    end

end

figure
hold on
x=[12:150];
patch([x flip(x)], [alarm' zeros(size(x))], 'y', 'EdgeColor','none','FaceAlpha',0.4)
plot([1:150],tsMI.data,'-k','linewidth',1.5)
plot([12:150],moransmax,'--r','linewidth',1.2)
xticks([1:4:150])
xticklabels([1982:2019])
ylabel('moran''s I')
legend('','morans''s I','mean + std','Location','northwest')
ylim([0.05,0.35])
%exportgraphics(gcf,'trendtest_France.png','Resolution',600)

figure
subplot(2,1,1)
hold on
plot([12:150],kendalltaus,'-k','linewidth',1.2)
plot([1,150],[0.9,0.9],'--r','linewidth',1.2)
xlim([1,150])
ylabel('kendall trend')

subplot(2,1,2)
hold on
plot([12:150],kendallp,'-k','linewidth',1.5)
plot([1,150],[0.05,0.05],'--r','linewidth',1.2)
xlim([1,150])
ylabel('kendall p')

%% sensitivity test ws
wss=[8,12,16];
figure
for j = 1:3
    ws=wss(j);
    moransmax=zeros(150-ws+1,1);
    kendalltaus=zeros(150-ws+1,1);
    kendallp=zeros(150-ws+1,1);
    alarm=zeros(150-ws+1,1);
    
    for i=1:150-ws+1
        
        moranswindow=tsMI.data(i:i+ws-1);
        moransmax(i)=mean(moranswindow)+1*std(moranswindow); 
        [rho,p]=corr([1:ws]',moranswindow,'Type','Kendall');
        kendalltaus(i)=rho;
        kendallp(i)=p;
        if kendalltaus(i)>0.9 && kendallp(i)<0.05 && moransmax(i)<moranswindow(ws)
            alarm(i)=1;
        end
    
    end
    
    subplot(3,1,j)
    hold on
    x=[ws:150];
    patch([x flip(x)], [alarm' zeros(size(x))], 'y', 'EdgeColor','none','FaceAlpha',0.4)
    plot([1:150],tsMI.data,'-k','linewidth',1.5)
    plot([ws:150],moransmax,'--r','linewidth',1.2)
    xticks([1:4:150])
    xticklabels([])
    ylabel('moran''s I')
    if j==1
        legend('','morans''s I','mean + std','Location','northwest')
    end
    ylim([0.05,0.35])
    title('window size = '+string(ws/4)+' years')
end
xticklabels([1982:2019])

%exportgraphics(gcf,'trendtest_France_sens_ws.png','Resolution',600)


%% sensitivity test rho
minKs=[0.8,0.9,0.95];
ws=12;
figure
for j = 1:3
    minK=minKs(j);
    moransmax=zeros(150-ws+1,1);
    kendalltaus=zeros(150-ws+1,1);
    kendallp=zeros(150-ws+1,1);
    alarm=zeros(150-ws+1,1);
    
    for i=1:150-ws+1
        
        moranswindow=tsMI.data(i:i+ws-1);
        moransmax(i)=mean(moranswindow)+1*std(moranswindow); 
        [rho,p]=corr([1:ws]',moranswindow,'Type','Kendall');
        kendalltaus(i)=rho;
        kendallp(i)=p;
        if kendalltaus(i)>minK && kendallp(i)<0.05 && moransmax(i)<moranswindow(ws)
            alarm(i)=1;
        end
    
    end
    
    subplot(3,1,j)
    hold on
    x=[ws:150];
    patch([x flip(x)], [alarm' zeros(size(x))], 'y', 'EdgeColor','none','FaceAlpha',0.4)
    plot([1:150],tsMI.data,'-k','linewidth',1.5)
    plot([ws:150],moransmax,'--r','linewidth',1.2)
    xticks([1:4:150])
    xticklabels([])
    ylabel('moran''s I')
    if j==1
        legend('','morans''s I','mean + std','Location','northwest')
    end
    ylim([0.05,0.35])
    title('minimum Kendall tau correlation = '+string(minK))
end
xticklabels([1982:2019])

%exportgraphics(gcf,'trendtest_France_sens_kendalltau.png','Resolution',600)

%% sensitivity test std
stdns=[1,1.5,2];
figure
for j = 1:3
    stdn=stdns(j);
    moransmax=zeros(150-ws+1,1);
    kendalltaus=zeros(150-ws+1,1);
    kendallp=zeros(150-ws+1,1);
    alarm=zeros(150-ws+1,1);
    
    for i=1:150-ws+1
        
        moranswindow=tsMI.data(i:i+ws-1);
        moransmax(i)=mean(moranswindow)+stdn*std(moranswindow); 
        [rho,p]=corr([1:ws]',moranswindow,'Type','Kendall');
        kendalltaus(i)=rho;
        kendallp(i)=p;
        if kendalltaus(i)>0.9 && kendallp(i)<0.05 && moransmax(i)<moranswindow(ws)
            alarm(i)=1;
        end
    
    end
    
    subplot(3,1,j)
    hold on
    x=[ws:150];
    patch([x flip(x)], [alarm' zeros(size(x))], 'y', 'EdgeColor','none','FaceAlpha',0.4)
    plot([1:150],tsMI.data,'-k','linewidth',1.5)
    plot([ws:150],moransmax,'--r','linewidth',1.2)
    xticks([1:4:150])
    xticklabels([])
    ylabel('moran''s I')
    if j==1
        legend('','morans''s I','mean + std','Location','northwest')
    end
    ylim([0.05,0.35])
    title('standard deviation multiplied by '+string(stdn))
end
xticklabels([1982:2019])

%exportgraphics(gcf,'trendtest_France_sens_std.png','Resolution',600)

%% aggregate data
aggrdata=sum(data');

ws=8;

AR1=zeros(150-ws+1,1);
vars=zeros(150-ws+1,1);

for i=1:150-ws-1
    
    datawindow=aggrdata(i:i+ws-1);
    cor=corrcoef(datawindow(1:end-1),datawindow(2:end));
    AR1(i)=cor(2,1);
    vars(i)=var(datawindow);

end

crisisypl5 = [600 600 1000 1000];
crisisypl6 = [-0.5 -0.5 1 1];
crisisypl7 = [0.0 0.0 10000 10000];

figure
subplot(3,1,1)
hold on
plot([1:150],aggrdata,'-k','linewidth',1.5)
ylabel('unemployment')
patch(crisis1xa,crisisypl5,'b','FaceAlpha',0.1)
patch(crisis2xa,crisisypl5,'b','FaceAlpha',0.1)
patch(crisis3xa,crisisypl5,'b','FaceAlpha',0.1)
patch(crisis4xa,crisisypl5,'b','FaceAlpha',0.1)
xticks([1:4:150])
xticklabels([])
text(0.5,980,'A','FontWeight', 'Bold')

subplot(3,1,2)
hold on
plot([ws:150],AR1,'-k','linewidth',1.5)
ylabel('autocorrelation')
patch(crisis1xa,crisisypl6,'b','FaceAlpha',0.1)
patch(crisis2xa,crisisypl6,'b','FaceAlpha',0.1)
patch(crisis3xa,crisisypl6,'b','FaceAlpha',0.1)
patch(crisis4xa,crisisypl6,'b','FaceAlpha',0.1)
xticks([1:4:150])
xticklabels([])
text(0.5,0.92,'B','FontWeight', 'Bold')

subplot(3,1,3)
hold on
plot([ws:150],vars,'-k','linewidth',1.5)
ylabel('variance')
patch(crisis1xa,crisisypl7,'b','FaceAlpha',0.1)
patch(crisis2xa,crisisypl7,'b','FaceAlpha',0.1)
patch(crisis3xa,crisisypl7,'b','FaceAlpha',0.1)
patch(crisis4xa,crisisypl7,'b','FaceAlpha',0.1)
xticks([1:4:150])
xticklabels([1982:2019])
text(0.5,9300,'C','FontWeight', 'Bold')

% exportgraphics(gcf,'univariate_analysis_France.png','Resolution',600)




