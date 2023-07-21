clear all;
close all;

addpath(genpath('/Users/umberto/surfdrive/Shared/IAS_paper/repo/IAS_project/ephys analysis/helper functions'));

addpath(genpath('/Users/umberto/Documents/GitRepositories/steinmetz-et-al-2019'));

load('/Users/umberto/surfdrive/Shared/INFLOW-PSTHs/PSTHoutput/dataPSTH.mat')

colors=[0.5 0 0;
    1 0 0;
    0 0.5 0;
    0 1 0];
t=(-50:149)*10;

smoothWin=0.025;
dt=0.01;

smWin = myGaussWin(smoothWin, 1/dt);
smWin(1:round(numel(smWin)/2)) = 0;
smWin = smWin./sum(smWin);

compute_sig=1;
alpha=0.025;

subplot(3,1,1);
PSTH_avg_norm=[];
PSTH_err_norm=[];
for i=1:4
    bsl{i}=mean(dataPSTH.V1.visual.PSTHmat{i}(:,1:99),2);
    PSTHmat_norm=dataPSTH.V1.visual.PSTHmat{i}-repmat(bsl{i},1,250);
    for j=1:length(bsl{i})
        PSTHmat_norm(j,:)=conv(PSTHmat_norm(j,:),smWin,'same');
    end
    if compute_sig
%         for j=1:size(PSTHmat_norm,2)
%             %[h(i,j),p(i,j)]=ttest(PSTHmat_norm(:,j));
%             [p(i,j),h(i,j)]=signrank(PSTHmat_norm(:,j),mean(PSTHmat_norm(:,51:99),2));
%         end
        dataPSTH.V1.visual.PSTHmat_norm{i}=PSTHmat_norm;
    end
    PSTH_avg_norm(:,i)=mean(PSTHmat_norm);
    PSTH_err_norm(:,i)=std(PSTHmat_norm)./sqrt(length(bsl{i}));
    plot(t,PSTH_avg_norm(51:end,i),'Color',colors(i,:),'LineWidth',2);
    hold on;
end
if compute_sig
    for j=1:size(PSTHmat_norm,2)
        %[hs(1,j),ps(1,j)]=ttest2(dataPSTH.V1.visual.PSTHmat_norm{1}(:,j),dataPSTH.V1.visual.PSTHmat_norm{2}(:,j));
        %[hs(2,j),ps(2,j)]=ttest2(dataPSTH.V1.visual.PSTHmat_norm{3}(:,j),dataPSTH.V1.visual.PSTHmat_norm{4}(:,j));
        %[hr(1,j),pr(1,j)]=ttest2(dataPSTH.V1.visual.PSTHmat_norm{1}(:,j),dataPSTH.V1.visual.PSTHmat_norm{3}(:,j));
        %[hr(2,j),pr(2,j)]=ttest2(dataPSTH.V1.visual.PSTHmat_norm{2}(:,j),dataPSTH.V1.visual.PSTHmat_norm{4}(:,j));
        [hs(1,j),ps(1,j)]=realdata_bootstrap(dataPSTH.V1.visual.PSTHmat_norm{1}(:,j),dataPSTH.V1.visual.PSTHmat_norm{2}(:,j),1000,alpha);
        [hs(2,j),ps(2,j)]=realdata_bootstrap(dataPSTH.V1.visual.PSTHmat_norm{3}(:,j),dataPSTH.V1.visual.PSTHmat_norm{4}(:,j),1000,alpha);
        [hr(1,j),pr(1,j)]=realdata_bootstrap(dataPSTH.V1.visual.PSTHmat_norm{1}(:,j),dataPSTH.V1.visual.PSTHmat_norm{3}(:,j),1000,alpha);
        [hr(2,j),pr(2,j)]=realdata_bootstrap(dataPSTH.V1.visual.PSTHmat_norm{2}(:,j),dataPSTH.V1.visual.PSTHmat_norm{4}(:,j),1000,alpha);
    end
    hs=fdr_bh(ps(:,51:end),0.05,'pdep');
    hr=fdr_bh(pr(:,51:end),0.05,'pdep');
    hs_V1=(hs(1,:) & hs(2,:));
    hr_V1=(hr(1,:) & hr(2,:));
end
for i=1:4
    shadedErrorBar(t,PSTH_avg_norm(51:end,i),PSTH_err_norm(51:end,i),{'-','LineWidth',2,'Color',colors(i,:)},1);
    hold on
end
set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
title('V1','FontSize',20);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
legend(dataPSTH.V1.visual.selectionLegend);

subplot(3,1,2);
PSTH_avg_norm=[];
PSTH_err_norm=[];
for i=1:4
    bsl{i}=mean(dataPSTH.PPC.visual.PSTHmat{i}(:,1:99),2);
    PSTHmat_norm=dataPSTH.PPC.visual.PSTHmat{i}-repmat(bsl{i},1,250);
    for j=1:length(bsl{i})
        PSTHmat_norm(j,:)=conv(PSTHmat_norm(j,:),smWin,'same');
    end
    if compute_sig
%         for j=1:size(PSTHmat_norm,2)
%             %[h(i,j),p(i,j)]=ttest(PSTHmat_norm(:,j));
%             [p(i,j),h(i,j)]=signrank(PSTHmat_norm(:,j),mean(PSTHmat_norm(:,51:99),2));
%         end
        dataPSTH.PPC.visual.PSTHmat_norm{i}=PSTHmat_norm;
    end
    PSTH_avg_norm(:,i)=mean(PSTHmat_norm);
    PSTH_err_norm(:,i)=std(PSTHmat_norm)./sqrt(length(bsl{i}));
    plot(t,PSTH_avg_norm(51:end,i),'Color',colors(i,:),'LineWidth',2);
    hold on;
end
if compute_sig
    for j=1:size(PSTHmat_norm,2)
        %[hs(1,j),ps(1,j)]=ttest2(dataPSTH.PPC.visual.PSTHmat_norm{1}(:,j),dataPSTH.PPC.visual.PSTHmat_norm{2}(:,j));
        %[hs(2,j),ps(2,j)]=ttest2(dataPSTH.PPC.visual.PSTHmat_norm{3}(:,j),dataPSTH.PPC.visual.PSTHmat_norm{4}(:,j));
        %[hr(1,j),pr(1,j)]=ttest2(dataPSTH.PPC.visual.PSTHmat_norm{1}(:,j),dataPSTH.PPC.visual.PSTHmat_norm{3}(:,j));
        %[hr(2,j),pr(2,j)]=ttest2(dataPSTH.PPC.visual.PSTHmat_norm{2}(:,j),dataPSTH.PPC.visual.PSTHmat_norm{4}(:,j));
        [hs(1,j),ps(1,j)]=realdata_bootstrap(dataPSTH.PPC.visual.PSTHmat_norm{1}(:,j),dataPSTH.PPC.visual.PSTHmat_norm{2}(:,j),1000,alpha);
        [hs(2,j),ps(2,j)]=realdata_bootstrap(dataPSTH.PPC.visual.PSTHmat_norm{3}(:,j),dataPSTH.PPC.visual.PSTHmat_norm{4}(:,j),1000,alpha);
        [hr(1,j),pr(1,j)]=realdata_bootstrap(dataPSTH.PPC.visual.PSTHmat_norm{1}(:,j),dataPSTH.PPC.visual.PSTHmat_norm{3}(:,j),1000,alpha);
        [hr(2,j),pr(2,j)]=realdata_bootstrap(dataPSTH.PPC.visual.PSTHmat_norm{2}(:,j),dataPSTH.PPC.visual.PSTHmat_norm{4}(:,j),1000,alpha);
    end
    hs=fdr_bh(ps(:,51:end),0.05,'pdep');
    hr=fdr_bh(pr(:,51:end),0.05,'pdep');
    hs_PPC=(hs(1,:) & hs(2,:));
    hr_PPC=(hr(1,:) & hr(2,:));
end
for i=1:4
    shadedErrorBar(t,PSTH_avg_norm(51:end,i),PSTH_err_norm(51:end,i),{'-','LineWidth',2,'Color',colors(i,:)},1);
    hold on
end
set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
title('PPC','FontSize',20);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
legend(dataPSTH.PPC.visual.selectionLegend);

subplot(3,1,3);
PSTH_avg_norm=[];
PSTH_err_norm=[];
for i=1:4
    bsl{i}=mean(dataPSTH.CG1.visual.PSTHmat{i}(:,1:99),2);
    PSTHmat_norm=dataPSTH.CG1.visual.PSTHmat{i}-repmat(bsl{i},1,250);
    for j=1:length(bsl{i})
        PSTHmat_norm(j,:)=conv(PSTHmat_norm(j,:),smWin,'same');
    end
    if compute_sig
%         for j=1:size(PSTHmat_norm,2)
%             %[h(i,j),p(i,j)]=ttest(PSTHmat_norm(:,j));
%             [p(i,j),h(i,j)]=signrank(PSTHmat_norm(:,j),mean(PSTHmat_norm(:,51:99),2));
%         end
        dataPSTH.CG1.visual.PSTHmat_norm{i}=PSTHmat_norm;
    end
    PSTH_avg_norm(:,i)=mean(PSTHmat_norm);
    PSTH_err_norm(:,i)=std(PSTHmat_norm)./sqrt(length(bsl{i}));
    plot(t,PSTH_avg_norm(51:end,i),'Color',colors(i,:),'LineWidth',2);
    hold on;
end
if compute_sig
    for j=1:size(PSTHmat_norm,2)
        %[hs(1,j),ps(1,j)]=ttest2(dataPSTH.CG1.visual.PSTHmat_norm{1}(:,j),dataPSTH.CG1.visual.PSTHmat_norm{2}(:,j));
        %[hs(2,j),ps(2,j)]=ttest2(dataPSTH.CG1.visual.PSTHmat_norm{3}(:,j),dataPSTH.CG1.visual.PSTHmat_norm{4}(:,j));
        %[hr(1,j),pr(1,j)]=ttest2(dataPSTH.CG1.visual.PSTHmat_norm{1}(:,j),dataPSTH.CG1.visual.PSTHmat_norm{3}(:,j));
        %[hr(2,j),pr(2,j)]=ttest2(dataPSTH.CG1.visual.PSTHmat_norm{2}(:,j),dataPSTH.CG1.visual.PSTHmat_norm{4}(:,j));
        [hs(1,j),ps(1,j)]=realdata_bootstrap(dataPSTH.CG1.visual.PSTHmat_norm{1}(:,j),dataPSTH.CG1.visual.PSTHmat_norm{2}(:,j),1000,alpha);
        [hs(2,j),ps(2,j)]=realdata_bootstrap(dataPSTH.CG1.visual.PSTHmat_norm{3}(:,j),dataPSTH.CG1.visual.PSTHmat_norm{4}(:,j),1000,alpha);
        [hr(1,j),pr(1,j)]=realdata_bootstrap(dataPSTH.CG1.visual.PSTHmat_norm{1}(:,j),dataPSTH.CG1.visual.PSTHmat_norm{3}(:,j),1000,alpha);
        [hr(2,j),pr(2,j)]=realdata_bootstrap(dataPSTH.CG1.visual.PSTHmat_norm{2}(:,j),dataPSTH.CG1.visual.PSTHmat_norm{4}(:,j),1000,alpha);
    end
    hs=fdr_bh(ps(:,51:end),0.05,'pdep');
    hr=fdr_bh(pr(:,51:end),0.05,'pdep');
    hs_CG1=(hs(1,:) & hs(2,:));
    hr_CG1=(hr(1,:) & hr(2,:));
end
for i=1:4
    shadedErrorBar(t,PSTH_avg_norm(51:end,i),PSTH_err_norm(51:end,i),{'-','LineWidth',2,'Color',colors(i,:)},1);
    hold on
end
set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
title('ACC','FontSize',20);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
legend(dataPSTH.CG1.visual.selectionLegend);

%% Plot significance
if compute_sig
    figure;
    hs=[hs_V1;hs_PPC;hs_CG1];
    hr=[hr_V1;hr_PPC;hr_CG1];
    hV1=[hs_V1;hr_V1];
    hPPC=[hs_PPC;hr_PPC];
    hCG1=[hs_CG1;hr_CG1];

    subplot(611);
    bar(t,hs_V1,1,'k');
    set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
    title('V1 stimulus coding')
    subplot(613);
    bar(t,hs_PPC,1,'k');
    set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
    title('PPC stimulus coding')
    subplot(615);
    bar(t,hs_CG1,1,'k');
    set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
    title('CG1 stimulus coding')

    subplot(612);
    bar(t,hr_V1,1,'k');
    set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
    title('V1 report coding')
    subplot(614);
    bar(t,hr_PPC,1,'k');
    set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
    title('PPC report coding')
    subplot(616);
    bar(t,hr_CG1,1,'k');
    set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
    title('CG1 report coding')
end