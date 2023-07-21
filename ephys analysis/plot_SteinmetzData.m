clear all;
close all;

addpath(genpath('/Users/umberto/surfdrive/Shared/IAS_paper/repo/IAS_project/ephys analysis/helper functions')); %---UO

addpath(genpath('/Users/umberto/Documents/GitRepositories/steinmetz-et-al-2019'));

mainFolder='/Users/umberto/Documents/GitRepositories/steinmetz-et-al-2019/Data/allData';

area_list={'VISp','VISa','ACA ','MOs '};

compute_sig=1;
alpha=0.025;

for ak=1:length(area_list)
    
    nArea=area_list{ak};
    dt=0.01;
    window=[-0.6 1.6];
    smoothWin=0.05;
    
    sessions=dir(mainFolder);
    sessions=sessions(4:end);
    
    psth_all=[];
    % all_areas={};
    % id=1;
    for i=1:length(sessions)
        s=loadSession(fullfile(mainFolder,sessions(i).name));
        %     for k=1:size(s.channels.brainLocation.allen_ontology,1)
        %         all_areas{id}=s.channels.brainLocation.allen_ontology(k,:);
        %         id=id+1;
        %     end
        
        for nProbe=1:size(s.probes.sitePositions,2)
            
            [psth,abins]=eventRasters_2ndbump(s,nProbe,nArea,window,dt,smoothWin);
            
            if ~isempty(abins)
                bins=abins;
            end
            
            if isempty(psth_all)
                psth_all=psth;
            else
                psth_all=cat(2,psth_all,psth);
            end
            
        end
    end
    
    jt=smoothWin/dt*2;
    psth_all=psth_all(:,:,jt+1:length(bins)-jt);
    bins=bins(jt+1:end-jt);
    
    %Remove bsl
    app=find(bins<0);
    bsl_all=nanmean(psth_all(:,:,app),3);
    
    % bsl_std=nanstd(psth_all(:,:,app),[],3);
    % bsl_all=repmat(bsl_all,1,1,length(bins));
    % sig_th=repmat(2*bsl_std,1,1,length(bins));
    
    psth_all=psth_all-bsl_all;
    
   
    if compute_sig
        for j=1:size(psth_all,3)
            %[hs(1,j),ps(1,j)]=ttest2(dataPSTH.PPC.visual.PSTHmat_norm{1}(:,j),dataPSTH.PPC.visual.PSTHmat_norm{2}(:,j));
            %[hs(2,j),ps(2,j)]=ttest2(dataPSTH.PPC.visual.PSTHmat_norm{3}(:,j),dataPSTH.PPC.visual.PSTHmat_norm{4}(:,j));
            %[hr(1,j),pr(1,j)]=ttest2(dataPSTH.PPC.visual.PSTHmat_norm{1}(:,j),dataPSTH.PPC.visual.PSTHmat_norm{3}(:,j));
            %[hr(2,j),pr(2,j)]=ttest2(dataPSTH.PPC.visual.PSTHmat_norm{2}(:,j),dataPSTH.PPC.visual.PSTHmat_norm{4}(:,j));
            [hs(1,j),ps(1,j)]=realdata_bootstrap(psth_all(3,:,j),psth_all(1,:,j),1000,alpha);
            [hs(2,j),ps(2,j)]=realdata_bootstrap(psth_all(6,:,j),psth_all(4,:,j),1000,alpha);
            [hr(1,j),pr(1,j)]=realdata_bootstrap(psth_all(3,:,j),psth_all(6,:,j),1000,alpha);
            [hr(2,j),pr(2,j)]=realdata_bootstrap(psth_all(1,:,j),psth_all(4,:,j),1000,alpha);
        end
        hs=fdr_bh(ps,0.05,'pdep');
        hr=fdr_bh(pr,0.05,'pdep');
        hs_all{ak}=(hs(1,:) & hs(2,:));
        hr_all{ak}=(hr(1,:) & hr(2,:));
    end


    psth_mean=squeeze(nanmean(psth_all,2));
    psth_sem=squeeze(nanstd(psth_all,[],2))./sqrt(size(psth_all,2));
    
    
    colors=[1 0 0;
        0.75 0 0;
        0.5 0 0;
        0 1 0;
        0 0.75 0;
        0 0.5 0;
        0.5 0.5 0.5;
        0 0 0];
    
    %figure;
    subplot(length(area_list),1,ak);
    for i=1:size(psth_mean,1)
        plot(bins*1000,psth_mean(i,:),'Color',colors(i,:),'LineWidth',2);
        hold on;
    end
    for i=1:size(psth_mean,1)
        shadedErrorBar(bins*1000,psth_mean(i,:),psth_sem(i,:),{'-','LineWidth',2,'Color',colors(i,:)},1);
        hold on
    end
    legend('Hit, contrast 25%','Hit, contrast 50%','Hit, contrast 100%','Miss, contrast 25%','Miss, contrast 50%','Miss, contrast 100%');%,'Correct rejection','False Alarm');
    
    xlabel('Time (ms)');
    ylabel('Firing rate (Hz)');
    title(nArea,'FontSize',20);
    set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
    %grid;
    
end

%% Plot significance
if compute_sig
    figure;

    subplot(811);
    bar(bins*1000,hs_all{1},1,'k');
    set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
    title('V1 stimulus coding')
    subplot(813);
    bar(bins*1000,hs_all{2},1,'k');
    set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
    title('PPC stimulus coding')
    subplot(815);
    bar(bins*1000,hs_all{3},1,'k');
    set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
    title('CG1 stimulus coding')
    subplot(817);
    bar(bins*1000,hs_all{4},1,'k');
    set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
    title('PFC stimulus coding')

    subplot(812);
    bar(bins*1000,hr_all{1},1,'k');
    set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
    title('V1 report coding')
    subplot(814);
    bar(bins*1000,hr_all{2},1,'k');
    set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
    title('PPC report coding')
    subplot(816);
    bar(bins*1000,hr_all{3},1,'k');
    set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
    title('CG1 report coding')
    subplot(818);
    bar(bins*1000,hr_all{4},1,'k');
    set(gca,'FontSize',16,'XTick',-500:100:1500,'XGrid', 'on');
    title('PFC report coding')
end