function [psth_all,bins]=eventRasters_2ndbump(s, nProbe, nArea, window, dt, smoothWin)
% Wrapper for evRastersGUI from /cortex-lab/spikes
%
% Inputs:
%   - s - is a struct created with loadSession
%   - nProbe - is an index of which probe to include (1-indexed, i.e. cannot
%   be zero)
%
% Example usage:
% >> s = loadSession('Muller_2017-01-07');
% >> eventRasters(s)

if nargin<2
    nProbe = 1; 
end

psth_all=[];
bins=[];

s.clusters.areaName=s.channels.brainLocation.allen_ontology(s.clusters.peakChannel,:);
appAreaID=zeros(size(s.clusters.areaName,1),1);
for i=1:size(s.clusters.areaName,1)
    if strcmp(s.clusters.areaName(i,:),nArea)
        appAreaID(i)=1;
    end
end

inclCID = find(s.clusters.probes==nProbe-1 & appAreaID)-1; 

if isempty(inclCID)
    return
end

inclSpikes = ismember(s.spikes.clusters, inclCID);

st = s.spikes.times(inclSpikes); 
clu = s.spikes.clusters(inclSpikes); 

lickTimes = s.licks.times;

contrastLeft = s.trials.visualStim_contrastLeft;
contrastRight = s.trials.visualStim_contrastRight;
feedback = s.trials.feedbackType;
choice = s.trials.response_choice;
%choice(choice==0) = 3; choice(choice==1) = 2; choice(choice==-1) = 1;

cweA = table(contrastLeft, contrastRight, feedback, choice); 

% appTrials{1}=find(cweA.contrastRight==0.25 & cweA.feedback==1 & choice==-1);
% appTrials{2}=find(cweA.contrastRight==0.5 & cweA.feedback==1 & choice==-1);
% appTrials{3}=find(cweA.contrastRight==1 & cweA.feedback==1 & choice==-1);
% appTrials{4}=find(cweA.contrastRight==0.25 & cweA.feedback==-1 & choice==0);
% appTrials{5}=find(cweA.contrastRight==0.5 & cweA.feedback==-1 & choice==0);
% appTrials{6}=find(cweA.contrastRight==1 & cweA.feedback==-1 & choice==0);
% appTrials{7}=find(cweA.contrastRight==0 & cweA.feedback==1 & choice==0);
% appTrials{8}=find(cweA.contrastRight==0 & cweA.feedback==-1 & choice==-1);


appTrials{1}=find(abs(cweA.contrastRight-cweA.contrastLeft)==0.25 & cweA.feedback==1);
appTrials{2}=find(abs(cweA.contrastRight-cweA.contrastLeft)==0.5 & cweA.feedback==1);
appTrials{3}=find(abs(cweA.contrastRight-cweA.contrastLeft)==1 & cweA.feedback==1);
appTrials{4}=find(abs(cweA.contrastRight-cweA.contrastLeft)==0.25 & cweA.feedback==-1 & choice==0);
appTrials{5}=find(abs(cweA.contrastRight-cweA.contrastLeft)==0.5 & cweA.feedback==-1 & choice==0);
appTrials{6}=find(abs(cweA.contrastRight-cweA.contrastLeft)==1 & cweA.feedback==-1 & choice==0);
%appTrials{7}=find(abs(cweA.contrastRight-cweA.contrastLeft)==0 & cweA.feedback==1 & choice==0);
%appTrials{8}=find(abs(cweA.contrastRight-cweA.contrastLeft)==0 & cweA.feedback==-1 & choice~=0);

for k=1:length(appTrials)
    
    stimOn = s.trials.visualStim_times(appTrials{k});
    beeps = s.trials.goCue_times(appTrials{k});
    feedbackTime = s.trials.feedback_times(appTrials{k});
    
    cwtA = table(stimOn, beeps, feedbackTime);
    
    unClu=unique(clu);
    for i=1:length(unClu)
        
        app=find(clu==unClu(i));
        try
            
            [psth, bins, rasterX, rasterY, spikeCounts] = psthRasterAndCounts(st(app), stimOn, window, dt);
            %         plot(bins,psth);
            %         hold on;
            
            if smoothWin>0
                smWin = myGaussWin(smoothWin, 1/dt);
                smWin(1:round(numel(smWin)/2)) = 0;
                smWin = smWin./sum(smWin);
                
                psth=conv(psth,smWin,'same');
            end
            
            psth_all(k,i,:)=psth/(dt*length(stimOn));
            
        catch error
            
            psth_all(k,i,:)=nan;
        end
    end
    
end