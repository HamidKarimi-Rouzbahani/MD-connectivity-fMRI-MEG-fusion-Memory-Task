dclc;
clear all;
close all;
desired_trigger=1; % 1 = stim; 2 = cue; 3 = probe

potentialGabor=[1 2 3 4 5;6 7 8 9 10;11 12 13 14 15;16 17 18 19 20;21 22 23 24 25]; % rows = freq; columns = orient; 

for Subj_num=2:11
    %     clearvars -except OnsetResponse Subject NoC RDM_correct_4 RDM_error_4 NoC RDM_correct_8 RDM_error_8 NoC RDM_correct_16 RDM_error_16
    %% load data
    if desired_trigger==1
        str_name='Stim';
    elseif desired_trigger==2
        str_name='Cue';
    elseif desired_trigger==3
        str_name='Prb';
    end
    load(['MEG_',sprintf('P%d',Subj_num) '_',str_name,'_aligned.mat']);
    % Behavioural data columns:[
    % 1:'NumOfTrial,2:NumOfFixatedTrial,3:Fixated(1/0),4:LeftStim(1-25),5:RightStim(1-25),6:Original_Cue(1-4),7:ProbeDifferent(1/0),8:DefaultCueColors(1/0),9:CueColours(1-4),
    % 10:SpatialFreqLeft(1-5),11:OrientLeft(1-5),12:PhaseLeft(1-5),13:SpatialFreqRight(1-5),14:OrientRight(1-5),15:PhaseRight(1-5),16:ProbeSpatialFreqRand(1-5),17:ProbeOrientRand(1-5),18:ProbePhaseRand(1-5),19:RespondedDifferent?(1/0),20:RespondedCorrectly?(1/0),
    % 21:ReactionTime(s),22:TrialDuration(s),23:ExperimentTime(s),24:FixationDotDuration(s),25:StimDuration(s),26:PreCueDuration(s),27:CueDuration(s),28:PostCueDuration(s),29:FeedbackDuration(s),30:TriggerCodeStim,
    % 31:TriggerCodeCue,32:TriggerCodeProbe']
    
    % 4 and 5 stim; 6 cue; Stim/Memory: stim left: freq (10) orient (11) phase (12); stim right: freq (13) orient (14) phase (15)

%      conds=nchoosek([1:2],2);
         conds=nchoosek([1:5],2);

    accuracy=nan(size(conds,1),size(Data{1,1},2));
    for cond=1:size(conds,1)
        %% Stimulus 10:15; left and right freq, orient and phase       
        indCat1=find(Behavioural(:,14)==conds(cond,1));
        indCat2=find(Behavioural(:,14)==conds(cond,2) );          

        %% Cues (6 column) 1-4 %Left SF %Left Orient %Right SF %Right Orient       
        % cue sides
%         indCat1=find(Behavioural(:,6)<3);
%         indCat2=find(Behavioural(:,6)>2);
        
        % cue feature
%         indCat1=find(Behavioural(:,6)==1 | Behavioural(:,6)==3);
%         indCat2=find(Behavioural(:,6)==2 | Behavioural(:,6)==4);

        %% Memory %cues (6 columns) (1-4) Left SF %Left Orient %Right SF %Right Orient and Stim (10:15) 
%         indCat1=find(Behavioural(:,6)==4 & Behavioural(:,14)==conds(cond,1)); %
%         indCat2=find(Behavioural(:,6)==4 & Behavioural(:,14)==conds(cond,2)); %

       %% 
        
        clearvars ClassA ClassB
        for i=1:length(indCat1)
            ClassA(:,:,i)=Data{1,indCat1(i)};
        end
        for i=1:length(indCat2)
            ClassB(:,:,i)=Data{1,indCat2(i)};
        end
        
        time_sample=0;
        for time=1:5:size(ClassA,2)
            time_sample=time_sample+1;
            ClassA_time=squeeze(ClassA(:,time,:));
            ClassB_time=squeeze(ClassB(:,time,:));
            Classifier_Model = fitcsvm([ClassA_time';ClassB_time'],[ones(size(ClassA_time,2),1);zeros(size(ClassB_time,2),1)],'KernelFunction','linear','KernelScale','auto');
            accuracy(cond,time_sample)=1-kfoldLoss(crossval(Classifier_Model));
            [Subj_num cond time_sample]
        end
    end
        save(['Decoding_',sprintf('P%d',Subj_num) '_',str_name,'_aligned_stm_right_ornt.mat'],'accuracy','time_sample');
end
ccc
%% plotting stim/memory
clc;clear all;
conds=[3 4 7]; % Done: left freq (not very good;all and correct) and ornt (good); right freq done (not good), orientation not good
smoothing=10;   % Done: left ornt and freq (good); right ornt
desired_trigger=2; % 1 = stim; 2 = cue; 3 = probe
cor_all=''; % '_cor' for correct and '' for all trials

if desired_trigger==1
    str_name='Stim';
    stim_mem='stm'; % 'stm'
    trigger=20;
    data_time_samples=[-500:25:1500];
elseif desired_trigger==2
    str_name='Cue';
    stim_mem='mem'; %'mem', 'cue'
    trigger=40;
     data_time_samples=[-500:25:5500];   
elseif desired_trigger==3
    str_name='Prb';
end

s=0;
for Subj_num=2:11
    s=s+1;
    load(['Decoding_',sprintf('P%d',Subj_num) '_',str_name,'_aligned_',stim_mem,'_left_ornt',cor_all,'.mat']);
    accuracies(s,:)=smooth(nanmean(accuracy(conds,1:time_sample),1),smoothing);
end
% plot(nanmean(accuracies))
shadedErrorBar(data_time_samples,nanmean(accuracies),(nanstd(accuracies)./sqrt(10)),{'color',[0.1 0.1 0.8],'LineWidth',3,'lineStyle','-'},1);
% ccc
 hold on;
% load(['Decoding_',sprintf('P%d',Subj_num) '_',str_name,'_aligned_',stim_mem,'_left_ornt',cor_all,'.mat']);
% plot(smooth(nanmean(accuracy(conds,1:time_sample),1),smoothing))
% load(['Decoding_',sprintf('P%d',Subj_num) '_',str_name,'_aligned_',stim_mem,'_right_freq',cor_all,'.mat']);
% plot(smooth(nanmean(accuracy(conds,1:time_sample),1),smoothing))
% load(['Decoding_',sprintf('P%d',Subj_num) '_',str_name,'_aligned_',stim_mem,'_right_ornt',cor_all,'.mat']);
% plot(smooth(nanmean(accuracy(conds,1:time_sample),1),smoothing))
line([1 time_sample],[0.5 0.5])
line([trigger trigger],[0.4 .6])
legend LeftOrient LeftFreq RightOrient RightFreq
ccc
%% plotting cuecolor/side/feature
conds=1;
smoothing=5;
cor_all=''; % '_cor' for correct and '' for all trials
str_name='Cue';
stim_mem='cue_feat'; %'cue', 'side', 'feat

trigger=40;
s=0;
for Subj_num=2:11
    s=s+1;
    load(['Decoding_',sprintf('P%d',Subj_num) '_',str_name,'_aligned_',stim_mem,cor_all,'.mat']);
    accuracies(s,:)=smooth(nanmean(accuracy(conds,1:time_sample),1),smoothing);
end
plot(nanmean(accuracies))
hold on;
line([1 time_sample],[0.5 0.5])
line([trigger trigger],[0.25 .75])
%% plotting probe
conds=1:6;
smoothing=5;
cor_all='_cor'; % '_cor' for correct and '' for all trials
str_name='Prb';
stim_mem='probe';
trigger=20;
load(['Decoding_',sprintf('P%d',Subj_num) '_',str_name,'_aligned_',stim_mem,'_lf.mat']);
plot(smooth(nanmean(accuracy(conds,1:time_sample),1),smoothing))
hold on;
line([1 time_sample],[0.5 0.5])
line([trigger trigger],[0.25 .75])