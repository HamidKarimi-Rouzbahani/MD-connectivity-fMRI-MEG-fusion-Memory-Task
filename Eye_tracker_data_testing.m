clc;
clear all;
close all;
% load('Order_of_blocks_for_Psych.mat');
% Within_trial_or_whole_block= 0;     % 0=whole block; 1=within trial
% TP(1),TN(2),FP(3),FN(4)

Subjects={'1'};

% Eye_blink_counter_Correct_Missed=nan(1,30,2);
% Eye_saccade_counter_Correct_Missed=nan(1,30,2);
% Eye_fixation_counter_Correct_Missed=nan(1,30,2);
% Eye_fix_dur_counter_Correct_Missed=nan(1,30,2);
% Corrects_Missed=nan(1,30,2);

for Subject=[1]
    
    %     Attended_unattended=double(reshape(Cued_color_block==green_red_color,[64 30])');
    %     Correct=double(reshape(Response_to_dots(:,3)<3,[64 30])');
    %     Missed=double(reshape(Response_to_dots(:,3)>3,[64 30])');
    %     Eye_blink_counter_Correct_Missed(Subject,:,:)=0;
    %     Eye_saccade_counter_Correct_Missed(Subject,:,:)=0;
    %     Eye_fixation_counter_Correct_Missed(Subject,:,:)=0;
    %     Eye_fix_dur_counter_Correct_Missed(Subject,:,:)=0;
    %     Corrects_Missed(Subject,:,:)=0;
    
    
    eyedata=Edf2Mat([Subjects{Subject},'.edf']);
    ccc
    stim=0;
    for count=1:length(eyedata.Events.Messages.info)
        if strcmp(eyedata.Events.Messages.info{1,count},'Fixation Start')
            % Gabor on; Gabor off; Response different/same; Isfixed= 1;
            % TRIAL_RESULT 0/1
            stim=stim+1;
            Trial_start(stim)=eyedata.Events.Messages.time(count);
        end
    end
    clearvars -except desired_trial_fixations Trial_start Eye_fixation_counter_Correct_Missed Eye_fix_dur_counter_Correct_Missed Corrects_Missed Correct Missed Eye_saccade_counter_Correct_Missed Eye_blink_counter_Correct_Missed Attended_unattended Active_monitoring Within_trial_or_whole_block Stim_appeared Trial_started eyedata Block_to_analyze Subject Subjects Block_to_analyze Subject Eye_saccade_counter Eye_blink_counter
    
    pre_stim=100;
    pst_stim=10000;
    
    % Removing artifactual data
    % Efix
    % Eblink
    % Esacc
%     desired_trial_fixations=zeros(length(-pre_stim:pst_stim));
    for i=1:length(Trial_start)  
        Trials_data(i,1,:)=eyedata.Samples.posX(find(eyedata.Samples.time==Trial_start(i))-pre_stim:find(eyedata.Samples.time==Trial_start(i))+pst_stim,1);
        Trials_data(i,2,:)=eyedata.Samples.posY(find(eyedata.Samples.time==Trial_start(i))-pre_stim:find(eyedata.Samples.time==Trial_start(i))+pst_stim,1);
    
%     blink_duration=mean(eyedata.Events.Eblink.duration(eyedata.Events.Eblink.duration<1000))*2.5;
%     for i=1:length(eyedata.Events.Eblink.start)
%         if find(eyedata.Samples.time==eyedata.Events.Eblink.start(i))>blink_duration
%             desired_trial_fixations(find(eyedata.Samples.time==eyedata.Events.Eblink.start(i))-blink_duration:find(eyedata.Samples.time==eyedata.Events.Eblink.end(i))+blink_duration)=0;
%         end
%     end
%     
%     saccade_duration=mean(eyedata.Events.Esacc.duration(eyedata.Events.Esacc.duration<1000))*2.5;
%     for i=1:length(eyedata.Events.Esacc.start)
%         if find(eyedata.Samples.time==eyedata.Events.Esacc.start(i))>saccade_duration
%             desired_trial_fixations(find(eyedata.Samples.time==eyedata.Events.Esacc.start(i))-saccade_duration:find(eyedata.Samples.time==eyedata.Events.Esacc.end(i))+saccade_duration)=0;
%         end
%     end
    end

    ccc
    [Subject]
    clearvars eyedata
end
