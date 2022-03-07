clc;
clear all;
close all;
Subjects={'1'};

for Subject=[1]  
    
    eyedata=Edf2Mat([Subjects{Subject},'.edf']);
    ccc
    stim=0;
    for count=1:length(eyedata.Events.Messages.info)
        if strcmp(eyedata.Events.Messages.info{1,count},'Fixation Start')
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
    for i=1:length(Trial_start)  
        Trials_data(i,1,:)=eyedata.Samples.posX(find(eyedata.Samples.time==Trial_start(i))-pre_stim:find(eyedata.Samples.time==Trial_start(i))+pst_stim,1);
        Trials_data(i,2,:)=eyedata.Samples.posY(find(eyedata.Samples.time==Trial_start(i))-pre_stim:find(eyedata.Samples.time==Trial_start(i))+pst_stim,1);

    end

    ccc
    [Subject]
    clearvars eyedata
end
