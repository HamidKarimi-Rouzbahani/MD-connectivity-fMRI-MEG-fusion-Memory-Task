clc;
close all;
clear all;

ft_defaults;

toolboxdir = '/group/woolgar-lab/projects/Hamid/Projects/Vigilance/MEG/Analyses';
addpath([toolboxdir '/fieldtrip-lite-20190410']);
dataset_address='/megdata/cbu/mdconnectivity/';


%% start with an empty data struct and load MEG file
for Subj_num=[8]
    clearvars -except toolboxdir dataset_address Subj_num
    cfg = [];
    % Subj_num=10;%2,3,4,5,6,7,8,9,10,11
    
    desired_trigger=1; % 1 = stim; 2 = cue; 3 = probe
    
    photo_diode_based_triggers=0; % 1=photo_diode_based; 0=other_triggers
    
    if desired_trigger==1
        prestim=500;   % pre
        posstim=1500;   % post
        trigger_span=100;
    elseif desired_trigger==2
        prestim=500;   % pre
        posstim=5500;   % post
        trigger_span=500;
    elseif desired_trigger==3
        prestim=500;   % pre
        posstim=2000;   % post
        trigger_span=1000;
    end
    baseline_span=1:100;
    
    if Subj_num==1 % Dorian (photo-diodes stopped half-way through); triggers are inaccurate
        subject_code='0175';
        date='210921';
        beh_file='WorkMem_1_2021-09-21_13-27-42.csv';
    elseif Subj_num==2 % Jade
        subject_code='0180';
        date='210929';
        beh_file='WorkMem_2_2021-09-29_09-05-41.csv';
    elseif Subj_num==3 % Claire
        subject_code='0185';
        date='211006';
        beh_file='WorkMem_3_2021-10-06_12-24-42.csv';
    elseif Subj_num==4 % Hannah
        subject_code='0177';
        date='210922';
        beh_file='WorkMem_4_2021-09-22_09-43-57.csv';
    elseif Subj_num==5 % Nicola
        subject_code='0181';
        date='210930';
        beh_file='WorkMem_3_2021-09-30_13-12-34.csv';  %incorrectly saved as 3
    elseif Subj_num==6 % Bhavna; missed 2 trials on run 8
        subject_code='0193';
        date='211012';
        beh_file='WorkMem_8_2021-10-12_09-51-00.csv';  %incorrectly saved as 8
    elseif Subj_num==7 % Marco
        subject_code='0189';
        date='211008';
        beh_file='WorkMem_7_2021-10-08_09-32-45.csv';
    elseif Subj_num==8 % Runhao
        subject_code='0195';
        date='211013';
        beh_file='WorkMem_8_2021-10-13_09-10-16.csv';
    elseif Subj_num==9 % Jenessa
        subject_code='0191';
        date='211011';
        beh_file='WorkMem_9_2021-10-11_12-47-15.csv';
    elseif Subj_num==10 % Maria; might have missed trials on block 1 (but triggers are enough)
        subject_code='0194';
        date='211012';
        beh_file='WorkMem_9_2021-10-12_13-16-31.csv';  %incorrectly saved as 9
    elseif Subj_num==11 % Margreet
        subject_code='0190';
        date='211011';
        beh_file='WorkMem_11_2021-10-11_09-57-52.csv';
    elseif Subj_num==12 % Ruiyan
        subject_code='212';
        date='211216';
        beh_file='.csv';
    elseif Subj_num==13 % Andromachi; missed trials on block 7
        subject_code='211';
        date='211214';
        beh_file='.csv';
    end
    Behavioural=csvread(beh_file,2);
    % Columns: 'NumOfTrial,NumOfFixatedTrial,Fixated?(1/0),LeftStim(1-25),RightStim(1-25),Original_Cue(1-4),ProbeDifferent?(1/0),CuedStimAccoringToPrevTrial(1-4),DefaultCueColours?(1/0),CueColor(1-4),OrientLeft(1-5),PhaseLeft(1-5),SpatialFreqRight(1-5),OrientRight(1-5),PhaseRight(1-5),ProbeSpatialFreqRand(1-5),ProbeOrientRand(1-5),ProbePhaseRand(1-5),RespondedDifferent?(1/0),RespondedCorrectly?(1/0),ReactionTime(s),TrialDuration(s),ExperimentTime(s),FixationDotDuration(s),StimDuration(s),PreCueDuration(s),CueDuration(s),PostCueDuration(s),FeedbackDuration(s),TriggerCodeStim,TriggerCodeCue,TriggerCodeProbe
    
    Data=[];
    for run_num=1:8
        clearvars triggers trigger_start_code ft ft_downsampled cfg data triggers triggers_code triggers_sent
        %% tackling incostincencies in saving the data
        if Subj_num==1
            run_str=sprintf('%.2d',run_num);
        else
            run_str=num2str(run_num);
        end
        if Subj_num==1 || Subj_num==4
            trigger_shift=-12;
        else
            trigger_shift=0;
        end
        if Subj_num==9 || (Subj_num==10 && run_num~=5) || Subj_num==12
            Run_char='Run';
        else
            Run_char='run';
        end
        %% Loading the data
        cfg.dataset = fullfile([dataset_address,'meg21_',subject_code,'/',date,'/',Run_char,run_str,'_raw.fif']);
        data = ft_preprocessing(cfg);
        
        %% Triggers
        % EOG  = 1 and 2
        % Data = 3 to 308
        % MISC = 309 to 320
        % Triggers = 321 to 328 (STI001 to STI008)
        % triggers 1:125: Stimulus appearance is 1:25*5 (125 codes) based on the phase of left stim
        % triggers 126:129: Cues in order about left SF, left Orient, right SF, right Orient
        % triggers 131:140: Probes 130+(1/2)*(1:5)
        % Button presses = 317 to 324 (STI009 to STI0016)
        % Photodiode trigger = 330; (STI010)
        % Combined triggers and button = 337 (STI101)
        % 338 = Sys201
        if photo_diode_based_triggers==1
            photo_diode_channel=330+trigger_shift;
            triggers_code = +(data.trial{1}(photo_diode_channel,:)>4);
        else
            triggerchannels=[321:328]+trigger_shift;
            all_events_decimal_codes=bi2de([(data.trial{1}(triggerchannels,:)>4)'],'right-msb');
            if desired_trigger==1 % stimulus-aligned
                triggers_code = all_events_decimal_codes>0 & all_events_decimal_codes<126;
            elseif desired_trigger==2 % cue_aligned
                triggers_code = all_events_decimal_codes>125 & all_events_decimal_codes<130;
            elseif desired_trigger==3 % probe_aligned
                triggers_code = all_events_decimal_codes>130 & all_events_decimal_codes<141;
            end
        end
        
        % Now we want to find the start of each trigger signal.
        trigger_start_code(1)=0;
        for count = 2:length(triggers_code)
            %          if triggers_code(count)>triggers_code(count-1)
            if triggers_code(count)>triggers_code(count-1) && triggers_code(count+trigger_span*0.8)>triggers_code(count-1)
                trigger_start_code(count)=1;
            else
                trigger_start_code(count)=0;
            end
        end
        % In triggers_sent we have all the positions a trigger started.
        triggers_sent = find(trigger_start_code);
        
        if photo_diode_based_triggers==1
            c=0;
            for i=desired_trigger:3:length(triggers_sent)
                c=c+1;
                triggers(c)=triggers_sent(i);
            end
        else
            triggers=triggers_sent;
        end
        if length(triggers)~=50
            ccc
        end
        %% Define/epoch trials
        cfg=[];
        cfg.trl = [triggers-prestim;triggers+posstim;repmat(-prestim,1,length(triggers))]';
        
        % ft contains the cut data. One piece = one trial
        ft = ft_redefinetrial(cfg,data);
        
        %% removing unneeded channels
        for trl=1:length(ft.trial)
            ft.trial{1,trl}([1 2 309:end],:)=[];
        end
        c=0;
        for ch=3:308
            c=c+1;
            label2{c,1}=ft.label{ch,1};
        end
        ft.label=label2;
        ccc
        %% Filtering
        % figure;
        % subplot(2,1,1)
        % for t=1:10;plot(squeeze(nanmean(ft.trial{1,t}(10,:),1)));hold on;end
        % subplot(2,1,2)
        % x=ft.trial{1,10}(10,:);
        % periodogram(x,rectwin(length(x)),length(x))
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfreq   = [0.03 200];
        cfg.bpfilttype   = 'brickwall';
        cfg.bsfilter = 'yes';
        cfg.bsfreq   =[48 52];
        ft = ft_preprocessing(cfg,ft);
        
        cfg = [];
        cfg.resamplefs = 200;
        ft_downsampled = ft_resampledata(cfg,ft);
        ccc
        % figure;
        % subplot(2,1,1)
        % for t=1:10;plot(squeeze(nanmean(ft_downsampled.trial{1,t}(10,:),1)));hold on;end
        % subplot(2,1,2)
        % x=ft_downsampled.trial{1,10}(10,:);
        % periodogram(x,rectwin(length(x)),length(x))
        % ccc
        %% Removing the baseline
        for tt=1:length(ft_downsampled.trial)
            for ch=1:size(ft_downsampled.trial{1,1},1)
                ft_downsampled.trial{1,tt}(ch,:)=squeeze(ft_downsampled.trial{1,tt}(ch,:))-repmat(nanmean(squeeze(ft_downsampled.trial{1,tt}(ch,baseline_span)),2),[1 size(ft_downsampled.trial{1,tt}(ch,:),2)]);
            end
        end
        
        Data=horzcat(Data,ft_downsampled.trial);
        [run_num]
    end
    %% Save data
    if desired_trigger==1
        str_name='Stim';
    elseif desired_trigger==2
        str_name='Cue';
    elseif desired_trigger==3
        str_name='Prb';
    end
    % Remove missed MEG data fron behavioural file too
    if Subj_num==6
        Behavioural([351 352],:)=[];
    end
    save(['MEG_',sprintf('P%d',Subj_num) '_',str_name,'_aligned.mat'],'Data','Behavioural');
    [Subj_num]
end


