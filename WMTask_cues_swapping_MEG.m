function WMTask_cues_swapping_MEG
% written by Romain Quentin
% September 2015
% Screen('Preference', 'SkipSyncTests', 1);

% Experimenter='Romain Quentin';

% Revised by Hamid Karimi-Rouzbahani, Aug 2021
Experimenter='Hamid Karimi-Rouzbahani';
HideCursor;

%%
% Make sure this is running on OpenGL Psychtoolbox:
AssertOpenGL;

%% INPUTS

promptID='Subject ID (press enter for default): ';
defaultID='test';
promptEnt='Training before MEG ? (1:yes, 0:no, press enter for default (no))';
defaultEnt=false;


subjectID = input([promptID '[' defaultID '] '],'s');
doTraining = input([promptEnt '['  num2str(defaultEnt) '] ']);

if isempty(subjectID)
    subjectID=defaultID;
end
if isempty(doTraining)
    doTraining=defaultEnt;
end

if doTraining
    useOfEyelink=0; % 1 with Eyelink, 0 without eyelink
    whichComputer = 5; %(1 for MacbookRetina 15.4, 2 for testDesktop in office, 3 for MEG, 4 for Dell laptop, 5 for eye-tracking lab)
    useTrigger=0;
    button_box=0;
    photo_diode=0;
else
    useOfEyelink=1; % 1 with Eyelink, 0 without eyelink
    whichComputer = 3; %(1 for MacbookRetina 15.4, 2 for testDesktop in office, 3 for MEG)
    MegScale=120/71; %participants are seated 130 cm from the screen in the MEG, this number need to be apply to
    %all distance number when inside the MEG to keep one degree of visual field.
    useTrigger=1;
    button_box=1;
    photo_diode=1;
end


%% KEY INPUTS

KbName('UnifyKeyNames');
space=KbName('space');
escape = KbName('ESCAPE');
yes = KbName('y');%for redo the eyelink calibration
if button_box==0
    similar = KbName('j');
    different = KbName('k');
else
    similar = 'RB';
    similar_button_indx = 2;
    different = 'RY';
    different_button_indx = 3;
end
photo_diode_trigger='LY';
%% Make beep
beep=MakeBeep(500,0.1);

%% PREPARE RESPONSE DOCUMENT

path='C:\Users\test135user\Documents\Hamid_MD2021\Working_memory\PTB_task\data\';
cd C:\Users\test135user\Documents\Hamid_MD2021\Working_memory\PTB_task\data
addpath('C:\Users\test135user\Documents\Hamid_MD2021\Working_memory')

mkdir(subjectID)
cd(subjectID)
dt = datestr(now,'yyyy-mm-dd_HH-MM-SS');
UserID = num2str(subjectID);
subjectID = UserID;


if doTraining
    name = ['WorkMem_' num2str(subjectID) '_training_' dt '.csv'];
else
    name = ['WorkMem_' num2str(subjectID) '_' dt '.csv'];
end

fp=fopen(num2str(name), 'w');

% create file header
fprintf(fp,'Subject name (ID):, %s, \r\n', subjectID);

fprintf(fp,'NumOfTrial,NumOfFixatedTrial,Fixated?(1/0),LeftStim(1-25),RightStim(1-25),Original_Cue(1-4),ProbeDifferent?(1/0),CuedStimAccoringToPrevTrial(1-4),DefaultCueColours?(1/0),CueColor(1-4),OrientLeft(1-5),PhaseLeft(1-5),SpatialFreqRight(1-5),OrientRight(1-5),PhaseRight(1-5),ProbeSpatialFreqRand(1-5),ProbeOrientRand(1-5),ProbePhaseRand(1-5),RespondedDifferent?(1/0),RespondedCorrectly?(1/0),ReactionTime(s),TrialDuration(s),ExperimentTime(s),FixationDotDuration(s),StimDuration(s),PreCueDuration(s),CueDuration(s),PostCueDuration(s),FeedbackDuration(s),TriggerCodeStim,TriggerCodeCue,TriggerCodeProbe\r\n' );

%% Save the script in data folder
FileNameAndLocation=mfilename('fullpath');
FileName=mfilename;
newbackup=sprintf('%sbackup.m',strcat(path, num2str(subjectID), '\', FileName));
currentfile=strcat(FileNameAndLocation, '.m');
copyfile(currentfile,newbackup);

%% define name of EL file
% trial defaults (EyeLink)
dummymode=0;
if ~IsOctave
    commandwindow
else
    more off;
end

if IsOctave
    edfFile = 'DEMO';
else
    edfFile = subjectID;
end

trigger=1;
%% SCREEN AND STIMULI, FIXATION PARAMETERS

rGray = 127;
gGray = 127;
bGray = 127;

if whichComputer ==1
    widthCm=33.17; %heightCm=20.73;
elseif whichComputer ==2
    widthCm=33.5773; %heightCm=26.8618;
elseif whichComputer ==3 % for MEG experiment
    widthCm=57;
elseif whichComputer ==4 % Dell (36.5 * 22.9 cm) ; MSI = 34.4 cm width
    widthCm=36.5;
elseif whichComputer ==5 % Dell (37.5 * 30 cm)
    widthCm=37.5;
end

screens = max(Screen('Screens'));
maxPriorityLevel = MaxPriority(screens);
Priority(maxPriorityLevel);

[widthPx,heightPx]=WindowSize(screens);
pxPerCM = widthPx/widthCm;

if doTraining %no MEG
    wrectStim = 19;
    sizeDot=0.5;
    sizeCue=0.8;
    stimExcentricity = 8.5; %in cm
    fix_dist = 2.5; %distance acceptable for eye position from the center during the fixation in cm
else
    wrectStim = 19*MegScale;
    sizeDot=0.5*MegScale;
    sizeCue=0.8*MegScale;
    stimExcentricity = 8.5*MegScale;
    fix_dist = 2.5*MegScale; %distance acceptable for eye position from the center during the fixation in cm
end
wpixelrectStim = (wrectStim*pxPerCM);
rectStim = [0,0,wpixelrectStim,wpixelrectStim];
pixelStimExcentricity = (stimExcentricity*pxPerCM);
pixelsizeDot=(sizeDot*pxPerCM);
pixelsizeCue = (sizeCue*pxPerCM);
fix_distPX = fix_dist*pxPerCM;

%% OPEN GRAPHICS AND STIM POSITION
[myscreen, rect] = Screen('OpenWindow', screens, [rGray gGray bGray]);
[staticX,staticY] = RectCenter(rect);

% Center of the position of each stimuli

rectStim_central = CenterRectOnPoint(rectStim,staticX,staticY);


rectStim_resp = CenterRectOnPoint(rectStim,staticX,staticY);

%% TIMING OF THE TASKS
takingaccountscreenrefresh=0.005;
durationRestingScreen = 0.8 - takingaccountscreenrefresh;
durationCrossBefRand = 0.35 - takingaccountscreenrefresh; %we add rand*0.1 at each trial so the mean of duration dot is 0.4
durationGabor = 0.1 - takingaccountscreenrefresh;
durationPreCueBefRand = 1.5 - takingaccountscreenrefresh; %we add rand*0.1 at each trial so the mean of duration dot is 1
durationCue = 0.5 - takingaccountscreenrefresh;
durationPostCueBefRand = 5.0 - takingaccountscreenrefresh; %we add rand*0.1 at each trial so the mean of duration dot is 2
durationProbe = 2.0 - takingaccountscreenrefresh; %we add rand*0.1 at each trial so the mean of duration dot is 2
durationFeedback = 0.1 - takingaccountscreenrefresh;
%% TRIAL MATRIX
%Create a matrix with 2400 trials. First column:left stimulus (which one of
%25 possibles gabors), Second column=right stimulus
if doTraining
    numberTrialPerCond=2; %number of trials, total will be 25*2N*2 (minimum 2 for the two conditions change or no change!!!) 1 leads to 100 trials
else
    numberTrialPerCond=4; %number of trials, total will be 25*2N*2 (minimum 2 for the two conditions change or no change!!!) 1 leads to 100 trials
end
Num_Trls_per_Blk=25;

A(1:48)=1;
B(1:48)=2;
C(1:48)=3;
D(1:48)=4;
E(1:48)=5;
F(1:48)=6;
G(1:48)=7;
H(1:48)=8;
I(1:48)=9;
J(1:48)=10;
K(1:48)=11;
L(1:48)=12;
M(1:48)=13;
N(1:48)=14;
O(1:48)=15;
P(1:48)=16;
Q(1:48)=17;
R(1:48)=18;
S(1:48)=19;
T(1:48)=20;
U(1:48)=21;
V(1:48)=22;
W(1:48)=23;
X(1:48)=24;
Y(1:48)=25;

stimLeft = [A B C D E F G H I J K L M N O P Q R S T U V W X Y];
stimLeft=stimLeft';

stimRight=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25];

stimRight2=zeros(25,24);
for i=1:25
    stimRight2(i,:)=stimRight(stimRight~=i);
end

stimRight=(repmat(stimRight2,1,2))';
stimRight=stimRight(:);

%type of Cue

AA(1:24)=1; %Left SF
BB(1:24)=2; %Left Orient
CC(1:24)=3; %Right SF
DD(1:24)=4; %Right Orient

cueLeft2 = [AA BB]';
cueLeft = repmat(cueLeft2,25,1);

cueRight2 =[CC DD]';
cueRight = repmat(cueRight2,25,1);
cueMatrix=[cueLeft;cueRight];
stimMatrix_cueLeft=[stimLeft,stimRight];
stimMatrix_cueRight=[stimRight,stimLeft];

stimMatrix=[stimMatrix_cueLeft;stimMatrix_cueRight];

SMatrix=[stimMatrix,cueMatrix]; %matrix with all trials (2400 trials)

%if not all trials
t=[zeros(1,24-numberTrialPerCond),ones(1,numberTrialPerCond)]';
t2=zeros(24,100);
for i=1:100
    t2(:,i)=t(randperm(24),:);
end
t2=t2(:);

isUse = any(t2(:,1)==1,2);
SMatrix = SMatrix(isUse,:);
isChange=repmat([0 1],1,25*2*numberTrialPerCond)';
SMatrix=[SMatrix,isChange];

trialsArray = SMatrix(randperm(length(SMatrix)),:);
% Columns 1 and 2 left and right gabors (each 1-25)
% Column 3: 1-4 show cues of left SF, left Orient, right SF, right Orient
% Column 4: 0 shows same (no change) and 1 shows different (change)

% Initialization of the cueing system based on prev trial
totalTrials=length(SMatrix);

%% Gabor parameters
typeOfGabor=1;
contrast=1;
widthOfGrid=wpixelrectStim;
zoom=1;

if doTraining
    lambdaBase=pxPerCM;
else
    lambdaBase=pxPerCM*MegScale;
end

extensionOfGab = widthOfGrid/zoom;
halfExtensionOfGab = extensionOfGab / 2;
widthArray = (-halfExtensionOfGab) : halfExtensionOfGab;

% two-dimensional square grid
[x1 y1] = meshgrid(widthArray, widthArray);

%Calculate the tukeyfunction
alpha=0.3;
n=sqrt(((x1 .^ 2) + (y1 .^ 2)));
tukeyFunctionMaskmatrix=zeros(size(widthOfGrid));

for i=1:length(n)
    for j=1:length(n)
        if abs(n(i,j)) < (length(n)/2)*(1-alpha/2)
            tukeyFunctionMaskmatrix(i,j)=1;
        elseif abs(n(i,j)) >= (length(n)/2)*(1-alpha/2) && abs(n(i,j)) <= (length(n)/2)
            tukeyFunctionMaskmatrix(i,j)=1/2*(1+cos(pi*(2*n(i,j)/(alpha*length(n)/2)-2/alpha+1)));
        elseif abs(n(i,j)) > length(n)/2
            tukeyFunctionMaskmatrix(i,j)=0;
        end
    end
end

%calculate the signal decrease at the center
centerDecreasematrix=zeros(size(n));
for j=1:length(n)/2
    if j < length(n)/2*(1-alpha/2)
        centerDecreasematrix(:,j)=1;
    elseif j >= length(n)/2*(1-alpha/2)
        centerDecreasematrix(:,j)=1/2*(1+cos(pi*(2*j/(alpha*length(n)/2)-2/alpha+1)));
    end
end
for j=ceil(length(n)/2):length(n)
    j2=length(n)-(j-1);
    if j2 < length(n)/2*(1-alpha/2)
        centerDecreasematrix(:,j)=1;
    elseif j2  >= length(n)/2*(1-alpha/2)
        centerDecreasematrix(:,j)=1/2*(1+cos(pi*(2*j2/(alpha*length(n)/2)-2/alpha+1)));
    end
end

tukey_final= tukeyFunctionMaskmatrix.*centerDecreasematrix;


whichSF=[1 1.5 2.25 3.375 5.0625];
whichOrient=[-72 -36 0 36 72];


potentialGabor=[1 2 3 4 5;6 7 8 9 10;11 12 13 14 15;16 17 18 19 20;21 22 23 24 25];
phaseMatrix=[0 0.2 0.4 0.6 0.8];
potentialPhase=[1 2 3 4 5];

%% Serial Port
if useTrigger ==1 || button_box==1
    %     portSpec = FindSerialPort([], 1);
    %     [handle, errmsg] = IOPort('OpenSerialPort', portSpec);
    %
    %% set up MEG
    MEG = MEGSynchClass;
    MEG.SendTrigger(0); % make sure all triggers are off
end

%% Show examples to the participant /si demandé au début du script
DrawFormattedText(myscreen, 'Press buttons to check stimuli with different spatial frequencies', 'center', 'center', [0 0 0]);
Screen('Flip', myscreen);
WaitSecs(0.4);

if button_box
    MEG.WaitForButtonPress([],similar_button_indx);         % Wait for a button to be pressed
else
    KbWait;
end

ornt=36;
phs=1;
g=0;
for sf=whichSF
    g=g+1;
    gabors(:,:,:,g)=makeGabor(typeOfGabor, widthPx, widthCm, contrast, widthOfGrid, zoom, ornt,lambdaBase/sf,phs,tukey_final,0);
    gabor_final_tex = Screen('MakeTexture', myscreen, gabors(:,:,:,g));
    Screen('DrawTexture', myscreen, gabor_final_tex,[],rectStim_central);
    Screen('Flip', myscreen );
    WaitSecs(0.2);
    if button_box
        MEG.WaitForButtonPress([],similar_button_indx);         % Wait for any button to be pressed
    else
        KbWait;
    end
end


DrawFormattedText(myscreen, 'Now press buttons to check stimuli in different orientations', 'center', 'center', [0 0 0]);
Screen('Flip', myscreen);
WaitSecs(0.4);
if button_box
    MEG.WaitForButtonPress([],similar_button_indx);         % Wait for any button to be pressed
else
    KbWait;
end
sf=2.25;
for ornt=whichOrient
    g=g+1;
    gabors(:,:,:,g)=makeGabor(typeOfGabor, widthPx, widthCm, contrast, widthOfGrid, zoom, ornt,lambdaBase/sf,phs,tukey_final,0);
    gabor_final_tex = Screen('MakeTexture', myscreen, gabors(:,:,:,g));
    Screen('DrawTexture', myscreen, gabor_final_tex,[],rectStim_central);
    Screen('Flip', myscreen );
    WaitSecs(0.2);
    if button_box
        MEG.WaitForButtonPress([],similar_button_indx);         % Wait for any button to be pressed
    else
        KbWait;
    end
end

% prepare matrix for collecting trial
data = zeros(totalTrials,32);

%% EyeLink set-ups

if useOfEyelink
    
    % initializations
    el=EyelinkInitDefaults(myscreen);
    
    el.backgroundcolour = [rGray gGray bGray];
    el.calibrationtargetcolour =  0;
    EyelinkUpdateDefaults(el);
    
    %connection with eyetracker, opening file
    if ~EyelinkInit(dummymode)
        fprintf('Eyelink Init aborted.\n');
        cleanup(useOfEyelink, edfFile, trigger);
        return;
    end
    i = Eyelink('Openfile', edfFile);
    if i~=0
        fprintf('Cannot create EDF file ''%s'' ', edfFile);
        cleanup(useOfEyelink, edfFile, trigger);
        return;
    end
    if Eyelink('IsConnected')~=1 && ~dummymode
        cleanup(useOfEyelink, edfFile, trigger);
        return;
    end
    
    %configure eye tracker
    Eyelink('command', 'add_file_preamble_text ''Recorded by %s''', Experimenter);
    % This command is crucial to map the gaze positions from the tracker to
    % screen pixel positions to determine fixation
    Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, widthPx-1, heightPx-1);
    Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, widthPx-1, heightPx-1);
    Eyelink('command', 'calibration_type = HV5 '); %Eyelink('command', 'calibration_type = HV9');
    Eyelink('command', 'generate_default_targets = YES');
    % set parser (conservative saccade thresholds)
    Eyelink('command', 'saccade_velocity_threshold = 35');
    Eyelink('command', 'saccade_acceleration_threshold = 9500');
    % set EDF file contents
    % retrieve tracker version and tracker software version
    [v,vs] = Eyelink('GetTrackerVersion');
    fprintf('Running experiment on a ''%s'' tracker.\n', vs );
    vsn = regexp(vs,'\d','match');
    if v ==3 && str2double(vsn{1}) == 4 % if EL 1000 and tracker version 4.xx
        % remote mode possible add HTARGET ( head target)
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT,HTARGET');
        % set link data (used for gaze cursor)
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT,HTARGET');
    else
        Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT');
        % set link data (used for gaze cursor)
        Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
    end
    % allow to use the big button on the eyelink gamepad to accept the
    % calibration/drift correction target
    Eyelink('command', 'button_function 5 "accept_target_fixation"');
    
    %enter Eyetracker camera setup mode, calibration and validation
    EyelinkDoTrackerSetup(el);
    
end



%% START and RUN experiment
NbTrial=1;
SuccflTrial=1; % a trial in which fixation was correct and the subject used cue from prev trial
default_cue_colors=1;
%Press any key to proceed
WaitSecs(0.2);
% if button_box
%     MEG.WaitForButtonPress([],similar_button_indx);         % Wait for any button to be pressed
% else
%     KbWait;
% end
Screen('TextSize', myscreen, 30)
if button_box
    DrawFormattedText(myscreen, 'Press the left button to see the rules when ready!', 'center', 'center', [0 0 0]);
else
    DrawFormattedText(myscreen, 'Press "j" to see the rules when ready!', 'center', 'center', [0 0 0]);
end
WaitSecs(0.2);
Screen('Flip', myscreen)

% % press j or space to start
if button_box
    MEG.WaitForButtonPress([],similar_button_indx);         % Wait for any button to be pressed
else
    [~, ~, keyCode] = KbCheck(-3);
    while keyCode(similar) ==0
        [~, ~, keyCode] = KbCheck(-3);
    end
end


block=1;
% randomize cue colours for each subject
Cue_colors=[255 0 255;
    255 255 0;
    0 0 255;
    0 255 255]; % violet, yellow, blue, cyan
if isempty(str2num(subjectID))==0
    rand_cues=circshift([1:4],str2num(subjectID));
    Cue_colors=Cue_colors(rand_cues,:);
end
if photo_diode==1
    photo_diode_on_color=[255 255 255];
    photo_diode_off_color=[0 0 0];
    photoDiodeRectHeight=100; %height (in pixel) of the blanck rectangle at the top of the screen for the photoresistor
    photoDiodeRectwidth= 100; %width (in pixel) of the blanck rectangle at the top of the screen for the photoresistor
    photoDiodeXshift=50;
    photoDiodeYshift=50;
end
ListenChar(2);
while SuccflTrial <= totalTrials
    offset_x = round(staticX/4);
    offset_y = -round(staticY/2);
    length_line=100;
    gapX=round(staticX/9);
    offset_yp=round(staticY/16);
    VertSteps=round(staticY/4);
    x_offset_ins1=round(staticX/19);
    y_offset_ins2=round(staticY/12);
    
    
    
    if SuccflTrial <= totalTrials/2
        if (mod(SuccflTrial,Num_Trls_per_Blk)==1)
            DrawFormattedText(myscreen, 'Cue set 1', staticX-x_offset_ins1*2,y_offset_ins2, [255 255 255]);
            
            %cue left SF
            Screen('FillRect',myscreen,Cue_colors(1,:),[staticX-offset_x,staticY+offset_y,staticX-offset_x+length_line,staticY+offset_y+length_line]');
            DrawFormattedText(myscreen, 'Left Frequency', staticX-offset_x+length_line+gapX,staticY+offset_y+offset_yp*2, [80 80 80]);
            offset_y=offset_y+VertSteps;
            
            %cue left Orient
            Screen('FillRect',myscreen,Cue_colors(2,:),[staticX-offset_x,staticY+offset_y,staticX-offset_x+length_line,staticY+offset_y+length_line]');
            DrawFormattedText(myscreen, 'Left Orientation', staticX-offset_x+length_line+gapX,staticY+offset_y+offset_yp*2, [80 80 80]);
            offset_y=offset_y+VertSteps;
            
            %cue right SF
            Screen('FillRect',myscreen,Cue_colors(3,:),[staticX-offset_x,staticY+offset_y,staticX-offset_x+length_line,staticY+offset_y+length_line]');
            DrawFormattedText(myscreen, 'Right Frequency', staticX-offset_x+length_line+gapX,staticY+offset_y+offset_yp*2, [80 80 80]);
            offset_y=offset_y+VertSteps;
            
            %cue right Orient
            Screen('FillRect',myscreen,Cue_colors(4,:),[staticX-offset_x,staticY+offset_y,staticX-offset_x+length_line,staticY+offset_y+length_line]');
            DrawFormattedText(myscreen, 'Right Orientation', staticX-offset_x+length_line+gapX,staticY+offset_y+offset_yp*2, [80 80 80]);
            
            if button_box
                DrawFormattedText(myscreen, 'Press the left button twice when you are finished to start the block!', staticX-x_offset_ins1*10,staticY+offset_y+offset_yp*6, [255 255 255]);
            else
                DrawFormattedText(myscreen, 'Press "j" twice when you are finished to start the block!', staticX-x_offset_ins1*8,staticY+offset_y+offset_yp*6, [255 255 255]);
            end
            
            if photo_diode ==1
                Screen('FillRect',myscreen,photo_diode_off_color,[photoDiodeXshift photoDiodeYshift photoDiodeRectwidth photoDiodeRectHeight]); %show the photoresistor blank rectangle
            end
            
            Screen('Flip', myscreen );
            WaitSecs(0.2);
            if button_box
                MEG.WaitForButtonPress([],similar_button_indx);         % Wait for any button to be pressed
            else
                KbWait;
            end
            WaitSecs(1);
            if button_box
                MEG.WaitForButtonPress([],similar_button_indx);         % Wait for any button to be pressed
            else
                KbWait;
            end
            default_cue_colors=1;
        end
    elseif SuccflTrial > totalTrials/2
        if (mod(SuccflTrial,Num_Trls_per_Blk)==1)
            DrawFormattedText(myscreen, 'Cue set 2', staticX-x_offset_ins1*2,y_offset_ins2, [255 255 255]);
            if SuccflTrial == totalTrials/2+1
                DrawFormattedText(myscreen, 'Cues have changed !!!', staticX-x_offset_ins1*6,y_offset_ins2*3, [255 255 255]);
            end
            
            %cue left SF
            Screen('FillRect',myscreen,Cue_colors(4,:),[staticX-offset_x,staticY+offset_y,staticX-offset_x+length_line,staticY+offset_y+length_line]');
            DrawFormattedText(myscreen, 'Left Frequency', staticX-offset_x+length_line+gapX,staticY+offset_y+offset_yp*2, [80 80 80]);
            offset_y=offset_y+VertSteps;
            
            %cue left Orient
            Screen('FillRect',myscreen,Cue_colors(3,:),[staticX-offset_x,staticY+offset_y,staticX-offset_x+length_line,staticY+offset_y+length_line]');
            DrawFormattedText(myscreen, 'Left Orientation', staticX-offset_x+length_line+gapX,staticY+offset_y+offset_yp*2, [80 80 80]);
            offset_y=offset_y+VertSteps;
            
            %cue right SF
            Screen('FillRect',myscreen,Cue_colors(2,:),[staticX-offset_x,staticY+offset_y,staticX-offset_x+length_line,staticY+offset_y+length_line]');
            DrawFormattedText(myscreen, 'Right Frequency', staticX-offset_x+length_line+gapX,staticY+offset_y+offset_yp*2, [80 80 80]);
            offset_y=offset_y+VertSteps;
            
            %cue right Orient
            Screen('FillRect',myscreen,Cue_colors(1,:),[staticX-offset_x,staticY+offset_y,staticX-offset_x+length_line,staticY+offset_y+length_line]');
            DrawFormattedText(myscreen, 'Right Orientation', staticX-offset_x+length_line+gapX,staticY+offset_y+offset_yp*2, [80 80 80]);
            
            if button_box
                DrawFormattedText(myscreen, 'Press the left button twice when you are finished to start the block!', staticX-x_offset_ins1*10,staticY+offset_y+offset_yp*6, [255 255 255]);
            else
                DrawFormattedText(myscreen, 'Press "j" twice when you are finished to start the block!', staticX-x_offset_ins1*8,staticY+offset_y+offset_yp*6, [255 255 255]);
            end
            if photo_diode ==1
                Screen('FillRect',myscreen,photo_diode_off_color,[photoDiodeXshift photoDiodeYshift photoDiodeRectwidth photoDiodeRectHeight]); %show the photoresistor blank rectangle
            end
            Screen('Flip', myscreen );
            WaitSecs(0.2);
            if button_box
                MEG.WaitForButtonPress([],similar_button_indx);         % Wait for any button to be pressed
            else
                KbWait;
            end
            WaitSecs(1);
            if button_box
                MEG.WaitForButtonPress([],similar_button_indx);         % Wait for any button to be pressed
            else
                KbWait;
            end
            default_cue_colors=0;
        end
    end
    
    if useTrigger ==1
        MEG.SendTrigger(0); % reset triggers
    end
    tic
    disp( ['TrialCounter: ' num2str(NbTrial)]);
    disp( ['SuccflTrialCounter: ' num2str(SuccflTrial)]);
    disp( ['UnfixedTrialCounter: ' num2str(NbTrial-SuccflTrial)]);
    %% PRE-START THE TRIAL FOR EYELINK
    if useOfEyelink==1
        % let eyelink knows which trial
        WaitSecs(0.05);
        Eyelink('Message', 'TRIALNB %d',NbTrial);
        % This supplies the title at the bottom of the eyetracker display
        Eyelink('command', 'record_status_message "Correct to tried trials= %d/%d, Block= %d"', SuccflTrial, NbTrial, block);
        Eyelink('Command', 'set_idle_mode');
        Eyelink('Command', 'clear_screen %d', 0);
        Eyelink('Command', 'set_idle_mode');
        WaitSecs(0.05);
    end
    
    % initialize the acquisition of online eye position
    eyeCounter = 1;
    eyePoints = zeros(10000,3);
    if useOfEyelink==1
        Eyelink('StartRecording'); % start recording eyelink
        eye_used = Eyelink('EyeAvailable'); % 0 (left), 1 (right) or 2 (both)
        if eye_used == 2
            eye_used = 1;
        end
    end
    if useOfEyelink==1 && dummymode==0
        error=Eyelink('CheckRecording');
        if(error~=0)
            break;
        end
    end
    
    %% resting screen
    if photo_diode ==1
        Screen('FillRect',myscreen,photo_diode_off_color,[photoDiodeXshift photoDiodeYshift photoDiodeRectwidth photoDiodeRectHeight]); %show the photoresistor blank rectangle
    end
    start=Screen('Flip', myscreen);
    
    %Cross
    Screen('BlendFunction', myscreen, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    Screen('DrawDots', myscreen,[staticX staticY],pixelsizeDot,[80 80 80],[],1);
    if photo_diode ==1
        Screen('FillRect',myscreen,photo_diode_off_color,[photoDiodeXshift photoDiodeYshift photoDiodeRectwidth photoDiodeRectHeight]); %show the photoresistor blank rectangle
    end
    
    onset_fixcross = Screen('Flip', myscreen,start+durationRestingScreen);
    
    fixcrosstart=tic;
    
    
    if useOfEyelink==1
        Eyelink('Message', 'Fixation Start');
    end
    
    [sfLeft,orientLeft]=find(potentialGabor==(trialsArray(SuccflTrial,1)));
    [sfRight,orientRight]=find(potentialGabor==(trialsArray(SuccflTrial,2)));
    
    phaseNbLeft=randsample(potentialPhase,1);
    phaseNbRight=randsample(potentialPhase,1);
    phaseLeft=phaseMatrix(phaseNbLeft);
    phaseRight=phaseMatrix(phaseNbRight);
    gaborLeft=makeGabor(typeOfGabor, widthPx, widthCm, contrast, widthOfGrid, zoom, whichOrient(orientLeft),lambdaBase/(whichSF(sfLeft)),phaseLeft,tukey_final,0);
    gaborRight=makeGabor(typeOfGabor, widthPx, widthCm, contrast, widthOfGrid, zoom, whichOrient(orientRight),lambdaBase/(whichSF(sfRight)),phaseRight,tukey_final,0);
    
    gabor_final=zeros(length(gaborLeft),length(gaborLeft),3);
    gabor_final(:,ceil((length(gabor_final)/2)):end,:)=gaborRight(:,ceil(length(gabor_final)/2):end,:);
    gabor_final(:,1:ceil((length(gabor_final)/2)),:)=gaborLeft(:,1:ceil(length(gabor_final)/2),:);
    
    gabor_final_tex = Screen('MakeTexture', myscreen, gabor_final);
    
    Screen('DrawTexture', myscreen, gabor_final_tex,[],rectStim_central);
    
    Screen('DrawDots',myscreen,[staticX staticY],pixelsizeDot,[80 80 80],[],1);
    
    if photo_diode ==1
        Screen('FillRect',myscreen,photo_diode_on_color,[photoDiodeXshift photoDiodeYshift photoDiodeRectwidth photoDiodeRectHeight]); %show the photoresistor blank rectangle
    end
    
    durationCross=durationCrossBefRand+rand*0.1;
    onset_gabor = Screen('Flip', myscreen,onset_fixcross+durationCross);
    fixcrossTime=toc(fixcrosstart);
    gaborstart=tic;
    
    triggGabor=trialsArray(SuccflTrial,1)*phaseNbLeft;
    if useTrigger==1
        MEG.SendTrigger(triggGabor);
    end
    
    if useOfEyelink==1
        Eyelink('Message', 'Gabor ON');
    end
    
    % check eye position during gabor presentation
    ttestgabor=onset_gabor;
    while ttestgabor < onset_gabor+durationGabor-0.005
        if useOfEyelink==1
            if Eyelink('NewFloatSampleAvailable') > 0
                evt = Eyelink( 'NewestFloatSample');
                eyePoints(eyeCounter,1) = GetSecs-onset_gabor;
                eyePoints(eyeCounter,2) = evt.gx(eye_used+1);
                eyePoints(eyeCounter,3) = evt.gy(eye_used+1);
                eyeCounter = eyeCounter +1;
            end
        else
            eyePoints(eyeCounter,1) = GetSecs-onset_gabor;
            eyePoints(eyeCounter,2)=1;
            eyePoints(eyeCounter,3)=1;
            eyeCounter = eyeCounter +1;
        end
        ttestgabor=GetSecs;
    end
    
    Screen('DrawDots',myscreen,[staticX staticY],pixelsizeDot,[80 80 80],[],1);
    if photo_diode ==1
        Screen('FillRect',myscreen,photo_diode_off_color,[photoDiodeXshift photoDiodeYshift photoDiodeRectwidth photoDiodeRectHeight]); %show the photoresistor blank rectangle
    end
    offset_gabor = Screen('Flip', myscreen,onset_gabor+durationGabor);
    
    gaborTime=toc(gaborstart);
    precuestart=tic;
    
    if useTrigger==1
        MEG.SendTrigger(0); % reset triggers
    end
    
    if useOfEyelink==1
        Eyelink('Message', 'Gabor OFF');
    end
    
    durationPreCue=durationPreCueBefRand+rand*0.1;
    
    if SuccflTrial<=totalTrials/2
        if trialsArray(SuccflTrial,3)== 1 %cue left SF
            Screen('FillRect',myscreen,Cue_colors(1,:),round([(staticX-round(pixelsizeDot/2)),staticY-round(pixelsizeDot/2),(staticX+round(pixelsizeDot/2)),staticY+round(pixelsizeDot/2)]'));
            Cue_color_used=1;
        elseif trialsArray(SuccflTrial,3)== 2 %cue left orient
            Screen('FillRect',myscreen,Cue_colors(2,:),round([(staticX-round(pixelsizeDot/2)),staticY-round(pixelsizeDot/2),(staticX+round(pixelsizeDot/2)),staticY+round(pixelsizeDot/2)]'));
            Cue_color_used=2;
        elseif trialsArray(SuccflTrial,3)== 3 %cue right SF
            Screen('FillRect',myscreen,Cue_colors(3,:),round([(staticX-round(pixelsizeDot/2)),staticY-round(pixelsizeDot/2),(staticX+round(pixelsizeDot/2)),staticY+round(pixelsizeDot/2)]'));
            Cue_color_used=3;
        elseif trialsArray(SuccflTrial,3)== 4 %cue right orient
            Screen('FillRect',myscreen,Cue_colors(4,:),round([(staticX-round(pixelsizeDot/2)),staticY-round(pixelsizeDot/2),(staticX+round(pixelsizeDot/2)),staticY+round(pixelsizeDot/2)]'));
            Cue_color_used=4;
        end
    elseif SuccflTrial>totalTrials/2
        if trialsArray(SuccflTrial,3)== 1 %cue left SF
            Screen('FillRect',myscreen,Cue_colors(4,:),round([(staticX-round(pixelsizeDot/2)),staticY-round(pixelsizeDot/2),(staticX+round(pixelsizeDot/2)),staticY+round(pixelsizeDot/2)]'));
            Cue_color_used=4;
        elseif trialsArray(SuccflTrial,3)== 2 %cue left orient
            Screen('FillRect',myscreen,Cue_colors(3,:),round([(staticX-round(pixelsizeDot/2)),staticY-round(pixelsizeDot/2),(staticX+round(pixelsizeDot/2)),staticY+round(pixelsizeDot/2)]'));
            Cue_color_used=3;
        elseif trialsArray(SuccflTrial,3)== 3 %cue right SF
            Screen('FillRect',myscreen,Cue_colors(2,:),round([(staticX-round(pixelsizeDot/2)),staticY-round(pixelsizeDot/2),(staticX+round(pixelsizeDot/2)),staticY+round(pixelsizeDot/2)]'));
            Cue_color_used=2;
        elseif trialsArray(SuccflTrial,3)== 4 %cue right orient
            Screen('FillRect',myscreen,Cue_colors(1,:),round([(staticX-round(pixelsizeDot/2)),staticY-round(pixelsizeDot/2),(staticX+round(pixelsizeDot/2)),staticY+round(pixelsizeDot/2)]'));
            Cue_color_used=1;
        end
    end
    
    if photo_diode ==1
        Screen('FillRect',myscreen,photo_diode_on_color,[photoDiodeXshift photoDiodeYshift photoDiodeRectwidth photoDiodeRectHeight]); %show the photoresistor blank rectangle
    end
    
    onset_cue= Screen('Flip', myscreen,offset_gabor+durationPreCue);
    precueTime=toc(precuestart);
    cuestart=tic;
    
    triggCue=trialsArray(SuccflTrial,3)+125;
    if useTrigger==1
        MEG.SendTrigger(triggCue);
    end
    
    Screen('DrawDots',myscreen,[staticX staticY],pixelsizeDot,[80 80 80],[],1);
    if photo_diode ==1
        Screen('FillRect',myscreen,photo_diode_off_color,[photoDiodeXshift photoDiodeYshift photoDiodeRectwidth photoDiodeRectHeight]); %show the photoresistor blank rectangle
    end
    offset_cue = Screen('Flip', myscreen,onset_cue+durationCue);
    
    postcuestart=tic;
    cueTime=toc(cuestart);
    
    if useTrigger==1
        MEG.SendTrigger(0); % reset triggers
    end
    
    durationPostCue=durationPostCueBefRand+rand*0.1;
    
    randomOrient=randi(5,1);
    randomSF=randi(5,1);
    
    if trialsArray(SuccflTrial,3)== 1 || trialsArray(SuccflTrial,3)== 2 %Cue on left
        while randomOrient==orientLeft %orientLeft
            randomOrient=randi(5,1);
        end
        while randomSF==sfLeft % sfLeft
            randomSF=randi(5,1);
        end
    elseif trialsArray(SuccflTrial,3)== 3 || trialsArray(SuccflTrial,3)== 4 %Cue on Right
        while randomOrient==orientRight % orientRight
            randomOrient=randi(5,1);
        end
        while randomSF==sfRight % sfRight
            randomSF=randi(5,1);
        end
    end
    
    phaseNbResp=randsample(potentialPhase,1);
    phaseResp=phaseMatrix(phaseNbResp);
    
    if trialsArray(SuccflTrial,3)== 1 %cue on left SF
        if trialsArray(SuccflTrial,4)==0 % same
            gaborResp=makeGabor(typeOfGabor, widthPx, widthCm, contrast, widthOfGrid, zoom, whichOrient(randomOrient),lambdaBase/(whichSF(sfLeft)),phaseResp,tukey_final,0);
        elseif trialsArray(SuccflTrial,4)==1 % different
            gaborResp=makeGabor(typeOfGabor, widthPx, widthCm, contrast, widthOfGrid, zoom, whichOrient(randomOrient),lambdaBase/(whichSF(randomSF)),phaseResp,tukey_final,0);
        end
    elseif trialsArray(SuccflTrial,3)== 2 %cue on left Orient
        if trialsArray(SuccflTrial,4)==0 % same
            gaborResp=makeGabor(typeOfGabor, widthPx, widthCm, contrast, widthOfGrid, zoom, whichOrient(orientLeft),lambdaBase/(whichSF(randomSF)),phaseResp,tukey_final,0);
        elseif trialsArray(SuccflTrial,4)==1 % different
            gaborResp=makeGabor(typeOfGabor, widthPx, widthCm, contrast, widthOfGrid, zoom, whichOrient(randomOrient),lambdaBase/(whichSF(randomSF)),phaseResp,tukey_final,0);
        end
    elseif trialsArray(SuccflTrial,3)== 3 %cue on right SF
        if trialsArray(SuccflTrial,4)==0 % same
            gaborResp=makeGabor(typeOfGabor, widthPx, widthCm, contrast, widthOfGrid, zoom, whichOrient(randomOrient),lambdaBase/(whichSF(sfRight)),phaseResp,tukey_final,0);
        elseif trialsArray(SuccflTrial,4)==1 % different
            gaborResp=makeGabor(typeOfGabor, widthPx, widthCm, contrast, widthOfGrid, zoom, whichOrient(randomOrient),lambdaBase/(whichSF(randomSF)),phaseResp,tukey_final,0);
        end
    elseif trialsArray(SuccflTrial,3)== 4 %cue on right Orient
        if trialsArray(SuccflTrial,4)==0 % same
            gaborResp=makeGabor(typeOfGabor, widthPx, widthCm, contrast, widthOfGrid, zoom, whichOrient(orientRight),lambdaBase/(whichSF(randomSF)),phaseResp,tukey_final,0);
        elseif trialsArray(SuccflTrial,4)==1 % different
            gaborResp=makeGabor(typeOfGabor, widthPx, widthCm, contrast, widthOfGrid, zoom, whichOrient(randomOrient),lambdaBase/(whichSF(randomSF)),phaseResp,tukey_final,0);
        end
    end
    
    
    gaborResptex = Screen('MakeTexture', myscreen, gaborResp);
    Screen('DrawTexture', myscreen, gaborResptex,[],rectStim_resp);
    
    if photo_diode ==1
        Screen('FillRect',myscreen,photo_diode_on_color,[photoDiodeXshift photoDiodeYshift photoDiodeRectwidth photoDiodeRectHeight]); %show the photoresistor blank rectangle
    end
    
    Screen('DrawDots',myscreen,[staticX staticY],pixelsizeDot,[80 80 80],[],1);
    
    onset_probe=Screen('Flip', myscreen,offset_cue+durationPostCue);
    
    postcueTime=toc(postcuestart);
    
    triggProbe=130+(trialsArray(SuccflTrial,4)+1)*phaseNbResp;
    
    if useTrigger==1
        MEG.SendTrigger(triggProbe);
    end
    
    
    if button_box
        
        MEG.WaitForButtonPress(durationProbe,[similar_button_indx different_button_indx]);         % Wait for two buttons
        pressed_button = MEG.LastButtonPress; % record the key pressed
        [~, secs, ~] = KbCheck;
        if sum(strcmp(similar,pressed_button))
            if useOfEyelink==1
                Eyelink('Message', 'Response same');
            end
            response=0; %same
        elseif sum(strcmp(different,pressed_button))
            if useOfEyelink==1
                Eyelink('Message', 'Response different');
            end
            response=1; %different
        else
            if useOfEyelink==1
                Eyelink('Message', 'Missed');
            end
            response=nan; %miss
        end
        
    else
        
        [~, secs, keyCode] = KbCheck;
        while keyCode(similar) == 0 && keyCode(different) == 0 && keyCode(escape) == 0 && secs < onset_probe+durationProbe
            [~, secs, keyCode] = KbCheck;
            
            if keyCode(escape)
                ListenChar(0);
                Screen('CloseAll');
                sca
            elseif keyCode(similar)
                if useOfEyelink==1
                    Eyelink('Message', 'Response same');
                end
                response=0; %same
            elseif keyCode(different)
                if useOfEyelink==1
                    Eyelink('Message', 'Response different');
                end
                response=1; %different
            else
                if useOfEyelink==1
                    Eyelink('Message', 'Missed');
                end
                response=nan; %miss
            end
        end
        
    end
    
    if useTrigger==1
        MEG.SendTrigger(0); % reset triggers
    end
    
    if (trialsArray(SuccflTrial,4)== 0 && response == 0) || (trialsArray(SuccflTrial,4)== 1 && response == 1)
        isCorrect=1;
    elseif (trialsArray(SuccflTrial,4)== 1 && response == 0) || (trialsArray(SuccflTrial,4)== 0 && response == 1)
        isCorrect=0;
    elseif isnan(response)
        isCorrect=nan;
    end
    
    %% is Fixed?
    if useOfEyelink==1
        isFixed = isInCircle(eyePoints, fix_distPX, widthPx, heightPx);
    else
        isFixed = 1;
    end
    
    
    if useOfEyelink==1
        WaitSecs(0.001);
        Eyelink('Message', 'IsFixed= %d', isFixed);
        Eyelink('Message', 'TRIAL_RESULT 0'); % for ending a trial in data viewer
    end
    
    if isFixed==0
        if button_box==1
            DrawFormattedText(myscreen, 'You did not fixate the fixation dot, press the left button to continue', 'center', 'center', [0 0 0]);
        else
            DrawFormattedText(myscreen, 'You did not fixate the fixation dot, press "j" to continue', 'center', 'center', [0 0 0]);
        end
        if photo_diode ==1
            Screen('FillRect',myscreen,photo_diode_off_color,[photoDiodeXshift photoDiodeYshift photoDiodeRectwidth photoDiodeRectHeight]); %show the photoresistor blank rectangle
        end
        Screen('Flip', myscreen);
        WaitSecs(0.4);
        
        if button_box==1
            MEG.WaitForButtonPress([],similar_button_indx);         % Wait for any button to be pressed
        else
            [~, ~, keyCode] = KbCheck(-3);
            while keyCode(similar) == 0
                [~, ~, keyCode] = KbCheck(-3);
            end
        end
    end
    
    %% Feedback
    if isFixed == 1
        if isCorrect==1 % correct
            Screen('DrawDots',myscreen,[staticX staticY],pixelsizeDot,[0 255 0],[],1);
        elseif isCorrect==0 % incorrect
            Screen('DrawDots',myscreen,[staticX staticY],pixelsizeDot,[255 0 0],[],1);
        elseif isnan(isCorrect) % miss
            DrawFormattedText(myscreen, 'Missed !', 'center', 'center', [255 0 0]);
        end
    end
    if photo_diode ==1
        Screen('FillRect',myscreen,photo_diode_off_color,[photoDiodeXshift photoDiodeYshift photoDiodeRectwidth photoDiodeRectHeight]); %show the photoresistor blank rectangle
    end
    %onset feedback
    onset_feedback=Screen('Flip', myscreen);
    feedbackstart=tic;
    
    % get reaction time
    rt= secs- onset_probe;
    
    if photo_diode ==1
        Screen('FillRect',myscreen,photo_diode_off_color,[photoDiodeXshift photoDiodeYshift photoDiodeRectwidth photoDiodeRectHeight]); %show the photoresistor blank rectangle
    end
    %offset feedback
    Screen('Flip',myscreen,onset_feedback+durationFeedback);
    feedbackTime=toc(feedbackstart);
    
    trialTime=toc;
    
    runningTime=GetSecs;
    data(NbTrial,1)  = NbTrial;
    data(NbTrial,2) = SuccflTrial;
    data(NbTrial,3) = isFixed;
    data(NbTrial,4:7)  = trialsArray(SuccflTrial,:); % left and right stimuli,
    data(NbTrial,8)  = default_cue_colors; % default_cue_color
    data(NbTrial,9)  = Cue_color_used;
    data(NbTrial,10) = sfLeft;
    data(NbTrial,11) = orientLeft;
    data(NbTrial,12) = phaseLeft;
    data(NbTrial,13) = sfRight;
    data(NbTrial,14) = orientRight;
    data(NbTrial,15) = phaseRight;
    data(NbTrial,16)  = randomSF;
    data(NbTrial,17)  = randomOrient;
    data(NbTrial,18)  = phaseResp;
    data(NbTrial,19)  = response;
    data(NbTrial,20)  = isCorrect;
    data(NbTrial,21)  = rt;
    data(NbTrial,22)  = trialTime;
    data(NbTrial,23)  = runningTime;
    data(NbTrial,24)  = fixcrossTime;
    data(NbTrial,25)  = gaborTime;
    data(NbTrial,26)  = precueTime;
    data(NbTrial,27)  = cueTime;
    data(NbTrial,28)  = postcueTime;
    data(NbTrial,29)  = feedbackTime;
    data(NbTrial,30)  = triggGabor;
    data(NbTrial,31)  = triggCue;
    data(NbTrial,32)  = triggProbe;
    data_discription{1}={'NumOfTrial','NumOfFixatedTrial','Fixated?(1,0)','LeftStim(1:25)','RightStim(1:25)',...
        'Original_Cue(1:4)','ProbeDifferent?(1,0)','CuedStimAccoringToPrevTrial(1:4)','DefaultCueColours?(1,0)','CueColor(1-4)','OrientLeft(1:5)',...
        'PhaseLeft(1:5)','SpatialFreqRight(1:5)','OrientRight(1:5)','PhaseRight(1:5)','ProbeSpatialFreqRand(1:5)',...
        'ProbeOrientRand(1:5)','ProbePhaseRand(1:5)','RespondedDifferent?(1,0)','RespondedCorrectly?(1,0)',...
        'ReactionTime(s)','TrialDuration(s)','ExperimentTime(s)','FixationDotDuration(s)','StimDuration(s)',...
        'PreCueDuration(s)','CueDuration(s)','PostCueDuration(s)','FeedbackDuration(s)','TriggerCodeStim','TriggerCodeCue','TriggerCodeProbe'};
    fprintf(fp, '%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%f,%d,%d,%f,%d,%d,%f,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d\r\n', data(NbTrial,:));
    
    if isFixed
        if mod(SuccflTrial,Num_Trls_per_Blk)==0
            Snd('play',beep)
            
            corResponse = round(nanmean(data((NbTrial-(Num_Trls_per_Blk-1)):NbTrial,20))*100);
            disp( ['Performance',num2str(corResponse),' %: ']);
            
            if doTraining
                DrawFormattedText(myscreen, 'Congratulation! End of a block', 'center', (staticY - 25), [0 0 0]);
                DrawFormattedText(myscreen, 'Please wait for instructions.', 'center', (staticY + 25), [0 0 0]);
                DrawFormattedText(myscreen, ['Accuracy = ' num2str(corResponse) '%'], 'center', (staticY + 200), [0 0 0]);
            else
                DrawFormattedText(myscreen, 'Congratulation! End of a block', 'center', (staticY - 25), [0 0 0]);
                DrawFormattedText(myscreen, 'Localising head... Please wait for instructions ', 'center', (staticY + 100), [0 0 0]);
                DrawFormattedText(myscreen, ['Accuracy = ' num2str(corResponse) '%'], 'center', (staticY + 200), [0 0 0]);
            end
            if photo_diode ==1
                Screen('FillRect',myscreen,photo_diode_off_color,[photoDiodeXshift photoDiodeYshift photoDiodeRectwidth photoDiodeRectHeight]); %show the photoresistor blank rectangle
            end
            Screen('Flip', myscreen);
            if button_box==1
                MEG.WaitForButtonPress([],similar_button_indx);         % Wait for any button to be pressed
            else
                [~, ~, keyCode] = KbCheck(-3);
                while keyCode(similar) == 0
                    [~, ~, keyCode] = KbCheck(-3);
                end
            end
            block=block+1;
        end
    end
    
    %% reshuffling the rest of undone trials
    NbTrial = NbTrial +1;
    if isFixed
        SuccflTrial=SuccflTrial+1;
    else
        newtrialsArray=trialsArray(SuccflTrial:totalTrials,:);
        newtrialsArray=newtrialsArray(randperm(totalTrials-(SuccflTrial-1)),:);
        trialsArray(SuccflTrial:totalTrials,:)=newtrialsArray;
    end
end
DrawFormattedText(myscreen, 'Congratulations! You have completed this task ', 'center', 'center', [0 0 0]);
if photo_diode ==1
    Screen('FillRect',myscreen,photo_diode_off_color,[photoDiodeXshift photoDiodeYshift photoDiodeRectwidth photoDiodeRectHeight]); %show the photoresistor blank rectangle
end
Screen('Flip', myscreen);
WaitSecs(3);

%% CLOSING

Screen('CloseAll');
if useOfEyelink==1
    cleanup(useOfEyelink, edfFile, trigger);
end
if useTrigger ==1
    MEG.SendTrigger(0); % reset triggers
    MEG.delete;% stop MEG from limiting button presses
end

end

function cleanup(useOfEyelink, edfFile, ~)

if useOfEyelink==1
    Eyelink('Command' , 'set_idle_mode');
    WaitSecs(0.5);
    Eyelink('CloseFile');
    %download data file
    try
        fprintf('Receiving data file ''%s''\n', edfFile );
        status=Eyelink('ReceiveFile');
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
        end
        if 2==exist(edfFile, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
        end
    catch
        fprintf('Problem receiving data file ''%s''\n', edfFile );
    end
    Eyelink('Shutdown');
end

Screen('CloseAll');

ListenChar(0);
ShowCursor;
end