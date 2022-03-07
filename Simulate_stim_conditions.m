clc;
clear all;
close all;
%% TRIAL MATRIX
%Create a matrix with 2400 trials. First column:left stimulus (which one of
%25 possibles gabors), Second column=right stimulus
    numberTrialPerCond=4; %number of trials, total will be 25*2N*2 (minimum 2 for the two conditions change or no change!!!) 1 leads to 100 trials
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

