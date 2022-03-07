%% aligned to stim
clc;
clear all;
conds=[1:10]; % Done: left freq (not very good;all and correct) and ornt (good); right freq done (not good), orientation not good
smoothing=3;   % Done: left ornt and freq (good); right ornt
str_name='Stim';
stim_mem='stm'; % 'stm'
trigger=0;
data_time_samples=[-500:25:1500];
xticks=[-500 0 500 1000 1500];
time_sample=81;
figure;
gca = axes('Position',[0.11 0.15 0.775 0.5]);
s=0;
for Subj_num=2:11
    load(['Decoding_',sprintf('P%d',Subj_num) '_',str_name,'_aligned_',stim_mem,'_left_ornt.mat']);
    s=s+1;
    accuraciesTask(s,:)=smooth(nanmean(accuracy(conds,1:time_sample),1),smoothing);
    load(['Decoding_',sprintf('P%d',Subj_num) '_',str_name,'_aligned_',stim_mem,'_left_freq.mat']);
    s=s+1;
    accuraciesTask(s,:)=smooth(nanmean(accuracy(conds,1:time_sample),1),smoothing);
end
accuraciesTask=accuraciesTask*100;
Stimulus=shadedErrorBar(data_time_samples,nanmean(accuraciesTask),(nanstd(accuraciesTask)./sqrt(20)),{'color',[0.1 0.1 0.1],'LineWidth',3,'lineStyle','-'},1);
ylim_min=48;
ylim_max=62;
line([min(data_time_samples) max(data_time_samples)],[50 50],'LineWidth',1,'Color','k','LineStyle',':');
line([trigger trigger],[ylim_min ylim_max],'LineWidth',1,'Color','k','LineStyle',':')
set(gca,'FontSize',14,'LineWidth',1,'XTick',...
    xticks,'XTickLabel',...
    {xticks},'YTick',...
    [linspace(ylim_min,ylim_max,8)],'YTickLabel',{[linspace(ylim_min,ylim_max,8)]},...
    'XMinorTick','on','YMinorTick','off','ycolor','k','tickdir','out','xcolor','k');
ylim([ylim_min ylim_max])
xlim([min(xticks) max(xticks)])
box off
xlabel('Time relative to stimulus onset (ms)')
ylabel('Decoding accuracy (%)')
legend([Stimulus.mainLine],{'Stimulus content'},'box','off')

%% cue aligned
clc;
clear all;
conds=[1:10]; % Done: left freq (not very good;all and correct) and ornt (good); right freq done (not good), orientation not good
smoothing=5;   % Done: left ornt and freq (good); right ornt
time_sample=241;
str_name='Cue';
trigger=0;
data_time_samples=[-500:25:5500];
xticks=[-500 0 1000 2000 3000 4000 5000];
figure;
gca = axes('Position',[0.11 0.15 0.775 0.5]);

s=0;
ss=0;
for Subj_num=2:11
    load(['Decoding_',sprintf('P%d',Subj_num) '_',str_name,'_aligned_','cue_feat.mat']);
    s=s+1;
    accuraciesTask(s,:)=smooth(nanmean(accuracy(1,1:time_sample),1),smoothing);
    load(['Decoding_',sprintf('P%d',Subj_num) '_',str_name,'_aligned_','cue_side.mat']);
    s=s+1;
    accuraciesTask(s,:)=smooth(nanmean(accuracy(1,1:time_sample),1),smoothing);
    
%     load(['Decoding_',sprintf('P%d',Subj_num) '_',str_name,'_aligned_','mem_left_freq.mat']);
%     ss=ss+1;
%     accuraciesMem(ss,:)=smooth(nanmean(accuracy(conds,1:time_sample),1),smoothing);
    load(['Decoding_',sprintf('P%d',Subj_num) '_',str_name,'_aligned_','mem_left_ornt.mat']);
    ss=ss+1;
    accuraciesMem(ss,:)=smooth(nanmean(accuracy(conds,1:time_sample),1),smoothing); 
%     load(['Decoding_',sprintf('P%d',Subj_num) '_',str_name,'_aligned_','mem_right_freq.mat']);
%     ss=ss+1;
%     accuraciesMem(ss,:)=smooth(nanmean(accuracy(conds,1:time_sample),1),smoothing);
%     load(['Decoding_',sprintf('P%d',Subj_num) '_',str_name,'_aligned_','mem_right_ornt.mat']);
%     ss=ss+1;
%     accuraciesMem(ss,:)=smooth(nanmean(accuracy(conds,1:time_sample),1),smoothing); 
end
accuraciesTask=accuraciesTask*100;
accuraciesMem=accuraciesMem*100;
Task=shadedErrorBar(data_time_samples,nanmean(accuraciesTask),(nanstd(accuraciesTask)./sqrt(20)),{'color',[0.1 0.1 0.8],'LineWidth',3,'lineStyle','-'},1);
hold on;
Memory=shadedErrorBar(data_time_samples,nanmean(accuraciesMem),(nanstd(accuraciesMem)./sqrt(20)),{'color',[0.8 0.1 0.1],'LineWidth',3,'lineStyle','-'},1);
ylim_min=48;
ylim_max=56;
line([min(data_time_samples) max(data_time_samples)],[50 50],'LineWidth',1,'Color','k','LineStyle',':');
line([trigger trigger],[ylim_min ylim_max],'LineWidth',1,'Color','k','LineStyle',':')
set(gca,'FontSize',14,'LineWidth',1,'XTick',...
    xticks,'XTickLabel',...
    {xticks},'YTick',...
    [linspace(ylim_min,ylim_max,5)],'YTickLabel',{[linspace(ylim_min,ylim_max,5)]},...
    'XMinorTick','on','YMinorTick','off','ycolor','k','tickdir','out','xcolor','k');
ylim([ylim_min ylim_max])
xlim([min(xticks) max(xticks)])
box off
xlabel('Time relative to cue onset (ms)')
ylabel('Decoding accuracy (%)')
legend([Task.mainLine, Memory.mainLine],{'Task','Memory content'},'box','off')
