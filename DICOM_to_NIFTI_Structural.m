clc;
clear all;
% close all;


Main_analysis_directory='/group/woolgar-lab/projects/Hamid/Projects/MDconnectivity';

% addpath(fullfile([Main_analysis_directory,'/spm12/']));
% addpath(fmri_data_address,'/Analyses/spm12')
subject=25;
Run_folder=['/mridata/cbu/CBU220082_MR21003/20220202_094435/Series005_CBU_MPRAGE_32chn'];
Directory_for_saving=[Main_analysis_directory,'/MRI_data_analysis/'];
h=[];
for run=1:1
    clearvars dirs_chosen hdrs
    dirs=dir(Run_folder);
    c=0;
    ind=[];
    for itm=3:length(dirs)
        if strcmp(dirs(itm).name(end-3:end),'.dcm')
            c=c+1;
            ind=vertcat(ind,1);
            matlabbatch{1, 1}.spm.util.import.dicom.data{c,1}=[Run_folder,'/',dirs(itm).name];
            h(1).hdrs{c}=spm_dicom_headers([Run_folder,'/',dirs(itm).name],0);
        end
    end
    matlabbatch{1, 1}.spm.util.import.dicom.root='flat';

    %     if run<length(Run_folder)
    %         mkdir ([Directory_for_saving, sprintf('%s%0.3d', 'S', subject)], ['run',num2str(run)]);
    %         matlabbatch{1, run}.spm.util.import.dicom.outdir={[Directory_for_saving,sprintf('%s%0.3d', 'S', subject),'/run',num2str(run)]};
    %     else
    mkdir ([Directory_for_saving, sprintf('%s%0.3d', 'S', subject)], ['structurals']);
    matlabbatch{1, 1}.spm.util.import.dicom.outdir={[Directory_for_saving,sprintf('%s%0.3d', 'S', subject),'/structurals']};
    %     end
    
    matlabbatch{1, 1}.spm.util.import.dicom.protfilter='.*';
    matlabbatch{1, 1}.spm.util.import.dicom.convopts.format='nii';
    matlabbatch{1, 1}.spm.util.import.dicom.convopts.meta=0;
    matlabbatch{1, 1}.spm.util.import.dicom.convopts.icedims=0;
end
csub = sprintf('%s%0.3d', 'S', subject); %'S01'; % ID number for the subject we're going to analyse
csub = csub(1:4);

save(fullfile(Directory_for_saving,num2str(csub,'S%03d'),[num2str(csub,'S%03d'),'_Header_Files.mat']),'h');
spm_jobman('run',matlabbatch);
for run=1:1
    for c=1:length(h(run).hdrs)
        h(run).hdrs{1,c}=h(run).hdrs{1,c}{1,1};
    end
end
save(fullfile(Directory_for_saving,num2str(csub,'S%03d'),[num2str(csub,'S%03d'),'_Header_Files.mat']),'h');
done=1;






