function ppp_inftod_fMRI_bids(fvol, svol, mout_dir, disc1, sigth, mvmt, pvec, temp)
% Preprocesses infant fMRI data
% For questions: Theodore_Turesky@gse.harvard.edu

% Inputs:

%   fvol     - 4D fMRI nifti
%   svol     - skull stripped anatomical nifti
%   mout_dir - desired output directory
%   disc1    - number of fMRI volumes to discard from beginning
%   sigth    - threshold: percent change mean global signal
%   mvmt     - threshold: inter-volume movement
%   pvec     - binary vector indicating desired processing steps
%              [T1 preproc, EPI preproc, ART, first-level stats]
%   temp     - template (e.g., UNC infant-1yr-withCerebellum.nii)




% ----------------- No need to change the parameters below -----------------
disc2 = 0; % #fMRI volumes to discard from end unless already discarded from run
[p_ppp, ~] = fileparts(which('ppp_inftod_fMRI_bids.m'));
dimTemp = fullfile(p_ppp, 'dimChangeTemp.mat');
realignTemp = fullfile(p_ppp, '01realignTemp.mat');
smoothTemp = fullfile(p_ppp, '05smoothTemp.mat');
ref = fullfile(p_ppp, 'BETA_Subject001_Condition001_Source011.nii');


% Build output structure
full_dirs = regexp(fvol,'/','split');
out_dir = fullfile(mout_dir, full_dirs{1, end-4});
f_dir = fullfile(out_dir, full_dirs{1, end-3}, full_dirs{1, end-2}, full_dirs{1, end-1});
s_dir = fullfile(out_dir, full_dirs{1, end-3}, 'anat');
an_dir = fullfile(out_dir, 'analysis');
eval(sprintf('!mkdir -p -- %s %s %s', f_dir, s_dir, an_dir));

% system(cmd);
% clearvars cmd

% Build files
mask = fullfile(s_dir, 'brainmask.nii.gz');
n4_img = fullfile(s_dir, 'N4_t1.nii');


% -- T1 Processing --
% Generate brain mask
if pvec(1,1) == 1
    
cmd = ['fslmaths ' svol ' -bin ' mask];

disp(['Running FSL Command: ' cmd]);
system(cmd);
clearvars cmd

% N4 Bias Correction
cmd = ['N4BiasFieldCorrection -d 3 -i ' svol ' -x ' mask ' -o ' n4_img ...
    ' -s 2 -c [100x100x100x100,0.0000000001] -b [200]'];

disp(['Running ANTs Command: ' cmd]);
system(cmd);
clearvars cmd

msk = 1;

if ~isfile(n4_img)
    cmd = ['N4BiasFieldCorrection -d 3 -i ' svol ' -o ' n4_img ...
    ' -s 2 -c [100x100x100x100,0.0000000001] -b [200]'];

disp(['Running ANTs Command: ' cmd]);
system(cmd);
clearvars cmd

msk = 0;

end

disp(['msk = ' msk]);

% Estimate structure-template registration
cmd = ['antsRegistrationSyN.sh -d 3 -f ' temp ' -m ' n4_img ...
    ' -o ' fullfile(s_dir, 'unc')];

disp(['Running ANTs Command: ' cmd]);
system(cmd);
clearvars cmd

end


% -- EPI Processing --
if pvec(1,2) == 1

% 4D to 3D convert
spm('defaults','fmri');
spm_jobman('initcfg');

load(dimTemp);
matlabbatch{1,1}.spm.util.split.outdir{1,1} = f_dir;
matlabbatch{1,1}.spm.util.split.vol{1,1} = fvol;

spm_jobman('run',matlabbatch);
clearvars matlabbatch

% Realign images
load(realignTemp);
fvols = dir(fullfile(f_dir, 'task-*.nii'));
n_fvols = size(fvols,1);

for i = (1+disc1):(n_fvols-disc2)
    matlabbatch{1,1}.spm.spatial.realign.estwrite.data{1,1}{i-disc1,1} = fullfile(f_dir, fvols(i).name);
end

spm_jobman('run',matlabbatch);
clearvars matlabbatch i

% Estimate structure-function coregistration
if msk
    cmd = ['antsRegistrationSyN.sh -d 3 -f ' n4_img ' -m ' fullfile(f_dir, 'meantask-*.nii') ' -x ' mask ...
        ' -o ' fullfile(s_dir, 'AF') ' -t a']; % made rigid for BB064: -t r
else
    cmd = ['antsRegistrationSyN.sh -d 3 -f ' n4_img ' -m ' fullfile(f_dir, 'meantask-*.nii')  ...
        ' -o ' fullfile(s_dir, 'AF') ' -t a']; % made rigid for BB064: -t r
end
disp(['Running ANTs Command: ' cmd]);
system(cmd);
clearvars cmd

% Apply transforms
rfvols = dir(fullfile(f_dir, 'rtask*.nii'));

for i = 1:size(rfvols,1)
    
    cmd = ['antsApplyTransforms -d 3 -i ' fullfile(f_dir, rfvols(i).name) ...
        ' -r ' ref ...
        ' -t ' fullfile(s_dir, 'unc1Warp.nii.gz') ...
        ' -t ' fullfile(s_dir, 'unc0GenericAffine.mat') ...
        ' -t ' fullfile(s_dir, 'AF0GenericAffine.mat') ...
        ' -o ' fullfile(f_dir, ['w' rfvols(i).name]) ];

    disp(['Running ANTs Command: ' cmd]);
    system(cmd);
    clearvars cmd

end

% Smooth images
load(smoothTemp);

for i = 1:size(rfvols,1)
    matlabbatch{1,1}.spm.spatial.smooth.data{i,1} = fullfile(f_dir, ['w' rfvols(i).name]);
end

spm_jobman('run',matlabbatch);
clearvars matlabbatch i

end

% -- Motion Artifact Detection --
if pvec(1,3) == 1
    
cfmiartrepair_general(20, 2, sigth, mvmt, f_dir, 'swrtask-') % change 2 to 3 to remove global mean signal artifacts
saveas(gcf, fullfile(f_dir, 'art.jpg'));
close(gcf)

end

% -- First-Level Stats --
if pvec(1,4) == 1
    
load(statsTemp);

matlabbatch{1,1}.spm.stats.fmri_spec.dir{1,1} = f_dir;

for i = 1:size(sfsm_vols,1)
    matlabbatch{1,1}.spm.stats.fmri_spec.sess(1).scans{i,1} = fullfile(f_dir, ['sw' rfvols(i).name]);
end

matlabbatch{1,1}.spm.stats.fmri_spec.sess(1).multi_reg{1,1} = fullfile(f_dir, ['multiReg_' sigth '_' mvmt '.mat']);

spm_jobman('run',matlabbatch);
clearvars matlabbatch i

end
