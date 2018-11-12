% SCRIPT: neu3ca_rt
%--------------------------------------------------------------------------
% Copyright (C) Neu3CA Research Group, Eindhoven University of Technology
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 3
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,USA.
%
% Author: Stephan Heunis, <j.s.heunis@tue.nl>, 2018

%--------------------------------------------------------------------------
% DEFINITION
%--------------------------------------------------------------------------

% The neu3ca_rt script peforms preprocessing and simulates real-time fMRI
% analysis by running through all volumes of an fMRI dataset and performing
% real-time processing steps.

% The script requires the following setup/installs:
%   -   Matlab 2014a or later
%   -   SPM12: http://www.fil.ion.ucl.ac.uk/spm/
%   -   Experimental data (e.g. https://openneuro.org/datasets/ds000157)
%   -   File/directory locations (see below in initialisation section)

%%

%--------------------------------------------------------------------------
% INITIALISATION (INPUT REQUIRED)
%--------------------------------------------------------------------------
% STEP 1 - Setup/download data with known location and directory structure
% In this example, we use open data from OpenNeuro. See: https://openneuro.org/datasets/ds000157
% STEP 2 - Set required variables:
% Specify Matlab directory
matlab_dir             =   '/Users/jheunis/Documents/MATLAB';
% Specify SPM installation directory
spm_dir             =   '/Users/jheunis/Documents/MATLAB/spm12';
% Specify parent directory that contains all data
data_dir            =   '/Users/jheunis/Desktop/neu3carttest';
% Specify specific subject directory
sub_dir = [data_dir filesep 'sub-opennft']; % just an example
% Specify functional and structural filenames (these are all based on the example data)
% functional4D_fnz     =   [sub_dir filesep 'func' filesep 'sub-' sprintf('%02d',sub) '_task-passiveimageviewing_bold.nii.gz'];
% gunzip(functional4D_fnz);
% functional4D_fn     = [sub_dir filesep 'func' filesep 'sub-' sprintf('%02d',sub) '_task-passiveimageviewing_bold.nii'];
functional0_fn      =   [sub_dir filesep 'template_func.nii'];
functional1_fn      =   [sub_dir filesep 'fanon-0007-00001-000001-01.nii'];
roiL_fn = [sub_dir filesep 'lROI_1.nii'];
roiR_fn = [sub_dir filesep 'rROI_2.nii'];
% structural_fnz       =   [sub_dir filesep 'anat' filesep 'sub-' sprintf('%02d',sub) '_T1w.nii.gz'];
% gunzip(structural_fnz);
structural_fn = [sub_dir filesep 'structScan_PSC.nii'];

% Experiment and processing details (here, using example data from OpenNeuro)
Nt                =   155;
TR = 1.92;
timing_units = 'secs';
task_onsets = [0; 40.1; 77.2; 111.3; 143.3; 179.4; 218.5; 251.5; 299.6; 334.7; 374.8; 411.9; 445.9; 478.0; 514.1; 553.2];
task_durations = [24.1000; 24.06; 24.07; 24.06; 24.06; 24.07; 24.04; 24.06; 24.07; 24.10; 24.06; 24.06; 24.09; 24.09; 24.06; 24.07];
smoothing_kernel    = [6 6 6];
voxel_size = [3 3 3.75];
Nskip = 5;

%%

%--------------------------------------------------------------------------
% PREPROCESSING
%--------------------------------------------------------------------------
[d, fn, ext] = fileparts(structural_fn);
if exist([d filesep 'rc1' fn ext], 'file')
    % Just load file/variable names, don't redo preprocessing
    disp('Preprocessing already done - loading variables')
    preproc_data = struct;
    [d, fn, ext] = fileparts(structural_fn);
    preproc_data.forward_transformation = [d filesep 'y_' fn ext];
    preproc_data.inverse_transformation = [d filesep 'iy_' fn ext];
    preproc_data.gm_fn = [d filesep 'c1' fn ext];
    preproc_data.wm_fn = [d filesep 'c2' fn ext];
    preproc_data.csf_fn = [d filesep 'c3' fn ext];
    preproc_data.bone_fn = [d filesep 'c4' fn ext];
    preproc_data.soft_fn = [d filesep 'c5' fn ext];
    preproc_data.air_fn = [d filesep 'c6' fn ext];
    preproc_data.rstructural_fn = [d filesep 'r' fn ext];
    preproc_data.rgm_fn = [d filesep 'rc1' fn ext];
    preproc_data.rwm_fn = [d filesep 'rc2' fn ext];
    preproc_data.rcsf_fn = [d filesep 'rc3' fn ext];
    preproc_data.rbone_fn = [d filesep 'rc4' fn ext];
    preproc_data.rsoft_fn = [d filesep 'rc5' fn ext];
    preproc_data.rair_fn = [d filesep 'rc6' fn ext];
else
    preproc_data = neu3ca_rt_preRtPreProc(functional0_fn, structural_fn, spm_dir);
end

%%

%--------------------------------------------------------------------------
% DATA INITIALIZATION
%--------------------------------------------------------------------------
% Volume dimensions, and reference image
funcref_spm = spm_vol(functional0_fn);
funcref_3D  = spm_read_vols(funcref_spm);
[Ni, Nj, Nk] = size(funcref_3D);
N_vox = Ni*Nj*Nk;

% Whole brain masking
[GM_img_bin, WM_img_bin, CSF_img_bin] = neu3ca_rt_getSegments(preproc_data.rgm_fn, preproc_data.rwm_fn, preproc_data.rcsf_fn, 0.5);
I_GM = find(GM_img_bin);
I_WM = find(WM_img_bin);
I_CSF = find(CSF_img_bin);
mask_3D = GM_img_bin | WM_img_bin | CSF_img_bin;
I_mask = find(mask_3D);
N_maskvox = numel(I_mask);

% ROI masking
roiL_img = spm_read_vols(spm_vol(roiL_fn));
roiR_img = spm_read_vols(spm_vol(roiR_fn));
I_roiL = find(roiL_img(:)>0);
I_roiR = find(roiR_img(:)>0);

% Real-time realignment and reslicing parameter initialisation
% This step is based on the real-time implementation of realignment and
% reslicing for the OpenNFT toolbox: https://github.com/OpenNFT/OpenNFT
flagsSpmRealign = struct('quality',.9,'fwhm',5,'sep',4,...
    'interp',4,'wrap',[0 0 0],'rtm',0,'PW','','lkp',1:6);
flagsSpmReslice = struct('quality',.9,'fwhm',5,'sep',4,...
    'interp',4,'wrap',[0 0 0],'mask',1,'mean',0,'which', 2);
dicomInfoVox   = sqrt(sum((funcref_spm.mat(1:3,1:3)).^2));
fwhm = smoothing_kernel ./ dicomInfoVox;
A0=[];x1=[];x2=[];x3=[];wt=[];deg=[];b=[];
R(1,1).mat = funcref_spm.mat;
R(1,1).dim = funcref_spm.dim;
R(1,1).Vol = funcref_3D;
N_skip = 0;

% Run SPM 1st level design in order to get design matrix
% sess_params = struct;
% sess_params.timing_units = timing_units;
% sess_params.timing_RT = TR;
% sess_params.cond_name = 'Task';
% sess_params.cond_onset = task_onsets;
% sess_params.cond_duration = task_durations;
% if ~exist([sub_dir filesep 'SPM.mat'], 'file')
%     cd(sub_dir)
%     spm_specify1stlevel_jsh(data_dir, functional4D_fn, '', sess_params)
% end
Ndyn = Nt-Nskip;
load([sub_dir filesep 'SPM_PSC_Int_155.mat']);
convolved_task_design = SPM.xX.X(:,2); % convolved task time course regressor
drift_regressors = SPM.xX.K.X0; % cosine basis set for drift regressors
X_design = [convolved_task_design drift_regressors ones(Nt-Nskip,1)]; % design matrix, including task + drift + constat regressors

% Predefine some matrices/structures
F = zeros(Ni*Nj*Nk,Ndyn);
F_dyn_denoised = zeros(Ni*Nj*Nk, Ndyn);
TS = nan(2, Ndyn);
rF = F;
srF = F;
MP = zeros(Ndyn,6);
T = zeros(Ndyn,8);


%%
%--------------------------------------------------------------------------

% REAL-TIME ANALYSIS
%--------------------------------------------------------------------------

% Create figure for real-time visualisation, if required
% fig = figure;

for i = (Nskip+1):Nt
    % STEP 1: LOAD CURRENT VOLUME
    % Using SPM
    tic;
    disp(num2str(i))
    % Set filename of expected dynamic image
    dynamic_fn = [sub_dir filesep 'fanon-0007-' sprintf('%.5d', i) '-' sprintf('%.6d', i) '-01.nii'];
    % Load dynamic data into matrix
    f_spm = spm_vol(dynamic_fn);
    f = spm_read_vols(f_spm);
    F(:,i) = f(:);
    t1=toc;
    
    % STEP 2 + 3: REALIGN AND RESLICE TO REFERENCE VOLUME
    % Using OpenNFT functionality
    tic;
    % Realign to reference volume
    R(2,1).mat = f_spm.mat;
    R(2,1).dim = f_spm.dim;
    R(2,1).Vol = f;
    [R, A0, x1, x2, x3, wt, deg, b, nrIter] = spm_realign_rt(R, flagsSpmRealign, i, N_skip + 1, A0, x1, x2, x3, wt, deg, b);
    % Get motion correction parameters
    tmpMCParam = spm_imatrix(R(2,1).mat / R(1,1).mat);
    if (i == N_skip + 1)
        offsetMCParam = tmpMCParam(1:6);
    end
    MP(i,:) = tmpMCParam(1:6) - offsetMCParam;
    t2=toc; 
    tic;
    % Reslice to reference image grid
    rf = spm_reslice_rt(R, flagsSpmReslice);
    rF(:,i) = rf(:);
    t3=toc;
    tic;
    
    % STEP 4: SMOOTH REALIGNED VOLUME
    % Using OpenNFT functionality and SPM
    srf = zeros(Ni, Nj, Nk);
    gKernel = smoothing_kernel ./ dicomInfoVox;
    spm_smooth(rf, srf, gKernel);
    srF(:,i) = srf(:);
    t4=toc;
    
    % STEP 5: CUMULATIVE GLM DENOISING
    % Detrend masked data using a cumulative GLM. Mask uses combination of
    % GM, WM and CSF masks.
    tic;
    x = X_design(1:i, :);
    beta_func = x\srF(I_mask,1:i)'; % func = X*beta + e ==> beta = X\func ==> func_detrended = mp - X(i)*beta(i)
    f_dyn_denoised = srF(I_mask,1:i)' - x(:, 2:(end-1))*beta_func(2:(end-1), :); % remove effects of all regressors except constant and task
    f_dyn_denoised = f_dyn_denoised';
    F_dyn_denoised(I_mask,i) = f_dyn_denoised(:,i);
    t5=toc;
    
    % STEP 6: CUSTOM CODE
    tic;
    % add your custom code here
    TS(1,i) = mean(F_dyn_denoised(I_roiL,i));
    TS(2,i) = mean(F_dyn_denoised(I_roiR,i));
    t6=toc;
    
    % STEP 7: VISUALISATIONS  
    tic;
    % add your visualisation code here
    drawnow;
    t7=toc;
    
    % STEP 8: SUMMARIZE TIMING  
    t8 = t1+t2+t3+t4+t5+t6+t7;
    T(i,:) = [t1, t2, t3, t4, t5, t6, t7, t8];
    
end