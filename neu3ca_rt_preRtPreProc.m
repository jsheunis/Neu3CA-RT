% FUNCTION: neu3ca_rt_PreRtPreProc
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

% Function to complete pre-real-time preprocessing of structural and
% functional data from a single subject. Steps include coregistering
% structural image to initial functional image, segmenting the coregistered
% structural image into tissue types, and reslicing the segments to the
% functional resolution image grid. Makes use of spm12 batch routines.
% If spm12 batch parameters are not explicitly set, defaults are assumed.
%
% INPUT:
% funcional0_fn     - filename of initial pre-real-time 3D functional scan
% structural_fn     - filename of T1-weighted structural scan
% spm_dir           - SPM12 directory
%
% OUTPUT:
% output            - structure with filenames and data

%--------------------------------------------------------------------------

function output = neu3ca_rt_preRtPreProc(functional0_fn, structural_fn, spm_dir)

output = struct;

% STEP 1 -- Coregister structural image to first dynamic image (estimate)
disp('1 - Coregistering structural to functional image space...');
coreg_estimate = struct;
% Ref
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.ref = {functional0_fn};
% Source
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.source = {structural_fn};
% Other
% coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.other = {};
% Eoptions
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
% Run
cfg_util('run',coreg_estimate.matlabbatch);
disp('done');

% STEP 2 -- Segmentation of coregistered structural image into GM, WM, CSF, etc
% (with implicit warping to MNI space, saving forward and inverse transformations)
disp('2 - Segmenting coregistered structural image into GM, WM, CSF, etc...');
segmentation = struct;
% Channel
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.vols = {structural_fn};
% Tissue
for t = 1:6
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).tpm = {[spm_dir filesep 'tpm' filesep 'TPM.nii,' num2str(t)]};
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).ngaus = t-1;
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).native = [1 0];
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).warped = [0 0];
end
segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
% Warp
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.write=[1 1];
% Run
cfg_util('run',segmentation.matlabbatch);
% Save filenames
[d, fn, ext] = fileparts(structural_fn);
output.forward_transformation = [d filesep 'y_' fn ext];
output.inverse_transformation = [d filesep 'iy_' fn ext];
output.gm_fn = [d filesep 'c1' fn ext];
output.wm_fn = [d filesep 'c2' fn ext];
output.csf_fn = [d filesep 'c3' fn ext];
output.bone_fn = [d filesep 'c4' fn ext];
output.soft_fn = [d filesep 'c5' fn ext];
output.air_fn = [d filesep 'c6' fn ext];
disp('done');


% STEP 3 -- Reslice all to functional-resolution image grid
disp('3 - Reslice all generated images to functional-resolution image grid');
reslice = struct;
% Ref
reslice.matlabbatch{1}.spm.spatial.coreg.write.ref = {functional0_fn};
% Source
source_fns = {};
source_fns{1} = structural_fn;
[d, fn, ext] = fileparts(structural_fn);
for i = 2:7
    source_fns{i} = [d filesep 'c' num2str(i-1) fn ext];
end
reslice.matlabbatch{1}.spm.spatial.coreg.write.source = source_fns';
% Roptions
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
% Run
cfg_util('run',reslice.matlabbatch);
% Save filenames
[d, fn, ext] = fileparts(structural_fn);
output.rstructural_fn = [d filesep 'r' fn ext];
output.rgm_fn = [d filesep 'rc1' fn ext];
output.rwm_fn = [d filesep 'rc2' fn ext];
output.rcsf_fn = [d filesep 'rc3' fn ext];
output.rbone_fn = [d filesep 'rc4' fn ext];
output.rsoft_fn = [d filesep 'rc5' fn ext];
output.rair_fn = [d filesep 'rc6' fn ext];
disp('done');
