% FUNCTION: neu3ca_rt_getSegments
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

% This function constructs 3D binary images for GM, WM and CSF based on the
% relative value of GM/WM/CSF tissue probability maps per voxel. It assumes
% that the image filename parameters are for images that are in the same
% space (i.e. they match in size, voxel by voxel). If a threshold is
% specified (0<= threshold <=1), the images are first thresholded before
% the binary masks are calculated. For no threshold, specify zero.
%
% INPUT:
% gm_fn     - filename of gray matter image output from SPM segment routine
% wm_fn     - filename of white matter image output from SPM segment routine
% csf_fn    - filename of CSF image output from SPM segment routine
% threshold - value between 0 and 1 for thresholding GM/WM/CSF segment
%             images before calculation of the binary images
% 
% OUTPUT: 
% Three filenames for binary 3D images
%--------------------------------------------------------------------------
function [GM_img_bin, WM_img_bin, CSF_img_bin] = neu3ca_rt_getSegments(gm_fn, wm_fn, csf_fn, threshold)

GM_spm = spm_vol(gm_fn);
WM_spm = spm_vol(wm_fn);
CSF_spm = spm_vol(csf_fn);

GM_img = spm_read_vols(GM_spm);
WM_img = spm_read_vols(WM_spm);
CSF_img = spm_read_vols(CSF_spm);

if threshold ~= 0
    GM_img_thresh = GM_img;
    WM_img_thresh = WM_img;
    CSF_img_thresh = CSF_img;
    GM_img_thresh(GM_img < threshold) = 0;
    WM_img_thresh(WM_img < threshold) = 0;
    CSF_img_thresh(CSF_img < threshold) = 0;
    GM_img_bin = (GM_img_thresh >= WM_img_thresh) & (GM_img_thresh >= CSF_img_thresh) & (GM_img_thresh ~= 0);
    WM_img_bin = (WM_img_thresh > GM_img_thresh) & (WM_img_thresh >= CSF_img_thresh) & (WM_img_thresh ~= 0);
    CSF_img_bin = (CSF_img_thresh > GM_img_thresh) & (CSF_img_thresh > WM_img_thresh) & (CSF_img_thresh ~= 0);
else
    GM_img_bin = (GM_img >= WM_img) & (GM_img >= CSF_img) & (GM_img ~= 0);
    WM_img_bin = (WM_img > GM_img) & (WM_img >= CSF_img) & (WM_img ~= 0);
    CSF_img_bin = (CSF_img > GM_img) & (CSF_img > WM_img) & (CSF_img ~= 0);
end


