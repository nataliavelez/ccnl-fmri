function bic = ccnl_bic_bids(EXPT,model,masks, subjects)

% Compute Bayesian information criterion (BIC) for a model.
%
% USAGE: bic = ccnl_bic(EXPT,model,mask,[subjects])
%
% INPUTS:
%   EXPT - experiment structure
%   model - model label (str)
%   masks - cell array of ROI files
%   subjects (optional) - which subjects to analyze (default all subjects)
%
% OUTPUTS:
%   bic - [nSubjects x 1] vector of BIC values
%
% Sam Gershman, Oct 2016
% Edited by Natalia Velez, May 2022

if nargin < 4; subjects = 1:length(EXPT.subject); end

% % load mask
% if ischar(mask)
%     Vmask = spm_vol(mask);
%     Ymask = spm_read_vols(Vmask);
%     mask = Ymask~=0 & ~isnan(Ymask);
% end

bic = nan(length(subjects),1);
for s = 1:length(subjects)
    subj = subjects(s);

    % Load mask
    mask_path = masks{s};
    if ischar(mask_path)
        mask = load_mask(mask_path);
        modeldir = fullfile(EXPT.subject(s).modeldir, 'func', ...
            ['task-teaching_model-' model]);
        bic(s) = model_evidence(mask, modeldir);
    else
        warning('No ROI found, skipping sub-%02d', subj);
    end    
end
end

function mask = load_mask(path)
    if contains(path, '.gz')
        in_file = strrep(path, '.gz', '');
        if ~exist(in_file, 'file')
            warning('Mask path is compressed, unzipping: %s', path);
            gunzip(path);
        end
    else
        in_file = path;
    end

    Vmask = spm_vol(in_file);
    Ymask = spm_read_vols(Vmask);
    mask = Ymask~=0 & ~isnan(Ymask);
end
    
function sub_bic = model_evidence(mask, modeldir)
    % Load model results
    load(fullfile(modeldir,'SPM.mat'), 'SPM');

    [N,K] = size(SPM.xX.X); % number of datapoints and parameters
    V = spm_vol(fullfile(modeldir,'ResMS.nii'));    % residual variance image
    Y = spm_read_vols(V);
    
    % Compute model evidences
    if all(size(Y) == size(mask))
        sub_bic = N*nansum(log(Y(mask))) + K*log(N);
    else
        [mask_cor1, mask_cor2, mask_cor3] = ind2sub(size(mask),find(mask==1));
        res_cor = mni2cor(cor2mni([mask_cor1 mask_cor2 mask_cor3], Vmask.mat),V.mat);
        res_inds = sub2ind(size(Y),res_cor(:,1),res_cor(:,2),res_cor(:,3));
        disp(sum(isnan(Y(res_inds))))
        sub_bic = N*nansum(log(Y(res_inds))) + K*log(N);
    end
    
end