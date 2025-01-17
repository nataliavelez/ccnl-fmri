function [C, H, T, P, all_subject_ceof, Behavioral, Neural] = ccnl_rsa(EXPT, rsa_idx, roi_masks, subjects, corr_type, distance_metric)

    % RSA for given ROIs. Also see ccnl_rsa_searchlight.m
    % Requires Kriegeskorte's RSA toolbox: http://www.mrc-cbu.cam.ac.uk/methods-and-resources/toolboxes/license/ (Nili et al., 2014)
    %
    % USAGE:
    %   [C, H, T, P, all_subject_coef, Behavioral, Neural] = ccnl_rsa(EXPT, rsa_idx, roi_masks, [corr_type], [distance_metric])    
    %
    % EXAMPLE:
    %   [C, ~, T, P] = ccnl_rsa(exploration_expt(), 1, 'masks/hippocampus.nii');
    %
    % INPUT:
    %   EXPT - experiment structure
    %   rsa_idx - which RSA to use 
    %   roi_masks - mask name or cell array of mask names. Could pass 3D masks instead.
    %   subjects - (optional) list of subjects
    %   corr_type - (optional) type of (rank) correlation between RDMs (default: Spearman)
    %   distance_metric - (optional) neural distance metric (default: cosine)
    %
    % OUTPUT:
    %   C - [nROIs x nModels] matrix of (rank) correlations (averaged across subjects)
    %   H - [nROIs x nModels] matrix of hypothesis outcomes
    %   T - [nROIs x nModels] matrix of t-statistics
    %   P - [nROIs x nModels] matrix of p-values
    %   all_subject_ceof - [nROIs x nModels x nSubjects] matrix of (rank) correlations
    %   Behavioral - struct array with behavioral RDMs (see ccnl_behavioral_rdms.m)
    %   Neural - struct array with neural RDMs (see ccnl_roi_rdms.m)
    %
    % Momchil Tomov, Oct 2018

    if ~exist('subjects', 'var')
        subjects = 1:length(EXPT.subject);
    end

    if ~exist('corr_type', 'var')
        corr_type = 'Spearman';
    end

    if ~exist('distance_metric', 'var') 
        distance_metric = 'cosine';
    end

    % create rsa folder if none exists
    if ~isdir(EXPT.rsadir); mkdir(EXPT.rsadir); end
    rsadir = fullfile(EXPT.rsadir,['rsa',num2str(rsa_idx)]);
    if ~isdir(rsadir); mkdir(rsadir); end

    % get behavioral (model) RDMs
    [Behavioral, control] = ccnl_behavioral_rdms(EXPT, rsa_idx, subjects);

    % get searchlight RDMs
    [Neural] = ccnl_roi_rdms(EXPT, rsa_idx, roi_masks, subjects, false, distance_metric);

    % compute second-order correlations (similarity match)
    [C, H, T, P, all_subject_coef] = ccnl_match_rdms(Neural, Behavioral, control, corr_type);

