function ccnl_rsa_searchlight(EXPT, rsa_idx, inds, radius, subbatch_size, subjects, corr_type, distance_metric)

    % Run searchlight RSA on a set of voxels. Useful to run in 
    % parallel in batches of e.g. 10000 voxels. Also see ccnl_rsa_view.m
    % Requires Kriegeskorte's RSA toolbox: http://www.mrc-cbu.cam.ac.uk/methods-and-resources/toolboxes/license/ (Nili et al., 2014)
    %
    % USAGE:
    %   ccnl_rsa_searchlight(EXPT, rsa_idx, inds, radius, [, subbatch_size, subjects, corr_type, distance_metric, radius])    
    %
    % EXAMPLE:
    %   ccnl_rsa_searchlight(exploration_expt(), 1, 1:10000);
    %
    % INPUT:
    %   EXPT - experiment structure
    %   rsa_idx - which RSA to use 
    %   inds - the batch: voxel indices in mask to use for searchlight centers
    %   radius - searchlight radius in voxels
    %   subbatch_size - (optional) partition inds in subbatches of this size. Higher values risk OOM, lower values reload the betas too often and make things slow. Defaults to 1000
    %   subjects - (optional) list of subjects
    %
    % Momchil Tomov, Oct 2018

    rng('shuffle'); % for parallelization

    if ~exist('subjects', 'var')
        subjects = 1:length(EXPT.subject);
    end

    if ~exist('subbatch_size', 'var')
        subbatch_size = 1000;
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

    % Need to split batch into sub-batches of voxels.
    % For each batch, extract the neural RDMs and correlate with the
    % behavioral RDMs.
    % This is necessary because we cannot fit all neural RDMs in memory (if subbatch_size = length(inds)).
    % On the other hand, if we compute RDMs 1 by 1 (subbatch_size = 1), we'll have to reload the betas every time -> it's SUPER slow. So this is a compromise.
    % Note that while in theory we could just use fewer voxels in inds (and have subbatch_size = length(inds)), then we would have to run too many jobs -> the cluster gods will be upset.
    %
    all_inds = sort(inds);
    for s = 1:subbatch_size:length(all_inds)
        e = min(length(all_inds), s + subbatch_size - 1);
        inds = all_inds(s:e);

        fprintf('\n------- subbatch %d-%d ------\n\n', s, e);

        % gen filename
        filename = sprintf('searchlight_%d-%d_c=%s_d=%s_%s.mat', min(inds), max(inds), corr_type, distance_metric, random_string());
        disp(fullfile(rsadir, filename));

        % get searchlight RDMs
        [Neural, cor] = ccnl_searchlight_rdms(EXPT, rsa_idx, inds, radius, subjects, distance_metric);

        % compute second-order correlations (similarity match)
        [C, H, T, P, all_subject_coef] = ccnl_match_rdms(Neural, Behavioral, control, corr_type);

        % save output 
        save(fullfile(rsadir, filename), 'cor', 'C', 'H', 'T', 'P', 'all_subject_coef', 'inds', 'rsa_idx', 'corr_type', 'distance_metric', '-v7.3');
    end

end


