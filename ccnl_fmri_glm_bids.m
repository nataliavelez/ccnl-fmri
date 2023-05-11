function ccnl_fmri_glm_bids(EXPT,task,model,subjects,fake)
    
    % Estimate subject-level GLM on BIDS-formatted dataset
    % Compared to ccnl_fmri_glm, this script expects EXPT to have a
    % additional inputs, enabling us to run files generated through other
    % processing pipelines:
    %
    % task: task label
    % outdir: where to save outputs
    % functionals: paths to functional files
    % structural: path to normalized structural file
    % model: name of GLM
    % mask: path to mask file
    
    %
    % If fake = true, only creates the SPM.mat files, without actually 
    % running the GLM. Use for quick sanity checks
    %
    % USAGE: ccnl_fmri_glm_bids(EXPT,model,[subjects])
    
    cdir = pwd;
    if nargin < 3; subjects = 1:length(EXPT.subject); end
    if nargin < 4; fake = false; end

    % create models folder if none exists
    outdir = EXPT.modeldir;
    if ~isfolder(outdir)
        mkdir(outdir)
    end
    
    % if ~isdir(fullfile(EXPT.modeldir,['model',num2str(model)])); mkdir(fullfile(EXPT.modeldir,['model',num2str(model)])); end
    
    % generic design specification
    def = spm_get_defaults;
    job.timing.RT = EXPT.TR;
    job.timing.units = 'secs';
    job.timing.fmri_t = def.stats.fmri.t;
    job.timing.fmri_t0 = def.stats.fmri.t0;
    job.volt = 1;
    job.fact = [];
    job.cvi = def.stats.fmri.cvi;
    job.global = 'None';
    job.mthresh = 0.4;
    
    % hrf specification; use canonical hrf by default
    if isfield(EXPT,'bases')
        job.bases = EXPT.bases;
    else
        job.bases.hrf.derivs = [0 0];
    end
    
    job0 = job;
    
    % subject-specific design specification
    for subj = subjects
        job = job0;
        
        S = EXPT.subject(subj);
        sub_no = regexp(S.inputdir, '(?<=sub-)[0-9]+', 'match');
        sub_no = str2double(sub_no{1});
        
        % Experiment directory & structural files
        eventdir = strrep(S.inputdir, 'model_inputs', 'model_events');
        outdir = sprintf(fullfile(S.modeldir, 'func', 'task-%s_model-%s'), ...
            task, model);
        job.dir{1} = outdir;
        job.mask{1} = fullfile(S.inputdir, S.mask);
        
        for i = 1:length(S.functional)
            run_no = regexp(S.functional{i}, '(?<=run-)[0-9]+', 'match');
            run_no = str2double(run_no{1});

            multi = EXPT.create_multi(model,sub_no,run_no);
                
            multi_file = sprintf(fullfile(S.inputdir, 'func', ...
                'sub-%02d_task-%s_run-%02d_model-%s_desc-multi.mat'), ...
                sub_no, task, run_no, model);
            save(multi_file,'-struct','multi');

            job.sess(i).hpf = def.stats.fmri.hpf;   % high-pass filter
            job.sess(i).scans{1} = fullfile(S.inputdir, S.functional{i});
            job.sess(i).multi{1} = multi_file;
            job.sess(i).cond = struct('name',{},'onset',{},'duration',{},'tmod',{},'pmod',{},'orth',{});
            job.sess(i).regress = [];
            job.sess(i).multi_reg{1} = spm_file(fullfile(S.inputdir, S.motion{i}));  % motion regressors from realignment
        end
        
        if ~isfolder(outdir)
            mkdir(outdir)
        end
        cd(outdir);

%         % drop runs with missing event files
%         keep_runs = ~cellfun(@isempty, {job.sess.scans});
%         job.sess = job.sess(keep_runs);
        job = spm_run_fmri_spec(job);

        if ~fake
            load(job.spmmat{1});
            spm_spm(SPM);
        end

        cd(cdir);
    end
end
