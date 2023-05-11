function ccnl_fmri_con_bids(EXPT,task,model,contrasts,subjects)
    
    % Construct contrasts and perform group-level analysis.
    %
    % USAGE: ccnl_fmri_con(EXPT,model,contrasts,[subjects])

    fprintf('Computing contrasts for task %s, GLM %s\n', task, model);
    
    if nargin < 4; subjects = 1:length(EXPT.subject); end
    cdir = pwd;
    C = containers.Map({'+', '-'}, [1, -1]);
    spm('defaults','fmri');
    spm_jobman('initcfg');
    if isstr(contrasts); contrasts = {contrasts}; end
    
    %% Construct contrasts
    for subj = subjects
        S = EXPT.subject(subj);
        modeldir = sprintf(fullfile(S.modeldir, 'func', 'task-%s_model-%s'), ...
            task, model);
        load(fullfile(modeldir,'SPM.mat'), 'SPM');
        
        matlabbatch = [];
        matlabbatch{1}.spm.stats.con.spmmat{1} = fullfile(modeldir,'SPM.mat');
        matlabbatch{1}.spm.stats.con.delete = 1;
        
        convec = zeros(size(SPM.xX.name));
        for j = 1:length(contrasts)
            con = regexp(contrasts{j},'[-+]','split');
            con = con(~cellfun('isempty', con)); % remove 1st element
            sgn = regexp(contrasts{j},'[-+]','match');
            %sgn = [{'+'}, sgn];
            contrasts{j}
            con
            sgn
            N = zeros(1,length(con));
            for c = 1:length(con)
                con{c} = strtrim(con{c});
                sgn{c} = strtrim(sgn{c});
                ix = logical(zeros(1,length(SPM.xX.name)));
                found = false;
                for i = 1:length(SPM.xX.name)
                    if isequal(strfind(SPM.xX.name{i},[con{c},'*']), 1) || ~isempty(strfind(SPM.xX.name{i},[' ',con{c},'*'])) || ~isempty(strfind(SPM.xX.name{i},[' ',con{c},'^']))  || ~isempty(strfind(SPM.xX.name{i},['x',con{c},'^']))
                        convec(j,i) = C(sgn{c});
                        N(c) = N(c) + 1;
                        ix(i) = 1;
                        found = true;
                    end
                end
                assert(found, ['Cannot find regressor ', con{c}, ' -- typo maybe?']);
%                convec(j,ix) = convec(j,ix)/N(c);
            end
            convec

% %             <<<< DEBUG: Fix unbalanced contrasts - NV
%             pos = convec(j,:) > 0;
%             pos_regressors = SPM.xX.name(pos);
%             pos_runs = cellfun(@(r) str2double(regexp(r, '(?<=Sn\()([0-9]+)', ...
%                 'match', 'once')), pos_regressors);
% 
%             neg = convec(j,:) < 0;
%             neg_regressors = SPM.xX.name(neg);
%             neg_runs = cellfun(@(r) str2double(regexp(r, '(?<=Sn\()([0-9]+)', ...
%                 'match', 'once')), neg_regressors);
% 
%             in_both = intersect(pos_runs, neg_runs);
%             reg_runs = cellfun(@(r) str2double(regexp(r, '(?<=Sn\()([0-9]+)', ...
%                 'match', 'once')), SPM.xX.name);
%             valid_runs = ismember(reg_runs, in_both);
% 
%             convec_filtered = zeros(size(convec));
%             convec_filtered(j,pos & valid_runs) = 1/sum(pos(valid_runs));
%             convec_filtered(j,neg & valid_runs) = -1/sum(neg(valid_runs));
%             convec = convec_filtered;
            % >>>>
            matlabbatch{1}.spm.stats.con.consess{j}.tcon.name = contrasts{j};
            matlabbatch{1}.spm.stats.con.consess{j}.tcon.convec = convec(j,:);
            matlabbatch{1}.spm.stats.con.consess{j}.tcon.sessrep = 'none';
        end
        
        spm_jobman('run',matlabbatch);
    end
    
    cd(cdir);
    
    %% save contrasts
    modeldir0 = fullfile(EXPT.modeldir, 'group', ['task-' task '_model-' model]);
    if ~isfolder(modeldir0)
        disp(['Saving outputs in new folder: ' modeldir0]);
        mkdir(modeldir0)
    end
    save(fullfile(modeldir0,'contrasts'),'contrasts','convec');
    
    %% Group-level analysis
    for j = 1:length(contrasts)
        modeldir = sprintf(fullfile(modeldir0, 'con-%04d'), j);
        if isfolder(modeldir)
            disp(['Old results found, clearing directory: ' modeldir]);
            rmdir(modeldir,'s')
        end
        mkdir(modeldir);
        cd(modeldir)
        job.dir{1} = modeldir;
        job.multi_cov = []; job.globalc.g_omit = 1;
        job.globalm.glonorm = 1; job.globalm.gmsca.gmsca_no = 1;
        job.masking.tm.tm_none = 1; job.masking.em{1} = [];
        
        % covariates - not working yet!!
        if isfield(EXPT,'cov')
            if ~isfield(EXPT.cov,'iCC')
                for i = 1:length(EXPT.cov); EXPT.cov(i).iCC = 1; end
            end
            if ~isfield(EXPT.cov,'iCFI')
                for i = 1:length(EXPT.cov); EXPT.cov(i).iCFI = 1; end
            end
            job.cov = EXPT.cov;
        else
            job.cov = [];
        end
        
        con = sprintf('con_%04d.nii',j);
        for s = 1:length(subjects)
            job.des.t1.scans{s} = fullfile(EXPT.subject(s).modeldir, 'func', ...
                ['task-' task '_model-' model], con);
        end
        
        out = spm_run_factorial_design(job);
        load(out.spmmat{1});
        SPM = spm_spm(SPM);
        
        % write contrasts and t-maps
        matlabbatch = [];
        matlabbatch{1}.spm.stats.con.spmmat{1} = fullfile(modeldir,'SPM.mat');
        matlabbatch{1}.spm.stats.con.delete = 1;
        for i = 1:length(SPM.xX.name)
            convec = zeros(size(SPM.xX.name)); convec(i) = 1;
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = [contrasts{j},', ',SPM.xX.name{i}];
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = convec';
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        end
        
        spm_jobman('run',matlabbatch);
        
        cd(cdir);
    end
