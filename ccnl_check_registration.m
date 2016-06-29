function ccnl_check_registration(EXPT,subj)
    
    % Check co-registration of structural and functional images. For the
    % functional image, this uses the mean of the first run following
    % normalization.
    %
    % USAGE: check_registration(EXPT,subj)
	
    S = EXPT.subject(subj);
    P{1} = fullfile(fileparts(S.functional{1}),'mean.nii');
    P{2} = fullfile(fileparts(S.functional{1}),'wBrain.nii');
    spm_check_registration(char(P));