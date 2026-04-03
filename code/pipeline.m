bidsFolder = '.';
outputdir = [bidsFolder, filesep, 'derivatives', filesep, 'clean32'];
studyName = 'GER_Dortmund_Eyes-open-EEG-followup-EEGLAB';
STUDYfile = [outputdir, filesep, studyName, '.study'];
%%
% ~42 GiB in total, remember to enable RAM save in pop_editioptions() first!
if ~exist(STUDYfile, 'file')
    filename = [bidsFolder, filesep, 'participants.tsv'];
    readtableOptions = {
        'FileType', 'text';
        'Delimiter', '\t';
        }';
    T = readtable(filename, readtableOptions{:});
    % import multi-session subjects only
    T = T(strcmp(T.session1, 'yes') & strcmp(T.session2, 'yes'), :);
    subjects = T.participant_id';
    importbidsOptions = {
        'studyName', studyName;
        'subjects', subjects;
        'bidschanloc', 'on';
        'outputdir', outputdir;
        'bidstask', 'EyesOpen';
        }';
    [STUDY, ALLEEG] = pop_importbids(bidsFolder, importbidsOptions{:});
    % drop unused STUDY.datasetinfo fields
    unusedFields = {'session1', 'late_ses1', 'session2', 'late_ses2'};
    STUDY.datasetinfo = rmfield(STUDY.datasetinfo, unusedFields);
else
    loadstudyOptions = {
        'filename', [studyName, '.study'];
        'filepath', outputdir;
        }';
    [STUDY, ALLEEG] = pop_loadstudy(loadstudyOptions{:});
end
%%
% ~4 hours run for the first time.
datasetinfo = STUDY.datasetinfo; % to bypass parfor limitation
finalsetname = 'EDF file preprocessed'; % to resume process
parfor i = 1:numel(ALLEEG)
    % get dataset info
    EEG = ALLEEG(1, i);
    if strcmp(EEG.setname, finalsetname) % skip if setname match the final
        continue
    end
    % load dataset
    loadsetOptions = {
        'filename', EEG.filename;
        'filepath', EEG.filepath;
        'loadmode', 'all';
        }';
    EEG = pop_loadset(loadsetOptions{:});
    % fill event field on missing
    if ~isfield(EEG.event, 'type')
        for j = 1:numel(EEG.event)
            EEG.event(1, j).type = 'empty'; % MNE-Python rely on this field
        end
    end
    % set condition field
    bidsEntities = split(EEG.filename, '_');
    isacq = startsWith(bidsEntities, 'acq-');
    acq = bidsEntities(isacq);
    acqtmp = split(acq, '-');
    EEG.condition = char(acqtmp(2));
    datasetinfo(1, i).condition = char(acqtmp(2)); % update for STUDY.datasetinfo
    % resample to 2 * nyquist
    EEG = pop_resample(EEG, 2 * 250);
    % filter dataset
    EEG = pop_eegfiltnew(EEG, 'locutoff', 1);
    EEG = pop_cleanline(EEG, 'LineFrequencies', [50, 100]);
    % reject bad channels
    chaneditOptions = {
        'insert', {3, 'labels', 'AFz', 'type', 'EEG', 'datachan', 1};
        'lookup', 'standard_1005.elc';
        }'; % channel AFz is required by a specific EEG device
    EEG.urchanlocs = pop_chanedit(EEG.chanlocs, chaneditOptions{:});
    cleanOptions = {
        'BurstCriterion', 'off';
        'WindowCriterion', 'off';
        'Highpass', 'off';
        'BurstRejection', 'off';
        }'; % DO NOT USE ASR
    EEG = pop_clean_rawdata(EEG, cleanOptions{:});
    % average reref channels
    indedit = find(strcmp('FC4', {EEG.chanlocs.labels}));
    if isempty(indedit) % in case channel FC4 is also rejected as bad
        indedit = 45;
    end
    EEG.data = [EEG.data(1:indedit, :); zeros(1, EEG.pnts); EEG.data(indedit + 1:end, :)];
    EEG.nbchan = EEG.nbchan + 1;
    chaneditOptions = {
        'insert', {indedit + 1, 'labels', 'FCz', 'type', 'EEG', 'datachan', 1};
        'lookup', 'standard_1005.elc';
        }'; % channel FCz is the original online reference
    EEG = pop_chanedit(EEG, chaneditOptions{:});
    EEG = pop_reref(EEG, []);
    % run ICA & reconstruct channels
    runicaOptions = {
        'icatype', 'runica';
        'pca', EEG.nbchan - 1;
        }';
    EEG = pop_runica(EEG, runicaOptions{:});
    EEG = pop_iclabel(EEG, 'default');
    EEG = pop_interp(EEG, EEG.urchanlocs);
    % clear ICA fields as select & reorder channels will cause mismatch
    icaFields = {'icaact', 'icawinv', 'icasphere', 'icaweights', 'icachansind'};
    for field = icaFields
        EEG.(char(field)) = [];
    end
    % select time/channel subset and reorder channels
    selectChannel = {
        'Fp1','Fp2','AFz', ...
        'F7','F3','Fz','F4','F8', ...
        'FC3','FCz','FC4', ...
        'FT7','FT8', ...
        'C3','Cz','C4', ...
        'T7','T8', ...
        'CP3','CP4', ...
        'TP7','TP8', ...
        'P7','P3','Pz','P4','P8', ...
        'PO3','PO4', ...
        'O1','Oz','O2', ...
        };
    EEG = pop_select(EEG, 'channel', selectChannel);
    [~, indord] = ismember(selectChannel, {EEG.chanlocs.labels});
    EEG.data = EEG.data(indord, :);
    EEG.urchanlocs = EEG.chanlocs;
    EEG.chanlocs = EEG.chanlocs(:, indord);
    % save dataset
    savesetOptions = {
        'filename', EEG.filename;
        'filepath', EEG.filepath;
        'check', 'on';
        'savemode', 'onefile';
        }';
    EEG.setname = finalsetname; % assign final setname here
    EEG = pop_saveset(EEG, savesetOptions{:});
    EEG.data = 'in set file'; % use magic string to fool lazy loading
    ALLEEG(1, i) = EEG;
end
ALLEEG = eeg_store(ALLEEG, ALLEEG, 1:numel(ALLEEG)); % fix EEG.saved state
STUDY.datasetinfo = datasetinfo; % assign updated datasetinfo here
STUDY = std_changroup(STUDY, ALLEEG); % update channel group
STUDY = std_checkset(STUDY, ALLEEG); % update condition group
STUDY = pop_savestudy(STUDY, ALLEEG, 'savemode', 'resave');