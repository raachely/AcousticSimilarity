function prepareDRsylls(birdID, syllName, params, varargin)
% calculates spectra and feature tables of all syllables across all
% sessions, and calculates distance matrices between syllables w/in each
% age
% 
% generates & saves allSpecs-[birdID].mat, which contains DRsylls, spectra,
% & featureTable variables. 
%
% call to DRcluster() generates distance measures that are later used for
% clustering. Distance measures are contained and saved in altClustDataAge
% ... .mat.

dbstop if error

if nargin < 2
    syllName = 'approvedSyllables';
end
if nargin < 3
    params = defaultParams;
end
params = processArgs(params, varargin{:});

birdPath = [pwd filesep 'data' filesep birdID];
if ~exist(birdPath, 'dir')
    error('Data for %s does not exist, aborting...', birdID);
end

saveFile = fullfile(birdPath, strcat('allSpecs-', birdID, '.mat'));
start = 1;

%% Compile syllable info across all sessions
rep = reportOnData(birdID, [], params, 'rejectNoNeuronSessions', false); %change to true to cluster only if a given session has associated spike file(s)
DRsylls = struct([]);
sessionNum = [];
for ii = start:numel(rep) %for each session
    iSession = rep(ii).sessionID;
    iMani = rep(ii).manifest;
    
    if ~findInManifest(iMani, syllName)
        warning('Variable %s not found assigned for session %s, skipping session...', ...
            syllName, iSession);
        continue;
    end
    
    sylls = loadFromManifest(iMani, syllName);
    [sylls.file] = deal([birdPath filesep iSession '.mat' ]); %'file' field = file w/ song info ('sAudio') from which each syllable came
    [sylls.age]  = deal(getAgeOfSession(iSession)); 
    %DRsylls = [DRsylls sylls];
    DRsylls = [DRsylls sylls']; %temporary change
    sessionNum = [sessionNum ii * ones(1,numel(sylls))]; %keeps track of which syllables came from which session 
end

% get sampling rate
ii = 1;
while ~findInManifest(rep(ii).manifest, 'metaStruct')
    ii = ii+1;
end
fs = readSamplingRate(rep(ii).manifest);
params.fs = fs;
params.fine.fs = fs;

% remove any syllables that are too short
isTooShort = (params.fine.windowSize / 1000 > [DRsylls.stop] - [DRsylls.start]);
DRsylls(isTooShort) = [];
sessionNum(isTooShort) = [];
fprintf('Removing %d syllables that are too short...\n', sum(isTooShort));

%% Get spectra and feature tables for all syllables
N = numel(DRsylls);
fieldsToKeep = {'AM','FM','pitchGoodness','wienerEntropy','fundamentalFreq','times'}; 
spectra = initEmptyStructArray(fieldsToKeep, N);
featureTable = cell(1,N);
progressbar(sprintf('Calculating spectra & features for regions (# = %d)',N));
for ii = start:N
    % get noisemask
    if ii==1 || sessionNum(ii-1) ~= sessionNum(ii)
        [nMExist, nMFile] = findInManifest(rep(sessionNum(ii)).manifest, 'noiseMask');
        if ~nMExist
            warning('Noise mask not available for session %s, please load or exit... ', ...
                rep(sessionNum(ii)).sessionID);
            keyboard
        else
            fprintf('Loading noise mask from %s...\n',nMFile);
            noiseMask = loadFromManifest(rep(sessionNum(ii)).manifest, 'noiseMask');
        end
        if ~exist('noiseMask', 'var'), error('Noise mask not loaded');
        end
    end
    
    cl = getClipAndProcess([],DRsylls(ii), params, 'noroll', ...
        'doFilterNoise', true, 'noiseFilter', noiseMask);
    pF = params.fine;
    pF.features = [pF.features 'pitchGoodness']; %RY changed to pitchGoodness, & modified getMTSpectrumStats to use 'pitchGoodness' instead of 'harmonicPitch' as call string
    tmpSpec = getMTSpectrumStats(cl, pF);
    featureTable{ii} = extractFeatures(tmpSpec); %121 fields of feature values for each syllable
    progressbar(ii/N);
    for jj = 1:numel(fieldsToKeep)
        spectra(ii).(fieldsToKeep{jj}) = tmpSpec.(fieldsToKeep{jj});
    end
end
featureTable = [featureTable{:}];

save(saveFile,'DRsylls','spectra','featureTable');
fprintf('Syllables, spectra, & features across all sessions for %s saved to %s.\n',birdID, saveFile);

%% Calculate distance matrices (by age)
% runs DRcluster on all syllables from within each age, resulting in
% altClustDataAge-[age]-[date/time stamp].mat, which contains distance
% measures between syllabes from that age.
fprintf('Starting distance matrix calculations:\n')
params.warpingCost = 1.0;

%make cluster folder if it doesn't yet exist for this bird
clustFolder = fullfile(pwd, 'data', strcat('cluster-', birdID), filesep);
if exist(clustFolder, 'dir') == 0
    mkdir(clustFolder) 
end

% loop through each age
ages = unique([DRsylls.age]);
startAge = 1; 
for ii = startAge:numel(ages) %for each age ..  
    thisAge = ages(ii);
    isThisAge = ([DRsylls.age]==thisAge);
    seld = find(isThisAge); %seld indexes into DRsylls; = rows that correspond to syllables from a particular age (isThisAge)
    
    timeFlag = datestr(clock, 'mm_dd_HH_MM');
    diary([clustFolder 'diary-' timeFlag '-age' num2str(thisAge) '.txt']);
    
    t1 = clock;
    [empMats, distMats, empDistrs] = DRcluster(DRsylls(seld), featureTable(seld), spectra(seld), params); %calculate distance matrices between syllables for syllables from this age
    fprintf('Time for total distance calculations: %0.2fs\n',etime(clock, t1));
    
    clusterFile = [clustFolder 'altClustDataAge-' num2str(ages(ii)) '-' timeFlag]; 
    save([clusterFile '.mat'], 'empMats', 'distMats', 'empDistrs', 'seld');
    
    diary off
end
fprintf('All done; allSpecs & cluster files for %s generated & saved.\n', birdID)
