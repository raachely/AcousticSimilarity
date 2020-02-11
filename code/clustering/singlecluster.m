function singlecluster(birdID, sylls, fs, noiseMask, params, varargin)
% "quick" cluster: calculates spectra & featureTable for sylls & clusters
% them. Does this for all entries in sylls. For now(?), assumes all
% syllables are from the same session (if fs & noiseMask aren't provided,
% uses 1 syllable's file to extract session ID, to find associated
% metaStruct (for sampling rate) & noiseMask files).

% Wrote to be called by segmentAndCluster to do preliminary
% clustering on auto-segmented syllables, but should work as standalone
% f(x) as well

% syllables should be typical syllable structure, + 'file' field. requires
% that associated meta-, sAudio, & noise mask files for these syllables
% have been generated

% specify the number of clusters you want w/ num_desired_clusters
% RY 3/2018

dbstop if error
%% Set-up
num_desired_clusters = 15; %based on results from DRcluster; must be within 4-25
[filepath, sessionID, ~] = fileparts(sylls(1).file); %use first syllable file

if nargin < 5
    params = defaultParams;
end
params = processArgs(params, varargin{:});

if nargin == 2
    %get sampling rate & noiseMask (even if they're already contained w/in params; easier to just have the variables here)
    rep = reportOnData(birdID, sessionID, params);
    
    fs = 1/getfield(loadFromManifest(rep.manifest,'metaStruct'),'interval');
    if ~isnan(params.fs) %if sampling rate was already fed in through params, check that it's correct
        assert(params.fs == fs, 'Sampling rates are different');
    end
    noiseMask = loadFromManifest(rep.manifest, 'noiseMask');
    if ~isempty(params.noiseFilter) %if noise mask was already fed in through params, check that it's correct
        assert(isequal(params.noiseFilter, noiseMask), 'Noise masks are different');
    end
end
params.fs = fs;
params.noiseFilter = noiseMask;

% remove any syllables that are too short
isTooShort = (params.fine.windowSize / 1000 > [sylls.stop] - [sylls.start]);
sylls(isTooShort) = [];
%fprintf('Removing %d syllables that are too short...\n', sum(isTooShort));

%% Get spectra & feature tables
N = numel(sylls);
fieldsToKeep = {'AM','FM','pitchGoodness','wienerEntropy','fundamentalFreq','times'};
spectra = initEmptyStructArray(fieldsToKeep, N);
featureTable = cell(1,N);
progressbar(sprintf('Calculating spectra & features for syllables (# = %d)',N));
for ii = 1:N
    cl = getClipAndProcess([],sylls(ii), params, 'noroll', ...
        'doFilterNoise', true, 'noiseFilter', noiseMask);
    params.fine.fs = fs;
    pF = params.fine;
    pF.features = [pF.features 'pitchGoodness'];
    tmpSpec = getMTSpectrumStats(cl, pF);
    featureTable{ii} = extractFeatures(tmpSpec);
    progressbar(ii/N);
    for jj = 1:numel(fieldsToKeep)
        spectra(ii).(fieldsToKeep{jj}) = tmpSpec.(fieldsToKeep{jj});
    end
end
featureTable = [featureTable{:}];

%% Cluster syllables
params = processArgs(defaultParams, 'warpingCost', 1.0);
t1 = clock;
[~, distMatrices, ~] = DRcluster(sylls, featureTable, spectra, params); %get distances
nClusters = 4:25;
linktree = linkage(distMatrices.cosim, 'complete'); %do clustering
clusterIdxs = cluster(linktree,'maxclust',nClusters);
fprintf('Time for total clustering: %0.2fs\n',etime(clock, t1));

%% Plot clusters
tempfoldername = strcat('TEMP-unref_syll_clusts-',sessionID);
mkdir(fullfile(filepath, tempfoldername));
tempfolder = fullfile(filepath, tempfoldername);

takenIdxs = clusterIdxs(:,(num_desired_clusters-3)); %columns in clusterIdxs correspond to each max cluster threshold, 4-25
nTypes = max(takenIdxs);
for jj = 1:nTypes
    fig = figure;
    theseSylls = sylls(takenIdxs==jj);
    fprintf('Plotting trained clusters for syllable %d (#=%d)...\n',jj, sum(takenIdxs==jj));
    
    mosaicDRSpec(theseSylls, params, 'dgram.minContrast', 1e-10, ...
        'preroll', 3, 'postroll', 3, 'maxMosaicLength', 5.5);
    set(fig, 'Name', sprintf('Syllable #%d',jj));
        
    fig_name = strcat('unrefinedSyllableCluster-c',num2str(jj));
    figFileName = fullfile(tempfolder, fig_name);
    savefig(fig, figFileName);
    %close(fig);
end
fprintf('Figures for clusters of unrefined syllables saved in (%s). \n', tempfolder);

