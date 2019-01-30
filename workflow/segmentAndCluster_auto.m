% Script for segmenting juvenile & tutor syllables from multiple sessions.
% Modified from segmentAndCluster.m to auto-segment and save syllables (no
% manual refinement) from multiple sessions; meant to be a quick & dirty
% rough pass at segmentation
%
% Need Spike2-exported .mat file that contains all field info for the song
% channel, named [session ID]-spike2.mat
%
% 11/2018

dbstop if error
clear variables

[files, matPath, ~] = uigetfile('*.mat', 'Select all -spike2 .mat files to analyze','C:\Users\rachelyu\Documents\MATLAB\', 'MultiSelect', 'on');
files = cellstr(files); %ensures that FileNames is a cell array (ex. if selected only 1 file)

for current_file = 1:length(files) %for each file/session
    matFile = files{current_file};
    fprintf('Working on file: %s.\n', matFile);
    [~, ID, ~] = fileparts(matFile);
    split = strsplit(ID, '-');
    sessionID = split{1};
    birdID = strtok(sessionID, '_');
    
    % check for existing params for this bird
    syllparams_filename = fullfile(matPath, strcat('syllparams-', birdID, '.mat'));
    
    %% Get song & bout onset/offset info from Spike2-exported files
    loaded = load(fullfile(matPath, matFile));
    fld = fieldnames(loaded);
    songStruct = loaded.(fld{1});
    fs = 1/songStruct.interval;
    
    %sAudio & metaStruct aren't used in this script, but are needed later in other scripts
    sAudio = songStruct.values;
    songfile = fullfile(matPath, strcat(sessionID, '.mat'));
    save(songfile, 'sAudio', '-v7.3'); %save audio only (song voltage values) into [sessionID].mat file; needs to be saved in v7.3 to enable partial loading later
    
    metaStruct = rmfield(songStruct,'values'); %metaStruct is everything in Spike2-exported structure except values field
    metaFile = fullfile(matPath, strcat('meta-', sessionID, '.mat'));
    save(metaFile, 'metaStruct');
    
    fprintf('Saved audio & meta files. Loaded spike audio %s.\n', sessionID);
    
    %Get bout onset/offset times from Spike2-exported text file:
    returnFil = fullfile(matPath, strcat('motifReturn-', sessionID, '.txt'));
    returnFil = strrep(returnFil, '.mat','');
    returnFid = fopen(returnFil);
    
    fgetl(returnFid); fgetl(returnFid); %read through header & dead line
    if string(fgetl(returnFid)) == '"LOW"' %if motifs are on level channel (+ reads through that line)
        times = textscan(returnFid,'%f');
        times = times{1};
        fclose(returnFid);
        
        returnStarts = times(1:2:end);
        returnStops = times(2:2:end);
    else
        frewind(returnFid) %start from beginning
        fgetl(returnFid); fgetl(returnFid); %read through header & dead line
        returnStarts = textscan(returnFid,'%f'); returnStarts = returnStarts{1}; %reads through first block of numbers (starts)
        
        fgetl(returnFid); fgetl(returnFid); %reads through next 2 lines (header & dead line for offsets channel)
        returnStops = textscan(returnFid,'%f'); returnStops = returnStops{1}; %reads through second block of numbers (stops)
        fclose(returnFid);
    end
    
    % check that starts & stops are matched
    returnSIndices = [zeros(numel(returnStarts),1); ones(numel(returnStops), 1)]; %zeros for return starts, 1s for return stops
    [fusedTimes,sIdxs] = sort([returnStarts; returnStops]); %sort starts & stops in 1 single vector
    returnSIndices = returnSIndices(sIdxs);
    if any(diff(returnSIndices) == 0) % returnSIndices should be alternating 1/0; i.e., diff's should = 1
        isBadIndices = [diff(returnSIndices) == 0; sIdxs(end)==0];
        fusedTimes(isBadIndices)
    end
    assert(numel(returnStarts) == numel(returnStops), 'ERROR: Check your bout on/offsets!');
    assert(all(returnStarts < returnStops) && all(returnStarts(2:end) > returnStops(1:end-1)), 'ERROR: Check your bout on/offsets!');
    
    manualMotifs = eventFromTimes(returnStarts, returnStops, fs); %bout times
    
    %% Set parameters
    params = defaultParams; %to use alias, ex. params = progressArgs(defaultParams, 'Gy242')
    params.fs = fs;
    
    if exist(syllparams_filename, 'file') == 0
        % plot clip of song to help  determine appropriate minPower threshold:
        sample = getClip(manualMotifs(1),  songStruct); %sample first motif
        sample = highPassSample(sample,  params); %use high-passed version when  getting spectrum stats (requires defined fs)
        params.coarse.fs = fs; %need to  add this field if using a specified  parameter set for spectrogram plotting
        spec = getMTSpectrumStats(sample,  params.coarse);
        plotAllFigures(spec,  manualMotifs(1), params, 'optGraphs',  {'waveform','totalPower'});
        title('Select minPower threshold  on RMS Power plot')
        [~, minPower] = ginput(1);
        close
        clear sample
        params.syllable.minPower = minPower;
    elseif exist(syllparams_filename, 'file') == 2
        load(syllparams_filename); %loads syllparams variable
        params.syllable = syllparams(end); %load most recent syllable params
        params.syllable.session = sessionID; %need this for later, when concatenating syll params w/ previous sessions
    end
    %% Find noise
    notNoise = manualMotifs;
    approved = false;
    while ~approved
        candidateNoise = autodetectNoise(songStruct, notNoise, params);
        fprintf('Is this clip pure noise? Mark y if yes...\n');
        approved = markRegions(songStruct,candidateNoise);
        if ~approved
            notNoise = [notNoise; eventFromTimes(0,candidateNoise.stop,fs)];
        end
    end
    noiseMask = noiseAnalysis(songStruct, candidateNoise);
    params.noiseFilter = noiseMask;
    save(fullfile(matPath, strcat('noiseMask-', sessionID, '.mat')), 'noiseMask');
    
    %% Automatically parse juvenile motifs into syllables
    ROIs = manualMotifs;
    autoSyllables = ...
        parseRegionsIntoSyllables(songStruct, ROIs, params,...
        'editSpecType', 'inter', 'inter.freqBands', linspace(1,10240,params.inter.NfreqBands), ...
        'preroll', 30, 'doFilterNoise',true,...
        'noiseFilter', noiseMask,'nps.reduction',-18, ...
        'dgram.minContrast', 1e-9,'minCenterFreq', 800,...
        'plot', false, 'pause', false);
    [autoSyllables.file] = deal(songfile); %'file' field = file w/ song info ('sAudio') from which each syllable came
    
    %% Check auto-segmenting params
    ROIs = manualMotifs;
    [~,subROIs] = checkRegions(autoSyllables,fs);
    
    hh = figure;
    seg_check = 5; %# of motifs to check before asking about segmentation quality
    ii = 1;
    restart = 1;
    
%     click = -1;
%     while click == -1
%         sound([sin(1:.6:400), sin(1:.6:400), sin(1:.6:400)])
%         click = getkeywait(5);
%     end
    
    while restart == 1
        thisMotif = ROIs(ii);
        fprintf('Checking motif %d/%d...\n', ii, numel(ROIs));
        plotAndAdjust(songStruct,subROIs, thisMotif, params, ...
            'editSpecType', 'inter', 'inter.freqBands', linspace(1,10240,params.inter.NfreqBands), ...
            'dgram.minContrast',1e-9, 'nps.reduction', -12, ... %try modifying these parameters if spectrogram doesn't look good
            'doFilterNoise',true, 'noiseFilter', noiseMask, ...
            'optGraphs',{'waveform','deriv'});
        
        if ~ishandle(hh) %bring back figure if accidentally close handle
            hh = figure;
            continue
        end
        
        if ii == seg_check
            check = input('Does syllable segmentation seem OK? Y to continue; N to quit & adjust parameters. ','s');
            if check == 'Y' || check == 'y' || isempty(check) %if auto-segmentation is ok
                restart = 0; %break out of while loop
                continue
            elseif check == 'N' || check == 'n' %if auto-segmentation wasn't good
                [autoSyllables, params] = resegment(songStruct, ROIs, noiseMask, params); %re-auto-segment syllables w/ adjusted params
                [autoSyllables.file] = deal(songfile);
                save(fullfile(matPath, strcat('syllables-', sessionID, '.mat')), 'manualMotifs', 'syllables')
                fprintf('Re-saved auto-segmented syllables + manualMotifs, w/ adjusted (improved) syllable segmentation params (%s). \n', sessionID)            
                [~,subROIs] = checkRegions(autoSyllables,fs); %re-run checkRegions on 'new' auto-segmented syllables (w/ modified params)
                ii = 1; %re-start segmentation check
                continue
            end
        end
        ii = ii+1;
    end    
    close(hh);
    autoSyllables = subROIs;
    
    %quick checks that syll start & stop times are in order
    assert(issorted(cell2mat({autoSyllables.start})) && issorted(cell2mat({autoSyllables.stop})), 'Syllable start/stop times are not in ascending order..');
    
    % save segmented syllables + syllable segmentation params
    save(fullfile(matPath, strcat('autosyllables-', sessionID, '.mat')), 'manualMotifs', 'autoSyllables')
    if exist(syllparams_filename, 'file') == 0
        syllparams = params.syllable;
        syllparams.session = sessionID;
        save(syllparams_filename, 'syllparams');
    elseif exist(syllparams_filename, 'file') == 2
        syllparams = [syllparams params.syllable]; %(syllparams should already be loaded; expanding it here)
        save(syllparams_filename, 'syllparams');
    end
    disp('Finished auto syllable segmenting; saved. Getting baseline periods now ..');
    
    %% Find baseline periods
    params.verbose = false;
    baselinePeriodsS = constructBaselineSilence(songStruct, 2, autoSyllables, params);
    [baselinePeriodsS.type] = deal('baseline');
    baselinePeriodFile = fullfile(matPath, strcat('baselinePeriods_autosylls-', sessionID, '.mat'));
    save(baselinePeriodFile,'baselinePeriodsS');
    
    fprintf('Baseline stuff calculated.\nDone with auto (unrefined) syllable segmentation for (%s); moving on ...\n', sessionID);
    clearvars -except matPath files current_file
end
Fprintf('ALLLLL DONEEEE');

%% f(x) to change parameters if auto syllable segmentation was no good
function [redosyllables, redoparams] = resegment(songStruct, ROIs, noiseMask, oldparams)

oldparams.syllable %display old params
fprintf(['borderRise, comboLength, & minLength change syll boundaries:\n', ...
    'greater borderRise = tighter boundaries.\nsylls closer than comboLength are linked.\n\n', ...
    'Other params change syll inclusion; should only need to adjust minPower, possibly flatFactor:\n', ...
    'flatFactor = factor at which to call minima/maxima the same; larger = reject more extrema.\n', ...
    'maybe change if sylls are being neglected (flatFactor too high) or cut in half (too low).']);

syll_params = fieldnames(oldparams.syllable);
[idx, cancel] = listdlg('PromptString', 'Which params do you want to adjust?', ...
    'ListString', syll_params, 'CancelString', 'Default Params', 'SelectionMode', 'multiple');

if cancel == 0 %if default params is selected
    defaultparams = processArgs(defaultParams);
    oldparams.syllable = defaultparams.syllable;
elseif cancel == 1
    for j = 1:length(idx) %for each param that was selected to be changed, prompt user to enter new value
        prompt = sprintf('Enter the new value for %s: ', syll_params{idx(j)});
        newval = input(prompt);
        oldparams.syllable.(syll_params{idx(j)}) = newval; %change in params
    end
end
redoparams = oldparams; %renaming here just for clarity
redoparams.syllable %display new params
disp('Check new parameters ^. Press any key to continue on to re-segmenting syllables');
pause

redosyllables = ... %re-do auto syllable segmentation w/ new params
    parseRegionsIntoSyllables(songStruct, ROIs, redoparams,...
    'editSpecType', 'inter', 'inter.freqBands', linspace(1,10240,redoparams.inter.NfreqBands), ...
    'preroll', 30, 'doFilterNoise',true,...
    'noiseFilter', noiseMask,'nps.reduction',-18, ...
    'dgram.minContrast', 1e-9,'minCenterFreq', 800,...
    'plot', false, 'pause', false);
end