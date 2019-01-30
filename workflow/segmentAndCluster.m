% Script for segmenting juvenile & tutor syllables from 1 session. RY
% edited segmentAndCluster.m from Bottjer lab GitHub
% 02/2018

dbstop if error
clear variables
%% Prep: file info
% [session-ID]-spike2.mat file exported from Spike2 contains all field
% information for the song channel. [session ID].mat contains only the
% audio information (.values field); all other information is stored in
% metaStruct-[session ID].mat
[matFile, matPath] = uigetfile('*.mat','Please choose the song Spike2 file','data'); %spike2 files should be named [sessionID]-spike2.mat
fprintf('Working on file: %s.\n', matFile);
[~, ID, ~] = fileparts(matFile);
split = strsplit(ID, '-');
sessionID = split{1};
birdID = strtok(sessionID, '_');

%% Juvenile syllable segmentation
% first check to see if segmenting has already been started
paused_sylls_filename = fullfile(matPath, strcat('paused_unrefinedsyllables-', sessionID, '.mat'));
syllparams_filename = fullfile(matPath, strcat('syllparams-', birdID, '.mat'));

if exist(paused_sylls_filename, 'file') == 0 %if fresh file (not starting from partway through segmentation)
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
    syllables = ...
        parseRegionsIntoSyllables(songStruct, ROIs, params,...
        'editSpecType', 'inter', 'inter.freqBands', linspace(1,10240,params.inter.NfreqBands), ...
        'preroll', 30, 'doFilterNoise',true,...
        'noiseFilter', noiseMask,'nps.reduction',-18, ...
        'dgram.minContrast', 1e-9,'minCenterFreq', 800,...
        'plot', false, 'pause', false);
    [syllables.file] = deal(songfile); %'file' field = file w/ song info ('sAudio') from which each syllable came. this is only necessary if you're going to cluster these syllables below, but doesn't hurt to just add it
    save(fullfile(matPath, strcat('syllables-', sessionID, '.mat')), 'manualMotifs', 'syllables')
    fprintf('Saved auto-segmented syllables + manualMotifs (%s). \n', sessionID)
    
    %% (Optional) Plot motifs &/or do initial round of clustering on auto-segmented syllables
    % Attempts to get a sense of juv syllables before manually refining
%   % Plot motifs for visualiation. Code is fine, but commented out b/c I didn't think it was useful
%     h=figure;
%     for ii = 1:numel(ROIs) %plot each motif (no adjustment regions; just for visualization) 
%         thisMotif = ROIs(ii);
%         fprintf('Checking motif %d/%d...\n', ii, numel(ROIs));
%         plotAndAdjust(songStruct, [], thisMotif, params, ...
%             'editSpecType', 'inter', 'inter.freqBands', linspace(1,10240,params.inter.NfreqBands), ...
%             'dgram.minContrast',1e-9, 'nps.reduction', -12, ... %try modifying these parameters if spectrogram doesn't look good
%             'doFilterNoise',true, 'noiseFilter', noiseMask, ...
%             'optGraphs',{'waveform','deriv'});
%         if ~ishandle(h) %if user closes figure, break out of loop & continue on
%             break
%         end
%     end
%     fprintf('Moving on to clustering (unrefined syllables). \n');
    
    singlecluster(birdID, syllables, params, 'noiseFilter', noiseMask); %cluster (auto-segmented, unrefined) syllables
    fprintf('Paused after initial clustering of auto-segmented syllables. Press any key to continue onto manual refinement.')
    pause

elseif exist(paused_sylls_filename, 'file') == 2
    load(paused_sylls_filename);
end
%% Manually refine syllable segmentation
% this part can be quite time-intensive. script allows user to save & quit
% partway through here; script will then load the saved info & resume on
% next run

ROIs = manualMotifs;
[~,subROIs] = checkRegions(syllables,fs);
if exist(paused_sylls_filename, 'file') == 0 %if starting from scratch
    do_segcheck = 1; %do segmentation quality check
    tmpSyllables = initEvents;
elseif exist(paused_sylls_filename, 'file') == 2 %if starting from save point
    if exist('cont', 'var') == 0 %if didn't get to segmentation quality check, start over
        do_segcheck = 1; 
        tmpSyllables = initEvents;
        disp('Starting manual syllable segmentation/refinement..')
    else
        do_segcheck = 0;
        cont = jj;
        disp('Resuming from partway through manual syllable segmentation/refinement..')
    end
end

hh = figure;
if do_segcheck == 1 
    seg_check = 5; %# of motifs to check before asking about segmentation quality
    ii = 1;
    restart = 1;
    while restart == 1
        thisEv = ROIs(ii);
        fprintf('Editing %d/%d...\n', ii, numel(ROIs));
        revisedSylls = plotAndAdjust(songStruct, subROIs, thisEv, params, ...
            'editSpecType', 'fine', 'inter.freqBands', linspace(1,10240,params.inter.NfreqBands), ...
            'dgram.minContrast',1e-9, 'nps.reduction', -12, ... %try modifying these parameters if spectrogram doesn't look good
            'doFilterNoise',true, 'noiseFilter', noiseMask, ...
            'adjustLabels',true);
        
        if ~ishandle(hh) %if user closes figure (before quality check)
            fprintf('Quitting before auto-segmenting quality check (motif #%d)!\nWill save workspace if you quit now, but will re-start from 1st motif next time.\n',seg_check);
            quit = input('Are you sure you want to quit? Y to quit, N to return to refinement: ', 's');
            if quit == 'Y' || quit == 'y' 
                paused_sylls_filename = fullfile(matPath, strcat('paused_unrefinedsyllables-', sessionID, '.mat')); 
                disp('Saving workspace..')
                save(paused_sylls_filename); %save current workspace; script will load this next time
                disp('Workspace saved.')
                return
            elseif quit == 'N' || quit == 'n'
                hh = figure; %bring back figure handle & continue (script will resume on the same motif # that fig handle was closed on)
                continue
            end
        end
            
        if ii == seg_check 
            check = input('Does syllable segmentation seem OK? Y to continue; N to quit & adjust parameters. ','s');
            if check == 'Y' || check == 'y' || isempty(check) %if auto-segmentation is ok
                restart = 0; %break out of while loop & continue on w/ the rest of motifs
                cont = ii+1;
                continue
            elseif check == 'N' || check == 'n' %if auto-segmentation wasn't good 
                %close(hh);
                [syllables, params] = resegment(songStruct, ROIs, noiseMask, params); %re-auto-segment syllables w/ adjusted params
                [syllables.file] = deal(songfile);
                save(fullfile(matPath, strcat('syllables-', sessionID, '.mat')), 'manualMotifs', 'syllables')
                fprintf('Re-saved auto-segmented syllables + manualMotifs, w/ adjusted (improved) syllable segmentation params (%s). \n', sessionID)
    
                [~,subROIs] = checkRegions(syllables,fs); %re-run checkRegions on 'new' auto-segmented syllables (w/ modified params)
                ii = 1; %re-start segmentation check
                tmpSyllables = initEvents;
                continue
            end
        end
        tmpSyllables = [tmpSyllables revisedSylls'];
        ii = ii+1;
    end
end


% manually check/adjust syllable boundaries, after initial segmentation quality check
for jj = cont:numel(ROIs) %should start at 1 after seg_check
    thisEv = ROIs(jj);
    fprintf('Editing %d/%d...\n', jj, numel(ROIs));
    revisedSylls = plotAndAdjust(songStruct, subROIs, thisEv, params, ...
        'editSpecType', 'fine', 'inter.freqBands', linspace(1,10240,params.inter.NfreqBands), ...
        'dgram.minContrast',1e-9, 'nps.reduction', -12, ... %try modifying these parameters if spectrogram doesn't look good
        'doFilterNoise',true, 'noiseFilter', noiseMask, ...
        'adjustLabels',true);
    
    if ~ishandle(hh) %if user closes figure, check if it was because of a mistake with previous motif. Otherwise, break out of loop and save progress
        reason = input('Do you want to go back (1), or quit (2)? ');
        if reason == 1 %if made mistake and want to go back 1 motif
            lastmotifsylls = find(cell2mat({tmpSyllables.start}) >= ROIs(jj-1).start & cell2mat({tmpSyllables.stop}) <= ROIs(jj-1).stop); %delete tmp sylls from last motif
            tmpSyllables(lastmotifsylls) = [];
            hh = figure;
            for kk = jj-1:jj %(re-)segment sylls from last motif + current motif
                thisEv = ROIs(kk);
                fprintf('Editing %d/%d...\n', kk, numel(ROIs));
                revisedSylls = plotAndAdjust(songStruct, subROIs, thisEv, params, ...
                    'editSpecType', 'fine', 'inter.freqBands', linspace(1,10240,params.inter.NfreqBands), ...
                    'dgram.minContrast',1e-9, 'nps.reduction', -12, ... 
                    'doFilterNoise',true, 'noiseFilter', noiseMask, ...
                    'adjustLabels',true);
                tmpSyllables = [tmpSyllables revisedSylls']; 
            end
            continue %(skips directly to next jj)
        elseif reason == 2            
            fprintf('Paused syllable refinement at bout %d. \n',jj);
            paused_sylls_filename = fullfile(matPath, strcat('paused_unrefinedsyllables-', sessionID, '.mat'));
            disp('Saving workspace..')
            save(paused_sylls_filename); %save current workspace; script will load this next time
            disp('Workspace saved.')
            return
        end
    end        
    tmpSyllables = [tmpSyllables revisedSylls'];
end

%quick checks that (1) syll start & stop times are in order & (2) didn't miss any motifs
assert(issorted(cell2mat({tmpSyllables.start})) && issorted(cell2mat({tmpSyllables.stop})), 'Syllable start/stop times are not in ascending order..');

for i = 1:size(manualMotifs)
    try
        assert(~isempty(find(cell2mat({tmpSyllables.start}) >= manualMotifs(i).start & cell2mat({tmpSyllables.stop}) <= manualMotifs(i).stop, 1))); %check that each motif has marked syllable(s)
        
    catch
        fprintf('Motif %d is missing syllable labels! Marking now...\n', i); %if found motif w/ no syllables marked
        prev = find(cell2mat({tmpSyllables.start}) <= manualMotifs(i).start,1,'last'); %last syllable from previous motif
        next = find(cell2mat({tmpSyllables.start}) >= manualMotifs(i).stop,1); %first syllable from next motif
        assert(next == prev+1);
        
        thisEv = ROIs(i); %plot that motif & mark syllables
        revisedSylls = plotAndAdjust(songStruct, subROIs, thisEv, params, ...
            'editSpecType', 'fine', 'inter.freqBands', linspace(1,10240,params.inter.NfreqBands), ...
            'dgram.minContrast',1e-9, 'nps.reduction', -12, ... %try modifying these parameters if spectrogram doesn't look good
            'doFilterNoise',true, 'noiseFilter', noiseMask, ...
            'adjustLabels',true);
        
        tmpSyllables = [tmpSyllables(1:prev) revisedSylls' tmpSyllables(next:end)]; %add to tmpSyllables
    end
end

% save approved syllables + syllable segmentation params
close(hh);
approvedSyllables = tmpSyllables';
save(fullfile(matPath, strcat('approvedSyllables-', sessionID, '.mat')), 'approvedSyllables')
if exist(syllparams_filename, 'file') == 0
    syllparams = params.syllable;
    syllparams.session = sessionID;
    save(syllparams_filename, 'syllparams');
elseif exist(syllparams_filename, 'file') == 2
    syllparams = [syllparams params.syllable]; %(syllparams should already be loaded; expanding it here)
    save(syllparams_filename, 'syllparams');
end 
delete(paused_sylls_filename)
clear tmpSyllables
disp('Finished manual syllable refinement; saved. Getting baseline periods now ..');

%% Find baseline periods
baselinePeriodsS = constructBaselineSilence(songStruct, 2, approvedSyllables, params);
[baselinePeriodsS.type] = deal('baseline');
baselinePeriodFile = fullfile(matPath, strcat('baselinePeriods-', sessionID, '.mat'));
save(baselinePeriodFile,'baselinePeriodsS');
fprintf('Baseline stuff calculated.\n');

%% Parse tutor syllables & get features
tutFile = fullfile(matPath, strcat('tutor-', birdID, '.mat'));
if exist(tutFile, 'file') == 2 %skip if TUT song has already been segmented
    fprintf('TUT song already segmented for this bird (%s).', birdID)
else
    %select tutor file
    [tutwav, tutpath, ~] = uigetfile('*.wav', 'Select tutor song (.wav) file');
    tut_wavfile = fullfile(tutpath, tutwav);
    tutorStruct = loadWavFile(tut_wavfile);
    
    % set tutor params
    tutorParams = processArgs(defaultParams,'fs',1/tutorStruct.interval,'preroll',0,'postroll',0);
    
    % highPass and amplify
    wholeRegion = getWholeSongAsRegion(tutorStruct); %get start/stop times/indices for entire song (i.e. treat as single motif, rather than using motifReturns.txt + eventFromTimes)
    tutorStruct.values = 3*highPassSample(getClip(wholeRegion, tutorStruct), tutorParams); %high passs (default 400Hz) & amplify
    
    % (loosely) select tutor song boundaries - better for display if orig .wav file is very long
    disp('Select tutor motif.');
    tutorSong = plotAndAdjust(tutorStruct,[],wholeRegion, tutorParams, ...
        'editSpecType', 'fine', 'adjustLabels',true,'dgram.minContrast',1e-8,...
        'optGraphs',{'waveform','deriv'});
    
    % automatically parse syllables
    [tutorSylls, features] = parseRegionsIntoSyllables(tutorStruct, tutorSong,tutorParams,...
        'doFilterNoise',false,'syllable.comboLength',5);
    
    % manual refinement of syllable parsing; label TUT syllables here
    tutorSylls = plotAndAdjust(tutorStruct,tutorSylls,tutorSong, tutorParams, ...
        'editSpecType', 'fine', 'adjustLabels',true,'dgram.minContrast',1e-8);     
    %[tutorSylls.file] = deal(tutFile);
    
    % get features for TUT syllables
    [tSFeats,tSSpecs] = getFeatures(tutorStruct,tutorSylls,tutorParams,...
        'plot',false,'verbose',true,'playsample',false,'doFilterNoise',false);
    
    % Save TUT syllables & features
    %initially edited to take out tutorStruct; info that was in tutorStruct
    %was split into meta-tutor-[birdID].mat & tutor-[birdID].mat, in order
    %to be compatible with getClip for easy plotting later. Later decided
    %this was unnecessary - easier to just open up TUT file in SAP for
    %visualization
    
    %sAudio = tutorStruct.values;
    %metaStruct.interval = tutorStruct.interval;
    %metaStruct.length = tutorStruct.length;
    %tut_metaFile = fullfile(matPath, strcat('meta-', tutFile, '.mat'));
    %save(tut_metaFile, 'metaStruct');
    %save(tutFile, 'sAudio','-v7.3','tutorSylls','tSFeats','tSSpecs'); 
    
    save(tutFile, 'tutorStruct','tutorSylls','tSFeats','tSSpecs');
    disp('Done w/ TUT segmentation.')
end

fprintf('All done with juvenile (& tutor) syllable segmentation (%s)!\n', sessionID)


%% f(x) to change parameters if auto syllable segmentation was no good
function [redosyllables, redoparams] = resegment(songStruct, ROIs, noiseMask, oldparams)

oldparams.syllable %display old params
fprintf('borderRise, comboLength, & minLength change syll boundaries.\n')
fprintf('greater borderRise = tighter boundaries.\nsylls closer than comboLength are linked.\n') 
fprintf('other params change syll inclusion; should only need to adjust minPower.\n')

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