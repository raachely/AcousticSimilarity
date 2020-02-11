function database = reportOnData(birdIDs, sessions, params, varargin)
% reportOnData looks in the data folder and sees what is generated for what
% folder
% if birdIDs are provided, only looks at those birdIDs
% if sessions are provided, only looks at those sessions
%
% requisite data directory structure: current path must be in parent
% directory to folder titled "data"; "data" contains folder(s) titled by
% bird IDs, and each bird ID folder contains all related data files for
% that bird. Within each bird file, additional subfolders contain files
% titled *_times.mat ("spikeFiles") that are lists of spike timepoints for
% each neuron; in order to be found, these must be in separate subfolders
% within the bird folder - not n the parent bird folder itself.


% note: to clear old contents files, run clearSummaries.m
dataDir = [pwd filesep 'data']; %directory for "data" folder
relDataDir = 'data';
% find all the bird IDs from directory
if nargin < 1 || isempty(birdIDs)
    files = dir(dataDir);
    birdIDs = {files([files.isdir]).name};
    isID = ~cellfun('isempty',regexp(birdIDs, '^[A-Z][a-z]?\d{1,3}')); %could make tighter by using cap [????]
    birdIDs = birdIDs(isID);
end
if nargin < 2, sessions = ''; end
if nargin < 3 || isempty(params), params = defaultParams; end
if nargin > 3, params = processArgs(params,varargin{:}); end

if ~iscell(birdIDs)
    birdIDs = {birdIDs};
end

summaryFields = {'sessionID','manifest','spikeFiles'};

database = cell(1,numel(birdIDs));
for ii = 1:numel(birdIDs) % for each bird
    thisID = birdIDs{ii};
    thisBirdPath = [dataDir filesep thisID];
    relBirdPath  = [relDataDir filesep thisID];
    matFiles = dir([thisBirdPath filesep '*.mat']);
    matFiles = {matFiles.name};
    
    % look for spike-sorted neuron files - list of neuron spiking
    % timepoints; files all titled *_times.mat. These should all be within
    % subdirectories (folders) in each bird directory (folder)
    subDirs = dir(thisBirdPath); subDirs = {subDirs([subDirs.isdir]).name}; %get subdirectories
    subDirs(1:2)=[]; % get rid of the first two . / .. directories
    
    birdNeuronFiles = '';
    for jj = 1:numel(subDirs)
        thisSubDir = [thisBirdPath filesep subDirs{jj}]; % in each subdirectory..
        relSubDir = [relBirdPath filesep subDirs{jj}];
        possSpikeFiles = dir([thisSubDir filesep thisID '*times.mat']); % ..get all the neuron spike files
        if isempty(possSpikeFiles)
            continue
        end
        birdNeuronFiles = [birdNeuronFiles strcat([relSubDir filesep], {possSpikeFiles.name})];
        if params.rejectSpikeFiles
            filesToReject = ~cellfun('isempty', strfind(birdNeuronFiles,'REJECTME')); % OK to delete I think?
            birdNeuronFiles = birdNeuronFiles(~filesToReject);
        end
    end
    isNFClaimed = false(1,numel(birdNeuronFiles));
    
    % look for song channel files exported from spike2; all should be named
    % with birdID at front. RY modified b/c started including
    % spike2-exported files, which were getting mixed up. Subsequent report
    % generation depends on existence of sAudio file (had originally tried
    % modifying so that report building was dependent instead on existence
    % of spike2 file, but realized that was more trouble than it was worth.
    % Left in code, commented out)
    %isAudioFile = strncmp(matFiles, thisID, numel(thisID));
    sessionFiles = find(strncmp(matFiles, thisID, numel(thisID))); %index of files that start w/ bird ID; either spike2 or audio files
    isSpike2File = find(contains(matFiles, 'spike2')); %index of files that contain 'spike2' in filename
    %spike2Files = matFiles(isSpike2File); %spike2 files for this bird
    isAudioFile = setdiff(sessionFiles, isSpike2File); %audio file = file that doesn't have 'spike2' in the file name
    audioFiles = matFiles(isAudioFile);  %sAudio files for this bird
    
    %spike2Files = strrep(spike2Files, '.mat', '');
    audioFiles = strrep(audioFiles, '.mat',''); % get rid of the mat suffixes, in order to aid matching w/ intersect below
    if ~isempty(sessions) %if specific session is specified, get sAudio file for just that session
        audioFiles = intersect(audioFiles, sessions);
        %spike2Files = intersect(spike2Files, strcat(sessions,'-spike2'));
    end
    
    sessionRecords = initEvents(numel(audioFiles),summaryFields);
    %sessionRecords = initEvents(numel(spike2Files),summaryFields);
    
    for jj = 1:numel(audioFiles) % for each session
        thisSession = audioFiles{jj};
        %thisSession = erase(spike2Files{jj},'-spike2'); %loop through spike2Files instead of audioFiles, but clear '-spike2' appendage from session ID
        thisSpikeStem = strtok(thisSession, '.');
        
        isRelatedFile = ~cellfun('isempty',strfind(matFiles,thisSession));
        relatedFiles = matFiles(isRelatedFile); % includes the original spikeFile
        
        % RY commented out summary file generation:
        %         % summary file -> not in the data directory,
        %         % but its parent directory
        %         summaryDir = [dataDir filesep '..' filesep 'summaries'];
        %         if ~exist(summaryDir, 'dir'), mkdir(summaryDir); end
        %
        %         summaryFile = [summaryDir filesep ...
        %             'contents-' thisSession '.mat'];
        %         % get a manifest of all the files
        %
        %         % first, check to see if the contents are updated
        %         freshContents = false;
        %         if exist(summaryFile,'file')
        %             freshContents = true;
        %             contentUpdate = getModifiedStamp(summaryFile);
        %             for kk = 1:numel(relatedFiles)
        %                 if(getModifiedStamp([thisBirdPath filesep relatedFiles{kk}]) > contentUpdate)
        %                     freshContents = false;
        %                     break;
        %                 end
        %             end
        %
        %             if ~isempty(birdNeuronFiles)
        %                 % are we related to this particular session?
        %                 isRelatedNeuronFile = ~cellfun('isempty',strfind(birdNeuronFiles,thisSpikeStem));
        %                 relatedNeuronFiles = birdNeuronFiles(isRelatedNeuronFile);
        %
        %                 for kk = 1:numel(relatedNeuronFiles)
        %                     if(getModifiedStamp(relatedNeuronFiles{kk}) > contentUpdate)
        %                         freshContents = false;
        %                         break;
        %                     end
        %                 end
        %             end
        %         end
        %         % if the summary file is fresher than its related data files
        %         if freshContents
        %             % read off the manifest from the previous file, converting the one-field
        %             % structure to a struct.
        %             tmp = load(summaryFile);
        %             flds = fieldnames(tmp);
        %             sessionRecords(jj) = tmp.(flds{1});
        %
        %             if ~isempty(birdNeuronFiles)
        %                 isNFClaimed = isNFClaimed | cellfun(...
        %                     @(x) any(strcmpi(x,sessionRecords(jj).spikeFiles)), birdNeuronFiles);
        %             end
        %         else
        [~,idxShortest] = min(cellfun('length',relatedFiles));
        [~,contents.sessionID,~] = fileparts(relatedFiles{idxShortest}); %assumes shortest file name is session ID only ..
        
        relatedFiles = strcat([relBirdPath filesep], relatedFiles);
        contents.manifest = getManifest(relatedFiles);
        contents.spikeFiles = [];
        
        % find the neuron files that are linked with the session
        if ~isempty(birdNeuronFiles)
            isRelatedNeuronFile = ~cellfun('isempty',strfind(birdNeuronFiles,thisSpikeStem));
            
            isNFClaimed = isNFClaimed | isRelatedNeuronFile;
            
            contents.spikeFiles = birdNeuronFiles(isRelatedNeuronFile);
        end
        
        sessionRecords(jj) = contents;
        %            save(summaryFile,'contents');
        %         end
               
        
        if params.verbose && ~params.quiet
            spikeData = loadSpikeData(sessionRecords(jj).spikeFiles);
            nNeurons = numel(spikeData);
            fprintf('%s/%s with %d recorded variables, %d neurons...\n', thisID, thisSession, numel(sessionRecords(jj).manifest), nNeurons);
            for kk = 1:numel(sessionRecords(jj).manifest)
                fprintf('\t%s --> %s\n', sessionRecords(jj).manifest(kk).originalFile, ...
                    sessionRecords(jj).manifest(kk).name);
            end
        end
    end
    
    unclaimedSessionRecords = initEvents(0,summaryFields);
    if any(~isNFClaimed) && isempty(sessions)
        unclaimedNeuronFiles = birdNeuronFiles(~isNFClaimed);
        % get the session IDs from the neuron files
        sessionID = cell(1,numel(unclaimedNeuronFiles));
        for jj = 1:numel(unclaimedNeuronFiles)
            [~, birdFile] = fileparts(unclaimedNeuronFiles{jj});
            sessionID{jj} = birdFile(1:strfind(birdFile, '_ch')-1);
        end
        [uniqIDs,~,crossIdxs] = unique(sessionID);
        
        unclaimedSessionRecords = initEvents(numel(uniqIDs),summaryFields);
        
        for jj = 1:numel(uniqIDs)
            unclaimedSessionRecords(jj).sessionID = sessionID{jj};
            iNeuronFiles = unclaimedNeuronFiles(crossIdxs==jj);
            if params.verbose && ~params.quiet
                fprintf('UNPROCESSED: %s - %d attached neurons...\n', uniqIDs{jj}, ...
                    numel(loadSpikeData(iNeuronFiles)));
            end
            % let's create an entry for these files now
            unclaimedSessionRecords(jj).manifest=[];
            unclaimedSessionRecords(jj).spikeFiles = iNeuronFiles;
        end
    end
    database{ii} = [sessionRecords unclaimedSessionRecords];
    if params.rejectNoNeuronSessions
        sessionsWithNoNeurons = arrayfun(@(x) isempty(x.spikeFiles), database{ii});
        database{ii}(sessionsWithNoNeurons) = [];
    end
end
if numel(database) == 1
    database = database{1};
end


function manifest = getManifest(files)
manifest = initEvents(0, {'name','originalFile'});
for iFile = 1:numel(files)
    thisFile = files{iFile};
    
    % this who step is slow, so we presave manifests
    vars = who('-file', thisFile); %variable name contained w/in each related file
    fileVars = struct('name',vars);
    % make the path relative
    [fileVars.originalFile] = deal(thisFile);
    % we don't know how many variables each file has, hence no
    % preallocation
    manifest = [manifest; fileVars];
end

function timestamp = getModifiedStamp(filename)
if ~iscell(filename)
    timestamp = getfield(dir(filename),'datenum');
else
    timestamp = cellfun(@(x) getfield(dir(x),'datenum'), filename);
end