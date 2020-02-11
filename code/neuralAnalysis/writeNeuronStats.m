function writeNeuronStats(birdIDs)
% writes the average firing rates and response strengths for each day where 
% both neurons and clusters of syllables have be identified

% birdIDs must be a cell:
if ~iscell(birdIDs), birdIDs = {birdIDs}; end
    
nBirds = numel(birdIDs);
params = defaultParams;
%%
for hh = 1:nBirds
    thisBird = birdIDs{hh};
    sessionData = reportOnData(thisBird, '', [],'verbose',false);
    sessions = {sessionData.sessionID};
    
    ages = getAgeOfSession(sessions);
    uAges = unique(ages(~isnan(ages)));
    
    for gg = 1:numel(uAges)       
        thisAge = uAges(gg);        
        birdDir    = [pwd filesep 'data' filesep thisBird filesep];
        
        % load syllables and corresponding feature information
        [DRsylls, sForSylls] = loadAgeSylls(thisBird, thisAge);
        uSessions = unique(sForSylls);
        nSessions = numel(uSessions);

        % if we have acceptedLabels, then use those
        nTotalSylls = numel(DRsylls);
        try
            clusterIdxs = loadAcceptedLabels(thisBird, thisAge);
        catch excep
            if strcmp(excep.identifier,'loadAcceptedLabels:fileNotFound')               
                if nTotalSylls < 1000
                    fprintf(['Not enough syllables (%d) for labels for %s '...
                    'age %d, going to next age...\n'],nTotalSylls, thisBird, thisAge);
                    continue;
                else
                    warning('writeNeuronStats:NoLabels','%s should have labels at age %d (%d syllables)\n',...
                        thisBird, thisAge, nTotalSylls);
                end
            else
                rethrow(excep);
            end
        end
        
        for ii = 1:nSessions % loop through sessions
            thisSession = uSessions{ii};
            
            matchingRecord = sessionData(strcmp(thisSession, sessions));
            matchingSpikeFiles = matchingRecord.spikeFiles;
            %[spikes, nNeuronsPerFile, isMUA] = loadSpikeData(matchingSpikeFiles);
            [spikes, nNeuronsPerFile] = loadSpikeData(matchingSpikeFiles); %RY changed to ignore MUA stuff (for now?)
            if isempty(spikes) 
                fprintf('Session %s has no neurons associated, continuing to next session...\n', thisSession);
                continue;
            end
            
            % get core/shell identity
            isCoreUnit = false(numel(spikes),1);
            cumIndex = [0 cumsum(nNeuronsPerFile)];
            for jj = 1:numel(matchingSpikeFiles)
                unitsFromFile = (cumIndex(jj)+1):cumIndex(jj+1);
                
                % core/shell ness is indicated by the file folder
                isCoreUnit(unitsFromFile) = ~isempty(strfind(matchingSpikeFiles{jj}, 'core'));
            end
            % get the plasticness of the session
            %[~, thisSessionIsPlastic, thisPlasticScore] = getAgeOfSession(thisSession); %RY commented out to ignore plastic-ness stuff (for now?)
            
            isThisSession = strcmp(sForSylls, thisSession);
            syllables = DRsylls(isThisSession);
            clusterNum = clusterIdxs(isThisSession);
            
            % transfer the labels to the age-specific syllables
            foo = num2cell(clusterNum); [DRsylls(isThisSession).type] = foo{:};
            
            % treat NaNs as a separate group
            hasNanGroup = false;
            if any(isnan(clusterNum))
                hasNanGroup = true;
                clusterNum(isnan(clusterNum)) = max(clusterNum)+1;
            end
            
            % each syllable cluster will have its own column given by getRS
            
            neuronData = struct([]); % empty initialization is simpler
            clusterTypes = unique(clusterNum);
            for jj = 1:numel(clusterTypes) % loop through syllable types
                thisIdx = find(clusterNum == clusterTypes(jj));
                thisCluster = syllables(thisIdx);
                
                % note: in the case there are no syllables of this type, we
                % still want to generate the empty structure with getRS
                
                % internally name the syllables - this is important and should
                % correspond to getRS
                if ~isempty(thisCluster)
                    [thisCluster.type] = deal('syllable');
                end
                
                % RY moved baseline calculation period to segmentAndCluster.m, because it requires calling original song structure file 
                baselinePeriodFile = [birdDir 'baselinePeriods-' thisSession '.mat'];
                load(baselinePeriodFile);
%                 if ~exist(baselinePeriodFile, 'file')
%                     songStruct = loadFromManifest(sessionData.manifest,'songStruct');
%                     params.fs = 1/songStruct.interval;
%                     aS = loadFromManifest(sessionData.manifest, 'approvedSyllables');
%                     baselinePeriodsS = constructBaselineSilence(songStruct, 2, aS, params);
%                     %save(baselinePeriodFile,'baselinePeriods','baselinePeriodsS');
%                     save(baselinePeriodFile,'baselinePeriodsS');
%                 else
%                     load(baselinePeriodFile);
%                 end
                [baselinePeriodsS.age] = deal(thisAge);
                [baselinePeriodsS.file] = deal('whatever,dude'); % this doesn't matter
                [baselinePeriodsS.type] = deal('baseline');
                
                %adding 50ms premotor - fixed lag for determining significance
                syllsWPre = addPrePost(thisCluster', params, 'preroll', 50, 'postroll', 0); 
                
                %this put all neurons into neuronData for each syllable type, rewrite so each neuron gets responses to all syllable types
                %if there are no syllables in the cluster, this returns an
                %empty data 
                [neuronData{jj}, rawFiring] = getRS([syllsWPre; baselinePeriodsS], spikes, ... % 3/1 check this
                    'syllable', 'baseline', ...
                    'p_ttest', 'meanRS', 'meanRSI', 'standRS',...
                    'cvarFR',params, 'verbose',false);
                % first column is firing during syllable , second column is baseline firing
                [neuronData{jj}.rawRates] = rawFiring{:}; 
                foo = num2cell(isCoreUnit); [neuronData{jj}.isCore] = foo{:};
                %foo = num2cell(isMUA); [neuronData{jj}.isMUA] = foo{:}; %RY commented out to ignore MUA stuff (for now?)
                [neuronData{jj}.syllIndex] = deal(thisIdx);
                %[neuronData{jj}.isPlastic] = deal(thisSessionIsPlastic); %RY commented out to ignore plastic-ness stuff (for now?)
                %[neuronData{jj}.plasticScore] = deal(thisPlasticScore);
                
                % JMA adding in burst fraction
                BF = cell(length(spikes),1);
                for kk = 1:length(spikes)
                    [~, syllableBF] = getISIs(spikes{kk},syllsWPre);
                    BF{kk} = syllableBF;
                end
                [neuronData{jj}.burstFraction] = BF{:};
                
                % remember, we reassigned the max+1 value to be NaN and
                % flipped a flag
                if hasNanGroup && clusterTypes(jj) == max(clusterNum)
                    [neuronData{jj}.syllID] = deal(NaN);
                else
                    [neuronData{jj}.syllID] = deal(clusterTypes(jj));
                end                
            end
            %end syllable type loop            
           
            % compile all the neurons together
            neuronSyllableData = [neuronData{:}]; % rows are neurons, columns are syllables
            fil = [birdDir 'neuronSyllable-' thisSession '.mat'];
            fprintf('Writing %d neurons, %d syllables into file %s...\n', size(neuronSyllableData,1),...
                size(neuronSyllableData,2),fil);
            save(fil, 'neuronSyllableData');
        end %end session loop
    end % end ages loop
end % end birds loop
fprintf('Done!\n');
% notes from a while back:
    % Do neurons fire to one syllable type?
    % I want to look at PSTHs of neurons firing to each syllable type
    % Where are the significant peaks and troughs in the PSTHs
    
    % Do neurons fire more when there is comparatively more variability in a
    % Syllable compared to the other syllables of its type?
    % Do neurons fire more similarly temporally for more similar syllables?
    
    % Do the neurons fire randomly during a syllable type or in a time-aligned
    %manner? cross-correlations might suck i know
    
    
  