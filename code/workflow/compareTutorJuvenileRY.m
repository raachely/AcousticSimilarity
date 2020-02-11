function compareTutorJuvenileRY(birdID, age, allSpecs, clusterIdxs, normDists)
% get tutor song syllables distance to juvenile syllables, taking into
% account both local and global distances.

% allSpecs = loaded allSpecs-[birdID].mat file (structure containing DRsylls, etc)
% clusterIdxs = mtarix of (approved) cluster identifiers/labels
% normDists = loaded altClustDataAge ... .mat file (structure containing distMats, etc)

% RY modified compareTutorJuvenile.m from BottjerLab repo so that it could
% be called as function instead of separate script

fprintf('Calculating tutor-juv distance comparisons for bird %s, age %d...\n', birdID, age)

%% Load/assign variables
dataDir  = ['data' filesep            birdID filesep];

% juv syllable info
isAge = [allSpecs.DRsylls.age] == age;
juvSylls = allSpecs.DRsylls(isAge);
juvSpectra = allSpecs.spectra(isAge);
juvFeatures = allSpecs.featureTable(isAge);
% get the session ID for each syllable
juvSessions = cell(1,numel(juvSylls)); 
for ii = 1:numel(juvSylls) 
    [~,juvSessions{ii}] = fileparts(juvSylls(ii).file);
end
 
% tutor syllable info
tutorFil = fullfile(dataDir, strcat('tutor-', birdID, '.mat'));
if exist(tutorFil,'file') ~= 2
    error('compareTutorJuvenile:tutorFileNotFound',...
        'Tutor file %s for bird %s not found...',tutorFil, birdID);
end
load(tutorFil);

%% Calculate distances between tutorSpecs/tutorFeatures and juv DRspectra/DRfeatures 
% that have labels only within the correct matchings
% goes through the same steps as in DRcluster

% step 0: z-normalize the juvenile syllable statistics
fn = fieldnames(juvFeatures);
juvFeaturesTabled = cellfun(@(x) [juvFeatures.(x)]', fn', 'UniformOutput',false);
juvFeaturesTabled = [juvFeaturesTabled{:}];
[zNormedTable, featMu, featSigma] = zscore(juvFeaturesTabled);

% step 1: get local and global distances to each other syllable
nTutor = numel(tutorSylls);
nJuv = numel(juvSylls);
localDist = struct('raw', zeros(nJuv, nTutor), ...% initialization
                    'emp', zeros(nJuv, nTutor),...
                    'co', zeros(nJuv, nTutor));
globalDist = localDist; fusedDist = rmfield(localDist, 'emp'); 
 localEmp = squareform(normDists.empMats.warpedLocal);
globalEmp = squareform(normDists.empMats.global);
 fusedEmp = squareform(normDists.distMats.cosim);
progressbar('tutor syllables','juvenile syllables');
for ii = 1:nTutor
    tutorFeats = ((cellfun(@(x) tSFeats(ii).(x), fn') - featMu) ./ featSigma);
    tutorFeats(featSigma == 0) = 0;
    lenTut = tutorSylls(ii).stop - tutorSylls(ii).start;
    for jj = 1:nJuv
        progressbar([],jj/nJuv);
        lenJuv = juvSylls(jj).stop - juvSylls(jj).start;
        % local distance is normalized by the average length of syllables
        localDist.raw(jj, ii) = timeWarpedDistance(tSSpecs(ii), juvSpectra(jj)) / mean([lenTut, lenJuv]);
    end
    globalDist.raw(:, ii) = pdist2(tutorFeats, zNormedTable, 'euclidean');    
    localDist.emp(:, ii) = interp1(normDists.empDistrs.warpedLocal(2,:), normDists.empDistrs.warpedLocal(1,:), ...
        localDist.raw(:, ii));
    
    % no great way to extrapolate if distance is too far, so just divide by
    % the max for now
    inExcess = find(localDist.raw(:,ii) > normDists.empDistrs.warpedLocal(2,end)); 
    if numel(inExcess) > 0
        localDist.emp(inExcess,ii) = localDist.raw(inExcess,ii) / normDists.empDistrs.warpedLocal(2,end);        
    end
    globalDist.emp(:, ii) = interp1(normDists.empDistrs.global(2,:), normDists.empDistrs.global(1,:), ...
        globalDist.raw(:, ii));
    % no great way to extrapolate if distance is too far, so just divide by
    % the max for now
    inExcess = find(globalDist.raw(:,ii) > normDists.empDistrs.global(2,end)); 
    if numel(inExcess) > 0
        globalDist.emp(inExcess,ii) = globalDist.raw(inExcess,ii) / normDists.empDistrs.global(2,end);        
    end
    if any(isnan(globalDist.emp(:,ii))), keyboard; end;
    localDist.co(:, ii) = pdist2( localDist.emp(:, ii)', localEmp, 'correlation');
    globalDist.co(:, ii) = pdist2(globalDist.emp(:, ii)', globalEmp, 'correlation');
    progressbar(ii/nTutor, 0);
end

fusedDist.raw = sqrt(localDist.co .* globalDist.co);
fusedDist.co  = pdist2(fusedDist.raw', fusedEmp, 'correlation')'; % each tutor syllable against each juvenile syllable

%% For each juvenile syllable, average the distance for each type of the tutor syllable
% there are multiple instances of the tutor syllable
[tutorTypes,~,matchesTutor] = unique({tutorSylls.type});
nTutorTypes = numel(tutorTypes);
distanceToTutor = NaN(nTutorTypes, nJuv);

for ii = 1:nTutorTypes   
   distanceToTutor(ii, :) = mean(fusedDist.co(:, ii == matchesTutor),2)';
end

%% For each juvenile syllable, get different scores to tutor syllable:
% bestDistanceScore <- distance from individual juv syllable to closest tutor syll 
% distanceToConsensus <- distance from individual juv syllable to (tutor syll closest to all
% syllables of that juv syllable's cluster) 
% distanceToCentral <- distance from individual juv syllable to (tutor syll
% closest to the most central syllable of that juv syllable's cluster)
% distanceToHumanMatch <- distance to the matched tutor as described by the
% matchToTutor file

% best distance to any tutor syllable
[bestDistanceScore, bestInd] = min(distanceToTutor, [],1); %distanceToTutor is mxn matrix of distance scores, m = # tut syllable types, n = # syllables; min = i.e. for each juv syllable, which TUT syllable (index) is it closest to (min distance)?
bestTutorMatch = tutorTypes(bestInd); %'convert' from index to TUT syllable label (1->a, 2->b, etc)

 distanceToConsensus  = NaN(1, nJuv);
 distanceToCentral    = NaN(1, nJuv);
% distanceToHumanMatch = NaN(1, nJuv);
 
 consensusMatch = NaN(1, nJuv);
 centralMatch   = NaN(1, nJuv);
% humanMatch     = NaN(1, nJuv);
 
% %try to load human matches for a type
% typeMatch = [];
% matchFile = [dataDir 'matchToTutor-age' num2str(age) '.mat'];
% if exist(matchFile, 'file') == 2
%     load(matchFile); % loading typeMatch, cell array of chars;
%     % this is generated by assignClusterToTutor
% end

 if ~isempty(clusterIdxs)
     % distanceToCluster
     % first find the closest syllable to each cluster, by summing distances 
     nClusters = nanmax(clusterIdxs);
     centralOrder = findMostCentral(normDists.distMats.cosim, clusterIdxs); %centralOrder = list of some measure for each syllable, split/organized by juvenile syllable type
     
    for ii = 1:nClusters
        inThisCluster = find(clusterIdxs == ii);
        thisCentralIdx = centralOrder{ii}(1);
        
        distancesToThisCluster = sum(distanceToTutor(:,inThisCluster), 2); %sum of distance scores for all juv sylls in a cluster to each tut syll
        [~,closestByConsensus] = min(distancesToThisCluster); %pick lowest summed distance score - total distance score across juv sylls in the cluster is lowest for tut syll x
        distanceToConsensus(inThisCluster)  = distanceToTutor(closestByConsensus, inThisCluster);
        consensusMatch(inThisCluster) = closestByConsensus;
        
        %for each juv syllable type, select 1 rendition of that syllable type (thisCentralIdx; based on results of findMostCentral) to use in determining which TUT syllable that juv cluster is closest to        
        [~,closestByCentral] = min(distanceToTutor(:,thisCentralIdx)); %closestByCentral = index of TUT syllable; 1 would correspond to syllable a, 2 to b, etc. 
        distanceToCentral(inThisCluster)    = distanceToTutor(closestByCentral, inThisCluster);
        centralMatch(inThisCluster) = closestByCentral;
        
%         if ~isempty(typeMatch)
%             closestByHumanMatch = typeMatch{ii} - 'a' + 1;
%             distanceToHumanMatch(inThisCluster) = distanceToTutor(closestByHumanMatch, inThisCluster);
%             humanMatch(inThisCluster) = closestByHumanMatch;
%         end
     end
     % get the closest syllable to every cluster       
 end

%% split scores among sessions
[uJSessions,~, rIdx] = unique(juvSessions);
for ii = 1:numel(uJSessions)
    thisSession = uJSessions{ii};
    isInSession = (rIdx==ii);
    
    distToTutor = distanceToTutor(:,isInSession);
    bestDistScore = bestDistanceScore(isInSession);
    bestTutMatch = bestTutorMatch(isInSession);
    
     distToConsensus  = distanceToConsensus(isInSession); 
     distToCentral    = distanceToCentral(isInSession);
%     distToHumanMatch = distanceToHumanMatch(isInSession);
%     
     consMatch = consensusMatch(isInSession);
     centMatch = centralMatch(isInSession);
%     humMatch  = humanMatch(isInSession);
%     %the juvenile dimension is the last one in theis savefile
    saveFile = [dataDir 'tutorSyllableCompare-' thisSession '.mat'];
    fprintf('Saving session %s tut-syll comparison data to %s...\n', thisSession, saveFile);
%     save(saveFile, 'distToTutor', 'bestDistScore', 'bestTutMatch',...
%         'distToConsensus','distToCentral','distToHumanMatch',...
%         'consMatch', 'centMatch', 'humMatch');
    save(saveFile, 'distToTutor', 'bestDistScore', 'bestTutMatch', ...
        'distToConsensus', 'distToCentral', 'consMatch', 'centMatch');
end
