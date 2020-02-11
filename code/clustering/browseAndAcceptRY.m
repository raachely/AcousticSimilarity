function [labels, augmentedLabels] = browseAndAcceptRY(birdID)
% clustering + interactive refinement/inspection; RY modified/commented
% based on browseAndAccept.m from BottjerLab repo.
% 
% modified from BottjerLab browseAndAccept: generates & saves
% acceptedLabels-[birdID]-[age].mat, which contains cluster indices
% structure that has accepted, regrouped, and augmented labels. 

% also generates & saves session-specific cluster label files - i.e. splits
% acceptedLabels into session-specific labels. This was originally done
% with separate call to standalone splitAcceptedLabels.m from Bottjerlab
% repo

% also generates & saves labelstereotypy-[sessionID].mat, which are
% intra- and inter-cluster distances. This was originally done with
% additional call to standalone saveStereotypy.m from BottjerLab repo. 

% also generates tutorSyllableCompare-[sessionID].mat, which has similarity
% comparisons between juvenile & tutor syllables
% 
% 03/2018 RY

close all;
dbstop if error
%% Load cluster file, syllables, & feature tables. 
plotParams = processArgs(defaultParams,...
    'dgram.minContrast', 1e-12, 'doFilterNoise', false, 'preroll', 3, 'postroll', 3);

dataDir    = [pwd filesep 'data' filesep            birdID filesep]; %directory where allSpecs-[birdID].mat file is
clusterDir = [pwd filesep 'data' filesep 'cluster-' birdID filesep]; %directory where cluster file (altClustDataAge ... .mat) is found; only used to set starting directory in uigetfile
[clusterfile, pathName] = uigetfile([clusterDir '*.mat'],'Select the file to grab clusters from'); % choose the cluster file. age to cluster is dictated by the selected file
[~, clusterfilename, ~] = fileparts(clusterfile);

% parse the age - dictated by selected cluster file
[~, remain] = strtok(clusterfile,'-'); %"delete" altClustDataAge part
thisAge = sscanf(remain(2:end), '%d'); %convert remaining file name (age + date/time stamp) to integer
thisAge = thisAge(1); %take just first value, which is the age

% load syllables+feature table
fprintf('Loading spectra & distance matrix for bird %s (age %d)...', birdID, thisAge);
allSpecs = load(fullfile(dataDir, strcat('allSpecs-', birdID, '.mat'))); %load allSpecs (spectra isn't needed for clustering, but is needed in compareTutorJuvenile call
%subselect the syllables from just this age:
isAge = ([allSpecs.DRsylls.age] == thisAge);
ageDRsylls = allSpecs.DRsylls(isAge);
agefeatureTable = allSpecs.featureTable(isAge);

%load distance matrix (use cosim)
altClustFile = load(fullfile(pathName, clusterfile)); %only distMats is used for clustering, but other variables (empMats,etc) are used in compareTutorJuvenile call
dists = squareform(altClustFile.distMats.cosim); %replot distance values as square matrix dists; dists(i,j) gives distance between syllable i & j

fprintf(' done loading.\n');

% Obsoleted - option to pick the kind of distance measure to use for clustering
%if isstruct(clusterIdxs) 
%    fn = fieldnames(clusterIdxs);
%    seldClusterType = fn{listdlg('ListString', fn, 'SelectionMode', 'single','Name', 'Which type of clustering?')};
%    if isempty(seldClusterType), error('Did not select option, buggin'' out.'); end
%    cIdxs = clusterIdxs.(seldClusterType);
%else
%    cIdxs = clusterIdxs;
%    seldClusterType = 'cosim';
%end
%maxes = max(cIdxs,[],1);
%dists = squareform(distMats.(seldClusterType)); %uses 'cosim' clustering type by default. distMats.cosim is always 1xN vector; squareform turns this into NxN square, w/ values reflected across a diagonal of 0's
%clear distMats; %not sure why this is specifically cleared - it's not used later.. but whatever.

%% Pick which level of clustering you want:
% Do clustering based on desired end number of syllables; pre-calculates
% all possibilities across range of max values (maxes)
maxes = 4:30;
linktree = linkage(dists, 'complete');
cIdxs = cluster(linktree,'maxclust',maxes); %cluster using max cluster values ranging maxes (i.e. 4-30); end with nX27 matrix, n = # of syllables & each column corresponds to cluster indices if using 4, 5, 6, etc as max cluster value

%ask user how many clusters they want (within range of maxes)
nDefClusts = 12; %default number of clusters, used only in user display prompt
nClusts = str2double(inputdlg(...
    sprintf('How many clusters are you looking for? (%d-%d)',min(maxes),max(maxes)),...
    'Number of clusters', 1, {num2str(nDefClusts)}));

%get the (pre-calculated) clustering based on user input for max clusters
thisSetPtr = find(maxes==nClusts);
thisSetIdxs = cIdxs(:,thisSetPtr); %cluster indices (labels) for each syllable, if clustering using the selected max # of desired clusters (nClusts)

%% Plot/edit clusters
% Loops through the different clusters in order of decreasing size (i.e.
% start with cluster that has most syllables, down to least). First plots &
% saves all clusters, then brings them up one by one & asks user to accept
% or edit
counts = hist(thisSetIdxs,1:nClusts); %i.e. how many syllables are in cluster 1, 2, 3... nClusts
[~, clustlabel_sorted] = sort(counts, 'descend'); %clustlabel_sorted(1) = number of the cluster (label) that has the most syllables

% plot all clusters first (just for visualization)
fprintf('Plotting & saving %d clusters (most->least syllables) ... \n', nClusts);
for jj = 1:nClusts
    fig = figure(jj);
    mosaicDRSpec(ageDRsylls(thisSetIdxs == clustlabel_sorted(jj)), plotParams, 'maxMosaicLength', 10.0);
    set(fig, 'Name',sprintf('Cluster #%d (%d/%d)',clustlabel_sorted(jj),jj,nClusts));
    figFileName = fullfile(clusterDir, strcat(clusterfilename, '-c', num2str(jj), '.jpg')); %specify .jpg to override .tif default in saveCurrFigure()
    saveCurrFigure(figFileName);
    close(fig);
end
fprintf('Look through plotted clusters. Press any key to continue on to approving/editing each one.\n')
pause

% inspect each cluster
labels = -ones(size(thisSetIdxs)); %column of -1s, length = # of syllables. this gets updated with approved, possibly new labels for each syllable.
isFused = false(1,nClusts); %at start, no clusters have been merged
grpCtr = 1; %new count for number of clusters
fprintf('Inspecting clusters ...')
for ii = 1:nClusts %for each cluster (starting from the one that has the most syllables, down to least)
    if isFused(clustlabel_sorted(ii)) %if this cluster has already been merged with another, skip it
        continue
    end
    
    clf  %clear figure window
    mosaicDRSpec(ageDRsylls(thisSetIdxs == clustlabel_sorted(ii)), plotParams, 'maxMosaicLength', 7.0); %plot the cluster
    set(gcf,'Name',sprintf('Cluster #%d (%d/%d)',clustlabel_sorted(ii),ii,nClusts));
    
    retry = true;
    while retry
        retry = false;
        switch nm_questdlg('How to treat this cluster?', ...
                sprintf('Cluster %d, # = %d', clustlabel_sorted(ii), sum(thisSetIdxs == clustlabel_sorted(ii))),...
                'Accept', 'Reject', 'Edit', 'Accept')
            case 'Accept'
                labels(thisSetIdxs == clustlabel_sorted(ii)) = grpCtr; %if accept all syllables in this cluster, label the syllalbes in this cluster w/ current cluster count
                grpCtr = grpCtr + 1;
            case 'Reject'
                labels(thisSetIdxs == clustlabel_sorted(ii)) = -1; %if reject, labels for syllables in this cluster = -1
                continue;
            case 'Edit'
                % todo: report metrics of conformity/similarity?
                [newLabelGroups, retry] = editSub(thisSetIdxs, clustlabel_sorted(ii));
                % implement labeling groups
                for jj = 1:numel(newLabelGroups)
                    labels(newLabelGroups{jj}) = grpCtr; %after cluster editing, update syllable labels w/ current cluster count
                    grpCtr = grpCtr + 1;
                end
        end
    end
end

labels(labels == -1) = NaN; %rejected clusters

%% Alternative editing
% directly choose clusters to merge:
regroupedLabels = mergeClustersByHand(ageDRsylls, dists, labels);

% match labels to established groups:
% only matters if there are unlabeled syllables (NaNs; ex. if a cluster is
% rejected above) - loop through unlabeled syllables & look for matches
% w/in range of acceptable distance; otherwise augmentedLabels will be
% the same as regroupedLabels
augmentedLabels = semiLabelSyllables(ageDRsylls, dists, regroupedLabels);

%% Save to an acceptedLabels file
clusterIdxs = struct('accepted', labels, 'regrouped', regroupedLabels, 'augmented', augmentedLabels);
saveFileName = fullfile(dataDir, strcat('acceptedLabels-', birdID, '-age', num2str(thisAge), '.mat'));
save(saveFileName, 'clusterIdxs');
fprintf('Done w/ clustering. Cluster identifiers saved to %s. \n', saveFileName);

%% Split accepted labels into session-specific files
allSyllLabels = clusterIdxs.regrouped; %change cluster type here (.approved, .regrouped, or .augmented)

[syllSessions, ~, syllSessIdxs] = unique({ageDRsylls.file}); %unique sessions (ageDRsylls here only contains single age already)
for jj = 1:numel(syllSessions)
    [~,syllSessions{jj}] = fileparts(syllSessions{jj});
end

fprintf('Splitting labels (%s, age %d) into session-specific files ...\n', birdID, thisAge); 
for jj = 1:numel(syllSessions) %for each session 
    acceptedLabels = allSyllLabels(syllSessIdxs == jj); %get cluster labels for this session
    fileName = fullfile(dataDir, strcat('acceptedLabels-', syllSessions{jj}, '.mat'));
    save(fileName, 'acceptedLabels');
end

%% Find & save intra- & inter-cluster distances
%fprintf('Total syllables IDed / total: (%d/%d, age %d)\n', sum(~isnan(allSyllLabels)), numel(allSyllLabels), thisAge);
% if sum(~isnan(allSyllLabels)) == 0
%     fprintf('No IDed syllables, exiting...\n');
%     return
% end

fprintf('Calculating + saving intra- & inter- cluster distances...\n');
% find intra- & inter- cluster dists
clusterMat = altClustFile.distMats.cosim;
[intraCDists, interCDists] = clusterDistances(clusterMat, allSyllLabels);

%split into sessions 
for ii = 1:numel(syllSessions)
    intraClusterDists = intraCDists(syllSessIdxs==ii);
    interClusterDists = interCDists(syllSessIdxs==ii);
    saveFile = fullfile(dataDir, strcat('labelStereotypy-', syllSessions{ii}, '.mat'));
    save(saveFile, 'intraClusterDists', 'interClusterDists');    
end
fprintf('Calculated & saved session-specific acceptedLabels & labelStereotypy files.\n');

%% Compare juvenile syll distances to tutor sylls
fprintf('Last step: tutor-juvenile distance calculations (%s, age %d)\n', birdID, thisAge)
compareTutorJuvenileRY(birdID, thisAge, allSpecs, allSyllLabels, altClustFile); %saves tutorSyllableCompare file
fprintf('All done! Saved acceptedLabels, labelStereotypy, & tutorSyllableCompare files for %s, age %d. \n', birdID, thisAge)

%% Cluster editing functions
function [newLabelGrps, tryAgain] = editSub(labels, oldLabel) %if choose to edit a cluster
    nThisType = sum(labels == oldLabel); %number of syllalbes in this cluster (/ with this cluster label)
    ttl = sprintf('Cluster %d, # = %d', oldLabel, nThisType); %title of dialog box
    
    openFigures = [];
    
    tryAgain = false;
    newLabelGrps = {};
    switch nm_questdlg('How to edit this cluster?', ttl, 'Merge', 'Split', 'Sort by Ear', 'Merge') %'merge' is default
        case 'Merge'
            % follow the tree backwards
            mergedCluster = findMerge(thisSetPtr, oldLabel, cIdxs); %mergedCluster = number of cluster that this one is to be merged with
            
            % now presentIden and clustIden should fuse
            fprintf('Fusing candidate...');
            
            %RY commented out - couldn't find where this fig ever existed
%             clustSession = strrep(clusterfile(1:end-5), 'altClustDataAge-',''); %is this a mistake? clustSession = age-datestamp_timestamp, but final number of time stamp is cut off
%             figureName = sprintf('%sexpClusters-%s-%s-c%d.jpg', pathName, clustSession, ...
%                 seldClusterType, mergedCluster); %presumably should plot syllables in mergedCluster, but don't think this figure is ever created anywhere
%             
%             if exist(figureName,'file') == 2
%                 fOp = figure; openFigures = [openFigures fOp];
%                 imshow(figureName);
%                 set(gcf,'Name',[figureName ' - candidate cluster for merging']);
%             else
%             end

            switch nm_questdlg('Fuse?',ttl,'Yes','No, reject all', 'Retry', 'Yes')
                case 'Yes'
                    isFused(oldLabel) = true;
                    isFused(mergedCluster) = true;
                    newLabelGrps = [find(labels == oldLabel | labels == mergedCluster) newLabelGrps]; %if accept merge, then both syllables w/ this cluster label + w/ mergedCluster label are included in newLabelGrps
                case 'Accept original'
                    newLabelGrps = [find(labels == oldLabel) newLabelGrps];
                case 'No, reject'
                    % do nothing
                case 'Retry'
                    tryAgain = true;
            end
        case 'Split'
            % how do we split?
            splitFeatOptions = ['tree'; fieldnames(agefeatureTable)];
            % follow the tree forwards
            
            repeatSplit = true;
            while repeatSplit % loop around try-catchtry for clustering robustness
                repeatSplit = false;
                [splitOpts, hitOk] = listdlg('ListString', splitFeatOptions, ...
                    'Name', 'Which features to split?'); %select feature(s) on which to split, &/or on 'tree'
                
                if ~hitOk
                    tryAgain = true;
                    closeAll(openFigures);
                    return;
                end
                
                if any(splitOpts == 1) %if select 'tree'
                    [typeA, typeB] = findSplit(thisSetPtr, oldLabel, cIdxs);
                else
                    % try a two-component mixture of gaussian clusters
                    nFeats = numel(splitOpts);
                    stats = zeros(nThisType, nFeats);
                    omitFeature = false(1, nFeats);
                    varStat = zeros(1,nFeats);
                    for kk = 1:nFeats
                        stats(:,kk) = [agefeatureTable(labels == oldLabel).(splitFeatOptions{splitOpts(kk)})];
                    end
                    varStat = var(stats);
                    omitFeature = (varStat < eps(max(varStat))*nThisType);
                    
                    if all(omitFeature), warning('All features are constant...'); end;
                    if any(omitFeature)
                        omitStr = '';
                        omitList = splitOpts(omitFeature);
                        
                        % create the join string
                        for kk = 1:numel(omitList),
                            omitStr = strcat(omitStr , splitFeatOptions(omitList(kk)));
                            if kk < nFeats, omitStr = strcat(omitStr,', '); end
                        end
                        omitStr = omitStr{1};
                        fprintf('Removing features %s...\n', omitStr);
                        
                        stats(:,omitFeature) = [];
                        splitOpts(omitFeature) = [];
                        nFeats = numel(splitOpts);
                    end
                    
                    try
                        fitObj = gmdistribution.fit(stats,2,'Regularize',0.0001);
                        newIdxs = cluster(fitObj,stats);
                        
                        % do some plotting
                        fTab = figure('Name','Clustergram');
                        openFigures = [openFigures fTab];
                        if size(stats,2) == 1
                            nBins = 40; bins = zeros(1,nBins);
                            bins(2:end-1) = linspace(min(stats),max(stats),nBins-2);
                            bins(end) = 2 * bins(end-1) - bins(end-2);
                            bins(1) = 2 * bins(2) - bins(3);
                            bar(bins,histc(stats, bins),1);
                            xlim(bins([1 end]));
                            hold on;
                            plot(bins, pdf(fitObj,bins'),'r-', 'LineWidth', 2);
                            hold off;
                            xlabel(splitFeatOptions{splitOpts}, 'Interpreter','none'); ylabel('Count');
                        else
                            % pick the two most diagnostic features
                            dprime = zeros(1,nFeats);
                            for kk = 1:size(stats,2)
                                dMu = diff(fitObj.mu(:,kk));
                                sumSigma = sqrt(sum(fitObj.Sigma(kk,kk,:)));
                                dprime(kk) = dMu/sumSigma;
                            end
                            [~,bestDims] = sort(dprime, 'descend'); bestDims = bestDims(1:2);
                            plot(stats(:,bestDims(1)), stats(:,bestDims(2)),'k.');
                            xlabel(splitFeatOptions(splitOpts(bestDims(1))), 'Interpreter','none');
                            ylabel(splitFeatOptions(splitOpts(bestDims(2))), 'Interpreter','none');
                        end
                        
                        % split done
                        subset = find(labels == oldLabel);
                        typeA = subset(newIdxs == 1);
                        typeB = subset(newIdxs == 2);
                    catch err
                        repeatSplit = questdlg(['Error in clustering: [', err.message, ']; try again?'],...
                            'Clustering Error', ...
                            'Yes','No','Yes');
                        repeatSplit = strcmp('Yes', repeatSplit);
                    end
                    
                end
                
            end
            % give option to view/sing/approve blind
            switch nm_questdlg(sprintf('How to review the split (A = %d,B = %d)?', ...
                    numel(typeA), numel(typeB)), ...
                    ttl,'View','Listen','Continue w/o review', 'View')
                case 'View'
                    totalLenA = sum([ageDRsylls(typeA).stop] - [ageDRsylls(typeA).start]);
                    totalLenB = sum([ageDRsylls(typeB).stop] - [ageDRsylls(typeB).start]);
                    
                    defaultVal = min(totalLenA, 4.0);
                    tPrev = nm_inputdlg(...
                        sprintf('Number of seconds to preview (max %.1fs): ', totalLenA),...
                        'Cluster A', 1, {num2str(defaultVal)});
                    
                    if isempty(tPrev), tryAgain = true; return;
                    else tPrev = str2double(tPrev); end
                    
                    fprintf('Plotting split cluster A (# = %d)...\n', numel(typeA));
                    figA = figure;
                    openFigures = [openFigures figA];
                    
                    mosaicDRSpec(ageDRsylls(typeA), plotParams, 'maxMosaicLength', tPrev);
                    set(gcf,'Name',sprintf('Cluster A (# = %d)',numel(typeA)));
                    
                    defaultVal = min(totalLenB, tPrev);
                    tPrev = nm_inputdlg(...
                        sprintf('Number of seconds to preview (max %.1fs): ', totalLenB),...
                        'Cluster B', 1, {num2str(defaultVal)});
                    
                    if isempty(tPrev), tryAgain = true; return;
                    else tPrev = str2double(tPrev); end
                    
                    fprintf('Plotting split cluster B (# = %d)...\n', numel(typeB));
                    figB = figure;
                    openFigures = [openFigures figB];
                    
                    mosaicDRSpec(ageDRsylls(typeB), plotParams, 'maxMosaicLength', tPrev);
                    set(gcf,'Name',sprintf('Cluster B (# = %d)',numel(typeB)));
                case 'Listen'
                    % simply play the sounds: it's faster to load but
                    % slower to be complete
                    fprintf('Splitting candidate... branch A audio: \n');
                    markRegions([],ageDRsylls(typeA));
                    fprintf('Splitting candidate... branch B audio: \n');
                    markRegions([],ageDRsylls(typeB));
                    
                case 'Accept w/o review'
            end
            
            choices = {'Both','Branch A', 'Branch B','None'};
            [respStr, isAccept] = nm_listdlg('Name', 'Which branches to accept (none is an option)?', 'ListString',choices, ...
                'SelectionMode', 'Single','OKString','Accept',...
                'CancelString','Retry');
            
            if ~isAccept
                tryAgain = true;
                closeAll(openFigures);
                return;
            end
            respStr = choices{respStr};
            switch respStr
                case 'Both'
                    % will flatten these later
                    newLabelGrps = {typeA, typeB};
                case 'Branch A'
                    newLabelGrps = {typeA};
                case 'Branch B'
                    newLabelGrps = {typeB};
                case 'None'
                    newLabelGrps= {};
            end
        case 'Sort by Ear'
            fprintf('Split into no more than 10..., based on audio: \n');
            newIdxs = multiMark([], ageDRsylls(labels == oldLabel));
            if any(isnan(newIdxs))
                if strcmp('Start over',questdlg(['Hand labeling is stopped early, accept or start over?'],...
                        'Labeling stopped', ...
                        'Accept','Start over','Start over'))
                    tryAgain = true;
                    closeAll(openFigures);
                    return
                end
            end
            [~,~,rIdxs] = unique(newIdxs);
            
            newLabelGrps = cell(1,numel(rIdxs));
            for kk = 1:numel(rIdxs)
                newLabelGrps{kk} = find(rIdxs == kk);
            end
    end
    closeAll(openFigures);
end
end

function closeAll(figH)
for kk = 1:numel(figH)
    close(figH(kk));
end
end

function [otherCluster, mergeDepth] = findMerge(grpPtr, origClust, clusterLabels)
isMerged = false;

origPtr = grpPtr;
clustTrav = origClust;
while ~isMerged && grpPtr > 1
    ct = crosstab(clusterLabels(:,grpPtr), clusterLabels(:, grpPtr-1));
    % the row that contains the current cluster - does it
    % contain another cluster?
    dstClusters = find(ct(clustTrav,:)>0);
    if sum(ct(:,dstClusters) > 0) > 1
        otherCluster = sum(find(ct(:,dstClusters))) - clustTrav;
        isMerged = true;
    else
        clustTrav = dstClusters;
        grpPtr = grpPtr - 1;
    end
end
mergeDepth = grpPtr;

% follow the new cluster back to the current division? (we
% don't have to do this, it might split again)
for jj = grpPtr:origPtr-1
    ct = crosstab(clusterLabels(:,jj), clusterLabels(:,jj+1));
    [~, otherCluster] = max(ct(otherCluster,:));
end
end

function [clusterInds, splitInds, splitDepth] = findSplit(grpPtr, origClust, clusterLabels)
hasSplit = false;
clustTrav = origClust;

nClusterings = size(clusterLabels, 2);
while ~hasSplit && grpPtr < nClusterings
    ct = crosstab(clusterLabels(:,grpPtr), clusterLabels(:, grpPtr+1));
    % the row that contains the current cluster - does it
    % contain another cluster?
    dstClusters = find(ct(clustTrav,:) > 0);
    
    if numel(dstClusters) == 2
        clustTrav = dstClusters(1);
        otherCluster = dstClusters(2);
        hasSplit = true;
    else
        clustTrav = dstClusters;
    end
    grpPtr = grpPtr + 1;
end
if hasSplit
    clusterInds = find(clusterLabels(:,grpPtr) == clustTrav);
    splitInds = find(clusterLabels(:,grpPtr) == otherCluster);
    splitDepth = grpPtr;
else
    clusterInds = find(clusterLabels(:,grpPtr) == clustTrav);
    splitInds = [];
    splitDepth = nClusterings; % excepted
end
end