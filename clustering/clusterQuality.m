function [DB, d_i] = clusterQuality(clusterMat, clusterIds)
% DB returns the davies-bouldin index for the entire clustering scheme
% d_i returns the individual cluster qualities
if isvector(clusterMat)
    clusterMat = squareform(clusterMat);
end
% clusterIds must be numeric, not all NaN

nTypes = max(clusterIds);

centerIdx = zeros(1,nTypes);
sigma = zeros(1,nTypes);
for ii = 1:nTypes %loop through juv syll clusters
	thisType = find(clusterIds == ii);
	subMat = clusterMat(thisType, thisType); %get distance values for current syll cluster, i.e. distances between all syll pairs w/in this cluster

	% find the average distance of every element to its centroid
	%[sigma(ii), minInType] = min(mean(subMat,1));
    
    % find the RMS distance of every element to its centroid
    % NB: the centroid is defined as the element for which the rms distance 
    % from every other element is minimal
	[sigma(ii), minInType] = min(mean(subMat.^2,1)); %sigma = smallest mean square distance value
    sigma(ii) = sqrt(sigma(ii)); %sqrt -> sigma = smallest RMS distance value (root mean distance of centroid)
	centerIdx(ii) = thisType(minInType); %centerIdx = index of the syllable that had the smallest mean square distance (i.e. syllable that was, on average, closest to all other syllables in the cluster)
end

% get the distances between centroids
centD = clusterMat(centerIdx, centerIdx); %centD = NxN matrix of distances between cluster centroids (ex. column 1 = distance between centroid of cluster 1 & centroid of cluster 2, of cluster 3, etc ... cluster n)
	
% find the d_is
d_i = zeros(1,nTypes);
for ii = 1:nTypes
	notI = 1:nTypes; notI(ii) = [];
	d_i(ii) = mean((sigma(ii) + sigma(notI)) ./ centD(ii, notI));
end
DB = mean(d_i);
end
