function [neuronData, pairedFiring, pairedBF] = getRS(events, spikes, keyType1, keyType2, params, varargin)
%function [neuronData, pairedFiring] = getRS(events, spikes, keyType1, keyType2, ptype, mRSType, mRSIType, standType, cvarFRType, params, varargin)
%
% this compares the firing rates between events of key type 1 and events of
% key type 2
%
% spikes is a cell array of different neurons' spike times
% ptype is the field name of the probability (by t-test) that the syllables
% are different
% mRStype is the field name of the response strength (difference in firing)
% mRSItype is the field name of the response strength normalized by the sum
% standType is the field name of the response strength index normalized by
% the standard deviation
%
% warning: this method only takes string types
%
% pairedFiring is a cell array of the firing rate of each event and its paired baseline
%
% RY changed to hard-code ptype, MRSType, etc arguments. seemed unnecessary
% to have to feed in arguments for them every time just to use as labels
%

if nargin <= 4 || isempty(params) %RY changed from 9
    params = defaultParams;
end
if nargin > 4 %RY changed from 9
    params = processArgs(params, varargin{:});
end

nNeurons = numel(spikes);

%syllLengths = [events.stop] - [events.start]; Jenny commented out
isKey1 = strcmp(keyType1, {events.type});
isKey2 = strcmp(keyType2, {events.type});
key1 = events(isKey1);
key2 = events(isKey2);

% if there are no intervals that are labeled with the correct key
if sum(isKey1) < 1 || sum(isKey2) < 1
    %fieldNames = {['FR_' keyType1],           ['FR_' keyType2],             ptype, mRSType, mRSIType, standType, cvarFRType};
    fieldNames = {['FR_' keyType1],           ['FR_' keyType2],             'p_ttest', 'pBF_ttest', 'meanRS', 'meanRSI', 'standRS', 'cvarFR'}; %RY change to hard-code these field names
    kosherFieldNames = strrep(fieldNames, ' ', '_');
    neuronData = initEmptyStructArray(kosherFieldNames, nNeurons);
    for ii = 1:numel(kosherFieldNames)
        [neuronData.(kosherFieldNames{ii})] = deal(NaN);
    end
    pairedFiring = cell(1,nNeurons);
    return
end

meanRS = NaN(1,nNeurons);
meanRSI = NaN(1,nNeurons);
standRS = NaN(1,nNeurons);
cvarFR = NaN(1,nNeurons);
pbase = NaN(1,nNeurons);
matchBase = cell(1,nNeurons);
pairedFiring = cell(1,nNeurons);
matchBaseBF = cell(1,nNeurons); %for BF; RY added
pBF = NaN(1,nNeurons); %for BF; RY added
pairedBF = cell(1,nNeurons); %for BF; RY added

for ii = 1:nNeurons %for each neuron
    % count the spikes within each event
    [~,~, spikeRates1{ii}] = countSpikes(key1, spikes{ii});
    [~,~, spikeRates2{ii}] = countSpikes(key2, spikes{ii});
    % to do: add SEMs to output
    
    if all(spikeRates1{ii} == 0) && all(spikeRates2{ii} == 0) % if no spikes in either type of event
        pairedFiring{ii} = zeros(2,sum(isKey1)); % the default for each paired firing vector
        continue;
    end
    % we can't fill out the fields if there's no firing ever - there should
    % be a more specific place for this
     
    %[~,pbase(ii)]=ttest2(spikeRates1{ii}, spikeRates2{ii});
    
    %burst fraction for each event. RY added
    [~, BF1{ii}] = getISIs(spikes{ii}, key1);
    [~, BF2{ii}] = getISIs(spikes{ii}, key2);
    
    %Response strength (FR during event1 - FR during event2) - uses average
    %of nearest two baseline events (event2) as baseline (matchBase). This
    %average is matched with its corresponding key1 event in pairedFiring.
    %JMA wrote.
    FR = zeros(1,length(key1));
    RS = NaN(1,length(key1));   
    RSI = zeros(1,length(key1));
    matchBase{ii} = zeros(1,length(key1)); %the events of baseline (event 2) that match the original event
    matchBaseBF{ii} = zeros(1,length(key1)); %RY added for baseline burst fraction
    for jj = 1:length(key1) %for each key1 event
        startList = abs(vertcat(key2.start) - key1(jj).start); %baselines closest to event start
        stopList =  abs(vertcat(key2.start) - key1(jj).stop); %baselines closest to event stop
        [~, bothInd] = sort([startList; stopList]);
        bothInd(bothInd > length(startList)) = bothInd(bothInd > length(startList)) - length(startList);
        base1 = bothInd(1);
        if bothInd(2) == bothInd(1) %don't use the same baseline event twice
            base2 = bothInd(3);
        else
            base2 = bothInd(2);
        end
        matchBase{ii}(jj) = (spikeRates2{ii}(base1) + spikeRates2{ii}(base2)) / 2;
        matchBaseBF{ii}(jj) = (BF2{ii}(base1) + BF2{ii}(base2))/2;
        
        FR(jj) = spikeRates1{ii}(jj);
        RS(jj) = spikeRates1{ii}(jj) - matchBase{ii}(jj);
        RSI(jj) = (spikeRates1{ii}(jj) - matchBase{ii}(jj)) / (spikeRates1{ii}(jj) + matchBase{ii}(jj));
        if isnan(RSI(jj)), RSI(jj) = 0; end
    end
    pairedFiring{ii} = [spikeRates1{ii}; matchBase{ii}];
    pairedBF{ii} = [BF1{ii}'; matchBaseBF{ii}]; %RY added
    
    % t-test to distinguish spiking from baseline in events
    [~,pbase(ii)]=ttest(spikeRates1{ii}, matchBase{ii});
    [~,pBF(ii)] = ttest(BF1{ii}',matchBaseBF{ii}); %RY added
    
    %coeffient of variation of key1 event
    cvarFR(ii) = std(FR)/mean(FR); 
    
    % response strengths: raw (RS), normalized to sum of response (RSI), and
    % standardized to variance of response (standRS)
    meanRS(ii) = mean(RS);
    meanRSI(ii) = mean(RSI);
     
    covarSB = cov(spikeRates1{ii},matchBase{ii}); % get the covariance of the spikes, for standRS
    if numel(covarSB) == 1 % if no spikes ever occur in either baseline regions or ROI, then the covariance matrix will be scalar
        covarSB = 0;
    else
        covarSB = covarSB(2,1);
    end 
    standRS(ii) = sqrt(numel(key1)) * (mean(spikeRates1{ii}) - mean(matchBase{ii}))...
        /sqrt(var(spikeRates1{ii}) + var(matchBase{ii}) - 2 * covarSB); %JMA wrote; singing - baseline, in that order
    
    % the standardize response could be 0 or infinite if the variance is poorly behaved
    if isinf(standRS(ii)) || isnan(standRS(ii))
        standRS(ii) = 0.0;
    end
    
    % report
    if params.verbose
        fprintf('Neuron %d: ''%s'' rate = %0.3f Hz, matched ''%s'' rate = %0.3f Hz, p = %0.2f, RSI = %0.2f, standRS = %0.2f \n', ...
            ii, keyType1, mean(spikeRates1{ii}), keyType2, mean(matchBase{ii}), pbase(ii), meanRSI(ii), standRS(ii));
    end
end

%%% plating
% package data for output argument

% NB: the p-value here is unsigned
dataTable =   [cellfun(@mean,spikeRates1); cellfun(@mean, matchBase);   pbase;  pBF; meanRS; meanRSI; standRS; cvarFR];
%fieldNames = {['FR_' keyType1],           ['FR_' keyType2],             ptype, mRSType, mRSIType, standType, cvarFRType};
fieldNames = {['FR_' keyType1],           ['FR_' keyType2],             'p_ttest', 'pBF_ttest', 'meanRS', 'meanRSI', 'standRS', 'cvarFR'}; %RY change to hard-code these field names
kosherFieldNames = strrep(fieldNames, ' ', '_');
neuronData = initEmptyStructArray(kosherFieldNames, nNeurons);
for ii = 1:numel(kosherFieldNames)
    foo = num2cell(dataTable(ii,:));
    [neuronData.(kosherFieldNames{ii})] = foo{:};
end