# AcousticSimilarity
Segmenting syllables & calculating/comparing similarity between them 

# General Info
## Basics
What you'll need:  
1. Song channel exported from Spike2 to Matlab, named [sessionID]-spike2.mat
    * `sessionID` variable should be a unique identifier for each file that starts with the bird ID and has no dashes within it (underscores OK).
2. Song bout times exported from Spike2 to Matlab, named motifReturn-[sessionID].txt
3. Tutor song .wav  

To run:  
1. segmentAndCluster.m: segment juvenile & tutor syllables (doesn't cluster)
    * segmentAndCluster_auto does the same, but without manual refinement, + able to loop through multiple files
2. prepareDRsylls.m: calculate syllable spectra and cluster
    * original fullClusterByDayGy242.m was merged into this
    * calculates spectra & clusters syllables from all sessions within each day, so run this only after segmenting syllables from all sessions within an age. even more efficient to run segmentAndCluster on _all_ recordings from a bird, and then run prepareDRsylls just once to cluster all of them (within each age) at once
3. browseAndAccept.m: refine clustering
4. writeNeuronStats.m: 
5. compareTutorJuvenile.m: compare juv vs TUT syllables

## Parameters
`defaultParams.m` holds default parameters for a variety of functions   
`params = defaultParams` to load defaults into `params` variable. 
 
`processArgs.m` to load current parameters and change/overwrite individual parameters  
```
new_params = procesArgs(defaultParams, '[param_x]', [new_x_value], '[param_y]', [new_y_value]); %to specify new values for parameters x & y
this_fx_params = processArgs(params, varargin{:}); %to process multiple params fed into a function 
```

`defineAliases.m` can be used to load a set of parameters, rather than loading defaults and then individually specifying/changing values  

## Event structures
Regions of elapsed time - syllables, motifs, baseline periods, etc - are all considered "events"; these are always stored in the same event structure so that there is a consistent reference syntax for extracting information about a given event (time period): Nx1 structure where N = number of events, with 'type', 'start', 'stop', 'idxStart', 'idxStop', and 'file' as field names. Use `initEvents` to initialize an empty one.

## Potentially useful utilities
`reportOnData` returns a 1xN cell array, where N is the number of birds (bird ID folders in the data folder). Each cell contains a structure that lists the session ID, manifest (structure listing all the related files for each session), and spikeFiles (files w/ neuron spiking time points)

