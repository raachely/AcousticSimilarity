# AcousticSimilarity
Segmenting syllables & calculating/comparing similarity between them 

Modified from original BottjerLab repo (https://github.com/BottjerLab/Acoustic_Similarity), in attempts to consolidate/clarify analysis pipeline 

## Basics
What you'll need:  
1. Song channel exported from Spike2 to Matlab, named [sessionID]-spike2.mat
    * `sessionID` variable should be a unique identifier for each file that starts with the bird ID and has no dashes within it (underscores OK).
2. Song bout times exported from Spike2 to Matlab, named motifReturn-[sessionID].txt
3. Tutor song .wav  

Set up a root folder titled "data" with subfolders for each bird. The -spike2.mat and motifReturn- files should be located within the appropriate bird folder. If analyzing spike data, files containing spike times should be located within a subfolder within each bird folder and named [sessionID]_times.mat. Ex:  

* data  
   * Gy100  
      * Gy100_file1-spike2.mat
      * motifReturn-Gy100_file1.txt
      * Gy100_file2-spike2.mat
      * motifReturn-Gy100_file2.txt
      * spike files
         * Gy100_file1_times.mat
         * Gy100_file2_times.mat
         
   * Gy101
      * Gy101_file1-spike2.mat
      * motifReturn-Gy101_file1.txt  

Single spike files can contain spike times for multiple units from a single session. There are no restrictions on where the tutor song .wav file needs to be. Your current directory should be set to the 'data' folder when running these analyses.
These restrictions stem from  `reportOnData.m` and `loadSpikeDAta.m`, which are functions that find associated files/metadata for a given bird/session. Modify these accordingly if your data structure differs.

To run:  
1. segmentAndCluster.m: segment juvenile & tutor syllables (doesn't cluster)
    * segmentAndCluster_auto does the same, but without manual refinement, + able to loop through multiple files
2. prepareDRsylls.m: calculate syllable spectra and cluster
    * original fullClusterByDayGy242.m from BottjerLab repo was merged into this
    * calculates spectra & clusters syllables from all sessions within each day, so run this only after segmenting syllables from all sessions within an age. even more efficient to run segmentAndCluster on _all_ recordings from a bird, and then run prepareDRsylls just once to cluster all of them (within each age) at once
3. browseAndAccept.m: refine clustering & compare juvenile vs tutor syllables
4. writeNeuronStats.m: run basic calculations of spiking activity during syllable types/clusters

## Event structures
Regions of elapsed time - syllables, motifs, baseline periods, etc - are all considered "events"; these are always stored in the same event structure so that there is a consistent reference syntax for extracting information about a given event (time period): Nx1 structure where N = number of events, with 'type', 'start', 'stop', 'idxStart', 'idxStop', and 'file' as field names. Use `initEvents` to initialize an empty one.

## Parameters
`defaultParams.m` holds default parameters for a variety of functions   
`params = defaultParams` to load defaults into `params` variable. 
 
`processArgs.m` to load current parameters and change/overwrite individual parameters  
```
new_params = procesArgs(defaultParams, '[param_x]', [new_x_value], '[param_y]', [new_y_value]); %to specify new values for parameters x & y
this_fx_params = processArgs(params, varargin{:}); %to process multiple params fed into a function 
```

`defineAliases.m` can be used to load a set of parameters, rather than loading defaults and then individually specifying/changing values  
