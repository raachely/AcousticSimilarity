function songStruct = loadWavFile(filename)

[y, Fs] = audioread(filename,'native');
% populate some fields
songStruct.values = y;
songStruct.interval = 1/Fs;
songStruct.length = numel(y);
end