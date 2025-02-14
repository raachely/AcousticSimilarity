% JMA/JS script for writing objective cluster quality to birdSummaries.xlsx
% Evaluates cluster quality by day/age, for each bird in data folder
rep = reportOnData;
nBirds = numel(rep);

%% open Excel server
Excel = actxserver('Excel.Application');

% open the specific file (note: the file must exist)
xlFilename = [pwd filesep 'birdSummaries.xlsx'];
if ~exist(xlFilename,'file')
    error('File %s does not exist...', xlFilename);
end
%Workbook = invoke(Excel.Workbooks,'Open', xlName);
Workbook = Excel.Workbooks.Open(xlFilename);
%%
rowCount = 3; %starting row on excel sheet
for ii = 1:nBirds
    sessions = {rep{ii}.sessionID};
    birdID = strtok(sessions{1}, '_');
    ages = getAgeOfSession(sessions); %requires 'session records' sheet that lists age of each session
    uAges = unique(ages(~isnan(ages)));
      
    fprintf('Getting cluster quality for %s...\n',birdID);
    clusterQualities = examineClusterQuality(birdID, uAges);
    
    for jj = 1:numel(uAges)
        % write into excel file
        fprintf('\tWriting for age %d...\n', uAges(jj));
        thisAgeQuality = clusterQualities{jj};
        xlswrite1(xlFilename, {birdID, uAges(jj)}, 'Objective Cluster Quality', ...
            sprintf('A%d:B%d',rowCount,rowCount)); %write birdID & age to columns A & B
        if ~isempty(thisAgeQuality)
            beginCol = 'D'; endCol = char(beginCol + numel(thisAgeQuality) - 1);
            xlswrite1(xlFilename, thisAgeQuality, 'Objective Cluster Quality', ...
                sprintf('%s%d:%s%d', beginCol, rowCount, endCol, rowCount)); %write cluster quality values 
        end
        rowCount = rowCount + 1;
    end
end

%% wrap up
fprintf('Please review and save AFTER unpausing to accept changes to [%s] (paused)...\n', xlFilename);
Excel.visible = true;
pause;
%% close up shop

Workbook.Close;
Excel.Quit;