
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Authors: Joshua Welsh & Verity Ford
% Date: 2024-07-24
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rules to keep the code working
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Do not remove the AnimalID or Time Point columns from the Timepoint
% metadata spreadsheet
%
% 2) Column names will be graph names.
%
% 3) The spreadsheet has to be closed and saved in order for MATLAB to 
% process it
%
% 4) Do not edit sheet names or Database name (name of whole excel document) in excel
% close all hidden; clear; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Spreadsheet = 'Database.xlsx';
ProcessDatabase(Spreadsheet)


%%
function [] = ProcessDatabase(Spreadsheet)

% remove folder if it exists already to avoid confusion between
% different analyses
if isfolder('outputs')
    rmdir('outputs', 's');
end

mkdir('outputs')

ParseDatabase(Spreadsheet)
NormalizeDatabase()

normTypes = {'Raw','BS'};
plotTypes = {'Average','Animal'};

for plotType = plotTypes
    for normType = normTypes
        PlotData(plotType{1},normType{1})
    end
end

end
%%
function [] = ParseDatabase(Spreadsheet)
% Script for processing the database spreadsheet

% remove databases before processing
databases = {'Normalized database.mat','Processed database.mat'};
for i = 1:numel(databases)
    if ~isfile(databases{i})
        delete(databases{i})
    end
end

% import data
AnimalMeta = readtable(Spreadsheet,Sheet="Animal metadata");
TimepointMeta = readtable(Spreadsheet,Sheet="Timepoint metadata");
GroupMeta = readtable(Spreadsheet,Sheet="Group metadata");

% list all unique AnimalIDs in the Animal and Timepoint metadata sheets
AnimalID = unique(AnimalMeta.AnimalID,'stable');
TimepointAnimalIDs = unique(TimepointMeta.AnimalID,'stable');

% check all listed Group IDs within Animal metadata exists
uniqueGroups = unique(AnimalMeta.GroupID);
indexedGroups = uniqueGroups(~ismember(uniqueGroups, GroupMeta.GroupID));
if sum(ismember(uniqueGroups, GroupMeta.GroupID)) ~= numel(uniqueGroups)
    for i = 1:numel(indexedGroups)
        disp(['GroupID: ',num2str(indexedGroups(i)),' is not listed in the Animal metadata spreadsheet but appears in Timepoint metadata sheet'])
    end
    error('A GroupID listed in the Animal metadata does not exist in the Group metadata sheet')
end

% find AnimalIDs listed in Timepoint metadata that are not in AnimalID
% metadata.
indexedAnimalIDs = TimepointAnimalIDs(~ismember(TimepointAnimalIDs, AnimalID));
if ~isempty(indexedAnimalIDs)
    for i = 1:numel(indexedAnimalIDs)
        disp([indexedAnimalIDs{i},': AnimalID is not listed in the Animal metadata spreadsheet but appears in Timepoint metadata sheet'])
    end
    error('All listed Timepoint metadata AnimalIDs must be correctly indexed in the Animal metadata sheet in order to be correctly indexed')
end

% creating for loop - going through every unique animal ID in the Timepoint
% metadata
for i = 1:numel(AnimalID)

    % if using this animal then proceed, if not ignore
    if AnimalMeta.Include(i) == true && ismember(AnimalMeta.AnimalID(i), TimepointMeta.AnimalID)

        % obtain group ID name to use a variable name for structuring
        % database
        groupID = GroupMeta.GroupName{GroupMeta.GroupID==AnimalMeta.GroupID(i)};

        % check survival grouping
        if AnimalMeta.Survived(i) == true
            Survival = "SURVIVOR";
        elseif AnimalMeta.Survived(i) == false
            Survival = "NONSURVIVOR";
        else
            error(['Animal ID: ',AnimalID{i}, ' is not correctly assigned survival status in Animal metadata spreadsheet'])
        end

        % Group animals by their study ID
        StudyID = ['Study',num2str(AnimalMeta.StudyID(i))];

        %finding animals in database that fit unique ID
        tempind = strcmp(TimepointMeta.AnimalID, AnimalID{i});
        tempdatabase = TimepointMeta(tempind, :);

        % Nesting database under all groups
        database.(StudyID).(groupID).(Survival).(AnimalID{i}).rawData = tempdatabase;

        % baseline normalize database
        [database] = NormaliseDatabase(database, tempdatabase, StudyID, groupID, Survival, AnimalID, i);

    end
end

% save database
save('Processed database.mat',"AnimalMeta","TimepointMeta","database",'-mat')
end


function [database] = NormaliseDatabase(database, tempdatabase, StudyID, groupID, Survival, AnimalID, i)

% define variables to normalize
varIDs = tempdatabase.Properties.VariableNames;
removeIdx = ismember(varIDs,{'AnimalID','TimePoint'});
varIDs(removeIdx) = [];

%define timepoints
timepointID = tempdatabase.TimePoint;

%identify baseline
baselineind = timepointID == 0;

%check that baseline exists
if max(baselineind) ==1

    %go through each variable ID to normalize it
    for j = 1:numel(varIDs)
 
        % create a temporary variable of the column data
        tempVarData = tempdatabase.(varIDs{j});

        % check if variable is formatted as text or not
        if iscell(tempVarData)
            % convert text data to numeric formatting
            try
                checkEmpty = cellfun(@isempty, tempVarData, 'UniformOutput', true);
                tempVarData(checkEmpty) = {'NaN'};
                tempVarData = cell2mat(cellfun(@str2num, tempVarData, 'UniformOutput', false));
            catch
                error(['error occurred parsing variable: ', varIDs{j}])
            end
        end

        if isempty(tempVarData)
        else
            if numel(timepointID) == numel(tempVarData)
                database.(StudyID).(groupID).(Survival).(AnimalID{i}).(varIDs{j}).BS = tempVarData - tempVarData(baselineind);
                database.(StudyID).(groupID).(Survival).(AnimalID{i}).(varIDs{j}).Raw = tempVarData;
            else
                error(['Number of timepoints does not match data inputs for AnimalID: ', AnimalID{i}, ', variable:', varIDs{j}])
            end
        end
    end

else
    disp(['No_baseline_found: ',AnimalID{i}])
end
end

%%
function [] = NormalizeDatabase()

% import processed data
proc_data = load("Processed database.mat");

% Get the data types of each variable
variableTypes = varfun(@class, proc_data.TimepointMeta, 'OutputFormat', 'cell');

% obtain all unique variable names in the Timepoint metadata spreadsheet
defaultVariables = proc_data.TimepointMeta.Properties.VariableNames;

% remove the indexing columns: AnimalID, TimePoint
removeIdx = ismember(defaultVariables,{'AnimalID','TimePoint'});
variableTypes(removeIdx) = [];
defaultVariables(removeIdx) = [];

% create the default processed data output for every listed timepoint in
% the Timepoint metadata spreadsheet
defaultTimepoints(1,:) = sort(unique(proc_data.TimepointMeta.TimePoint),'ascend');

% create an empty cell array for normalized data output
defaultOutput = cell(numel(defaultVariables)+1,numel(defaultTimepoints)+1);
defaultOutput(1,:) = [{'Timepoint'},num2cell(defaultTimepoints)];

% define output types
normTypes = {'BS','Raw'};
statTypes = {'mean','SEM'};
normalized = [];

% obtain the unique study IDs
studyID = fieldnames(proc_data.database);

% iterate through each study ID
for i = 1:numel(studyID)

    exportPath = fullfile('outputs',studyID{i});

    % make dir
    if ~isfolder(exportPath)
        mkdir(exportPath)
    end

    % define output filename based on StudyID
    filename = fullfile(exportPath,'Summary statistics.xlsx');

    % obtain the groups within the study
    groupID = fieldnames(proc_data.database.(studyID{i}));

    % iterate through each group ID
    for j = 1:numel(groupID)

        % obtain the survival states within the group ID
        conditionID = fieldnames(proc_data.database.(studyID{i}).(groupID{j}));

        % iterate through each survival state
        for k = 1:numel(conditionID)

            % assign the path to the current database index to a variable
            % in order to make it easier for indexing downstream
            tempDatabase = proc_data.database.(studyID{i}).(groupID{j}).(conditionID{k});

            % obtain the animal IDs within each survival state
            tempanimals = fieldnames(tempDatabase);

            % reset the default output per condition
            [tempOutput] = resetOutput(defaultOutput);

            % iterate through every variable in the Timepoint metadata
            % spreadsheet to aggregate summary statistics
            for l = 1:numel(defaultVariables)

                % select the current variable being aggregated
                tempVariable = defaultVariables{l};

                % iterate through each animal to isolate the statistics
                % at each timepoint
                for n = 1:numel(tempanimals)

                    % create default values for all statistics so
                    % non-existent timepoints will be NaN
                    tempvarstat.BS(n,:) = nan(1,numel(defaultTimepoints));
                    tempvarstat.Raw(n,:) = nan(1,numel(defaultTimepoints));

                    % iterate through each timepoint in order to calculate
                    % statistics
                    for m = 1:numel(defaultTimepoints)

                        % determine if the selected timepoint is recorded
                        % for selected animal
                        temp_timepoint_ind = ismember(tempDatabase.(tempanimals{n}).rawData.TimePoint,defaultTimepoints(m));

                        % if timepoint is associated with animals and the
                        % variable is associated with the animal database
                        % proceed
                        if max(temp_timepoint_ind) == true & isfield(tempDatabase.(tempanimals{n}), tempVariable)
                            % obtain the stat for the timepoint
                            tempvarstat.BS(n,m) = tempDatabase.(tempanimals{n}).(tempVariable).BS(temp_timepoint_ind);
                            tempvarstat.Raw(n,m) = tempDatabase.(tempanimals{n}).(tempVariable).Raw(temp_timepoint_ind);
                        end
                    end
                   
                end

                % calculate the statistics
                [normalized] = calculateStats(tempvarstat, tempVariable, defaultTimepoints, studyID{i}, groupID{j}, conditionID{k}, normalized);

                % append cell array with processed stats for output
                [tempOutput] = writeStat2OutputArray(l+1, tempOutput, tempVariable, studyID{i}, groupID{j}, conditionID{k}, normalized);

            end

            % export the concatenated cell array of processed statistics to spreadsheet
            exportSpreadsheet(tempOutput, filename, groupID{j}, conditionID{k}, statTypes, normTypes)

        end
    end
end

save('Normalized database.mat', "normalized","normTypes",'-mat')

end

%%
function [] = exportSpreadsheet(tempOutput, filename, groupID, condition, statTypes, normTypes)

% iterate through combinations of stat types and groupings to define which
% sheet the processed statistics are outputted to
for o1 = 1:numel(statTypes)
    for o2 = 1:numel(normTypes)
        % define sheet name
        sheet = [groupID,', ', condition,', ',statTypes{o1},' (',normTypes{o2},')'];
        % write data to spreadsheet
        writecell(tempOutput.(normTypes{o2}).(statTypes{o1}), filename,'Sheet',sheet)
    end
end

end

%%
function [normalized] = calculateStats(tempvarstat, tempVariable, defaultTimepoints, studyID, groupID, condition, normalized)

% SEM function
SEMCalc = @(data) std(data,"omitmissing")./sqrt(sum(~isnan(data)));

% write stats to structure
normalized.(studyID).(groupID).(tempVariable).(condition).Timepoint = defaultTimepoints;
normalized.(studyID).(groupID).(tempVariable).(condition).mean = mean(tempvarstat.BS, "omitmissing");
normalized.(studyID).(groupID).(tempVariable).(condition).SEM = SEMCalc(tempvarstat.BS);
normalized.(studyID).(groupID).(tempVariable).(condition).meanRaw = mean(tempvarstat.Raw, "omitmissing");
normalized.(studyID).(groupID).(tempVariable).(condition).SEMRaw = SEMCalc(tempvarstat.Raw);
normalized.(studyID).(groupID).(tempVariable).(condition).animalBS = tempvarstat.BS;
normalized.(studyID).(groupID).(tempVariable).(condition).animalRaw = tempvarstat.Raw;

end

%%
function [tempOutput] = resetOutput(defaultOutput)

tempOutput.BS.mean = defaultOutput;
tempOutput.BS.SEM = defaultOutput;
tempOutput.Raw.mean = defaultOutput;
tempOutput.Raw.SEM = defaultOutput;

end

%%
function [tempOutput] = writeStat2OutputArray(i, tempOutput, tempVariable, studyID, groupID, condition, normalized)

% create temporary array to reduce length in downstream access
tempStats = normalized.(studyID).(groupID).(tempVariable).(condition);

% append cell array with processed stats for output
tempOutput.BS.mean(i,:) = [tempVariable, num2cell(tempStats.mean)];
tempOutput.BS.SEM(i,:)  = [tempVariable, num2cell(tempStats.SEM)];
tempOutput.Raw.mean(i,:) = [tempVariable, num2cell(tempStats.meanRaw)];
tempOutput.Raw.SEM(i,:) = [tempVariable, num2cell(tempStats.SEMRaw)];

end

%%
function [] = PlotData(plotType, normType)

close all hidden

% import normalized data
load("Normalized database.mat")

[fig,t] = createFigure();
ExportCount = 0;

% obtain the unique study IDs
studyID = fieldnames(normalized);


% iterate through each study ID
for i = 1:numel(studyID)

    % create output folder
    exportFolder = fullfile('outputs', studyID{i}, normType, plotType);

    % create output folder
    mkdir(exportFolder)

    % obtain the groups within the study
    groupID = fieldnames(normalized.(studyID{i}));

    % iterate through each group ID
    for j = 1:numel(groupID)

        varID = sort(fieldnames(normalized.(studyID{i}).(groupID{j})));

        for k = 1:numel(varID)

            % obtain the survival states within the group ID
            conditionID = fieldnames(normalized.(studyID{i}).(groupID{j}).(varID{k}));

            % iterate through each survival state
            for l = 1:numel(conditionID)

                % create an axis handle for first condition of selected
                % varID and reset plot handle
                if l == 1
                    [ax] = createAxis();
                    p = [];
                end

                % write the selection to a temporary variable for easier access
                % to plotting scripts
                tempData = normalized.(studyID{i}).(groupID{j}).(varID{k}).(conditionID{l});

                % obtain the summary data statistics
                [tempStats] = returnData(tempData, plotType, normType);
                tempStats.VarName = varID{k};
                tempStats.Condition = conditionID{l};

                % create plot
                [ax, p] = plotData(ax, p, l, tempStats, plotType);

            end

            % format the plot axis
            [ax] = formatAxis(ax, tempStats);

            if rem(k,8) == 0
                % create global legend
                createLegend(conditionID)
                % export figure
                [fig, ExportCount] = ExportFigure(ExportCount, fig, exportFolder);
            end
        end
    end
    if rem(k,8) ~= 0
        % create global legend
        createLegend(conditionID)
        % export figure
        [fig, ExportCount] = ExportFigure(ExportCount, fig, exportFolder);
    end
end



end

%%
function [tempStats] = returnData(tempData, plotType, normType)

tempStats.xData = tempData.Timepoint;

switch plotType
    case 'Average'
        switch normType
            case 'BS'
                tempStats.yData = tempData.mean;
                tempStats.yErr = tempData.SEM;
            case 'Raw'
                tempStats.yData = tempData.meanRaw;
                tempStats.yErr = tempData.SEMRaw;
        end
        % find missing values
        tempStats.nanInd = isnan(tempStats.yData);
    case 'Animal'
        switch normType
            case 'BS'
                tempStats.yData = tempData.animalBS;
            case 'Raw'
                tempStats.yData = tempData.animalRaw;
        end
end

end

%%
function [ax] = createAxis()

ax = nexttile();
set(ax, 'NextPlot', 'add')

end

%%
function [ax] = formatAxis(ax, stats)
maxX = ceil(max(stats.xData)/12)*12;

set(ax,...
    'XLim', [0 maxX], ...
    'XTick', 0:12:maxX, ...
    'XTickLabelRotation', 90, ...
    'FontSize', 14, ...
    'LineWidth', 2,...
    'TickDir','out');
xlabel(ax, 'Timepoint (hours)');
ylabel(ax, replace(stats.VarName, '_', ' '));
end

%%
function [ax, p] = plotData(ax, p, count, stats, plotType)

[col] = getPlotColor(count);

switch plotType
    case 'Average'
        % remove nan values so the lines connect each point
        xData = stats.xData(~stats.nanInd);
        yData = stats.yData(~stats.nanInd);

        yErrHigh = stats.yData((~stats.nanInd))+stats.yErr(~stats.nanInd);
        yErrLow = stats.yData((~stats.nanInd))-stats.yErr(~stats.nanInd);

        if ~isempty(xData)
            % plot stats
            p(count) = plot(ax, xData, yData,'o-', ...
                'linewidth',2,...
                'MarkerEdgeColor',col(count,:), ...
                'MarkerFaceColor',col(count,:), ...
                'Color',col(count,:), ...
                'displayname',replace(stats.Condition,'_','-'));

            % plot error bars
            for m = 1:numel(xData)
                line(ax, repmat(xData(m),1,2), ...
                    [yErrLow(m),yErrHigh(m)],...
                    'linewidth',2,'color',col(count,:))
            end
        end

    case 'Animal'
        for i = 1:size(stats.yData,1)
            xData = stats.xData;
            yData = stats.yData(i,:);

            stats.nanInd = isnan(yData);
            xData = xData(~stats.nanInd);
            yData = yData(~stats.nanInd);

            if ~isempty(xData)
                % plot stats
                p = plot(ax, xData, yData, ...
                    '-','linewidth',2,...
                    'Color',col(count,:));
            end
        end
end

end

function [col] = getPlotColor(count)

% define colors for survivor / non survivor
if count <= 2
    col = [0.5 0 0; 0 0 0];
else
    %if more than two groups in future create more colors
    col = lines(count);
end

end
%%
function [fig, ExportCount] = ExportFigure(ExportCount, fig, folder)

ExportCount = ExportCount +1;

print(fig,fullfile(folder, ['Output ', num2str(ExportCount)]), '-dpng','-r300')

close(fig)

[fig,t] = createFigure();

end

%%
function [fig,t] = createFigure()
fig = figure('Units','centimeters','Position',[0 0 21 29.7],'visible','off','Renderer','painters');
t = tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
end

%%
function [] = createLegend(conditionID)

[col] = getPlotColor(numel(conditionID));

for i = 1:size(col,1)
    p(i) = plot(NaN, NaN,'o-', ...
        'linewidth',2,...
        'MarkerEdgeColor',col(i,:), ...
        'MarkerFaceColor',col(i,:), ...
        'Color',col(i,:), ...
        'displayname',replace(conditionID{i},'_','-'));
end

% create global legend
lg = legend(p, 'location','northoutside','NumColumns',numel(conditionID),'box','off');
lg.Layout.Tile = 'north';
end