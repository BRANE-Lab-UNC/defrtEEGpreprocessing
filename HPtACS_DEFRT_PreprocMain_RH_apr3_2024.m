% Preprocessing pipeline for DEEFRT task for DEFRT EEG Study
% The following can be run as a function or as a script (i.e., for line-by-line
% running of a dataset through the pipeline). This script converts the raw EEG
% data from MFF format to Matlab usable data. The photodiode is parsed to
% identify relevant task events. Then the data is filtered (1 Hz highpass, 58 Hz
% lowpass) then downsampled (200 Hz). This data is then fed into artifact
% subspace reconstruction where bad channels are removed and bursty artifacts
% are corrected. Removed channels are then interpolated before average
% re-referencing. Next, ICA is computed. High probability artifacts are
% automatically labeled before a final manual review of all components. Any
% additional artifact components can be selected at this stage.
%
%   NOTE: If running this script as a loop, the loop will wait during the ICA
%   artifact selection step. Once the review is complete, hit enter in the
%   Command Window to progress to the next dataset.
%
% 6/22/23, CPW -- Updated DEFRT script to bring into alignment with GNG script
% 6/27/23, CPW -- Corrected ASR step to add 'WindowCriterion','off'.

%clear all;  % Clears variables, globals, functions, etc.
clear;      % Clears workspace
clc;        % Clears the command window.
%close all;  % Closes all open figures.
dbstop in HPtACS_DEFRT_PreprocMain_RH_apr3_2024.m at 447 %breakpoint on executable code

%% Defining Variables
runManualICA_Rejection = 0; % If 0, then skip. If 1, then do it.

%% Set up environment
% Get the full path of the currently executing script
currentScriptFullPath = mfilename('fullpath');

% Extract the directory of the current script
[currentScriptPath, ~, ~] = fileparts(currentScriptFullPath);

% Get the username in a platform-independent way
if ispc  % For Windows
    username = getenv('USERNAME');
else  % For macOS/Linux
    username = getenv('USER');
end

% Define the file name where the directory path will be stored
settingsFilePath = fullfile(currentScriptPath, 'MATLab Root', username, 'currentScriptDir.txt');

% Before attempting to open the file, ensure the directory exists
[settingsFileDir, ~, ~] = fileparts(settingsFilePath);

if ~exist(settingsFileDir, 'dir')
    mkdir(settingsFileDir); % Creates the directory and any necessary parent directories
end

% Check if the settings file exists
if exist(settingsFilePath, 'file') ~= 2
    % User input MATLab folder containing data
    % If the file does not exist, prompt the user for the MATLAB folder path
    MATLabDir = input('Enter the full path to your MATLAB folder: ', 's');
    
    % Open the file for writing
    fileId = fopen(settingsFilePath, 'wt');
    
    % Write the directory path to the file
    fprintf(fileId, '%s', MATLabDir);
    
    % Close the file
    fclose(fileId);
else
    % If the file exists, read the stored directory path
    fileId = fopen(settingsFilePath, 'rt');
    MATLabDir = fgets(fileId);
    fclose(fileId);
    
    % Remove newline character at the end if it exists
    MATLabDir = strtrim(MATLabDir);
end

% Assuming your desired root directory is up a few levels or in a specific
% relative path from the current script directory, adjust the path accordingly.
% For example, if the root directory is two levels up:
MATLabROOT = fullfile(MATLabDir, '/');
ROOT = fullfile(currentScriptPath, '/');

rawEEG_dir = [ROOT 'DEFRT_Raw_EEG/'];
if exist(rawEEG_dir,'dir')~=7
    error('EEG analysis folder does not exist');
end
PROC_EEG_DIR = [ROOT 'DEFRT_Proc_EEG/'];
if exist(PROC_EEG_DIR,'dir')~=7
    mkdir(PROC_EEG_DIR);
end

matches = dir(fullfile(MATLabROOT, 'eeglab*'));
dirMatches = matches([matches.isdir]);
if isempty(dirMatches)
    error('EEGLAB folder does not exist')
else
    % Select the first matching directory
    eeglabFolder = fullfile(MATLabROOT, dirMatches(1).name);
    
    % Add the found directory to the MATLAB path
    addpath(eeglabFolder);
end

defrtCodeFolder = [ROOT 'Code/'];
if exist(defrtCodeFolder,'dir')~=7
    error('DEFRT code folder does not exist');
end
addpath(defrtCodeFolder);

% Initialize EEGlab
eeglab; close all;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% DEFRT PREPROCESSING PIPELINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% DEEFRT TASK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify some default paramters
pdChan = 130;   % photodiode channel number
refChan = 129;  % reference channel number
nChan = 128;    % number of original EEG channels
newSrate = 200; % sampling rate for downsampled data

% Create subject list, 49 removed bc already processed. 50 removed
% because issue w/ data
SUBJECTS = {...
    '63'};
numSub = length(SUBJECTS);

% numMissingSub 1 by 2 matrix
% First column is subject ID
% Second column is a list of the missing triggers IN ORDER
PROBLEM_DATA = {
    '63', {'MISSING_Block2_Trial2_Trigger3', 'MISSING_Block3_Trial16_Trigger3', 'MISSING_Block6_Trial6_Trigger3', 'EXTRA_Block6_Trial21_Trigger1', 'EXTRA_Block6_Trial21_Trigger2'}
};
problem_SUB = PROBLEM_DATA(:,1);

%% For Loop for Data
%% MAIN LOOP --> (1) File Loading (2) Event Handling (3) Filtering (4) ASR (5) ICA
for iSub = 1:numSub

    %% Load and prepare EEG (1)
    subject = SUBJECTS{iSub};

    % Create a folder in the processed EEG directory for this participant
    SUB_PROC_EEG = [PROC_EEG_DIR 'DEFRT_' subject '_DEEFRT' '/'];
    if exist(SUB_PROC_EEG,'dir')~=7
        mkdir(SUB_PROC_EEG);
    end

    % Now check if their raw file exists in the raw directory
    dfrtFileFinder = dir([rawEEG_dir 'DEFRT_' subject '_EEG/' 'DEFRT_' subject '_DEEfRT_*.mff']);
    numPrimaryFiles = length(dfrtFileFinder);

    % some people have two files: the naming convention is to say "DEEfRT2"
    secondFile_dfrtFileFinder = dir([rawEEG_dir 'DEFRT_' subject '_DEEfRT2_*.mff']);
    numSecondFiles = length(secondFile_dfrtFileFinder);

    % New root file name for the processed data
    rootFilename = ['DEFRT_' subject '_DEEfRT'];

    % Check to see if participant is already completed
    preIcaFile = [SUB_PROC_EEG rootFilename '_preICA.mat'];
    postIcaFile = [SUB_PROC_EEG rootFilename '_postICA.mat'];
   % if exist(preIcaFile,'file') == 2
%        disp([rootFilename ' already preprocessed through ICA. Skipping for now. If this is a mistake, reprocess manually.'])
%    else
        disp(['Starting preprocessing steps on ' rootFilename '.'])

        % Check if an appropriate file was found
        if numPrimaryFiles == 0         % No file found
            error(['No DEEFRT file found for ' subject '. Check filename.'])
            
        elseif (numPrimaryFiles==1) && (numSecondFiles == 0)     % One file found (this is the expected behavior)
            
            % Raw DEfRT file
            rawDefrtFile = [rawEEG_dir 'DEFRT_' subject '_EEG/' dfrtFileFinder(1).name];
            % Load the dataset
            DFRT_EEG = mff_import(rawDefrtFile);
            
        elseif (numPrimaryFiles == 1) && (numSecondFiles == 1)
            
            % NOTE: So far no one has had multiple blocks in multiple files. All cases
            % have been that one is invalid for some reason. These cases should be added
            % as special cases above to ensure the correct input file is selected rather
            % that the code below which would load and append datasets. Uncomment the
            % lines below if that type of routine is ever needed.
            
            % Raw DEfRT file
            rawDefrtFile = [rawEEG_dir dfrtFileFinder(1).name];
            % Load the dataset
            DFRT_EEG1 = mff_import(rawDefrtFile);
            
            % Raw DEfRT file
            secondFile_rawDefrtFile = [rawEEG_dir secondFile_dfrtFileFinder(1).name];
            % Load the dataset
            DFRT_EEG2 = mff_import(secondFile_rawDefrtFile);
            DFRT_EEG = pop_mergeset(DFRT_EEG1,DFRT_EEG2);
            
        else
            error(['Uh oh! Something is weird about ' subject]);
        end
        
        % Exit the loop if the data is already preprocessed (the 'continue' command
        % will proceed to the next iteration of the loop, but will error running code
        % line-by-line (i.e., with F9), so skip this as needed.
        %  if skipSub; end %- KB commented out to eliminate error 10.25.23
        
        %% Begin preprocessing -- (2) Event handling
        % Center the data by removing the mean offset and linear trend of the data
        DFRT_EEG.data(1:128,:) = detrend(DFRT_EEG.data(1:128,:).').'; %KB 11/16/23 edit from DFRT.data(1:128,:) = 'detrend(DFRT.data(1:128,:))';
        
        % Extract the photodiode (PD) channel and remove it from the EEG data structure
        DFRT_EEG.photodiode.data = DFRT_EEG.data(130,:);      % Copy the photodiode channel into a separate field %error in this line KB 11/6/23: Index in position 1 exceeds array bounds (must not exceed 1)
        DFRT_EEG.photodiode.times = DFRT_EEG.times;           % Copy the time vector as well
        DFRT_EEG.data(pdChan,:) = []; % Remove the PD channel from the data matrix                   % Remove the PD channel from the data matrix
        DFRT_EEG.chanlocs(pdChan) = [];                   % Remove the channel entry in chanlocs
        DFRT_EEG.nbchan = refChan;                        % Reduce the EEG channels field to reflect the correct number
        
       

        %% Pulvinar Code
        % Identify the inflection points of the photodiode channel to find event
        % onsets
        [DFRT_EEG.pdevent,DFRT_EEG.photodiode.data] = parsePhotodiode_sept7_2023(DFRT_EEG.photodiode.data,DFRT_EEG.photodiode.times,true);
        DFRT_EEG.urevent = DFRT_EEG.event;      % Store existing event structure as backup
        DFRT_EEG.urpdevent = DFRT_EEG.pdevent;  % Store the initial PD event structure as backup
        
        % Plot initial view of all events -- This should show all the DIN and PD
        % events plotted over the trace of the photodiode channel.
        pdSamp = [DFRT_EEG.pdevent.latency];        % Create vector of photoddiode events

        % Create a large figure
        figure('Units','normalized','OuterPosition',[0 0 1 1]);
        
        % In the first subplot, we show all identified events in the file.
        %subplot(2,1,1)
        %plot(DFRT_EEG.times,         DFRT_EEG.photodiode.data,...
        %    DFRT_EEG.times(pdSamp), DFRT_EEG.photodiode.data(pdSamp), '*r',...
        %    DFRT_EEG.times(dinSamp),DFRT_EEG.photodiode.data(dinSamp),'*g');
%{
        plot(DFRT_EEG.times,         DFRT_EEG.photodiode.data,...
            DFRT_EEG.times(pdSamp), DFRT_EEG.photodiode.data(pdSamp), '*r');
        legend({'Photodiode Channel','Photodiode Events'});
        title('Overview of all events identified in the dataset')
%}
      
        % Merge event structures
        onsetEvts = [DFRT_EEG.pdevent.code] > 0; % Find times where the PD channel when high (i.e., event onset)
        % Time between each photodiode event
        ISIs = diff([DFRT_EEG.pdevent.latency]); % Create a vector of the time between events
        % Remove events with time difference around 1 second - because at
        % the start of each block we flash the photodiode on and off
        pdCheckEvts = logical([ISIs<1010 0]); % ISIs less than ~500 ms are from the pre-task event check routines
        % remove a single event between each block - when the screen
        % flashes on
        breakEvts = [ISIs>8500 0];           % ISIs greater than ~9 seconds are taken as inter-block periods
        % Delete photodiode events specified above ("break events" and
        % "check events")
        DFRT_EEG.pdevent = DFRT_EEG.pdevent(onsetEvts & ~(pdCheckEvts | breakEvts));  % Keep task-relevant events
        keepers = find(onsetEvts & ~(pdCheckEvts | breakEvts));
        latencies = [DFRT_EEG.pdevent.latency];
        print([SUB_PROC_EEG '/' rootFilename '_revisedeventParse.png'],'-dpng')

        %% OLD CODE FROM PREVIOUS PERSON
        % Handle special/odd cases
        if strcmpi(rootFilename,'DEFRT_22_DEEfRT') ...  % special case, extra PD event at the end for some reason
                || strcmpi(rootFilename,'DEFRT_24_DEEFRT') ...
                || strcmpi(rootFilename,'DEFRT_26_DEEFRT') ...
                || strcmpi(rootFilename,'DEFRT_28_DEEFRT') ...
                || strcmpi(rootFilename,'DEFRT_36_DEEFRT') ...
                || strcmpi(rootFilename,'DEFRT_08_DEEFRT') ...
                || length(DFRT_EEG.pdevent) == 379 % This exact number seems to indicate this occurrence.
            
            % Simple correction, remove last artifactual event
            DFRT_EEG.pdevent(end) = [];
            
        elseif strcmpi(rootFilename,'DEFRT_12_DEEFRT') % Special case
            % Recording started part-way through the first block, so those first trials
            % are removed
            DFRT_EEG.pdevent([1:12 253]) = []; % remove first partial block & an extra artifactual event b/w blocks 4 and 5
        end

        %% THE NEW CODE TO GET OUR EVENTS TIMED OUT PROPERLY
        % ADD EVENTS THAT ARE MISSING AND REMOVE EXTRA EVENTS
        % Structure to be editted
        pdevent = DFRT_EEG.pdevent;
        fieldNames = fieldnames(pdevent); % Get all field names of the pdevent struct array
        defaultValues = struct(); % Initialize a struct to hold default values for these fields
        for i = 1:length(fieldNames)
            defaultValues.(fieldNames{i}) = NaN; % Assuming NaN is an appropriate default value
        end

        % Does this participant have missing triggers?
        problem_subIdx = find(strcmpi(subject,problem_SUB));
        % If there are disrupted triggers for this person
        if ~isempty(problem_subIdx)
            problemTriggers = PROBLEM_DATA{problem_subIdx,2};
            numProblems = length(problemTriggers);
            for problemIdx = 1:numProblems
                thisProblem = problemTriggers{problemIdx};
                % problem parts
                problemParts = regexp(thisProblem,'_','split');
                % four parts to the problem
                % add or remove
                addOrRm_BOOL = problemParts{1};
                blockIdx_str = problemParts{2};
                blockIdx = str2double(blockIdx_str(6:end));
                trialIdx_str = problemParts{3};
                trialIdx = str2double(trialIdx_str(6:end));
                triggerIdx_str = problemParts{4};
                triggerIdx = str2double(triggerIdx_str(8:end));

                numTriggersInPractice = 18;
                numTriggersInTrial = 3;
                numTrials = 20;
                numTriggersInBlock = 3 * numTrials;

                problem_eventIdx = numTriggersInPractice + ...
                    (numTriggersInBlock * (blockIdx - 1)) + ...
                    ((trialIdx-1) * numTriggersInTrial) + ...
                    triggerIdx;

                if strcmpi(addOrRm_BOOL,'EXTRA')
                    % If we have one to many here, then remove it
                    pdevent(problem_eventIdx-1) = [];
                else
                    % If we have one missing, then we need to add it back in

                    % If we are missing the first trigger, then use
                    % intertrial time from the last trial to calculate the
                    % current trials first trigger
                    blockIdxStr = num2str(blockIdx);
                    base = ['D-EEfRT/DEFRT_sub0' subject '/DEFRT_sub0' subject '_' blockIdxStr];
                    files = dir([ROOT base '_*.mat']);
                    if ~isempty(files)
                        fullPath = fullfile(files(1).folder, files(1).name);
                        % Now you can load the file or do whatever is needed with fullPath
                        disp(['Accessing file: ' fullPath]);
                    else
                        disp('No files found matching the pattern.');
                    end
                    load(fullPath, 'trialInfo')
                    for i = 1:numel(trialInfo)
                        interTrialTimes(i) = trialInfo(i).iti;
                    end
                    missing_firstTrigger = interTrialTimes(trialIdx) + 3000;
                    % if we are missing the second trigger, then count
                    % forwards from the first trigger.
                    missing_secondTrigger = 4000;
                    % If we are missing the third trigger, then count
                    % forwards from the second trigger.
                    missing_thirdTrigger = 5500;

                    % Logic for handling missing triggers (as you've provided)
                    if triggerIdx == 1
                        previousEvent_time = pdevent(problem_eventIdx-1).latency;
                        new_eventTime = previousEvent_time + missing_firstTrigger;
                    elseif triggerIdx == 2
                        previousEvent_time = pdevent(problem_eventIdx-1).latency;
                        new_eventTime = previousEvent_time + missing_secondTrigger;
                    elseif triggerIdx == 3
                        previousEvent_time = pdevent(problem_eventIdx-1).latency;
                        new_eventTime = previousEvent_time + missing_thirdTrigger;
                    end
            
                    % Creating a new event structure for the new latency
                    newEvent = defaultValues; % Start with default values for all fields
                    newEvent.latency = new_eventTime; % Set the calculated latency
            
                    % Insert the new event into the pdevent struct array at the correct position
                    pdevent = [pdevent(1:problem_eventIdx-1), newEvent, pdevent(problem_eventIdx:end)];
                end
            end
        end
        % now that we made our edits put the event structure back in
        DFRT_EEG.pdevent = pdevent;

        % Loop through each entry and extract the numeric part
        for i = 1:numel(DFRT_EEG.pdevent)
            % Extract the numeric part from the description
            numericString = DFRT_EEG.pdevent(i).description(3:end);
            if ~isempty(numericString) % Make sure numericString is not empty
                % Convert the numeric string to a number
                numericValue = str2double(numericString);
                % Assign the numeric value to the keepers field
                DFRT_EEG.pdevent(i).keepers = numericValue;
            else
                % If numericString is empty, handle the case appropriately, perhaps by setting to NaN or 0
                DFRT_EEG.pdevent(i).keepers = NaN;
            end
        end

        % Now label the trials based on trial phase
        % The remaining events should only be onsets and will occur in 3s.
        trialOrder = {'decideON'; 'taskON'; 'rewardON'}; % Label for each trial phase
        if (mod(length(DFRT_EEG.pdevent)-18,60) == 0)   % After subtracting 3 events x 6 trials (i.e., practice),
            % do we have even blocks of 60 events?
            DFRT_EEG.pdevent(1:18) = [];                  % If so, remove the practice blocks
        end
        
        c = 1;                                      % Counter variable for loop
        if mod(length(DFRT_EEG.pdevent),60) == 0        % Check to ensure that we have even blocks of 60 (2 events x 20 trials)
            for trIdx = 1:length(DFRT_EEG.pdevent)        % Loop through each event
                DFRT_EEG.pdevent(trIdx).type = trialOrder{c}; % Label the event with the correct trial phase
                if mod(c,3)==0                          % Reset the counter variable c to if we've reached
                    c=1;                                  % The end fo the list.
                else
                    c=c+1;                                % Otherwise, count up.
                end
            end
        else                                        % If block structure is unclear
            ISIs2 = diff([DFRT_EEG.pdevent.latency]);     % Try to find the initial block start
            blockStarts = find(ISIs2>12000);          % and clip the practice blocks.
            firstEvent = find(diff(blockStarts)>18,1,'first');
            firstEvent = blockStarts(firstEvent)+1;
            DFRT_EEG.pdevent = DFRT_EEG.pdevent(firstEvent:end);
        end
        
        % Extract latencies from the structure array
        latencies = [DFRT_EEG.pdevent.latency];
        
        % Create a timeline from the first to the last latency
        % Assuming latencies are in some units of time (e.g., milliseconds)
        timeline = min(latencies):max(latencies);
        
        % Initialize an array to store the y-values, defaulting to 0
        y_values = zeros(size(timeline));
        
        % For each latency, find the corresponding index in the timeline and set it to 1
        for i = 1:length(latencies)
            % Find the index in the timeline where the latency occurs
            % The 'find' function returns an index, so the actual value in timeline at that index should match the latency
            idx = find(timeline == latencies(i));
            % Set the y-value to 1 at this index
            y_values(idx) = 1;
        end
        
        % Plot the timeline and the y-values
        plot(timeline, y_values, '-o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); % Blue line with blue circle markers
            xlabel('Time (milliseconds)');
            ylabel('Event Occurrence');
            title('Latency Events Over Time');
            ylim([-0.1 1.1]) 
        test = [];
            %{
        plot(DFRT_EEG.times,         DFRT_EEG.photodiode.data,...
            DFRT_EEG.times(pdSamp), DFRT_EEG.photodiode.data(pdSamp), '*r');
        legend({'Photodiode Channel','Photodiode Events'});
        title('Overview of all events identified in the dataset')
            %}
        
        % Update plot with new events
        %subplot(2,1,2); cla
        %eventSamp = round([DFRT_EEG.pdevent.latency]);
        %plot(DFRT_EEG.times,           DFRT_EEG.photodiode.data,...
            %DFRT_EEG.times(eventSamp),DFRT_EEG.photodiode.data(eventSamp),'*r');
        %legend({'Photodiode Channel','Retained Events'});
        %title('Selected events only, ensure no stray event labels remain')
        
        % Evalue the success of the event identifications. Red markers should only
        % occur on the rising edges.
        
        % Save figure for review
        %print([SUB_PROC_EEG '/' rootFilename '_eventParse.png'],'-dpng')
        
        
        
        %% Preprocess the data for ASR (3)
        % Apply high-pass filter and low pass separately, then resample
        disp('Beginning filtering and downsampling...'); tic;
        DFRT_EEG = pop_eegfiltnew(DFRT_EEG,'locutoff',1);  % Highpass 1 Hz
        DFRT_EEG = pop_eegfiltnew(DFRT_EEG,'hicutoff',58); % Lowpass 58 Hz (cuts line noise)
        DFRT_EEG = pop_resample(DFRT_EEG,newSrate);        % Resample to 200 Hz
        disp(['Filtering and resampling completed with processing time = ' num2str(toc/60,3) ' min.']);
        
        %% Run ASR and post-ASR processing (4)
        % Plot initial look at data to compare to ASR cleaned data;
        clf % clear previous figure
        channelVarZ = zscore(std(DFRT_EEG.data,[],2));  % Find high variance channels for plotting
        likelyBadChan = channelVarZ > 5;            % Identify channels where the variance over time is
        % >5 SD from the mean variance over channels
        % If any bad channels are found, plot in red
        if any(likelyBadChan)
            h0 = plot(DFRT_EEG.times,DFRT_EEG.data(likelyBadChan,:),'r','LineWidth',.5); hold on;
        end
        % Then plot all "good" channels in black. This is just to get a sense of
        % overall noise and the noise reduction post-ASR
        h1 = plot(DFRT_EEG.times,DFRT_EEG.data(~likelyBadChan,:),'k','LineWidth',.5); hold on;
        
        % Begin ASR. Default BurstCriterion = 5 has been too agressive in my previous
        % experience. Chang (2019) suggested higher criterion levels in the past up to
        % 20 or 30. I've found that values around 10-15 are usually decent at removing
        % artifacts while not removing too much real data. The MaxMem may be increased
        % if needed to prevent windowing artifacts (i.e., sharp edges around the
        % windows).
        disp('Beginning ASR preprocessing...'); tic
        DFRT_EEG.chanlocsOrig = DFRT_EEG.chanlocs;  % Storing for interpolation after ASR step
        DFRT_EEG = clean_artifacts(DFRT_EEG, ...                    % ASR FUNCTION
            'Highpass','off', ...        % Turn off filtering, we already did that
            'BurstCriterion',10, ...     % Use k=10 burst criterion
            'WindowCriterion','off', ... % no windowed data removal
            'MaxMem',2^13);              % Increases RAM available to the process
        disp(['... ASR concluded after ' num2str(toc/60,3) ' min.']);
        
        % Plot updated data (prior to rereferencing for comparability)
        % This should be reviewed for overall Gestalt. The blue traces should be +/-
        % 100 uV if everything went well.
        h2 = plot(DFRT_EEG.times,DFRT_EEG.data,'b','LineWidth',.5);
        title([rootFilename ': ASR Performance Evaluation'],'Interpreter','none');
        if any(likelyBadChan)
            legend([h1(1) h0(1) h2(1)],{'PreASR (Good Channels)'; 'PreASR (High Variance Channels)'; 'PostASR'});
        else
            legend([h1(1) h2(1)],{'PreASR'; 'PostASR'});
        end
        ylabel('Voltage (\muV)'); xlabel('Time (sec)');
        ylim([-1500 1500]);
        
        % Save output for review
        print([SUB_PROC_EEG '/' rootFilename '_ASRreview.png'],'-dpng')
        
        % Interpolate removed channels using original data channels (idx=1:128)
        DFRT_EEG.chanlocsASR = DFRT_EEG.chanlocs;
        DFRT_EEG = pop_interp(DFRT_EEG,DFRT_EEG.chanlocsOrig(1:nChan));
        DFRT_EEG.chanlocs = DFRT_EEG.chanlocsOrig;
        
        % Add a zero channel at 129 for the reference channel
        DFRT_EEG.data(refChan,:) = 0;
        DFRT_EEG.nbchan = refChan;
        
        % Re-reference the data to the average reference
        DFRT_EEG = pop_reref(DFRT_EEG,[]);
        
        %% Conduct ICA (5)
        % Compute the rank of the data to account for the reduction introduced by the
        % ASR and rereferencing steps.
        nComps = rank(cov(DFRT_EEG.data'));
        DFRT_EEG.rank = nComps;
        
        % Run the extended infomax ICA function (runica) to capture subgaussian
        % sources. EEG source distributions tend to be subgaussian, so this extended
        % option is optimal for EEG. The 'pca' flag ensures that we don't have faulty
        % components returned due to "overlearning" the dataset. This is necessary
        % because the ASR step reduces the actual rank and thus number of components
        % extracted from the data.
        disp('Beginning ICA step...'); tic
        DFRT_EEG = pop_runica(DFRT_EEG,'extended',1,'pca',nComps);
        disp(['... ICA concluded after ' num2str(toc/60,3) ' min.\n\nManually select components for removal.']);
        
        save(preIcaFile,'-struct','DFRT_EEG','-v7.3');
end
    
    % Now preprocessing is done - do we want to manually reject or not?
    if exist(preIcaFile,'file')==2 && exist(postIcaFile,'file')~=2 && runManualICA_Rejection
        
        %% Label component types with ICLabel
        % Use ICLabel to generate component classification estimates for review
        DFRT_EEG = load(preIcaFile);
        DFRT_EEG = iclabel(DFRT_EEG);     % Uses contains.m so this doesn't work with octave
        
        % Classification thresholds by type, adjust as needed. The left value is a
        % probability threshold (greater gets identified), the right value is logical
        % indicating "mark for removal" (1) or do nothing (0).
        compThresh = [0    0; ...        % brain
            0.7  1; ...        % muscle
            0.7  1; ...        % eye
            0.7  1; ...        % heart
            0.7  1; ...        % line noise
            0.7  0; ...        % channel noise
            0.7  0];           % other
        
        % Apply thresholds. Identified "bad components" are in DFRT.reject.gcompreject
        DFRT_EEG = pop_icflag(DFRT_EEG,compThresh);
        
        % Plot components for manual review and selection
        EEG = DFRT_EEG; % The eeglab plotting functions rely on the EEG structure as a
        % global variable, so store DFRT temporarily in EEG.
        EEG = pop_selectcomps(DFRT_EEG,1:nComps);  % Plot the component maps for selection
        disp('Evaluated ICs for artifact. High likelihood components already highlighted. Once done, press any key to continue.')
        
        % This step pause the loop so you can make selections. Click in the command
        % window and hit any key to continue.
        pause
        
        
        DFRT_EEG = EEG;       % Transfer the selections from EEG to DFRT
        DFRT_EEG.icaact = []; % Remove the icaact data to save diskspace
        
        %% Save outputs
        save([SUB_PROC_EEG rootFilename '_postICA'],'-struct','DFRT_EEG','-v7.3');
        close all
    end
%end
