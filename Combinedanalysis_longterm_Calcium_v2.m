%% 
%LZY script_20250626_update for long-term data
clear
clc


%read data
filename = 'F:\PhD Thesis\Long-term photometry\Cre line_Ca2+ GCamP\Processed_file\20241210_1224_Processed\All_rawdata_3minper60min.csv';
%'C:\Users\Wanglab\Desktop\Temp\20240906—1518-72h.
percentile_file = 'F:\PhD Thesis\Long-term photometry\Cre line_Ca2+ GCamP\Processed_file\Percentile\1-All_20241210_fromZT18-Percentile for fitting.csv';

data = csvread(filename,1,0);
percentile = csvread(percentile_file,1,0);

%%
cv_FOF_all_2fuse = [];
mean_FOF_all_2fuse = [];
cv_FOF_all_2d = [];
mean_FOF_all_2d = [];
all_update_signal = [];
all_update_reference = [];


for a = 1
%
j = 2*a;
m = 3*a;
all_percentile_value = percentile(:,j);
%set title as ZT time
start_value = percentile(1,j-1);
distance_to_ZT18 = 18 - start_value;
%Sampling frequency
Fs = 10;
%Order of polynomial fit
order  = 1;
%Set analyzed columns
signal_chnl = m-2;
reference_chnl = m-1;
time_chnl = m;

signal = data(:, signal_chnl);
signal = signal(signal ~= 0);
reference = data(:, reference_chnl);
reference = reference(reference ~= 0);
time = data(:,time_chnl);
time = time(time ~= 0);


%%
%数据拆分
[signal_blocks, reference_blocks, time_blocks] = splitDataByTime(signal, reference, time)

% Initialize cell arrays to store results
fitted_reference_blocks = {};
fit_coeffs_blocks = {};  
deltaFOF_blocks = {};
curr_signal_blocks = {};
curr_time_blocks = {};
modified_reference_blocks = {};
F470_O_F405_blocks = {};
percentile_blocks = [];
mean_470_blocks = [];
mean_405_blocks = [];
baseline_405_blocks = [];
baseline_470_blocks = [];
deltaFOF = [];
std_deltaF_F = [];
std_deltaFOF_2d = [];
std_deltaFOF_1d = [];
std_deltaFOF_1d_2fuse = [];
cv_FOF_2d = [];
mean_FOF_2d = [];
cv_FOF_2dfuse = [];
mean_FOF_2dfuse = [];
updated_signal = [];
update_reference = [];
detrend_baseline = [];
detrend_baseline_2fuse = [];

for i = (1 + distance_to_ZT18):(48 + distance_to_ZT18) %length(signal_blocks)
% Get the current block and remove the first 100 rows
  end_idx = min(1800, length(signal_blocks{i}));
  curr_signal = signal_blocks{i}(101:end_idx);
  curr_reference = reference_blocks{i}(101:end_idx);
  
  %Signal and reference vector
  updated_signal = [updated_signal; curr_signal];
  update_reference = [update_reference; curr_reference];
  mean_signal = mean(updated_signal);
  mean_reference = mean(update_reference);
  difference_mean = mean_signal - mean_reference;

end

% Process each block-take the 13-49 blocks for further analysis
for i = (1 + distance_to_ZT18):(48 + distance_to_ZT18) %length(signal_blocks)
% Get the current block and remove the first 100 rows
  end_idx = min(1800, length(signal_blocks{i}));
  curr_signal = signal_blocks{i}(101:end_idx);
  curr_reference = reference_blocks{i}(101:end_idx) + difference_mean;
  curr_time = time_blocks{i}(101:end_idx);  


curr_signal_blocks{end + 1} = curr_signal; %without the first 100rows
modified_reference_blocks{end + 1} = curr_reference;

%计算F470/F405和mean CV
F470_O_F405 = calculateFOF(curr_signal,  curr_reference);
mean_FOF_2d = [mean_FOF_2d; mean(F470_O_F405)];
cv_FOF_2d = [cv_FOF_2d; std(F470_O_F405)];%./ mean(F470_O_F405)
end

   


for c = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23]
    cv_FOF_2dfuse = [cv_FOF_2dfuse;(cv_FOF_2d(c) + cv_FOF_2d(c+1) + cv_FOF_2d(c+24) + cv_FOF_2d(c+25)) / 4 ];
    mean_FOF_2dfuse = [mean_FOF_2dfuse;(mean_FOF_2d(c) + mean_FOF_2d(c+1) + mean_FOF_2d(c+24) + mean_FOF_2d(c+25)) / 4 ];

end
%Store the detread_baseline and std_deltaFof each group

cv_FOF_all_2fuse(:,a) = cv_FOF_2dfuse;
mean_FOF_all_2fuse(:,a) = mean_FOF_2dfuse;
cv_FOF_all_2d(:,a) = cv_FOF_2d;
mean_FOF_all_2d(:,a) = mean_FOF_2d;
%plot(normalized_detrend_baseline,'g', 'LineWidth', 0.5);

end
%%  
%%Creat a matrix for curr_signal
curr_signal_matrix = [];
modified_reference_matrix = [];
curr_FOF_matrix = [];
gap = NaN(500,1);
for i= 1:48
curr_signal_matrix = [curr_signal_matrix; curr_signal_blocks{i};gap];
modified_reference_matrix = [modified_reference_matrix; modified_reference_blocks{i};gap];
%curr_FOF_matrix = [curr_FOF_matrix; F470_O_F405_blocks{i};gap];
end
%plot(curr_FOF_matrix,'g', 'LineWidth', 0.5);

%% 

function filtered_signal = BandPassFilter(signal)
    % 设计一个 Butterworth 带通滤波器
    f_low = 0.001; % 设置截止频率下限 0.015为采样率一半的0.015倍
    f_high = 0.2; % 设置截止频率上限
    [b, a] = butter(4, [f_low, f_high], 'bandpass');

    % Apply the filter to the signal
    filtered_signal = filtfilt(b, a, signal);
end

function smoothed_signal = smoothSignal(signal, span, method)
    smoothed_signal = smooth(signal, span, method);
end

function normalized_data= normalizeData(data)
    mean_val = mean(data);
    std_val = std(data);
    normalized_data = (data - mean_val) / std_val;
end

function zero_hundred_normalized_data= zero_hundred_normalize(data)
    
    zero_hundred_normalized_data = 100*(data - min(data)) / (max(data)-min(data));
end

function data = inverseNormalizeData(normalized_data, mean_val, std_val)
    data = normalized_data * std_val + mean_val;
end

function [fit_coeffs, fitted_reference] = fitReferenceToSignal(curr_reference, fitting_signal,fitting_reference, order)
    % Fit a linear model ( signal = a * reference + b)
    fit_coeffs = polyfit(fitting_reference, fitting_signal, order);
    % Generate the fitted signal
    fitted_reference = polyval(fit_coeffs, curr_reference);
end


function [signal_blocks, reference_blocks, time_blocks] = splitDataByTime(signal, reference, time)
    % Identify the indices where time resets
    time_reset_indices = find(diff(time) < 0) + 1;
    % Add start and end indices
    time_reset_indices = [1; time_reset_indices; length(time) + 1];

    % Initialize cell arrays to store the separated blocks
    signal_blocks = {};
    reference_blocks = {};
    time_blocks = {};

    % Split the data into blocks
    for i = 1:length(time_reset_indices) - 1
        start_idx = time_reset_indices(i);
        end_idx = time_reset_indices(i + 1) - 1;

        signal_blocks{end + 1} = signal(start_idx:end_idx);
        reference_blocks{end + 1} = reference(start_idx:end_idx);
        time_blocks{end + 1} = time(start_idx:end_idx);
    end
end

function [intensity_delta, z_score] = calculate(curr_signal, fitted_reference)
    intensity_delta = 100*(curr_signal - fitted_reference) ./ fitted_reference;
    average_delta = mean(intensity_delta);
    std_delta = std(intensity_delta);
    z_score = (intensity_delta - average_delta) ./ std_delta;
end
function [intensity_delta, z_score] = calculateCV(curr_signal, fitted_reference)
    intensity_delta = 100 * curr_signal ./ fitted_reference;
    average_delta = mean(intensity_delta);
    std_delta = std(intensity_delta);
    z_score = (intensity_delta - average_delta) ./ std_delta;
end

function deltaF = calculateFOF(curr_signal, fitted_reference)
    deltaF = curr_signal ./ fitted_reference;
end

function deltaF = calculateDF(curr_signal, fitted_reference)
    deltaF = curr_signal - fitted_reference;
end
% Function to detect events
function [count, event_indices] = detectEvents(z_score, median_val, std_val, Fs)
    % Set the threshold for detecting events
    threshold = max(median_val + 2 * std_val, 1.5 / 100); % 1.5% in decimal

    % Initialize the event count and indices
    count = 0;
    event_indices = [];

    % Define event duration in seconds and required number of samples for an event
    event_duration = 1; % Event duration in seconds
    min_samples = event_duration * Fs; % Minimum number of samples for an event

    % Identify points where the z_score exceeds the threshold
    above_threshold = z_score > threshold;

    % Find the starting and ending indices of events
    event_starts = find(diff([0; above_threshold; 0]) == 1);
    event_ends = find(diff([0; above_threshold; 0]) == -1) - 1;
    
    % Loop through each detected event and count if it exceeds the minimum duration
    for j = 1:length(event_starts)
        if (event_ends(j) - event_starts(j) + 1) >= min_samples
            count = count + 1;
            event_indices = [event_indices, event_starts(j):event_ends(j)];
        end
    end
end

% Function to integrate \Delta F/F values
function integrated_value = integrateDFOF(deltaF_F, Fs)
    positive_values = deltaF_F(deltaF_F > 0);
    integrated_value = trapz(positive_values) / Fs; % Integration over time
end