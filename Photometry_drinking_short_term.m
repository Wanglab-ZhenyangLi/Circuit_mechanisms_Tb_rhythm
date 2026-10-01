```matlab
%% 
clear all; clc;

% ---- Read data ----
filename = ''; 
data = csvread(filename,1,0);

%% ---- Channel selection ----
channel = 3;   % Please modify according to actual setup
if channel == 0
    signal_chnl = 1;   reference_chnl = 3;
elseif channel == 1
    signal_chnl = 4;   reference_chnl = 6;
elseif channel == 2
    signal_chnl = 7;   reference_chnl = 9;
elseif channel == 3
    signal_chnl = 10;  reference_chnl = 12;
end
signal = data(:, signal_chnl);
reference = data(:, reference_chnl);
time = data(:, end);

% ---- Trim first two hours ----
offset_sec = 0; duration_sec = 7200;
idx_trim = (time >= offset_sec) & (time <= offset_sec + duration_sec);
time = time(idx_trim);
signal = signal(idx_trim);
reference = reference(idx_trim);

%Check the whole picture of the signal and reference data 
figure;
subplot(2,1,1);
    plot(signal,'color','g', 'LineWidth', 1.5);
subplot(2,1,2);
    plot(reference,'b', 'LineWidth', 1.5);
    grid on;
    
% ---- Reference channel outlier repair (same as original code) ----
threshold = 0.2 * median(reference);
outlier_idx = (reference < threshold) | (reference == 0);
if any(outlier_idx)
    fprintf('Detected %d abnormal reference values, performing linear interpolation repair.\n', sum(outlier_idx));
    % Temporarily replace with NaN
    ref_corrected = reference;
    ref_corrected(outlier_idx) = NaN;
    % Fill by linear interpolation
    ref_corrected = fillmissing(ref_corrected, 'linear');
    % If NaN appears at beginning/end (extreme case), fill with nearest valid value
    ref_corrected = fillmissing(ref_corrected, 'nearest');
else
    ref_corrected = reference;
end

% Assign corrected reference signal to variable used later
reference = ref_corrected;

figure;
plot(reference,'b', 'LineWidth', 1.5);
grid on;

%% ---- Read drinking events ----
drink_file = '';
if channel == 0
    col_drink = 4;
elseif channel == 1
    col_drink = 3;
elseif channel == 2
    col_drink = 1;
elseif channel == 3
    col_drink = 2;
end
has_header_drink = false;

% ---- Read data ----
try
    if has_header_drink
        T_drink = readtable(drink_file);
        drink_raw = table2array(T_drink(:, col_drink));
    else
        drink_raw = readmatrix(drink_file);
        if size(drink_raw,2) >= col_drink
            drink_raw = drink_raw(:, col_drink);
        else
            error('Insufficient number of columns');
        end
    end
    drink_raw = drink_raw(:);
catch ME
    error('Failed to read drinking file: %s', ME.message);
end

% Time axis (500 Hz)
fs_drink = 500;
drink_time = (0:length(drink_raw)-1)' / fs_drink;

% Trim to same length as calcium signal
max_ca_time = max(time);
idx_drink = drink_time <= max_ca_time;
drink_time = drink_time(idx_drink);
drink_raw = drink_raw(idx_drink);
fprintf('Number of drinking data points = %d, time range %.2f s\n', length(drink_raw), drink_time(end));

% ---- Detect lick events (0→1 transition) ----
diff_drink = [0; diff(drink_raw)];
lick_onset_idx = find(diff_drink == 1);
lick_times = drink_time(lick_onset_idx);
fprintf('Number of raw lick events: %d\n', length(lick_times));

if isempty(lick_times)
    fprintf('No lick events detected, skipping drinking analysis\n');
    return;
end

%% ---- Analyze lick interval distribution (to help select merge threshold) ----
if length(lick_times) > 1
    intervals = diff(lick_times);
    fprintf('\n--- Lick interval statistics ---\n');
    fprintf('Interval mean: %.3f s\n', mean(intervals));
    fprintf('Interval median: %.3f s\n', median(intervals));
    fprintf('Interval std: %.3f s\n', std(intervals));
    fprintf('Interval quartiles (25%%-50%%-75%%): %.3f - %.3f - %.3f s\n', ...
        prctile(intervals,25), prctile(intervals,50), prctile(intervals,75));
    fprintf('Minimum interval: %.3f s\n', min(intervals));
    fprintf('Maximum interval: %.3f s\n', max(intervals));
    
    % Plot interval distribution
    figure('Name', 'Lick interval distribution');
    subplot(2,1,1);
    histogram(intervals, 'BinWidth', 0.02, 'FaceColor', [0.6 0.6 0.8], 'EdgeColor', 'none');
    hold on;
    [f, xi] = ksdensity(intervals);
    plot(xi, f * length(intervals) * 0.02, 'k-', 'LineWidth', 2);
    xlabel('Interval (s)');
    ylabel('Frequency');
    title('Lick interval distribution');
    xline(median(intervals), 'r--', sprintf('Median=%.2fs', median(intervals)), 'LineWidth', 1.5);
    xline(mean(intervals), 'g--', sprintf('Mean=%.2fs', mean(intervals)), 'LineWidth', 1.5);
    legend('Frequency', 'Kernel density', 'Median', 'Mean');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    subplot(2,1,2);
    boxplot(intervals, 'Orientation', 'horizontal', 'Symbol', '+');
    xlabel('Interval (s)');
    title('Interval box plot');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    
    fprintf('\n--- Threshold suggestions ---\n');
    fprintf('Suggested threshold 1 (0.5 s, commonly used): 0.50 s\n');
    fprintf('Suggested threshold 2 (2x median): %.2f s\n', 2*median(intervals));
    fprintf('Suggested threshold 3 (75%% quantile): %.2f s\n', prctile(intervals,75));
    fprintf('Please select an appropriate threshold based on the histogram valley.\n');
end

%% ---- Merge licks into bouts (threshold 15 s, can be adjusted based on above analysis) ----
bout_gap = 15;   % seconds
if length(lick_times) > 1
    inter_lick = diff(lick_times);
    bout_end_idx = find(inter_lick > bout_gap);
    bout_start_idx = [1; bout_end_idx + 1];
    bout_end_idx = [bout_end_idx; length(lick_times)];
    
    n_bouts = length(bout_start_idx);
    bout_onset = zeros(n_bouts, 1);
    bout_lick_count = zeros(n_bouts, 1);
    bout_duration = zeros(n_bouts, 1);
    
    for i = 1:n_bouts
        bout_licks = lick_times(bout_start_idx(i):bout_end_idx(i));
        bout_onset(i) = bout_licks(1);
        bout_lick_count(i) = length(bout_licks);
        bout_duration(i) = bout_licks(end) - bout_licks(1);
    end
    fprintf('Merged into %d lick bouts (interval threshold %.2f s)\n', n_bouts, bout_gap);
else
    n_bouts = 1;
    bout_onset = lick_times;
    bout_lick_count = 1;
    bout_duration = 0;
end

%% ---- Output bout length statistics ----
fprintf('\n--- Drinking bout length statistics ---\n');
fprintf('Total bouts: %d\n', n_bouts);
if n_bouts > 1
    fprintf('Duration (s) mean: %.3f, median: %.3f, std: %.3f\n', ...
        mean(bout_duration), median(bout_duration), std(bout_duration));
    fprintf('Duration quartiles (25%%-50%%-75%%): %.3f - %.3f - %.3f\n', ...
        prctile(bout_duration,25), prctile(bout_duration,50), prctile(bout_duration,75));
    fprintf('Duration min: %.3f, max: %.3f\n', min(bout_duration), max(bout_duration));
    fprintf('Lick count mean: %.2f, median: %.2f, std: %.2f\n', ...
        mean(bout_lick_count), median(bout_lick_count), std(bout_lick_count));
    fprintf('Lick count range: %d - %d\n', min(bout_lick_count), max(bout_lick_count));
    
    % Plot distributions
    figure('Name', 'Drinking bout duration and lick count');
    subplot(2,1,1);
    histogram(bout_duration, 'BinWidth', 0.1, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
    xlabel('Duration (s)');
    ylabel('Frequency');
    title('Drinking bout duration distribution');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    
    subplot(2,1,2);
    histogram(bout_lick_count, 'BinWidth', 1, 'FaceColor', [0.8 0.4 0.6], 'EdgeColor', 'none');
    xlabel('Lick count');
    ylabel('Frequency');
    title('Drinking bout lick count distribution');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
else
    fprintf('Only one bout, cannot compute distribution statistics.\n');
end

% ---- Analyze intervals between bouts ----
if n_bouts > 1
    bout_end = bout_onset + bout_duration;
    bout_intervals = bout_onset(2:end) - bout_end(1:end-1);
    fprintf('\n--- Drinking bout interval statistics ---\n');
    fprintf('Number of intervals: %d\n', length(bout_intervals));
    fprintf('Interval (s) mean: %.3f, median: %.3f, std: %.3f\n', ...
        mean(bout_intervals), median(bout_intervals), std(bout_intervals));
    fprintf('Interval quartiles (25%%-50%%-75%%): %.3f - %.3f - %.3f\n', ...
        prctile(bout_intervals,25), prctile(bout_intervals,50), prctile(bout_intervals,75));
    fprintf('Interval min: %.3f, max: %.3f\n', min(bout_intervals), max(bout_intervals));
    
    % Plot interval distribution histogram
    figure('Name', 'Drinking bout interval distribution');
    histogram(bout_intervals, 'BinWidth', 1, 'FaceColor', [0.6 0.6 0.8], 'EdgeColor', 'none');
    xlabel('Interval (s)');
    %xlim([0,30]);
    ylabel('Frequency');
    title('Drinking bout interval distribution');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
else
    fprintf('Fewer than 2 bouts, cannot compute intervals.\n');
end
%% ---- Event time-locking parameters (uniform time axis) ----
win_before = 30;   % seconds (first 10 s used for fitting)
win_after  = 60;   % seconds (last 10 s)
t_uniform = -win_before : 0.1 : win_after;   % 0.1 s resolution

% ---- Locate events on calcium signal time axis ----
event_indices = zeros(n_bouts, 1);
for i = 1:n_bouts
    [~, idx] = min(abs(time - bout_onset(i)));
    event_indices(i) = idx;
end

% ---- Filter events with complete windows and far from boundaries ----
valid = true(n_bouts, 1);
for i = 1:n_bouts
    t0 = time(event_indices(i));
    if (t0 - win_before < time(1)) || (t0 + win_after > time(end))
        valid(i) = false;
    end
end
event_indices = event_indices(valid);
bout_onset = bout_onset(valid);
n_events = length(event_indices);
fprintf('Number of valid bouts (complete window and far from boundary): %d\n', n_events);

if n_events == 0
    error('No valid bouts, please check data or adjust window parameters.');
end

%% ---- Perform independent reference fitting for each bout ----
event_dFF = zeros(n_events, length(t_uniform));

for i = 1:n_events
    t0 = time(event_indices(i));
    
    % 1. Fitting window: 10 s before event (excluding t0)
    fit_mask = (time >= t0 - win_before) & (time < t0);
    fit_idx = find(fit_mask);
    if length(fit_idx) < 5
        warning('Bout %d has insufficient fitting points, skipping', i);
        continue;
    end
    sig_fit = signal(fit_idx);
    ref_fit = reference(fit_idx);
    
    % 2. Linear fit (signal = a*ref + b)
    p = polyfit(ref_fit, sig_fit, 1);
    
    % 3. Entire event window (pre 10 to post 10 s)
    evt_mask = (time >= t0 - win_before) & (time <= t0 + win_after);
    evt_idx = find(evt_mask);
    t_rel = time(evt_idx) - t0;
    sig_evt = signal(evt_idx);
    ref_evt = reference(evt_idx);
    
    % 4. Use fit parameters to predict fitted_reference
    fitted_ref = polyval(p, ref_evt);
    
    % 5. Compute ΔF/F (%) = (signal - fitted_ref) / fitted_ref * 100
    dFF = (sig_evt - fitted_ref) ./ fitted_ref * 100;
    
    % 6. Interpolate to uniform time axis
    dFF_uniform = interp1(t_rel, dFF, t_uniform, 'linear', 'extrap');
    if any(isnan(dFF_uniform))
        dFF_uniform = fillmissing(dFF_uniform, 'nearest');
    end
    
    event_dFF(i, :) = dFF_uniform;
end

% Remove failed events (all NaN)
valid_events = ~isnan(event_dFF(:,1));
event_dFF = event_dFF(valid_events, :);
bout_onset = bout_onset(valid_events);
n_events = size(event_dFF, 1);
fprintf('Number of successfully processed bouts: %d\n', n_events);

if n_events == 0
    error('No successfully processed events, please check data.');
end

%% ---- Check first 10 events (4 per row) and export data to CSV ----
time_axis = t_uniform;
n_check = max(3, n_events);
if n_check > 0
    n_rows = ceil(n_check / 4);
    figure('Name', 'Check drinking bout events with behavior (interpolated)');
    
    % Preallocate data storage for export
    all_time = [];
    all_dFF = [];
    all_behavior = [];
    all_eventID = [];
    
    for i = 1:n_check
        subplot(n_rows, 4, i);
        
        yyaxis left;
        plot(time_axis, event_dFF(i,:), 'b-', 'LineWidth', 1);
        ylabel('ΔF/F (%)');
        %ylim([-10, 20]);
        idx_zero = find(time_axis == 0);
        if ~isempty(idx_zero)
            hold on;
            plot(0, event_dFF(i, idx_zero), 'ro', 'MarkerSize', 8);
            hold off;
        end
        
        yyaxis right;
        t0 = bout_onset(i);
        win_start = t0 - win_before;
        win_end = t0 + win_after;
        idx_beh = (drink_time >= win_start) & (drink_time <= win_end);
        t_beh = drink_time(idx_beh);
        beh_vals = drink_raw(idx_beh);
        if ~isempty(t_beh)
            time_abs = t0 + time_axis;
            beh_interp = interp1(t_beh, beh_vals, time_abs, 'nearest', 0);
        else
            beh_interp = zeros(size(time_axis));
        end
        stairs(time_axis, beh_interp, 'r-', 'LineWidth', 1.5);
        ylabel('Lick (0/1)');
        ylim([-0.1, 2.1]);
        
        xlabel('Time from onset (s)');
        title(sprintf('Bout %d at %.2fs', i, t0));
        xlim([-win_before, win_after]);
        grid off;
        ax = gca;
        ax.XAxis.TickDirection = 'out';
        ax.YAxis(1).TickDirection = 'out';
        ax.YAxis(2).TickDirection = 'out';
        ax.YAxis(1).Color = 'b';
        ax.YAxis(2).Color = 'r';
        
        % ---- Collect data for export ----
        all_time = [all_time; time_axis(:)];
        all_dFF = [all_dFF; event_dFF(i, :)'];
        all_behavior = [all_behavior; beh_interp(:)];
        all_eventID = [all_eventID; repmat(i, length(time_axis), 1)];
    end
    
    % ---- Export to CSV ----
    export_table = table(all_eventID, all_time, all_dFF, all_behavior, ...
        'VariableNames', {'EventID', 'TimeFromOnset_sec', 'dFF_percent', 'Behavior_01'});
    output_dir = "F:\PhD Thesis\Fiber photometry\Cre line_Ca2+ GCamP\Processed_file-short-term\Water_intake_simul_photometry\Sample_traces";
    export_filename = fullfile(output_dir, 'J944_licking_events_export.csv');
    writetable(export_table, export_filename);
    fprintf('Exported data for %d events to CSV file: %s\n', n_check, export_filename);
end

%% ---- Average response + SEM ----
mean_resp = mean(event_dFF, 1);
sem_resp = std(event_dFF, 0, 1) / sqrt(n_events);

figure('Name', 'Lick bout onset-triggered average (per-event fitting)');
upper = mean_resp + sem_resp;
lower = mean_resp - sem_resp;
fill([t_uniform, fliplr(t_uniform)], [upper, fliplr(lower)], 'c', ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;
plot(t_uniform, mean_resp, 'c-', 'LineWidth', 2);
xlabel('Time from bout onset (s)');
ylabel('ΔF/F (%)');
title(sprintf('Lick bout average (N=%d bouts, per-event ref fitting)', n_events));
xlim([-win_before, win_after]);
%ylim([-10, 20]);
grid off;
ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
hold off;

%% ---- Heatmap (sorted by peak) ----
if n_events > 1
    figure('Name', 'Lick bout heatmap (per-event fitting)');
    peak_vals = max(event_dFF, [], 2);
    [~, sort_idx] = sort(peak_vals, 'descend');
    sorted_data = event_dFF(sort_idx, :);
    imagesc(t_uniform, 1:n_events, sorted_data);
    axis xy;
    colormap('jet');
    colorbar;
    xlabel('Time from onset (s)');
    ylabel('Bout number (sorted by peak)');
    title('Lick bout calcium responses (heatmap)');
    %caxis([-5, 10]);
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    hold on;
    line([0, 0], [0.5, n_events+0.5], 'Color', 'w', 'LineWidth', 1.5, 'LineStyle', '--');
    hold off;
end

%% ---- Dual Y-axis figure: full-time-course calcium signal + drinking events ----
% First compute global ΔF/F₀ (if not already computed)
if ~exist('deltaF_F_global', 'var')
    fprintf('Computing global ΔF/F₀ for overview plot...\n');
    percentile_fit = 50;
    thresh_fit = prctile(signal, percentile_fit);
    idx_low_fit = signal <= thresh_fit;
    if sum(idx_low_fit) > 10
        p_global = polyfit(reference(idx_low_fit), signal(idx_low_fit), 1);
    else
        p_global = polyfit(reference, signal, 1);
    end
    fitted_ref_global = polyval(p_global, reference);
    deltaF_F_global = (signal - fitted_ref_global) ./ fitted_ref_global * 100;
end

% Interpolate drinking events to calcium signal time axis (nearest-neighbor)
drink_interp = interp1(drink_time, drink_raw, time, 'nearest', 0);

figure('Name', 'Calcium signal and drinking events over time');
yyaxis left;
plot(time / 60, deltaF_F_global, 'g-', 'LineWidth', 0.8);
ylabel('ΔF/F (%)');
%ylim([-10, 20]);
yyaxis right;
stem(time / 60, drink_interp, 'b-', 'LineWidth', 0.5, 'Marker', 'none');
ylabel('Drinking (0/1)');
ylim([-0.1, 1.1]);
xlabel('Time (minutes)');
title('Calcium signal and drinking events over 2 hours');
xlim([0, 120]);
xticks(0:30:120);
grid off;
ax = gca;
ax.XAxis.TickDirection = 'out';
ax.YAxis(1).TickDirection = 'out';
ax.YAxis(2).TickDirection = 'out';
```