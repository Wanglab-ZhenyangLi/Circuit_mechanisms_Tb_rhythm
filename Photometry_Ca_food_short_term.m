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

%% ---- Read eating events ----
eating_filename = '';
col_eat = 3; has_header_eat = true;
if has_header_eat
    T_eat = readtable(eating_filename);
    eat_raw = table2array(T_eat(:, col_eat));
else
    eat_raw = readmatrix(eating_filename);
    eat_raw = eat_raw(:, col_eat);
end
eat_raw = eat_raw(:);
fs_eat = 1;
eat_time = (0:length(eat_raw)-1)' / fs_eat;
% Trim to same length as calcium signal
max_ca_time = max(time);
idx_eat = eat_time <= max_ca_time;
eat_time = eat_time(idx_eat);
eat_raw = eat_raw(idx_eat);

% ---- Event detection (0→1 transition) ----
diff_eat = [0; diff(eat_raw)];
onset_idx = find(diff_eat == 1);
onset_time = eat_time(onset_idx);
fprintf('Detected %d eating onset events\n', length(onset_time));

%% ---- Merge consecutive events into bouts ----
if length(onset_time) > 1
    bout_threshold = 60;   % seconds, intervals <= this value are merged into the same bout
    bout_start_times = [];
    bout_end_times = [];
    bout_event_counts = [];
    bout_idx = 1;
    current_bout_start = onset_time(1);
    current_bout_end = onset_time(1);
    current_bout_count = 1;
    
    for i = 2:length(onset_time)
        if onset_time(i) - onset_time(i-1) <= bout_threshold
            % Belongs to same bout, update end time and count
            current_bout_end = onset_time(i);
            current_bout_count = current_bout_count + 1;
        else
            % End current bout, record
            bout_start_times(bout_idx) = current_bout_start;
            bout_end_times(bout_idx) = current_bout_end;
            bout_event_counts(bout_idx) = current_bout_count;
            bout_idx = bout_idx + 1;
            % Start new bout
            current_bout_start = onset_time(i);
            current_bout_end = onset_time(i);
            current_bout_count = 1;
        end
    end
    % Record last bout
    bout_start_times(bout_idx) = current_bout_start;
    bout_end_times(bout_idx) = current_bout_end;
    bout_event_counts(bout_idx) = current_bout_count;
    
    % Compute duration of each bout (seconds)
    bout_durations = bout_end_times - bout_start_times;
    
    % Update onset_time to bout start time
    onset_time = bout_start_times';
    
    fprintf('Number of bouts after merging: %d (threshold=%.0f s)\n', length(onset_time), bout_threshold);
else
    % Only one event, treat as a single bout
    bout_start_times = onset_time;
    bout_end_times = onset_time;
    bout_event_counts = ones(size(onset_time));
    bout_durations = zeros(size(onset_time));
    fprintf('Fewer than 2 events, each event treated as a separate bout.\n');
end

%% ---- Output bout length statistics ----
fprintf('\n--- Eating bout length statistics ---\n');
fprintf('Total bouts: %d\n', length(onset_time));
if length(bout_durations) > 1
    fprintf('Duration (s) mean: %.2f, median: %.2f, std: %.2f\n', ...
        mean(bout_durations), median(bout_durations), std(bout_durations));
    fprintf('Duration quartiles (25%%-50%%-75%%): %.2f - %.2f - %.2f\n', ...
        prctile(bout_durations,25), prctile(bout_durations,50), prctile(bout_durations,75));
    fprintf('Duration min: %.2f, max: %.2f\n', min(bout_durations), max(bout_durations));
    fprintf('Event count mean: %.2f, median: %.2f, std: %.2f\n', ...
        mean(bout_event_counts), median(bout_event_counts), std(bout_event_counts));
    fprintf('Event count range: %d - %d\n', min(bout_event_counts), max(bout_event_counts));
else
    fprintf('Only one bout, cannot compute distribution statistics.\n');
end

% ---- Plot bout length distributions (histograms) ----
if length(bout_durations) > 1
    figure('Name', 'Eating bout duration distribution');
    subplot(2,1,1);
    histogram(bout_durations, 'BinWidth', 10, 'FaceColor', [0.8 0.4 0.4], 'EdgeColor', 'none');
    xlabel('Duration (s)');
    ylabel('Frequency');
    title('Eating bout duration distribution');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    
    subplot(2,1,2);
    histogram(bout_event_counts, 'BinWidth', 1, 'FaceColor', [0.4 0.6 0.8], 'EdgeColor', 'none');
    xlabel('Event count (pellets)');
    ylabel('Frequency');
    title('Eating bout event count distribution');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
end

% ---- Analyze intervals between bouts ----
if length(bout_start_times) > 1
    bout_intervals = bout_start_times(2:end) - bout_end_times(1:end-1);
    fprintf('\n--- Eating bout interval statistics ---\n');
    fprintf('Number of intervals: %d\n', length(bout_intervals));
    fprintf('Interval (s) mean: %.2f, median: %.2f, std: %.2f\n', ...
        mean(bout_intervals), median(bout_intervals), std(bout_intervals));
    fprintf('Interval quartiles (25%%-50%%-75%%): %.2f - %.2f - %.2f\n', ...
        prctile(bout_intervals,25), prctile(bout_intervals,50), prctile(bout_intervals,75));
    fprintf('Interval min: %.2f, max: %.2f\n', min(bout_intervals), max(bout_intervals));
    
    % Plot interval distribution histogram
    figure('Name', 'Eating bout interval distribution');
    histogram(bout_intervals, 'BinWidth', 30, 'FaceColor', [0.6 0.6 0.8], 'EdgeColor', 'none');
    xlabel('Interval (s)');
    ylabel('Frequency');
    title('Eating bout interval distribution');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
else
    fprintf('Fewer than 2 bouts, cannot compute intervals.\n');
end


% ---- Recompute event indices (because onset_time updated to bout starts) ----
event_indices = zeros(length(onset_time),1);
for i = 1:length(onset_time)
    [~, idx] = min(abs(time - onset_time(i)));
    event_indices(i) = idx;
end

%% ---- Event time-locking parameters ----
win_before = 30;   % seconds
win_after  = 60;
t_uniform = -win_before : 0.1 : win_after;

% ---- Filter events with complete windows and far from boundaries ----
valid = true(size(event_indices));
for i = 1:length(event_indices)
    t0 = time(event_indices(i));
    if (t0 - win_before < time(1)) || (t0 + win_after > time(end))
        valid(i) = false;
    end
end
event_indices = event_indices(valid);
onset_time = onset_time(valid);
n_events = length(event_indices);
fprintf('Number of valid bouts (complete window and far from boundary): %d\n', n_events);

% ---- Perform reference fitting for each event ----
event_dFF = zeros(n_events, length(t_uniform));

for i = 1:n_events
    t0 = time(event_indices(i));
    
    % 1. Fitting window: 10 s before event (excluding t0)
    fit_mask = (time >= t0 - win_before) & (time < t0);
    fit_idx = find(fit_mask);
    if length(fit_idx) < 5
        warning('Event %d has insufficient fitting points, skipping', i);
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

% Remove failed events (rows all NaN)
valid_events = ~isnan(event_dFF(:,1));
event_dFF = event_dFF(valid_events, :);
onset_time = onset_time(valid_events);
n_events = size(event_dFF, 1);
fprintf('Number of successfully processed events: %d\n', n_events);

if n_events == 0
    error('No valid events, please check data.');
end

%% ---- Check first few events (plot) and export data to CSV ----
time_axis = t_uniform;
n_check = max(12, n_events);
if n_check > 0
    n_rows = ceil(n_check / 4);
    figure('Name', 'Check eating events with behavior (interpolated, reference fitted)');
    
    % Preallocate data storage for export
    all_time = [];
    all_dFF = [];
    all_behavior = [];
    all_eventID = [];
    
    for i = 1:n_check
        % ---- Left axis: Ca signal ΔF/F ----
        subplot(n_rows, 4, i);
        yyaxis left;
        plot(time_axis, event_dFF(i,:), 'b-', 'LineWidth', 1);
        ylabel('ΔF/F (%)');
        ylim([-10, 20]);
        idx_zero = find(time_axis == 0);
        if ~isempty(idx_zero)
            hold on;
            plot(0, event_dFF(i, idx_zero), 'ro', 'MarkerSize', 8);
            hold off;
        end
        
        % ---- Right axis: Eating marker (0/1) ----
        yyaxis right;
        t0 = onset_time(i);
        win_start = t0 - win_before;
        win_end = t0 + win_after;
        idx_beh = (eat_time >= win_start) & (eat_time <= win_end);
        t_beh = eat_time(idx_beh);
        beh_vals = eat_raw(idx_beh);
        if ~isempty(t_beh)
            time_abs = t0 + time_axis;
            beh_interp = interp1(t_beh, beh_vals, time_abs, 'nearest', 0);
        else
            beh_interp = zeros(size(time_axis));
        end
        stairs(time_axis, beh_interp, 'r-', 'LineWidth', 1.5);
        ylabel('Eating (0/1)');
        ylim([-0.1, 2.1]);
        
        % ---- Common settings ----
        xlabel('Time from onset (s)');
        title(sprintf('Event %d at t=%.2fs', i, t0));
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
    output_dir = "F:\PhD Thesis\Fiber photometry\Cre line_Ca2+ GCamP\Processed_file-short-term\Food_intake_simul_photometry\Sample traces";
    export_filename = fullfile(output_dir, 'J944_eating_events_export.csv');
    writetable(export_table, export_filename);
    fprintf('Exported data for %d events to CSV file: %s\n', n_check, export_filename);
end

%% ---- Average response and standard error ----
mean_resp = mean(event_dFF, 1);
sem_resp  = std(event_dFF, 0, 1) / sqrt(n_events);

figure('Name', 'Eating onset-triggered average (per-event ref fitting)');
upper = mean_resp + sem_resp;
lower = mean_resp - sem_resp;
fill([time_axis, fliplr(time_axis)], [upper, fliplr(lower)], 'm', ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;
plot(time_axis, mean_resp, 'm-', 'LineWidth', 2);
xlabel('Time from eating onset (s)');
ylabel('ΔF/F (%)');
title(sprintf('Eating onset average (N=%d events, per-event reference fitting)', n_events));
xlim([-win_before, win_after]);
%ylim([-10, 20]);
grid off;
ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
hold off;

% ---- Heatmap (sorted by peak) ----
if n_events > 1
    figure('Name', 'Eating heatmap (interpolated, reference fitted)');
    peak_vals = max(event_dFF, [], 2);
    [~, sort_idx] = sort(peak_vals, 'descend');
    sorted_data = event_dFF(sort_idx, :);
    imagesc(time_axis, 1:n_events, sorted_data);
    axis xy;
    colormap('jet');
    colorbar;
    xlabel('Time from onset (s)');
    ylabel('Event number (sorted by peak)');
    title('Eating calcium responses (heatmap, per-event reference fitting)');
    %caxis([-5, 10]);   % Note: if your MATLAB version does not support clim, use caxis
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    hold on;
    line([0, 0], [0.5, n_events+0.5], 'Color', 'w', 'LineWidth', 1.5, 'LineStyle', '--');
    hold off;
end

%% ---- Compute global ΔF/F₀ (for overview plot) ----
percentile_fit = 50;   % Use points below the 50th percentile as baseline
thresh_fit = prctile(signal, percentile_fit);
idx_low_fit = signal <= thresh_fit;
if sum(idx_low_fit) < 10
    warning('Insufficient baseline points, fitting with all data');
    p_global = polyfit(reference, signal, 1);
else
    p_global = polyfit(reference(idx_low_fit), signal(idx_low_fit), 1);
end
fitted_ref_global = polyval(p_global, reference);
deltaF_F = (signal - fitted_ref_global) ./ fitted_ref_global * 100;
fprintf('Global ΔF/F₀ computed.\n');

% ---- Plot full-time-course dual Y-axis figure ----
% This figure shows calcium signal and eating events over the entire two hours
figure('Name', 'Calcium signal and eating events over time');

% Left axis: Calcium signal ΔF/F
yyaxis left;
plot(time / 60, deltaF_F, 'g-', 'LineWidth', 0.8);
ylabel('ΔF/F (%)');
ylim([-10, 20]);  % Unified range

% Right axis: Eating marker (0/1)
yyaxis right;
% Interpolate eating data to calcium signal time points (nearest-neighbor interpolation, preserving 0/1)
eat_interp = interp1(eat_time, eat_raw, time, 'nearest', 0);
% Display eating events with vertical lines (draw a vertical line at each eating point)
stem(time / 60, eat_interp, 'b-', 'LineWidth', 0.5, 'Marker', 'none');
ylabel('Eating (0/1)');
ylim([-0.1, 2.1]);   % Keep 0/1 range

% Common settings
xlabel('Time (minutes)');
title('Calcium signal and eating events over 2 hours');
xlim([0, 120]);
xticks(0:30:120);
grid off;
ax = gca;
ax.XAxis.TickDirection = 'out';
ax.YAxis(1).TickDirection = 'out';
ax.YAxis(2).TickDirection = 'out';