```matlab
clear; close all; clc;

data = readtable("");

cell_names = data.Properties.VariableNames;

calcium_data = table2array(data);

if size(calcium_data, 1) < 150
    rows_to_use = size(calcium_data, 1);
else
    rows_to_use = 150;
end

calcium_data = calcium_data(1:rows_to_use, :);
%% 
[num_frames, num_cells] = size(calcium_data);

%% 
time_axis = (0:num_frames-1)';

corrected_data = zeros(size(calcium_data));  
baseline_traces = zeros(size(calcium_data)); 
fitted_params = cell(num_cells, 1);          
correction_factors = zeros(size(calcium_data)); 

for cell_idx = 1:num_cells
    
    raw_signal = calcium_data(:, cell_idx);
   
    if any(isnan(raw_signal))
        raw_signal = fillmissing(raw_signal, 'linear');
    end
    
    fit_signal = raw_signal; 
    
 ------------------------------------------------
  
    force_init_len = max(10, round(0.1 * num_frames));
    force_init_idx = 1:force_init_len;
    
    win_size = max(30, round(0.05 * num_frames));
    local_var = movvar(fit_signal, win_size);
    
    search_start = force_init_len + 1;
    seg_len = round(0.15 * num_frames);
    min_var_mean = inf;
    best_start = search_start;
    
    if (search_start + seg_len) <= num_frames
        for start_idx = search_start:10:(num_frames - seg_len + 1)
            end_idx = start_idx + seg_len - 1;
            current_var_mean = mean(local_var(start_idx:end_idx));
            current_signal_mean = mean(fit_signal(start_idx:end_idx));
            
            score = current_var_mean * (1 + 0.05*current_signal_mean);
            
            if score < min_var_mean
                min_var_mean = score;
                best_start = start_idx;
            end
        end
        stable_idx = best_start:(best_start + seg_len - 1);
    else
        stable_idx = [];
    end
    
    baseline_idx = union(force_init_idx, stable_idx)
   
    exp_model = @(p, t) p(1) * exp(-p(2) * t) + p(3);
    
    baseline_signal = fit_signal(baseline_idx);
    baseline_time = time_axis(baseline_idx) - time_axis(baseline_idx(1));
    
    init_c = quantile(baseline_signal, 0.25);
    
    if mean(baseline_signal(1:round(end/3))) > mean(baseline_signal(round(2*end/3):end))
        init_a = mean(baseline_signal(1:round(end/3))) - init_c;
    else
        init_a = (max(baseline_signal) - min(baseline_signal)) * 0.5;
    end
    
    % Estimate decay rate based on forced initial segment
    init_segment = fit_signal(force_init_idx);
    init_time_segment = time_axis(force_init_idx) - time_axis(1);
    
    if length(init_segment) >= 5
        p = polyfit(init_time_segment, init_segment, 1);
        slope = p(1);
        
        if slope < 0 && init_segment(1) > 0
            init_b = abs(slope) / (init_segment(1) - init_c + eps);
            init_b = min(max(init_b, 0.001), 0.2);
            fprintf('      -> Estimated decay rate from initial segment slope: b=%.4f\n', init_b);
        else
            init_b = 0.01;
        end
    else
        init_b = 0.02;
    end
    
    fprintf('      -> Dynamic initial values: a=%.2f, b=%.4f, c=%.2f\n', init_a, init_b, init_c);
    
    % Parameter bounds
    lb = [0, 1e-5, 0.5*init_c];
    ub = [3*init_a, 0.5, max(fit_signal)*1.5];
    
    options = optimset('Display', 'off', 'TolFun', 1e-8, 'TolX', 1e-8, ...
                       'MaxIter', 2000, 'MaxFunEvals', 3000);
    
    % Level 1: Exponential fit on full data --------------------------------------
    fprintf('      -> Attempting exponential fit on full data\n');
    try
        [params_full, resnorm_full, ~, exitflag_full] = lsqcurvefit(exp_model, ...
            [init_a, init_b, init_c], ...
            time_axis - time_axis(1), ...
            fit_signal, lb, ub, options);
        
        % Compute full fit quality metrics
        bleach_curve_full = exp_model(params_full, time_axis - time_axis(1));
        R2_full = 1 - (resnorm_full / sum((fit_signal - mean(fit_signal)).^2));
        global_ratio_range_full = max(bleach_curve_full ./ fit_signal) - min(bleach_curve_full ./ fit_signal);
        
        % Check conditions 1 and 3: R2<0.8 or global ratio range>1.0
        if R2_full >= 0.8 && global_ratio_range_full <= 1.0
            % Full fit successful
            fprintf('      -> Full exponential fit successful (R^2=%.3f, ratio range=%.2f)\n', R2_full, global_ratio_range_full);
            bleach_curve = bleach_curve_full;
            fitted_params{cell_idx}.method = 'exponential_full';
            fitted_params{cell_idx}.params = params_full;
            fitted_params{cell_idx}.R2 = R2_full;
            fitted_params{cell_idx}.exitflag = exitflag_full;
            
            % Perform correction
            baseline_traces(:, cell_idx) = bleach_curve;
            correction_factors(:, cell_idx) = bleach_curve(1) ./ bleach_curve;
            corrected_data(:, cell_idx) = raw_signal .* correction_factors(:, cell_idx);
            
            continue; % Skip subsequent steps, process next cell
        else
            fprintf('      -> Full exponential fit poor (R^2=%.2f, ratio range=%.2f)\n', R2_full, global_ratio_range_full);
            % Continue to try Level 2
        end
    catch ME
        fprintf('      -> Full exponential fit failed: %s\n', ME.message);
        % Continue to try Level 2
    end
    
   % Level 2: Initial segment + 10% quantile point exponential fit --------------------
    fprintf('      -> Attempting initial segment + 10%% quantile point exponential fit\n');
    
    % Compute 10% quantile points (lower envelope)
    envelope_window = max(5, min(31, round(num_frames * 0.07)));
    if mod(envelope_window, 2) == 0
        envelope_window = envelope_window + 1;
    end
    
    lower_envelope = zeros(num_frames, 1);
    half_win = floor(envelope_window / 2);
    
    for t = 1:num_frames
        start_idx = max(1, t - half_win);
        end_idx = min(num_frames, t + half_win);
        window_data = fit_signal(start_idx:end_idx);
        lower_envelope(t) = quantile(window_data, 0.10); % 10% quantile
    end
    
    % Select initial segment and lower envelope points as fit data
    combined_time = [time_axis(force_init_idx); time_axis];
    combined_signal = [fit_signal(force_init_idx); lower_envelope];
    
    % Re-estimate initial parameters for Level 2
    init_segment_2 = fit_signal(force_init_idx);
    init_time_2 = time_axis(force_init_idx) - time_axis(1);
    
    % Re-estimate parameters
    init_c_2 = quantile(init_segment_2, 0.25);
    init_a_2 = max(init_segment_2) - min(init_segment_2);
    
    if length(init_segment_2) >= 5
        p2 = polyfit(init_time_2, init_segment_2, 1);
        slope2 = p2(1);
        
        if slope2 < 0 && init_segment_2(1) > 0
            init_b_2 = abs(slope2) / (init_segment_2(1) - init_c_2 + eps);
            init_b_2 = min(max(init_b_2, 0.001), 0.2);
        else
            init_b_2 = 0.01;
        end
    else
        init_b_2 = 0.02;
    end
    
    fprintf('      -> Level 2 initial values: a=%.2f, b=%.4f, c=%.2f\n', init_a_2, init_b_2, init_c_2);
    
    % === Key fix: clear any leftover flags ===
    use_linear_fit = false; % Explicitly set to false first
    
    try
        [params_combined, resnorm_combined, ~, exitflag_combined] = lsqcurvefit(exp_model, ...
            [init_a_2, init_b_2, init_c_2], ...
            combined_time - combined_time(1), ...
            combined_signal, lb, ub, options);
        
        % Compute fit quality
        bleach_curve_combined = exp_model(params_combined, time_axis - time_axis(1));
        R2_combined = 1 - (resnorm_combined / sum((combined_signal - mean(combined_signal)).^2));
        global_ratio_range_combined = max(bleach_curve_combined ./ fit_signal) - min(bleach_curve_combined ./ fit_signal);
        
        % === Debug: print actual judgment values ===
        fprintf('      [Debug] Actual judgment values: R2_combined=%.3f, ratio range=%.3f\n', R2_combined, global_ratio_range_combined);
        
        % Check conditions 1 and 3: R2<0.7 or global ratio range>1.0
        if R2_combined >= 0.7 && global_ratio_range_combined <= 1.0
            % Level 2 fit successful
            fprintf('      -> Initial segment + 10%% quantile point exponential fit successful (R^2=%.3f, ratio range=%.2f)\n', R2_combined, global_ratio_range_combined);
            bleach_curve = bleach_curve_combined;
            fitted_params{cell_idx}.method = 'exponential_combined';
            fitted_params{cell_idx}.params = params_combined;
            fitted_params{cell_idx}.R2 = R2_combined;
            fitted_params{cell_idx}.exitflag = exitflag_combined;
            
            % === Key fix: explicitly set to prevent entering Level 3 ===
            use_linear_fit = false;
            
            % Perform correction and skip subsequent levels
            baseline_traces(:, cell_idx) = bleach_curve;
            correction_factors(:, cell_idx) = bleach_curve(1) ./ bleach_curve;
            corrected_data(:, cell_idx) = raw_signal .* correction_factors(:, cell_idx);
            
            continue; % === Key: jump directly to next cell, do not execute Level 3 ===
            
        else
            % Level 2 fit poor
            fprintf('      -> Initial segment + 10%% quantile point exponential fit poor (R^2=%.2f, ratio range=%.2f)\n', R2_combined, global_ratio_range_combined);
            use_linear_fit = true; % Allow entering Level 3
        end
        
    catch ME
        % Exponential fit process failed
        fprintf('      -> Initial segment + 10%% quantile point exponential fit failed: %s\n', ME.message);
        use_linear_fit = true; % Allow entering Level 3
    end
    % === End of Level 2 ===
    
    % Level 3: 10% quantile linear fit ----------------------------------------
    if exist('use_linear_fit', 'var') && use_linear_fit
        fprintf('      -> Executing fallback: linear fit based on 10%% quantile points\n');
        % Compute 5% quantile points (lower quantile, better excludes calcium transients)
        envelope_window = max(5, min(31, round(num_frames * 0.07)));
        if mod(envelope_window, 2) == 0
            envelope_window = envelope_window + 1;
        end
        
        lower_envelope = zeros(num_frames, 1);
        half_win = floor(envelope_window / 2);
        
        % Compute rolling 5% quantile to build a more conservative "lower envelope" curve
        for t = 1:num_frames
            start_idx = max(1, t - half_win);
            end_idx = min(num_frames, t + half_win);
            window_data = fit_signal(start_idx:end_idx);
            lower_envelope(t) = quantile(window_data, 0.01); % Changed to 5% quantile
        end
        % Prepare fit data: use points on the lower envelope
        X = time_axis - time_axis(1);
        Y = lower_envelope;
        
        % Use polyfit for linear fit, unconditionally accept its result
        % p_simple(1) is the slope, p_simple(2) is the intercept
        [p_simple, S] = polyfit(X, Y, 1);
        
        % Compute fitted values
        [Y_fit_simple, delta] = polyval(p_simple, X, S);
        
        % Compute R²
        residuals = Y - Y_fit_simple;
        ss_res_simple = sum(residuals .^ 2);
        ss_tot_simple = sum((Y - mean(Y)) .^ 2);
        
        if ss_tot_simple > eps
            R2_simple = 1 - (ss_res_simple / ss_tot_simple);
        else
            R2_simple = 0; % No variation in data, set R² to 0
        end
        
        final_slope = p_simple(1);
        final_intercept = p_simple(2);
        bleach_curve = Y_fit_simple;
        
        % Report fit results
        fprintf('         Linear fit complete (slope=%.2e, intercept=%.2f, R²=%.3f)\n', ...
                final_slope, final_intercept, R2_simple);
        
        % Store results
        fitted_params{cell_idx}.method = 'linear_10percentile_polyfit';
        fitted_params{cell_idx}.slope = final_slope;
        fitted_params{cell_idx}.intercept = final_intercept;
        fitted_params{cell_idx}.R2_linear = R2_simple;
        
        % Optional: light smoothing (commented out per "no smoothing" requirement)
        % bleach_curve = smoothdata(bleach_curve, 'movmean', 3);
        
    end
    % End of Level 3
        
    % Perform division correction ----------------------------------------------------
    baseline_traces(:, cell_idx) = bleach_curve;
    correction_factors(:, cell_idx) = bleach_curve(1) ./ bleach_curve;
    corrected_data(:, cell_idx) = raw_signal .* correction_factors(:, cell_idx);
    
end
fprintf('Bleach correction complete!\n');
%%
F_over_F0 = zeros(size(calcium_data));
for cell_idx = 1:num_cells
    raw = calcium_data(:, cell_idx);
    baseline = baseline_traces(:, cell_idx);
    
    % Compute F/F₀ = raw / baseline
    F_over_F0(:, cell_idx) = raw ./ baseline;
end
% Sort by mean F/F₀ in descending order
mean_F_over_F0 = mean(F_over_F0, 1);  % Compute mean of each column (cell)
[~, sort_idx] = sort(mean_F_over_F0, 'descend');  % Get sort indices
F_over_F0_sorted = F_over_F0(:, sort_idx);  % Reorder data by sort indices

%% 5. Calcium signal correction factor heatmap - Z-score normalization (Science journal standard)
% Create new figure window, set Science journal recommended dimensions
fig_width_cm = 17.8; % Science double-column width
fig_height_cm = 10; % Appropriate height
fig_width = fig_width_cm / 2.54 * 96; % Convert to pixels
fig_height = fig_height_cm / 2.54 * 96;

fig_heatmap = figure('Position', [100, 100, fig_width, fig_height], ...
    'Color', 'white', ...
    'Units', 'inches', ...
    'PaperUnits', 'inches', ...
    'PaperSize', [fig_width_cm/2.54, fig_height_cm/2.54]);

% Prepare heatmap data - each row is a cell
heatmap_data = deltaF_over_F0';

% Z-score normalization: subtract row mean, divide by row std
heatmap_data_zscore = zeros(size(heatmap_data));
row_mean = zeros(size(heatmap_data, 1), 1);
row_std = zeros(size(heatmap_data, 1), 1);

for i = 1:size(heatmap_data, 1)
    row_data = heatmap_data(i, :);
    row_mean(i) = mean(row_data);
    row_std(i) = std(row_data);
    
    if row_std(i) > 1e-10  % Avoid division by zero
        heatmap_data_zscore(i, :) = (row_data - row_mean(i)) / row_std(i);
    else
        heatmap_data_zscore(i, :) = zeros(size(row_data)); % If std is zero, set to 0
    end
end

% Sort by original mean correction factor (descending)
mean_correction = mean(heatmap_data, 2); % Use mean of original data for sorting
[mean_correction_sorted, sort_idx] = sort(mean_correction, 'descend');

% Reorder data by sort result
heatmap_data_sorted = heatmap_data_zscore(sort_idx, :);
row_mean_sorted = row_mean(sort_idx);
row_std_sorted = row_std(sort_idx);
cell_indices_sorted = 1:num_cells; % Original cell indices
cell_indices_sorted = cell_indices_sorted(sort_idx); % Sorted indices

% Create heatmap - use color map suitable for calcium signals
% For calcium signals, we typically use warm tones (red, orange, yellow) to highlight calcium transients
% Option 1: 'hot' - from black to red to yellow to white, ideal for calcium signals
% Option 2: 'parula' - MATLAB default, also suitable for calcium signals
% Option 3: 'inferno' or 'plasma' - from black to bright yellow, very prominent for transients
% Option 4: Custom red palette

% I recommend 'inferno' or 'plasma' as they highlight high values well
if exist('inferno', 'file') || exist('inferno.m', 'file')
    cmap = inferno(256); % inferno from black to bright yellow, very prominent
elseif exist('plasma', 'file') || exist('plasma.m', 'file')
    cmap = plasma(256); % plasma from purple to yellow, also very good
else
    % If these color maps are not available, use hot or create custom
    cmap = hot(256); % hot is built into MATLAB, from black to red to yellow to white
end

% Or create custom red palette optimized for calcium signals
% cmap = create_calcium_colormap();

imagesc(time_axis, 1:size(heatmap_data_sorted, 1), heatmap_data_sorted);
colormap(cmap);
cbar = colorbar;

% Science journal font settings
set(gca, 'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'normal');
set(cbar, 'FontName', 'Arial', 'FontSize', 8);

% Set axis labels
xlabel('Time (s)', 'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'bold');
ylabel('Cell index', 'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'bold');
cbar.Label.String = 'Z-score';
cbar.Label.FontName = 'Arial';
cbar.Label.FontSize = 9;
cbar.Label.FontWeight = 'bold';

% Set color range: calcium signals typically focus on positive values (calcium transients)
% For Z-score, negative values indicate below baseline, may not be the main focus
z_min = min(heatmap_data_sorted(:));
z_max = max(heatmap_data_sorted(:));

% To highlight calcium transients, we can focus the color range on positive values
% But keep some negative range to show baseline fluctuations
if z_max > 0
    % If data has positive values, set color range to [-1, max(3, z_max)]
    color_min = max(-1, min(-0.5, z_min)); % Smaller negative range
    color_max = max(3, z_max); % Ensure high Z-score calcium transients are visible
    caxis([color_min, color_max]);
    
    % Set colorbar ticks
    if color_max <= 5
        cbar_ticks = [color_min, 0, 1, 2, 3, color_max];
        cbar_ticklabels = arrayfun(@(x) sprintf('%.1f', x), cbar_ticks, 'UniformOutput', false);
    else
        cbar_ticks = [color_min, 0, 1, 2, 3, 5, color_max];
        cbar_ticklabels = {sprintf('%.1f', color_min), '0', '1', '2', '3', '5', sprintf('%.1f', color_max)};
    end
else
    % If no positive values, use full range
    caxis([z_min, z_max]);
    cbar_ticks = linspace(z_min, z_max, 5);
    cbar_ticklabels = arrayfun(@(x) sprintf('%.1f', x), cbar_ticks, 'UniformOutput', false);
end

% Set colorbar ticks
cbar.Ticks = cbar_ticks;
cbar.TickLabels = cbar_ticklabels;

% Set X-axis display
if length(time_axis) <= 20
    xticks(time_axis);
else
    % Show 6 time ticks
    num_ticks = min(6, length(time_axis));
    tick_indices = round(linspace(1, length(time_axis), num_ticks));
    xticks(time_axis(tick_indices));
    xticklabels(arrayfun(@(x) sprintf('%.0f', x), time_axis(tick_indices), 'UniformOutput', false));
end

% Set Y-axis display
num_cells_total = size(heatmap_data_sorted, 1);
if num_cells_total <= 30
    % Show all cell labels when cell count <= 30
    yticks(1:num_cells_total);
    yticklabels(cell_indices_sorted);
elseif num_cells_total <= 60
    % Show one label every 5 cells
    ytick_interval = 5;
    yticks(1:ytick_interval:num_cells_total);
    yticklabels(cell_indices_sorted(1:ytick_interval:end));
else
    % Show one label every 10 cells
    ytick_interval = 10;
    yticks(1:ytick_interval:num_cells_total);
    yticklabels(cell_indices_sorted(1:ytick_interval:end));
end

% Add title
title('Calcium signal correction factor Z-score heatmap (sorted by mean correction factor descending)', ...
    'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');

% Add grid lines
grid on;
set(gca, 'GridColor', [0.3, 0.3, 0.3], 'GridAlpha', 0.1, 'GridLineStyle', ':');

% Beautify figure
set(gca, 'Box', 'on', 'LineWidth', 0.5);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.008, 0.008]);

% Use vector graphics renderer
set(gcf, 'Renderer', 'painters');
%% 4. Result visualization
figure('Position', [200, 100, 1400, 900]);

% Show detailed analysis of first 4 cells as examples
num_examples = min(12, num_cells);
for i = 1:num_examples
    subplot(4, 4, i);
    
    % Raw signal vs fitted baseline
    plot(time_axis, calcium_data(:, i+24), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Raw signal');
    hold on;
    plot(time_axis, baseline_traces(:, i+24), 'r--', 'LineWidth', 2, 'DisplayName', 'Fitted baseline');
    
    title(sprintf('Cell %d: %s', i, cell_names{i+24}), 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Time (s)', 'FontSize', 9);
    ylabel('Fluorescence intensity', 'FontSize', 9);
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    
    % Display fit parameters
    if length(fitted_params{i}) >= 3
        if fitted_params{i}(3) == 0  % Linear fit
            text(0.05, 0.15, sprintf('Linear correction\nSlope: %.3e', fitted_params{i}(1)), ...
                'Units', 'normalized', 'FontSize', 8, 'BackgroundColor', 'white');
        else  % Exponential fit
            text(0.05, 0.15, sprintf('Exponential correction\nDecay constant: %.4f', fitted_params{i}(2)), ...
                'Units', 'normalized', 'FontSize', 8, 'BackgroundColor', 'white');
        end
    end
    
%     % Corrected signal (ΔF/F0)
%     subplot(4, 4, i+4);
%     plot(time_axis, corrected_data(:, i), 'g-', 'LineWidth', 1.5);
%     title(sprintf('Corrected ΔF/F0 (Cell %d)', i), 'FontSize', 10, 'FontWeight', 'bold');
%     xlabel('Time (s)', 'FontSize', 9);
%     ylabel('ΔF/F0', 'FontSize', 9);
%     grid on;
    
%     % Add zero line reference
%     hold on;
%     plot([time_axis(1), time_axis(end)], [0, 0], 'k--', 'LineWidth', 0.5);
end

% % Comparison of mean signal before and after correction
% subplot(4, 4, [11, 12, 15, 16]);
% % Compute mean signal
% mean_raw = mean(calcium_data, 2);
% mean_corrected = mean(corrected_data, 2);
% 
% yyaxis left;
% plot(time_axis, mean_raw, 'b-', 'LineWidth', 1.5);
% ylabel('Mean raw fluorescence (a.u.)', 'FontSize', 10);
% ylim([min(mean_raw)*0.9, max(mean_raw)*1.1]);
% 
% yyaxis right;
% plot(time_axis, mean_corrected, 'r-', 'LineWidth', 1.5);
% ylabel('Mean ΔF/F0', 'FontSize', 10);
% 
% title('Comparison of mean signal before and after correction', 'FontSize', 11, 'FontWeight', 'bold');
% xlabel('Time (s)', 'FontSize', 10);
% legend('Raw mean', 'Corrected mean ΔF/F0', 'Location', 'best', 'FontSize', 9);
% grid on;

% %% 5. Save results
% % Create result table
% result_table = array2table([time_axis, corrected_dfof], ...
%     'VariableNames', ['Time_sec', cell_names]);
% 
% % Save as CSV file
% output_filename = 'calcium_bleach_corrected_results.csv';
% writetable(result_table, output_filename);
% fprintf('Results saved to: %s\n', output_filename);
% 
% % Save fit parameters
% param_table = table();
% param_table.CellName = cell_names';
% for i = 1:num_cells
%     if length(fitted_params{i}) >= 3
%         param_table.Amp(i) = fitted_params{i}(1);
%         param_table.DecayRate(i) = fitted_params{i}(2);
%         param_table.Offset(i) = fitted_params{i}(3);
%     else
%         param_table.Amp(i) = NaN;
%         param_table.DecayRate(i) = NaN;
%         param_table.Offset(i) = NaN;
%     end
% end
% writetable(param_table, 'bleach_fitting_parameters.csv');
% fprintf('Fit parameters saved to: bleach_fitting_parameters.csv\n\n');

%% 
% Compute statistical metrics
num_cells = size(calcium_data, 2);

% Initialize statistics matrix
stats_table = table();
stats_table.CellName = cell_names';
stats_table.CellIndex = (1:num_cells)';
stats_table.Mean = zeros(num_cells, 1);
stats_table.Std = zeros(num_cells, 1);
stats_table.DynamicRange = zeros(num_cells, 1); % Max - Min
stats_table.CV = zeros(num_cells, 1); % Coefficient of variation = std / mean
stats_table.Kurtosis = zeros(num_cells, 1);
stats_table.Skewness = zeros(num_cells, 1); % Additional skewness
stats_table.Min = zeros(num_cells, 1);
stats_table.Max = zeros(num_cells, 1);

% Compute statistics for each cell
for i = 1:num_cells
    %signal = corrected_data(:, i);
    signal = F_over_F0(:, i);
    % Remove NaN values
    signal = signal(~isnan(signal));
    
    if ~isempty(signal)
        stats_table.Mean(i) = mean(signal);
        stats_table.Std(i) = std(signal);
        stats_table.DynamicRange(i) = max(signal) - min(signal);
        
        % Coefficient of variation (avoid division by zero)
        if stats_table.Mean(i) ~= 0
            stats_table.CV(i) = stats_table.Std(i) / abs(stats_table.Mean(i));
        else
            stats_table.CV(i) = NaN;
        end
        
        stats_table.Kurtosis(i) = kurtosis(signal);
        stats_table.Skewness(i) = skewness(signal);
    else
        stats_table.Mean(i) = NaN;
        stats_table.Std(i) = NaN;
        stats_table.DynamicRange(i) = NaN;
        stats_table.CV(i) = NaN;
        stats_table.Kurtosis(i) = NaN;
        stats_table.Skewness(i) = NaN;
    end
end
%%
% ========== Tau value extraction from skewed CV distribution ==========
% Assume you have cell CV values array: cell_cv_values

% Method 1: Inverse CDF method (most robust, recommended first choice)
sorted_cv = sort(stats_table.CV);
n = length(sorted_cv);
cdf = (1:n) / n;

% Find multiple characteristic decay points
tau_50 = interp1(cdf, sorted_cv, 0.50);  % Median
tau_63 = interp1(cdf, sorted_cv, 0.632); % Classic 1-1/e decay point
tau_80 = interp1(cdf, sorted_cv, 0.80);  % 80% decay point
tau_90 = interp1(cdf, sorted_cv, 0.90);  % 90% decay point

% Method 2: Double exponential mixture model fit (suitable for skewed distributions)
% Model: p(x) = w1*exp(-x/τ1) + w2*exp(-x/τ2), w1+w2=1
double_exp_pdf = @(p, x) p(1)*exp(-x/p(2)) + (1-p(1))*exp(-x/p(3));

% Prepare histogram data
[counts, bin_centers] = hist(stats_table.CV, 50);
bin_width = bin_centers(2) - bin_centers(1);
pdf_values = counts / (sum(counts) * bin_width);

% Initial parameter guesses
init_w1 = 0.7;  % Weight of first component
init_tau1 = mean(stats_table.CV) * 0.5;  % Fast decay component
init_tau2 = mean(stats_table.CV) * 2;    % Slow decay component

% Fit double exponential
try
    [params_double, resnorm] = lsqcurvefit(double_exp_pdf, ...
        [init_w1, init_tau1, init_tau2], ...
        bin_centers, pdf_values, ...
        [0.1, 0.001, 0.001], [0.9, Inf, Inf]); % Boundary constraints
    
    % Compute weighted average τ
    w1_fit = params_double(1);
    tau1_fit = params_double(2);
    tau2_fit = params_double(3);
    tau_weighted = w1_fit * tau1_fit + (1-w1_fit) * tau2_fit;
    
    fprintf('Double exponential fit results:\n');
    fprintf('  Fast component τ1 = %.4f (weight=%.2f)\n', tau1_fit, w1_fit);
    fprintf('  Slow component τ2 = %.4f (weight=%.2f)\n', tau2_fit, 1-w1_fit);
    fprintf('  Weighted average τ = %.4f\n', tau_weighted);
    
catch
    fprintf('Double exponential fit failed, using inverse CDF method\n');
    tau_weighted = tau_63; % Fall back to Method 1
end

% Method 3: Gamma distribution fit (designed for skewed positive distributions)
if exist('gamfit', 'file')
    [param_gam, ci_gam] = gamfit(stats_table.CV);
    shape_param = param_gam(1);  % Shape parameter k
    scale_param = param_gam(2);  % Scale parameter θ
    
    % Gamma distribution mean = k*θ, can be considered characteristic τ
    tau_gamma = shape_param * scale_param;
    
    fprintf('Gamma distribution fit:\n');
    fprintf('  Shape parameter k = %.4f, scale parameter θ = %.4f\n', shape_param, scale_param);
    fprintf('  Distribution mean τ = %.4f\n', tau_gamma);
end

% Method 4: Skewness-corrected τ (accounts for distribution skewness)
cv_skewness = skewness(stats_table.CV);
cv_mean = mean(stats_table.CV);
cv_median = median(stats_table.CV);

% Skewness correction formula: τ_skew = median * (1 + skewness correction factor)
if cv_skewness > 0
    % Right-skewed distribution: mean > median, use correction
    skew_correction = min(0.5, cv_skewness/5); % Limit correction magnitude
    tau_skew_adjusted = cv_median * (1 + skew_correction);
else
    % Left-skewed or symmetric, use median directly
    tau_skew_adjusted = cv_median;
end

% ========== Output all τ estimates ==========
fprintf('\n========== CV distribution τ value estimation summary ==========\n');
fprintf('Basic statistics:\n');
fprintf('  Mean = %.4f, Median = %.4f, Skewness = %.4f\n', cv_mean, cv_median, cv_skewness);
fprintf('Inverse CDF method:\n');
fprintf('  τ(50%%) = %.4f, τ(63%%) = %.4f\n', tau_50, tau_63);
fprintf('  τ(80%%) = %.4f, τ(90%%) = %.4f\n', tau_80, tau_90);
if exist('tau_weighted', 'var')
    fprintf('Double exponential model: weighted τ = %.4f\n', tau_weighted);
end
if exist('tau_gamma', 'var')
    fprintf('Gamma distribution: mean τ = %.4f\n', tau_gamma);
end
fprintf('Skewness-corrected: τ_skew = %.4f\n', tau_skew_adjusted);
%% Apply screening criterion: only consider CV >= 0.03   %dynamic range >= 100 and
stats_table.IsActive = stats_table.CV >= 0.035;
stats_table.ActivityType = cell(num_cells, 1);
for i = 1:num_cells
    if stats_table.IsActive(i)
        stats_table.ActivityType{i} = 'Active';
    else
        stats_table.ActivityType{i} = 'Inactive';
    end
end

% Get indices of active and inactive cells
active_indices = find(stats_table.IsActive);
inactive_indices = find(~stats_table.IsActive);

fprintf('\n=== Cell classification results ===\n');
fprintf('Number of active cells (dynamic range >= 100 and CV >= 0.02): %d\n', length(active_indices));
fprintf('Number of inactive cells: %d\n', length(inactive_indices));

% Display detailed information of active cells
if ~isempty(active_indices)
    fprintf('\n--- Active cell list ---\n');
    for i = 1:length(active_indices)
        idx = active_indices(i);
        fprintf('Cell %d: %s (dynamic range: %.2f, CV: %.4f)\n', ...
            idx, cell_names{idx}, stats_table.DynamicRange(idx), stats_table.CV(idx));
    end
end

% Display detailed information of inactive cells
if ~isempty(inactive_indices)
    fprintf('\n--- Inactive cell list ---\n');
    for i = 1:min(length(inactive_indices), 10) % Only show first 10 inactive cells
        idx = inactive_indices(i);
        fprintf('Cell %d: %s (dynamic range: %.2f, CV: %.4f)\n', ...
            idx, cell_names{idx}, stats_table.DynamicRange(idx), stats_table.CV(idx));
    end
    if length(inactive_indices) > 10
        fprintf('... and %d more inactive cells\n', length(inactive_indices) - 10);
    end
end

%% Select 5 cells from active and inactive cells respectively for visualization
% If a category has fewer than 5 cells, use all available cells
num_to_select = 5;

if length(active_indices) >= num_to_select
    selected_active = active_indices(randperm(length(active_indices), num_to_select));
else
    selected_active = active_indices;
    fprintf('\nNote: Only %d active cells, all used for visualization\n', length(active_indices));
end

if length(inactive_indices) >= num_to_select
    selected_inactive = inactive_indices(randperm(length(inactive_indices), num_to_select));
else
    selected_inactive = inactive_indices;
    fprintf('Note: Only %d inactive cells, all used for visualization\n', length(inactive_indices));
end

fprintf('\n=== Selected cells for visualization ===\n');
fprintf('Active cells: %s\n', mat2str(selected_active));
fprintf('Inactive cells: %s\n', mat2str(selected_inactive));

%% Visualize comparison of active vs inactive cells
figure('Position', [100, 100, 1400, 1000]);

% Subplot 1: Calcium signal traces of active cells
subplot(2, 3, 1);
if ~isempty(selected_active)
    colors_active = parula(length(selected_active));
    for i = 1:length(selected_active)
        idx = selected_active(i);
        plot(corrected_data(:, idx), 'Color', colors_active(i, :), 'LineWidth', 2);
        hold on;
    end
    xlabel('Time point');
    ylabel('Calcium signal intensity');
    title('Active cell calcium signal traces');
    % Create legend labels
    legend_labels_active = cell(length(selected_active), 1);
    for i = 1:length(selected_active)
        idx = selected_active(i);
        legend_labels_active{i} = sprintf('Cell %d', idx);
    end
    legend(legend_labels_active, 'Location', 'eastoutside');
    grid on;
else
    text(0.5, 0.5, 'No active cells', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title('Active cell calcium signal traces');
end

% Subplot 2: Calcium signal traces of inactive cells
subplot(2, 3, 2);
if ~isempty(selected_inactive)
    colors_inactive = parula(length(selected_inactive)); 
    for i = 1:length(selected_inactive)
        idx = selected_inactive(i);
        plot(corrected_data(:, idx), 'Color', colors_inactive(i, :), 'LineWidth', 2);
        hold on;
    end
    xlabel('Time point');
    ylabel('Calcium signal intensity');
    title('Inactive cell calcium signal traces');
    % Create legend labels
    legend_labels_inactive = cell(length(selected_inactive), 1);
    for i = 1:length(selected_inactive)
        idx = selected_inactive(i);
        legend_labels_inactive{i} = sprintf('Cell %d', idx);
    end
    legend(legend_labels_inactive, 'Location', 'eastoutside');
    grid on;
else
    text(0.5, 0.5, 'No inactive cells', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title('Inactive cell calcium signal traces');
end

% Subplot 3: Dynamic range distribution comparison
subplot(2, 3, 3);
if ~isempty(active_indices) && ~isempty(inactive_indices)
    % Create grouped data
    group = [ones(length(active_indices), 1); 2*ones(length(inactive_indices), 1)];
    data_dr = [stats_table.DynamicRange(active_indices); stats_table.DynamicRange(inactive_indices)];
    
    boxplot(data_dr, group, 'Labels', {'Active cells', 'Inactive cells'});
    ylabel('Dynamic range');
    title('Dynamic range distribution comparison');
    grid on;
else
    text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title('Dynamic range distribution comparison');
end

% Subplot 4: Coefficient of variation distribution comparison
subplot(2, 3, 4);
if ~isempty(active_indices) && ~isempty(inactive_indices)
    % Create grouped data
    group = [ones(length(active_indices), 1); 2*ones(length(inactive_indices), 1)];
    data_cv = [stats_table.CV(active_indices); stats_table.CV(inactive_indices)];
    
    boxplot(data_cv, group, 'Labels', {'Active cells', 'Inactive cells'});
    ylabel('Coefficient of variation (CV)');
    title('Coefficient of variation distribution comparison');
    grid on;
else
    text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title('Coefficient of variation distribution comparison');
end

% Subplot 5: Dynamic range vs coefficient of variation scatter plot
subplot(2, 3, 5);
if ~isempty(active_indices) && ~isempty(inactive_indices)
    scatter(stats_table.DynamicRange(active_indices), stats_table.CV(active_indices), ...
            50, 'g', 'filled', 'MarkerFaceAlpha', 0.7);
    hold on;
    scatter(stats_table.DynamicRange(inactive_indices), stats_table.CV(inactive_indices), ...
            50, 'r', 'filled', 'MarkerFaceAlpha', 0.7);
    
    % Add classification lines
    xline(100, '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    yline(0.02, '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    
    xlabel('Dynamic range');
    ylabel('Coefficient of variation (CV)');
    title('Dynamic range vs coefficient of variation');
    legend('Active cells', 'Inactive cells', 'Classification threshold', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title('Dynamic range vs coefficient of variation');
end

% Subplot 6: Statistical characteristics comparison of two cell types
subplot(2, 3, 6);
if ~isempty(active_indices) && ~isempty(inactive_indices)
    % Define statistical metrics to compare
    stat_names = {'Mean', 'Std', 'DynamicRange', 'CV', 'Kurtosis', 'Skewness'};
    stat_labels = {'Mean', 'Std', 'Dynamic range', 'CV', 'Kurtosis', 'Skewness'};
    
    % Compute mean statistics for two cell types
    active_means = zeros(1, length(stat_names));
    inactive_means = zeros(1, length(stat_names));
    
    for i = 1:length(stat_names)
        active_means(i) = mean(stats_table.(stat_names{i})(active_indices), 'omitnan');
        inactive_means(i) = mean(stats_table.(stat_names{i})(inactive_indices), 'omitnan');
    end
    
    % Create grouped bar chart - ensure X and Y lengths match
    X = 1:length(stat_names);
    bar_data = [active_means; inactive_means]';
    
    % Draw bar chart
    h = bar(X, bar_data);
    
    % Set colors
    h(1).FaceColor = [0.2, 0.8, 0.2]; % Active cells - green
    h(2).FaceColor = [0.8, 0.2, 0.2]; % Inactive cells - red
    
    xlabel('Statistical metric');
    ylabel('Mean value');
    title('Statistical characteristics comparison of two cell types');
    set(gca, 'XTickLabel', stat_labels, 'XTickLabelRotation', 45);
    legend('Active cells', 'Inactive cells', 'Location', 'best');
    grid on;
    
    % Add value labels
    for i = 1:length(X)
        % Value labels for active cells
        if ~isnan(active_means(i))
            text(X(i)-0.18, active_means(i)+max(active_means)/50, sprintf('%.2f', active_means(i)), ...
                'FontSize', 7, 'HorizontalAlignment', 'center');
        end
        
        % Value labels for inactive cells
        if ~isnan(inactive_means(i))
            text(X(i)+0.18, inactive_means(i)+max(inactive_means)/50, sprintf('%.2f', inactive_means(i)), ...
                'FontSize', 7, 'HorizontalAlignment', 'center');
        end
    end
else
    text(0.5, 0.5, 'Insufficient data', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title('Statistical characteristics comparison of two cell types');
end

sgtitle('Active vs inactive cell comparison analysis', 'FontSize', 14, 'FontWeight', 'bold');
%%
%%Visualization
% Method 1 as an example, you can modify these indices
selected_indices = [4,12,31]; % Modify to the cell indices you want to analyze
fprintf('\n=== Detailed statistical parameters of specified cells ===\n');
for i = 1:length(selected_indices)
    idx = selected_indices(i);
    fprintf('\n--- Cell %d: %s ---\n', idx, cell_names{idx});
    fprintf('Mean: %.4f\n', stats_table.Mean(idx));
    fprintf('Std: %.4f\n', stats_table.Std(idx));
    fprintf('Min: %.4f\n', stats_table.Min(idx));
    fprintf('Max: %.4f\n', stats_table.Max(idx));
    fprintf('Dynamic range: %.4f\n', stats_table.DynamicRange(idx));
    fprintf('Coefficient of variation (CV): %.4f\n', stats_table.CV(idx));
    fprintf('Kurtosis: %.4f\n', stats_table.Kurtosis(idx));
    fprintf('Skewness: %.4f\n', stats_table.Skewness(idx));
end
% Visualize signals and statistics of specified cells
figure('Position', [100, 100, 1400, 1000]);

% Subplot 1: Calcium signal traces of specified cells
subplot(2, 3, 1);
colors = hsv(length(selected_indices));
for i = 1:length(selected_indices)
    idx = selected_indices(i);
    plot(corrected_data(:, idx), 'Color', colors(i, :), 'LineWidth', 2);
    hold on;
end
xlabel('Time point');
ylabel('Calcium signal intensity');
title('Calcium signal traces of specified cells');
% Create legend labels - using indices and names
legend_labels = cell(length(selected_indices), 1);
for i = 1:length(selected_indices)
    idx = selected_indices(i);
    legend_labels{i} = sprintf('Cell %d: %s', idx, cell_names{idx});
end
legend(legend_labels, 'Location', 'eastoutside', 'Interpreter', 'none');
grid on;

% Subplot 2: Standard deviation comparison of specified cells
subplot(2, 3, 2);
bar(stats_table.Std(selected_indices), 'FaceColor', [0.2, 0.6, 0.8]);
% Use indices as x-axis labels
x_labels = cell(length(selected_indices), 1);
for i = 1:length(selected_indices)
    x_labels{i} = num2str(selected_indices(i));
end
set(gca, 'XTickLabel', x_labels);
xlabel('Cell index');
ylabel('Standard deviation');
title('Standard deviation comparison of specified cells');
grid on;

% Subplot 3: Dynamic range comparison of specified cells
subplot(2, 3, 3);
bar(stats_table.DynamicRange(selected_indices), 'FaceColor', [0.8, 0.4, 0.2]);
set(gca, 'XTickLabel', x_labels);
xlabel('Cell index');
ylabel('Dynamic range');
title('Dynamic range comparison of specified cells');
grid on;

% Subplot 4: Coefficient of variation comparison of specified cells
subplot(2, 3, 4);
bar(stats_table.CV(selected_indices), 'FaceColor', [0.4, 0.8, 0.4]);
set(gca, 'XTickLabel', x_labels);
xlabel('Cell index');
ylabel('Coefficient of variation (CV)');
title('Coefficient of variation comparison of specified cells');
grid on;

% Subplot 5: Kurtosis comparison of specified cells
subplot(2, 3, 5);
bar(stats_table.Kurtosis(selected_indices), 'FaceColor', [0.8, 0.2, 0.8]);
hold on;
plot(xlim, [3, 3], 'r--', 'LineWidth', 2); % Normal distribution reference line
set(gca, 'XTickLabel', x_labels);
xlabel('Cell index');
ylabel('Kurtosis');
title('Kurtosis comparison of specified cells (red line: normal distribution=3)');
grid on;

% Subplot 6: Skewness comparison of specified cells
subplot(2, 3, 6);
bar(stats_table.Skewness(selected_indices), 'FaceColor', [0.9, 0.7, 0.1]);
hold on;
plot(xlim, [0, 0], 'r--', 'LineWidth', 2); % Symmetric distribution reference line
set(gca, 'XTickLabel', x_labels);
xlabel('Cell index');
ylabel('Skewness');
title('Skewness comparison of specified cells (red line: symmetric distribution=0)');
grid on;

sgtitle('Calcium signal statistical analysis of specified cells', 'FontSize', 14, 'FontWeight', 'bold');

%% Visualize statistical distributions
figure('Position', [100, 100, 1400, 1000]);

% Subplot 1: Standard deviation distribution
subplot(2, 3, 1);
histogram(stats_table.Std, 20, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'black');
xlabel('Standard deviation');
ylabel('Number of cells');
title('Standard deviation distribution');
grid on;

% Subplot 2: Dynamic range distribution
subplot(2, 3, 2);
histogram(stats_table.DynamicRange, 20, 'FaceColor', [0.8, 0.4, 0.2], 'EdgeColor', 'black');
xlabel('Dynamic range');
ylabel('Number of cells');
title('Dynamic range distribution');
grid on;

% Subplot 3: Coefficient of variation distribution
subplot(2, 3, 3);
histogram(stats_table.CV, 20, 'FaceColor', [0.4, 0.8, 0.4], 'EdgeColor', 'black');
xlabel('Coefficient of variation (CV)');
ylabel('Number of cells');
title('Coefficient of variation distribution');
grid on;

% Subplot 4: Kurtosis distribution
subplot(2, 3, 4);
histogram(stats_table.Kurtosis, 20, 'FaceColor', [0.8, 0.2, 0.8], 'EdgeColor', 'black');
hold on;
% Mark normal distribution reference line
y_limits = ylim;
plot([3, 3], y_limits, 'r--', 'LineWidth', 2);
xlabel('Kurtosis');
ylabel('Number of cells');
title('Kurtosis distribution (red line: normal distribution=3)');
legend('Cell distribution', 'Normal distribution', 'Location', 'best');
grid on;

% Subplot 5: Skewness distribution
subplot(2, 3, 5);
histogram(stats_table.Skewness, 20, 'FaceColor', [0.9, 0.7, 0.1], 'EdgeColor', 'black');
hold on;
% Mark symmetric distribution reference line
y_limits = ylim;
plot([0, 0], y_limits, 'r--', 'LineWidth', 2);
xlabel('Skewness');
ylabel('Number of cells');
title('Skewness distribution (red line: symmetric distribution=0)');
legend('Cell distribution', 'Symmetric distribution', 'Location', 'best');
grid on;

% Subplot 6: Example plot of all signals
subplot(2, 3, 6);
% Randomly select a few representative signals to plot
if num_cells > 0
    if num_cells <= 5
        cells_to_plot = 1:num_cells;
    else
        cells_to_plot = randperm(num_cells, min(5, num_cells));
    end
    
    colors = lines(length(cells_to_plot));
    for i = 1:length(cells_to_plot)
        cell_idx = cells_to_plot(i);
        plot(corrected_data(:, cell_idx), 'Color', colors(i, :), 'LineWidth', 1.5);
        hold on;
    end
    xlabel('Time point');
    ylabel('Calcium signal intensity');
    title('Example calcium signal traces');
    legend(cell_names(cells_to_plot), 'Location', 'best', 'Interpreter', 'none');
    grid on;
end

sgtitle('Calcium signal statistical characteristics analysis', 'FontSize', 14, 'FontWeight', 'bold');

%%
% Get raw calcium signal data of active cells
active_cell_indices = find(stats_table.IsActive);
active_cell_names = cell_names(active_cell_indices);

active_calcium_data = calcium_data(:, active_cell_indices);
%% Compute Pearson correlation coefficient matrix of active cells
if ~isempty(active_calcium_data)
    fprintf('\n=== Computing Pearson correlation coefficient matrix of active cells ===\n');
    
    % Compute correlation coefficient matrix
    correlation_matrix = corr(active_calcium_data, 'Rows', 'pairwise');
    
    % Get number of active cells
    num_active_cells = size(active_calcium_data, 2);
    
    fprintf('Number of active cells: %d\n', num_active_cells);
    fprintf('Correlation matrix size: %d x %d\n', size(correlation_matrix));
    fprintf('Correlation coefficient range: [%.4f, %.4f]\n', min(correlation_matrix(:)), max(correlation_matrix(:)));
    
    % Extract upper triangle (excluding diagonal) for distribution analysis
    triu_indices = triu(true(size(correlation_matrix)), 1);
    correlation_values = correlation_matrix(triu_indices);
    
    fprintf('Correlation coefficient statistics:\n');
    fprintf('  Mean: %.4f\n', mean(correlation_values, 'omitnan'));
    fprintf('  Median: %.4f\n', median(correlation_values, 'omitnan'));
    fprintf('  Std: %.4f\n', std(correlation_values, 'omitnan'));
    
    % Compute proportion of significant correlations (e.g., |r| > 0.5)
    strong_positive = sum(correlation_values > 0.5) / length(correlation_values) * 100;
    strong_negative = sum(correlation_values < -0.5) / length(correlation_values) * 100;
    moderate_positive = sum(correlation_values > 0.3 & correlation_values <= 0.5) / length(correlation_values) * 100;
    moderate_negative = sum(correlation_values < -0.3 & correlation_values >= -0.5) / length(correlation_values) * 100;
    
    fprintf('Strong positive correlation (r > 0.5): %.2f%%\n', strong_positive);
    fprintf('Moderate positive correlation (0.3 < r <= 0.5): %.2f%%\n', moderate_positive);
    fprintf('Strong negative correlation (r < -0.5): %.2f%%\n', strong_negative);
    fprintf('Moderate negative correlation (-0.5 <= r < -0.3): %.2f%%\n', moderate_negative);
    
    %% Visualize correlation coefficient distribution
    figure('Position', [100, 100, 1200, 500]);
    
    % Subplot 1: Correlation coefficient distribution histogram
    subplot(1, 2, 1);
    histogram(correlation_values, 50, 'FaceColor', [0.3, 0.6, 0.9], 'EdgeColor', 'black');
    hold on;

    % Add reference lines - save handles first for legend
    y_limits = ylim;
    h_zero = plot([0, 0], y_limits, 'k--', 'LineWidth', 1.5);
    h_strong = plot([0.5, 0.5], y_limits, 'r--', 'LineWidth', 1);
    h_strong_neg = plot([-0.5, -0.5], y_limits, 'r--', 'LineWidth', 1);
    h_moderate = plot([0.3, 0.3], y_limits, 'g--', 'LineWidth', 1);
    h_moderate_neg = plot([-0.3, -0.3], y_limits, 'g--', 'LineWidth', 1);

    xlabel('Pearson correlation coefficient');
    ylabel('Frequency');
    title('Active cell correlation coefficient distribution');
    % Legend uses actual line handles to ensure color consistency
    legend([h_zero, h_strong, h_moderate], ...
        'Zero correlation', 'Strong correlation threshold (±0.5)', 'Moderate correlation threshold (±0.3)', 'Location', 'best');
    grid on;
    
    % Subplot 2: Correlation coefficient box plot
    subplot(1, 2, 2);
    boxplot(correlation_values, 'Orientation', 'horizontal');
    xlabel('Pearson correlation coefficient');
    title('Correlation coefficient box plot');
    grid on;
    
    sgtitle('Active cell Pearson correlation coefficient analysis', 'FontSize', 14, 'FontWeight', 'bold');
    
    %% Perform hierarchical clustering (get cluster_order)
    fprintf('\n=== Performing hierarchical clustering analysis ===\n');
    
    % Use hierarchical clustering to group correlation matrix
    % Compute distance matrix (1 - correlation coefficient)
    distance_matrix = 1 - correlation_matrix;
    
    % Initialize cluster_order
    cluster_order = 1:num_active_cells; % Default order
    
    if num_active_cells > 1
        % Convert distance matrix to vector format (pdist format)
        distance_vector = squareform(distance_matrix, 'tovector');
        
        % Perform hierarchical clustering
        linkage_tree = linkage(distance_vector, 'average');
        
        % Get cluster order - use dendrogram to obtain ordering
        [~, ~, cluster_order] = dendrogram(linkage_tree, 0);
        cluster_order = cluster_order';
        fprintf('Hierarchical clustering complete, obtained cell ordering.\n');
    else
        % If only one cell, use 1 directly
        cluster_order = 1;
        linkage_tree = [];
        fprintf('Only one active cell, no clustering needed.\n');
    end
    
    %% Create clustered heatmap
    figure('Position', [100, 100, max(800, num_active_cells*40), 800]);
    
    % Reorder correlation matrix and cell names according to cluster order
    sorted_correlation_matrix = correlation_matrix(cluster_order, cluster_order);
    sorted_cell_names = active_cell_names(cluster_order);
    
    % Create heatmap
    imagesc(sorted_correlation_matrix);
    colormap(jet); % Use jet color map, blue for negative correlation, red for positive correlation
    colorbar;
    caxis([-1, 1]); % Set color axis range from -1 to 1
    
    % Set axis labels - force display all labels
    xticks(1:num_active_cells);
    yticks(1:num_active_cells);
    xticklabels(sorted_cell_names);
    yticklabels(sorted_cell_names);
    
    % Dynamically adjust layout
    if num_active_cells > 20
        % Layout for many cells
        xtickangle(90);
        ytickangle(0);
        font_size = max(6, 10 - num_active_cells/20);
        set(gca, 'FontSize', font_size);
        
        % Adjust figure size to fit all labels
        fig_position = get(gcf, 'Position');
        set(gcf, 'Position', [100, 100, max(1000, num_active_cells*30), 800]);
    elseif num_active_cells > 10
        % Layout for moderate number of cells
        xtickangle(45);
        ytickangle(0);
        set(gca, 'FontSize', 9);
    else
        % Layout for few cells
        xtickangle(0);
        ytickangle(0);
        set(gca, 'FontSize', 10);
    end
    
    xlabel('Cell');
    ylabel('Cell');
    title('Active cell Pearson correlation coefficient heatmap (with clustering)');
    
    % Add grid lines for better observation
    hold on;
    for i = 1:num_active_cells
        plot([0.5, num_active_cells+0.5], [i+0.5, i+0.5], 'k-', 'LineWidth', 0.5);
        plot([i+0.5, i+0.5], [0.5, num_active_cells+0.5], 'k-', 'LineWidth', 0.5);
    end
    
    % Display dendrogram separately on the right side of heatmap
    if num_active_cells > 1
        dendro_axes = axes('Position', [0.85, 0.1, 0.12, 0.8]);
        dendrogram(linkage_tree, 0, 'Orientation', 'right');
        set(dendro_axes, 'XTick', [], 'YTick', [], 'YAxisLocation', 'right');
        title(dendro_axes, 'Cluster tree');
    end
    
    %% Create network visualization (if multiple cells)
    if num_active_cells > 1
        figure('Position', [100, 100, 1000, 800]);
        
        % Only show relatively strong correlations (e.g., |r| > 0.3)
        threshold = 0.3;
        strong_correlations = sorted_correlation_matrix;
        strong_correlations(abs(strong_correlations) < threshold) = 0;
        
        % Create graph object
        G = graph(strong_correlations, sorted_cell_names, 'upper');
        
        % Compute node sizes (based on dynamic range)
        node_sizes = zeros(length(sorted_cell_names), 1);
        for i = 1:length(sorted_cell_names)
            cell_idx = active_cell_indices(cluster_order(i));
            node_sizes(i) = 100 + stats_table.DynamicRange(cell_idx) / 10; % Scaling factor
        end
        
        % Compute node colors (based on mean calcium signal)
        node_colors = zeros(length(sorted_cell_names), 1);
        for i = 1:length(sorted_cell_names)
            cell_idx = active_cell_indices(cluster_order(i));
            node_colors(i) = stats_table.Mean(cell_idx);
        end
        
        % Draw network graph
        p = plot(G, 'Layout', 'force', 'UseGravity', true);
        
        % Set node properties
        p.MarkerSize = node_sizes / 20; % Further scaling
        p.NodeCData = node_colors;
        p.NodeLabel = sorted_cell_names;
        
        % Set edge properties
        edge_weights = G.Edges.Weight;
        positive_edges = edge_weights > 0;
        negative_edges = edge_weights < 0;
        
        % Red edges for positive correlation, blue edges for negative correlation
        p.EdgeColor = zeros(length(edge_weights), 3);
        p.EdgeColor(positive_edges, :) = repmat([1, 0, 0], sum(positive_edges), 1);
        p.EdgeColor(negative_edges, :) = repmat([0, 0, 1], sum(negative_edges), 1);
        
        % Set edge width based on correlation strength
        p.LineWidth = 2 * abs(edge_weights);
        
        colorbar;
        colormap(parula);
        title(sprintf('Active cell correlation network (|r| > %.1f)', threshold));
        xlabel('Node size: dynamic range | Node color: mean calcium signal | Edge color: red=positive correlation, blue=negative correlation');
    end
    
else
    fprintf('\nNo active cells, cannot compute correlation coefficient matrix.\n');
end