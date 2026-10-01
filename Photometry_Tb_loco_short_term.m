```matlab
%% 
clear all;
clc

%read data
filename = ''; 
data = csvread(filename,1,0);

%% 
%Set analyzed channel number
channel = 1;
%Order of polynomial fit
order  = 1;
%animal_id = 'J638';%for storage

if channel == 0 
        signal_chnl = 1;
        reference_chnl = 3;
end

if channel == 1
        signal_chnl = 4;
        reference_chnl = 6;
end
if channel == 2
        signal_chnl = 7;
        reference_chnl = 9;
end

if channel == 3
        signal_chnl = 10;
        reference_chnl = 12;
end
signal = data(:, signal_chnl);
reference = data(:, reference_chnl);
time = data(:,end);

% More accurate way: trim based on time values
% Define offset and duration
offset_sec = 0;
duration_sec = 7200;

% Trim
idx_trim = (time >= offset_sec) & (time <= offset_sec + duration_sec);
time = time(idx_trim);
signal = signal(idx_trim);
reference = reference(idx_trim);
fprintf('Number of points after data trimming = %d (time range %.2f - %.2f s)\n', length(time), time(1), time(end));

%Check the whole picture of the signal and reference data 
figure;
subplot(2,1,1);
    plot(signal,'color','g', 'LineWidth', 1.5);
subplot(2,1,2);
    plot(reference,'b', 'LineWidth', 1.5);
    grid on;
    
% ----- Outlier detection and repair -----
% Set threshold (here using 10% of median as example, adjust based on data)
threshold = 0.2 * median(reference);

% Find outliers (values below threshold, or equal to 0)
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
%
%signal = signal(101:end);
%reference = reference(101:end);

%% ----- Optimized fitting: use only low reference values -----
percentile_fit = 50;  % Use points below the 50th percentile (below median)
thresh_fit = prctile(signal, percentile_fit);
idx_low_fit = signal <= thresh_fit;

if sum(idx_low_fit) < 10
    warning('Insufficient low reference value points, using all data for fitting');
    fit_coeffs = polyfit(reference, signal, order);
else
    fit_coeffs = polyfit(reference(idx_low_fit), signal(idx_low_fit), order);
end
fitted_reference = polyval(fit_coeffs, reference);

% Modified: compute F/F0 instead of ΔF/F0
[F_F0, z_score] = calculate(signal, fitted_reference);

%F_F0 = smooth(F_F0);  % Optional smoothing

figure;
  subplot(2,1,1);
  hold on
    plot(signal,'color','g', 'LineWidth', 1.5);
    plot(fitted_reference,'color','black', 'LineWidth', 1.5);
  hold off
  subplot(2,1,2);
    plot(F_F0,'b', 'LineWidth', 1.5);
    %ylim([-10, 15]);  % No longer applicable, F/F0 is typically around 1
  grid on;
  

%% ===== Time-axis-based window analysis (original sampling rate) =====
% Read body temperature
combined_filename = 'F:\PhD Thesis\Fiber photometry\Cre line_Ca2+ GCamP\Processed_file-short-term\Tb-loco-simultaneous-recordings\1-J944-CH1-260712_time_120000-140000.csv';
data_all = csvread(combined_filename, 0, 0);   % No header row
temp_value = data_all(:, 1);   % First column: body temperature
loco_raw = data_all(:, 2);     % Second column: locomotion (total per 15 s)

N_temp = length(temp_value);
temp_time = (0:N_temp-1)' * 15;   % Temperature time axis, in seconds

% Trim temperature and locomotion to the same time window
idx_temp_trim = (temp_time >= offset_sec) & (temp_time <= offset_sec + duration_sec);
temp_time = temp_time(idx_temp_trim);
temp_value = temp_value(idx_temp_trim);
loco_raw = loco_raw(idx_temp_trim);   % Locomotion trimmed synchronously
fprintf('Number of points after temperature trimming = %d (time range %.2f - %.2f s)\n', length(temp_value), temp_time(1), temp_time(end));
fprintf('Number of points after locomotion trimming = %d\n', length(loco_raw));

% Ensure variable lengths match
if length(F_F0) ~= length(time)
    error('F_F0 and time lengths do not match, please check trimming');
end

% Window parameters (in seconds)
win_sec = 180;     % Window length 3 minutes
step_sec = 180;     % Step size 1 minute, windows overlap

% Generate window start times (from 0 until max time - window length)
max_start = max(time) - win_sec;
last_start = ceil(max_start / step_sec) * step_sec;
t_starts = offset_sec : step_sec : last_start;
n_windows = length(t_starts);

% Preallocate
ca_std = nan(n_windows, 1);
ca_mean = nan(n_windows, 1);
cv_window = nan(n_windows, 1);   % Added: CV for each window
temp_slope = nan(n_windows, 1);
temp_mean = nan(n_windows, 1);
loco_total = nan(n_windows, 1);  % Keep only locomotion total, remove loco_mean

for i = 1:n_windows
    t0 = t_starts(i);
    t1 = t0 + win_sec;
    
    % Calcium signal indices (based on time)
    idx_ca = (time >= t0) & (time < t1);
    if sum(idx_ca) < 5   % At least 5 points
        continue;
    end
    ca_win = F_F0(idx_ca);
    ca_std(i) = std(ca_win, 'omitnan');
    ca_mean(i) = mean(ca_win, 'omitnan');
    cv_window(i) = ca_std(i) / ca_mean(i);   % CV within window
    
    % Temperature indices (based on temp_time)
    idx_temp = (temp_time >= t0) & (temp_time < t1);
    if sum(idx_temp) < 2   % At least 2 points to compute slope
        continue;
    end
    temp_win = temp_value(idx_temp);
    temp_time_win = temp_time(idx_temp);
    % Linear regression for slope (°C/s)
    p = polyfit(temp_time_win, temp_win, 1);
    temp_slope(i) = p(1);
    temp_mean(i) = mean(temp_win, 'omitnan');
    
    % Locomotion indices (same time axis as temperature)
    idx_loco = (temp_time >= t0) & (temp_time < t1);
    if sum(idx_loco) < 2
        continue;
    end
    loco_win = loco_raw(idx_loco);
    loco_total(i) = sum(loco_win, 'omitnan');   % Compute total only
end

% Remove invalid windows (remove loco_mean related)
valid = ~isnan(ca_std) & ~isnan(cv_window) & ~isnan(temp_slope) & ~isnan(temp_mean) & ...
        ~isnan(loco_total);
ca_std = ca_std(valid);
ca_mean = ca_mean(valid);
cv_window = cv_window(valid);
temp_slope = temp_slope(valid);
temp_mean = temp_mean(valid);
loco_total = loco_total(valid);
t_starts = t_starts(valid);
fprintf('Number of valid windows = %d\n', length(ca_std));

% Correlation analysis
% 1. Temperature change rate vs Ca fluctuation (here Ca fluctuation uses ca_std, cv_window can also be used)
[rho_slope, p_slope] = corr(temp_slope, ca_std, 'Type', 'Spearman');
fprintf('\n--- Temperature change rate vs Ca fluctuation (std) ---\n');
fprintf('Spearman correlation coefficient = %.3f, p = %.4f\n', rho_slope, p_slope);
if p_slope < 0.05
    if rho_slope > 0, fprintf('Significant positive correlation (faster temperature change, larger Ca fluctuation)\n');
    else, fprintf('Significant negative correlation (faster temperature change, smaller Ca fluctuation)\n'); end
else
    fprintf('No significant correlation\n');
end

% 2. Mean temperature vs Ca fluctuation (std)
[rho_mean, p_mean] = corr(temp_mean, ca_std, 'Type', 'Spearman');
fprintf('\n--- Mean temperature vs Ca fluctuation (std) ---\n');
fprintf('Spearman correlation coefficient = %.3f, p = %.4f\n', rho_mean, p_mean);
if p_mean < 0.05
    if rho_mean < 0, fprintf('Significant negative correlation (lower temperature, larger Ca fluctuation)\n');
    else, fprintf('Significant positive correlation\n'); end
else
    fprintf('No significant correlation\n');
end

% 3. Mean temperature vs Ca mean (F/F0 mean)
[rho_mean_ca, p_mean_ca] = corr(temp_mean, ca_mean, 'Type', 'Spearman');
fprintf('\n--- Mean temperature vs Ca mean (F/F0) ---\n');
fprintf('Spearman correlation coefficient = %.3f, p = %.4f\n', rho_mean_ca, p_mean_ca);
if p_mean_ca < 0.05
    if rho_mean_ca < 0, fprintf('Significant negative correlation (lower temperature, larger Ca mean)\n');
    else, fprintf('Significant positive correlation\n'); end
else
    fprintf('No significant correlation\n');
end

% ---- Normalized coefficient of variation (CV_norm) vs body temperature correlation ----
% Compute global F/F0 mean (for normalization)
global_ca_mean = mean(ca_mean, 'omitnan');
if abs(global_ca_mean) < 0.01
    warning('Global F/F0 mean close to 0, CV_norm may be too large, please check data.');
end

% Compute CV_norm for each window
cv_norm = ca_std ./ abs(global_ca_mean);

% Correlation: CV_norm vs mean temperature
[rho_cv_norm, p_cv_norm] = corr(temp_mean, cv_norm, 'Type', 'Spearman');
fprintf('\n--- Normalized coefficient of variation (CV_norm) vs mean temperature ---\n');
fprintf('Spearman correlation coefficient = %.3f, p = %.4f\n', rho_cv_norm, p_cv_norm);
if p_cv_norm < 0.05
    if rho_cv_norm < 0
        fprintf('Significant negative correlation (lower temperature, larger CV_norm)\n');
    else
        fprintf('Significant positive correlation (higher temperature, larger CV_norm)\n');
    end
else
    fprintf('No significant correlation\n');
end

% Added: window CV (cv_window) vs mean temperature
[rho_cv_win_temp, p_cv_win_temp] = corr(temp_mean, cv_window, 'Type', 'Spearman');
fprintf('\n--- Window CV (std/mean) vs mean temperature ---\n');
fprintf('Spearman correlation coefficient = %.3f, p = %.4f\n', rho_cv_win_temp, p_cv_win_temp);
if p_cv_win_temp < 0.05
    if rho_cv_win_temp < 0
        fprintf('Significant negative correlation (lower temperature, larger CV)\n');
    else
        fprintf('Significant positive correlation (higher temperature, larger CV)\n');
    end
else
    fprintf('No significant correlation\n');
end

%% Correlation analysis (locomotion-related, using only loco_total)
if length(loco_total) > 5
    % Locomotion total vs Ca mean
    [rho_loco_total_mean, p_loco_total_mean] = corr(loco_total, ca_mean, 'Type', 'Spearman');
    % Locomotion total vs Ca fluctuation (cv_norm)
    [rho_loco_total_cvnorm, p_loco_total_cvnorm] = corr(loco_total, cv_norm, 'Type', 'Spearman');
    % Locomotion total vs window CV (cv_window)
    [rho_loco_total_cvwin, p_loco_total_cvwin] = corr(loco_total, cv_window, 'Type', 'Spearman');
    
    fprintf('\n--- Locomotion and calcium signal window correlation ---\n');
    fprintf('Locomotion total vs Ca mean: rho = %.3f, p = %.4f\n', rho_loco_total_mean, p_loco_total_mean);
    fprintf('Locomotion total vs CV_norm: rho = %.3f, p = %.4f\n', rho_loco_total_cvnorm, p_loco_total_cvnorm);
    fprintf('Locomotion total vs window CV: rho = %.3f, p = %.4f\n', rho_loco_total_cvwin, p_loco_total_cvwin);
else
    fprintf('Insufficient valid windows, skipping locomotion correlation analysis.\n');
end

%% Plot: locomotion vs Ca fluctuation scatter plot (if windows > 5)
if length(loco_total) > 5
    figure('Name', 'Locomotion total vs CV_norm');
    scatter(loco_total, cv_norm, 30, 'filled', 'MarkerFaceAlpha', 0.5);
    xlabel('Locomotion total (total activity in window)');
    ylabel('CV\_norm (std / |global mean|)');
    title('Locomotion total vs CV\_norm');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    coef = polyfit(loco_total, cv_norm, 1);
    x_fit = linspace(min(loco_total), max(loco_total), 100);
    y_fit = polyval(coef, x_fit);
    hold on; plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
    legend('Data points', 'Linear trend', 'Location', 'best'); hold off;

    figure('Name', 'Locomotion total vs Ca mean');
    scatter(loco_total, ca_mean, 30, 'filled', 'MarkerFaceAlpha', 0.5);
    xlabel('Locomotion total (total activity in window)');
    ylabel('Ca mean (F/F0)');
    title('Locomotion total vs Ca mean (F/F0)');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    coef2 = polyfit(loco_total, ca_mean, 1);
    x_fit2 = linspace(min(loco_total), max(loco_total), 100);
    y_fit2 = polyval(coef2, x_fit2);
    hold on; plot(x_fit2, y_fit2, 'k-', 'LineWidth', 2);
    legend('Data points', 'Linear trend', 'Location', 'best'); hold off;
end

% Plot locomotion and calcium signal over time (dual Y-axis)
figure('Name', 'Calcium signal and locomotion over time');
yyaxis left;
plot(time / 60, F_F0, 'g-', 'LineWidth', 0.8);
ylabel('F/F0');
% ylim automatically adjusted, no longer fixed
yyaxis right;
loco_interp = interp1(temp_time, loco_raw, time, 'nearest', 0);
plot(time / 60, loco_interp, 'b-', 'LineWidth', 0.8);
ylabel('Locomotion (activity)');
ylim([0, 300]);
xlabel('Time (minutes)');
title('Calcium signal (F/F0) and locomotion over 2 hours');
xlim([0, 120]);
xticks(0:30:120);
grid off;
ax = gca;
ax.XAxis.TickDirection = 'out';
ax.YAxis(1).TickDirection = 'out';
ax.YAxis(2).TickDirection = 'out';

% Added: separated F/F0 and locomotion total over time (stacked subplots)
figure('Name', 'Calcium signal and locomotion total separated');
subplot(2,1,1);
plot(time / 60, F_F0, 'g-', 'LineWidth', 0.8);
ylabel('F/F0');
title('Calcium signal (F/F0) over 2 hours');
grid off;
ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
xlim([0, 120]); xticks(0:30:120);
subplot(2,1,2);
loco_interp_total = interp1(temp_time, loco_raw, time, 'nearest', 0);
plot(time / 60, loco_interp_total, 'b-', 'LineWidth', 0.8);
ylabel('Locomotion total (activity)');
xlabel('Time (minutes)');
title('Locomotion total over 2 hours');
ylim([0, 50]);
grid off;
ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
xlim([0, 120]); xticks(0:30:120);

% Added: dual Y-axis plot of locomotion total and CV_norm
t_win_min = t_starts / 60;

figure('Name', 'Loco total vs CV_norm over time');
yyaxis right;
plot(t_win_min, loco_total, 'b-', 'LineWidth', 1.5);
ylabel('Locomotion total (sum within window)');
ax = gca; ax.YAxis(1).Color = 'b'; ax.YAxis(1).TickDirection = 'out';
yyaxis left;
plot(t_win_min, cv_norm, 'Color', [0 0.5 0], 'LineWidth', 1.5);
ylabel('CV\_norm (CaStd / |GlobalCaMean|)');
ax.YAxis(2).Color = [0 0.5 0]; ax.YAxis(2).TickDirection = 'out';
xlabel('Time (minutes)');
title('Activity total and CV\_norm over time');
xlim([0, 120]); xticks(0:30:120);
grid off;
ax.XAxis.TickDirection = 'out';
legend('Locomotion total', 'CV\_norm', 'Location', 'best');

%% ---- Plot CV_norm and body temperature over time (dual Y-axis) ----
figure('Name', 'CV_norm and temperature over time');
t_min = t_starts / 60;

yyaxis left;
plot(t_min, cv_norm, 'Color', [0 0.5 0], 'LineWidth', 1.5);
ylabel('CV\_norm (std / |global mean|)');
ax = gca;
ax.YAxis(1).Color = [0 0.5 0];
ax.YAxis(1).TickDirection = 'out';

yyaxis right;
plot(t_min, temp_mean, 'b-', 'LineWidth', 1.5);
ylabel('Temperature (°C)');
ylim([min(temp_mean)-0.5, max(temp_mean)+0.5]);
ax.YAxis(2).Color = 'b';
ax.YAxis(2).TickDirection = 'out';

xlabel('Time (minutes)');
title('CV\_norm and temperature over time');
xlim([0, max(t_min)]);
xticks(0:30:max(t_min));
legend('CV\_norm', 'Temperature', 'Location', 'best');
grid off;
ax.XAxis.TickDirection = 'out';

% Figure 2: Mean temperature vs Ca CV (scatter + linear trend line) using cv_norm
figure('Name', 'Mean temp vs CV_norm');
scatter(temp_mean, cv_norm, 30, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Mean temperature (°C)');
ylabel('CV\_norm');
title('Mean temperature within window vs CV\_norm');
grid off;
ax = gca;
ax.XAxis.TickDirection = 'out';
ax.YAxis.TickDirection = 'out';

coef3 = polyfit(temp_mean, cv_norm, 1);
x_fit3 = linspace(min(temp_mean), max(temp_mean), 100);
y_fit3 = polyval(coef3, x_fit3);
hold on;
plot(x_fit3, y_fit3, 'k-', 'LineWidth', 2);
legend('Data points', 'Linear trend', 'Location', 'best');

% Added: mean temperature vs window CV (cv_window)
figure('Name', 'Mean temp vs window CV');
scatter(temp_mean, cv_window, 30, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Mean temperature (°C)');
ylabel('Window CV (std/mean)');
title('Mean temperature within window vs window CV');
grid off;
ax = gca;
ax.XAxis.TickDirection = 'out';
ax.YAxis.TickDirection = 'out';
coef_cv = polyfit(temp_mean, cv_window, 1);
x_fit_cv = linspace(min(temp_mean), max(temp_mean), 100);
y_fit_cv = polyval(coef_cv, x_fit_cv);
hold on;
plot(x_fit_cv, y_fit_cv, 'k-', 'LineWidth', 2);
legend('Data points', 'Linear trend', 'Location', 'best');

%% Plot
% Figure 1: Temperature change rate vs Ca fluctuation (std)
figure('Name', 'Slope vs Ca std');
scatter(temp_slope, ca_std, 10, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Temperature change rate (°C/s)');
ylabel('Ca fluctuation (std F/F0)');
title('Temperature change rate within window vs Ca fluctuation (std)');
grid on;
coef = polyfit(temp_slope, ca_std, 1);
x_fit = linspace(min(temp_slope), max(temp_slope), 100);
y_fit = polyval(coef, x_fit);
hold on; plot(x_fit, y_fit, 'k--', 'LineWidth', 2);
legend('Data points', 'Linear trend');

% Figure 2: Mean temperature vs Ca fluctuation (std)
figure('Name', 'Mean temp vs Ca std');
scatter(temp_mean, ca_std, 30, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Mean temperature (°C)');
ylabel('Ca fluctuation (std F/F0)');
title('Mean temperature within window vs Ca fluctuation (std)');
grid off;
ax = gca;
ax.XAxis.TickDirection = 'out';
ax.YAxis.TickDirection = 'out';
coef2 = polyfit(temp_mean, ca_std, 1);
x_fit2 = linspace(min(temp_mean), max(temp_mean), 100);
y_fit2 = polyval(coef2, x_fit2);
hold on;
plot(x_fit2, y_fit2, 'k-', 'LineWidth', 2);
legend('Data points', 'Linear trend', 'Location', 'best');

% ===== Visualize Ca standard deviation and mean over time =====
t_min = t_starts / 60;

% Figure 1: Ca mean over time (green) and body temperature
figure('Name', 'Ca mean and temperature over time');
yyaxis left;
plot(t_min, ca_mean, 'g-', 'LineWidth', 1.5);
ylabel('Ca mean (F/F0)');
yyaxis right;
plot(t_min, temp_mean, 'b-', 'LineWidth', 1.5);
ylabel('Temperature (°C)');
ylim([34, 37]);
yticks(35:0.5:37);
xlabel('Time (minutes)');
title('Ca mean (F/F0) and temperature over time');
xlim([0, 180]);
xticks(0:30:180);
legend('Ca mean', 'Temperature', 'Location', 'best');
grid off;
ax2 = gca;
ax2.XAxis.TickDirection = 'out';
ax2.YAxis(1).TickDirection = 'out';
ax2.YAxis(2).TickDirection = 'out';

% Figure 2: Ca standard deviation (fluctuation) and temperature (dual Y-axis)
figure('Name', 'Ca std and temperature over time');
yyaxis left;
plot(t_min, ca_std, 'Color', [0 0.5 0], 'LineWidth', 1.5);
ylabel('Ca std (F/F0)');
ax = gca;
ax.YAxis(1).Color = [0 0.5 0];
ax.YAxis(1).TickDirection = 'out';
yyaxis right;
plot(t_min, temp_mean, 'b-', 'LineWidth', 1.5);
ylabel('Temperature (°C)');
ylim([34, 38]);
yticks(34:0.5:38);
ax.YAxis(2).Color = 'b';
ax.YAxis(2).TickDirection = 'out';
xlabel('Time (minutes)');
title('Ca fluctuation (std F/F0) and temperature over time');
xlim([0, 120]);
xticks(0:30:120);
legend('Ca std', 'Temperature', 'Location', 'best');
grid off;
ax.XAxis.TickDirection = 'out';

% Added: window CV over time and body temperature
figure('Name', 'Window CV and temperature over time');
yyaxis left;
plot(t_min, cv_window, 'Color', [0 0.5 0], 'LineWidth', 1.5);
ylabel('Window CV (std/mean)');
ax = gca;
ax.YAxis(1).Color = [0 0.5 0];
ax.YAxis(1).TickDirection = 'out';
yyaxis right;
plot(t_min, temp_mean, 'b-', 'LineWidth', 1.5);
ylabel('Temperature (°C)');
ylim([34, 38]);
yticks(34:0.5:38);
ax.YAxis(2).Color = 'b';
ax.YAxis(2).TickDirection = 'out';
xlabel('Time (minutes)');
title('Window CV and temperature over time');
xlim([0, 120]);
xticks(0:30:120);
legend('Window CV', 'Temperature', 'Location', 'best');
grid off;
ax.XAxis.TickDirection = 'out';

% Create new figure: raw signals
figure('Name', 'Raw signals over 2 hours');
ax1 = subplot(2,1,1);
plot(temp_time / 60, temp_value, 'b-', 'LineWidth', 1.5);
xlim([0, 180]);
ylabel('Temperature (°C)');
title('Body Temperature');
grid on;
ax2 = subplot(2,1,2);
plot(time / 60, F_F0, 'g-', 'LineWidth', 1);
xlim([0, 180]);
xlabel('Time (min)');
ylabel('F/F0');
title('Calcium signal (F/F0)');
grid on;
xticks(ax1, [0, 30, 60, 90, 120]);
xticks(ax2, [0, 30, 60, 90, 120]);
xlabel(ax1, '');

%% Function definitions
function [F_F0, z_score] = calculate(curr_signal, fitted_reference)
    % Compute F/F0
    F_F0 = curr_signal ./ fitted_reference;
    average_F_F0 = mean(F_F0);
    std_F_F0 = std(F_F0);
    z_score = (F_F0 - average_F_F0) ./ std_F_F0;
end

function [fit_coeffs, fitted_signal] = fitReferenceToSignal(curr_reference, fitting_signal, fitting_reference, order)
    fit_coeffs = polyfit(fitting_reference, fitting_signal, order);
    fitted_signal = polyval(fit_coeffs, curr_reference);
end