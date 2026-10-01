%% 
%LZY script_20260714
clear all;
clc

%read data
filename = 'F:\PhD Thesis\Fiber photometry\Cre line_Ca2+ GCamP\Processing_raw_data\Continuous-2-hr\1-20260712_120000.csv'; 
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

% 更准确的方式：根据time值截取
% 定义偏移和时长
offset_sec = 0;
duration_sec = 7200;

% 截取
idx_trim = (time >= offset_sec) & (time <= offset_sec + duration_sec);
time = time(idx_trim);
signal = signal(idx_trim);
reference = reference(idx_trim);
fprintf('数据截取后点数 = %d (时间范围 %.2f - %.2f 秒)\n', length(time), time(1), time(end));

%Check the whole picture of the signal and reference data 
figure;
subplot(2,1,1);
    plot(signal,'color','g', 'LineWidth', 1.5);
subplot(2,1,2);
    plot(reference,'b', 'LineWidth', 1.5);
    grid on;
    
% ----- 异常值检测与修复 -----
% 设定阈值（这里以中位数的10%为例，您可根据数据调整）
threshold = 0.2 * median(reference);

% 找出异常点（值小于阈值，或为0）
outlier_idx = (reference < threshold) | (reference == 0);

if any(outlier_idx)
    fprintf('检测到 %d 个异常参考值，将进行线性插值修复。\n', sum(outlier_idx));
    % 临时替换为NaN
    ref_corrected = reference;
    ref_corrected(outlier_idx) = NaN;
    % 线性插值填补
    ref_corrected = fillmissing(ref_corrected, 'linear');
    % 如果首尾出现NaN（极端情况），用最近的正常值填充
    ref_corrected = fillmissing(ref_corrected, 'nearest');
else
    ref_corrected = reference;
end

% 将修正后的参考信号赋给后续使用的变量
reference = ref_corrected;

figure;
plot(reference,'b', 'LineWidth', 1.5);
grid on;
%
%signal = signal(101:end);
%reference = reference(101:end);

%% ----- 优化拟合：仅使用低参考值 -----
percentile_fit = 50;  % 使用低于50%分位数的点（中位数以下）
thresh_fit = prctile(signal, percentile_fit);
idx_low_fit = signal <= thresh_fit;

if sum(idx_low_fit) < 10
    warning('低参考值点不足，改用全部数据拟合');
    fit_coeffs = polyfit(reference, signal, order);
else
    fit_coeffs = polyfit(reference(idx_low_fit), signal(idx_low_fit), order);
end
fitted_reference = polyval(fit_coeffs, reference);

% 修改：计算 F/F0 而不是 ΔF/F0
[F_F0, z_score] = calculate(signal, fitted_reference);

%F_F0 = smooth(F_F0);  % 可选平滑

figure;
  subplot(2,1,1);
  hold on
    plot(signal,'color','g', 'LineWidth', 1.5);
    plot(fitted_reference,'color','black', 'LineWidth', 1.5);
  hold off
  subplot(2,1,2);
    plot(F_F0,'b', 'LineWidth', 1.5);
    %ylim([-10, 15]);  % 不再适用，F/F0 通常在 1 附近
  grid on;
  

%% ===== 基于时间轴的窗口分析（原始采样率） =====
% 读取体温
combined_filename = 'F:\PhD Thesis\Fiber photometry\Cre line_Ca2+ GCamP\Processed_file-short-term\Tb-loco-simultaneous-recordings\1-J944-CH1-260712_time_120000-140000.csv';
data_all = csvread(combined_filename, 0, 0);   % 无标题行
temp_value = data_all(:, 1);   % 第一列：体温
loco_raw = data_all(:, 2);     % 第二列：活动量（每15秒总量）

N_temp = length(temp_value);
temp_time = (0:N_temp-1)' * 15;   % 体温时间轴，单位秒

% 截取体温和活动量到相同时间窗口
idx_temp_trim = (temp_time >= offset_sec) & (temp_time <= offset_sec + duration_sec);
temp_time = temp_time(idx_temp_trim);
temp_value = temp_value(idx_temp_trim);
loco_raw = loco_raw(idx_temp_trim);   % 活动量同步截取
fprintf('体温截取后点数 = %d (时间范围 %.2f - %.2f 秒)\n', length(temp_value), temp_time(1), temp_time(end));
fprintf('活动量截取后点数 = %d\n', length(loco_raw));

% 确保变量长度一致
if length(F_F0) ~= length(time)
    error('F_F0和time长度不一致，请检查截取');
end

% 窗口参数（单位：秒）
win_sec = 180;     % 窗口长度3分钟
step_sec = 180;     % 步长1分钟，窗口重叠

% 生成窗口起始时间（从0开始，直到最大时间-窗口长度）
max_start = max(time) - win_sec;
last_start = ceil(max_start / step_sec) * step_sec;
t_starts = offset_sec : step_sec : last_start;
n_windows = length(t_starts);

% 预分配
ca_std = nan(n_windows, 1);
ca_mean = nan(n_windows, 1);
cv_window = nan(n_windows, 1);   % 新增：每个窗口的 CV
temp_slope = nan(n_windows, 1);
temp_mean = nan(n_windows, 1);
loco_total = nan(n_windows, 1);  % 仅保留活动量总和，删除 loco_mean

for i = 1:n_windows
    t0 = t_starts(i);
    t1 = t0 + win_sec;
    
    % 钙信号索引（根据time）
    idx_ca = (time >= t0) & (time < t1);
    if sum(idx_ca) < 5   % 至少5个点
        continue;
    end
    ca_win = F_F0(idx_ca);
    ca_std(i) = std(ca_win, 'omitnan');
    ca_mean(i) = mean(ca_win, 'omitnan');
    cv_window(i) = ca_std(i) / ca_mean(i);   % 窗口内 CV
    
    % 体温索引（根据temp_time）
    idx_temp = (temp_time >= t0) & (temp_time < t1);
    if sum(idx_temp) < 2   % 至少2个点才可算斜率
        continue;
    end
    temp_win = temp_value(idx_temp);
    temp_time_win = temp_time(idx_temp);
    % 线性回归求斜率（℃/秒）
    p = polyfit(temp_time_win, temp_win, 1);
    temp_slope(i) = p(1);
    temp_mean(i) = mean(temp_win, 'omitnan');
    
    % 活动量索引（与体温时间轴一致）
    idx_loco = (temp_time >= t0) & (temp_time < t1);
    if sum(idx_loco) < 2
        continue;
    end
    loco_win = loco_raw(idx_loco);
    loco_total(i) = sum(loco_win, 'omitnan');   % 只计算总和
end

% 去除无效窗口（删除 loco_mean 相关）
valid = ~isnan(ca_std) & ~isnan(cv_window) & ~isnan(temp_slope) & ~isnan(temp_mean) & ...
        ~isnan(loco_total);
ca_std = ca_std(valid);
ca_mean = ca_mean(valid);
cv_window = cv_window(valid);
temp_slope = temp_slope(valid);
temp_mean = temp_mean(valid);
loco_total = loco_total(valid);
t_starts = t_starts(valid);
fprintf('有效窗口数 = %d\n', length(ca_std));

% 相关性分析
% 1. 温度变化率 vs Ca波动（这里 Ca波动 使用 ca_std，也可以使用 cv_window）
[rho_slope, p_slope] = corr(temp_slope, ca_std, 'Type', 'Spearman');
fprintf('\n--- 温度变化率 vs Ca波动 (std) ---\n');
fprintf('斯皮尔曼相关系数 = %.3f, p = %.4f\n', rho_slope, p_slope);
if p_slope < 0.05
    if rho_slope > 0, fprintf('显著正相关（温度变化越快，Ca波动越大）\n');
    else, fprintf('显著负相关（温度变化越快，Ca波动越小）\n'); end
else
    fprintf('无显著相关\n');
end

% 2. 平均温度 vs Ca波动 (std)
[rho_mean, p_mean] = corr(temp_mean, ca_std, 'Type', 'Spearman');
fprintf('\n--- 平均温度 vs Ca波动 (std) ---\n');
fprintf('斯皮尔曼相关系数 = %.3f, p = %.4f\n', rho_mean, p_mean);
if p_mean < 0.05
    if rho_mean < 0, fprintf('显著负相关（温度越低，Ca 波动越大）\n');
    else, fprintf('显著正相关\n'); end
else
    fprintf('无显著相关\n');
end

% 3. 平均温度 vs Ca mean (F/F0均值)
[rho_mean_ca, p_mean_ca] = corr(temp_mean, ca_mean, 'Type', 'Spearman');
fprintf('\n--- 平均温度 vs Ca mean (F/F0) ---\n');
fprintf('斯皮尔曼相关系数 = %.3f, p = %.4f\n', rho_mean_ca, p_mean_ca);
if p_mean_ca < 0.05
    if rho_mean_ca < 0, fprintf('显著负相关（温度越低，Ca mean越大）\n');
    else, fprintf('显著正相关\n'); end
else
    fprintf('无显著相关\n');
end

% ---- 归一化变异系数（CV_norm）与体温的相关性 ----
% 计算全局 F/F0 均值（用于归一化）
global_ca_mean = mean(ca_mean, 'omitnan');
if abs(global_ca_mean) < 0.01
    warning('全局 F/F0 均值接近 0，CV_norm 可能过大，请检查数据。');
end

% 计算每个窗口的 CV_norm
cv_norm = ca_std ./ abs(global_ca_mean);

% 相关性：CV_norm vs 平均温度
[rho_cv_norm, p_cv_norm] = corr(temp_mean, cv_norm, 'Type', 'Spearman');
fprintf('\n--- 归一化变异系数 (CV_norm) vs 平均温度 ---\n');
fprintf('斯皮尔曼相关系数 = %.3f, p = %.4f\n', rho_cv_norm, p_cv_norm);
if p_cv_norm < 0.05
    if rho_cv_norm < 0
        fprintf('显著负相关（温度越低，CV_norm 越大）\n');
    else
        fprintf('显著正相关（温度越高，CV_norm 越大）\n');
    end
else
    fprintf('无显著相关\n');
end

% 新增：窗口内 CV (cv_window) vs 平均温度
[rho_cv_win_temp, p_cv_win_temp] = corr(temp_mean, cv_window, 'Type', 'Spearman');
fprintf('\n--- 窗口内 CV (std/mean) vs 平均温度 ---\n');
fprintf('斯皮尔曼相关系数 = %.3f, p = %.4f\n', rho_cv_win_temp, p_cv_win_temp);
if p_cv_win_temp < 0.05
    if rho_cv_win_temp < 0
        fprintf('显著负相关（温度越低，CV 越大）\n');
    else
        fprintf('显著正相关（温度越高，CV 越大）\n');
    end
else
    fprintf('无显著相关\n');
end

%% 相关性分析（活动量相关，仅使用 loco_total）
if length(loco_total) > 5
    % 活动量总和 vs Ca均值
    [rho_loco_total_mean, p_loco_total_mean] = corr(loco_total, ca_mean, 'Type', 'Spearman');
    % 活动量总和 vs Ca波动 (cv_norm)
    [rho_loco_total_cvnorm, p_loco_total_cvnorm] = corr(loco_total, cv_norm, 'Type', 'Spearman');
    % 活动量总和 vs 窗口内 CV (cv_window)
    [rho_loco_total_cvwin, p_loco_total_cvwin] = corr(loco_total, cv_window, 'Type', 'Spearman');
    
    fprintf('\n--- 活动量与钙信号窗口相关性 ---\n');
    fprintf('活动量总和 vs Ca均值: rho = %.3f, p = %.4f\n', rho_loco_total_mean, p_loco_total_mean);
    fprintf('活动量总和 vs CV_norm: rho = %.3f, p = %.4f\n', rho_loco_total_cvnorm, p_loco_total_cvnorm);
    fprintf('活动量总和 vs 窗口内CV: rho = %.3f, p = %.4f\n', rho_loco_total_cvwin, p_loco_total_cvwin);
else
    fprintf('有效窗口数不足，跳过活动量相关性分析。\n');
end

%% 绘图：活动量 vs Ca波动 散点图（如果窗口数>5）
if length(loco_total) > 5
    figure('Name', 'Locomotion total vs CV_norm');
    scatter(loco_total, cv_norm, 30, 'filled', 'MarkerFaceAlpha', 0.5);
    xlabel('活动量总和 (窗口内总活动量)');
    ylabel('CV\_norm (std / |global mean|)');
    title('活动量总和 vs CV\_norm');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    coef = polyfit(loco_total, cv_norm, 1);
    x_fit = linspace(min(loco_total), max(loco_total), 100);
    y_fit = polyval(coef, x_fit);
    hold on; plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
    legend('数据点', '线性趋势', 'Location', 'best'); hold off;

    figure('Name', 'Locomotion total vs Ca mean');
    scatter(loco_total, ca_mean, 30, 'filled', 'MarkerFaceAlpha', 0.5);
    xlabel('活动量总和 (窗口内总活动量)');
    ylabel('Ca均值 (F/F0)');
    title('活动量总和 vs Ca均值 (F/F0)');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    coef2 = polyfit(loco_total, ca_mean, 1);
    x_fit2 = linspace(min(loco_total), max(loco_total), 100);
    y_fit2 = polyval(coef2, x_fit2);
    hold on; plot(x_fit2, y_fit2, 'k-', 'LineWidth', 2);
    legend('数据点', '线性趋势', 'Location', 'best'); hold off;
end

% 绘制活动量与钙信号随时间变化（双Y轴）
figure('Name', 'Calcium signal and locomotion over time');
yyaxis left;
plot(time / 60, F_F0, 'g-', 'LineWidth', 0.8);
ylabel('F/F0');
% ylim 自动调整，不再固定
yyaxis right;
loco_interp = interp1(temp_time, loco_raw, time, 'nearest', 0);
plot(time / 60, loco_interp, 'b-', 'LineWidth', 0.8);
ylabel('活动量 (activity)');
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

% 新增：分开的 F/F0 和活动总量随时间变化图（上下子图）
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
ylabel('活动总量 (activity)');
xlabel('Time (minutes)');
title('Locomotion total over 2 hours');
ylim([0, 50]);
grid off;
ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
xlim([0, 120]); xticks(0:30:120);

% 新增：活动总量与 CV_norm 的双Y轴图
t_win_min = t_starts / 60;

figure('Name', 'Loco total vs CV_norm over time');
yyaxis right;
plot(t_win_min, loco_total, 'b-', 'LineWidth', 1.5);
ylabel('活动总量 (窗口内总和)');
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
legend('活动总量', 'CV\_norm', 'Location', 'best');

%% ---- 绘制 CV_norm 与体温随时间变化的双Y轴图 ----
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

% 图2：平均温度 vs Ca CV (散点图 + 线性趋势线) 使用 cv_norm
figure('Name', 'Mean temp vs CV_norm');
scatter(temp_mean, cv_norm, 30, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('平均温度 (°C)');
ylabel('CV\_norm');
title('窗口内平均温度 vs CV\_norm');
grid off;
ax = gca;
ax.XAxis.TickDirection = 'out';
ax.YAxis.TickDirection = 'out';

coef3 = polyfit(temp_mean, cv_norm, 1);
x_fit3 = linspace(min(temp_mean), max(temp_mean), 100);
y_fit3 = polyval(coef3, x_fit3);
hold on;
plot(x_fit3, y_fit3, 'k-', 'LineWidth', 2);
legend('数据点', '线性趋势', 'Location', 'best');

% 新增：平均温度 vs 窗口内 CV (cv_window)
figure('Name', 'Mean temp vs window CV');
scatter(temp_mean, cv_window, 30, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('平均温度 (°C)');
ylabel('窗口内 CV (std/mean)');
title('窗口内平均温度 vs 窗口内 CV');
grid off;
ax = gca;
ax.XAxis.TickDirection = 'out';
ax.YAxis.TickDirection = 'out';
coef_cv = polyfit(temp_mean, cv_window, 1);
x_fit_cv = linspace(min(temp_mean), max(temp_mean), 100);
y_fit_cv = polyval(coef_cv, x_fit_cv);
hold on;
plot(x_fit_cv, y_fit_cv, 'k-', 'LineWidth', 2);
legend('数据点', '线性趋势', 'Location', 'best');

%% 绘图
% 图1：温度变化率 vs Ca波动 (std)
figure('Name', 'Slope vs Ca std');
scatter(temp_slope, ca_std, 10, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('温度变化率 (°C/s)');
ylabel('Ca波动 (std F/F0)');
title('窗口内温度变化率 vs Ca波动 (std)');
grid on;
coef = polyfit(temp_slope, ca_std, 1);
x_fit = linspace(min(temp_slope), max(temp_slope), 100);
y_fit = polyval(coef, x_fit);
hold on; plot(x_fit, y_fit, 'k--', 'LineWidth', 2);
legend('数据点', '线性趋势');

% 图2：平均温度 vs Ca波动 (std)
figure('Name', 'Mean temp vs Ca std');
scatter(temp_mean, ca_std, 30, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('平均温度 (°C)');
ylabel('Ca波动 (std F/F0)');
title('窗口内平均温度 vs Ca波动 (std)');
grid off;
ax = gca;
ax.XAxis.TickDirection = 'out';
ax.YAxis.TickDirection = 'out';
coef2 = polyfit(temp_mean, ca_std, 1);
x_fit2 = linspace(min(temp_mean), max(temp_mean), 100);
y_fit2 = polyval(coef2, x_fit2);
hold on;
plot(x_fit2, y_fit2, 'k-', 'LineWidth', 2);
legend('数据点', '线性趋势', 'Location', 'best');

% ===== 可视化 Ca 标准差和均值随时间变化 =====
t_min = t_starts / 60;

% 图1：Ca 均值随时间变化（绿色）与体温
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

% 图2：Ca 标准差（波动性）与温度（双Y轴）
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

% 新增：窗口内 CV 随时间变化与体温
figure('Name', 'Window CV and temperature over time');
yyaxis left;
plot(t_min, cv_window, 'Color', [0 0.5 0], 'LineWidth', 1.5);
ylabel('窗口内 CV (std/mean)');
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

% 创建新图：原始信号
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

% %% 导出窗口数据到CSV（可选，已注释）
% t_starts = t_starts(:);
% temp_mean = temp_mean(:);
% ca_mean = ca_mean(:);
% ca_std = ca_std(:);
% cv_window = cv_window(:);
% loco_total = loco_total(:);
% 
% 
% output_dir = 'F:\PhD Thesis\Fiber photometry\Cre line_Ca2+ GCamP\Processed_file-short-term\Output_Ca-CV_vs_Tb_loco';
% if ~exist(output_dir, 'dir')
%     mkdir(output_dir);
% end
% 
% T = table(t_starts, temp_mean, ca_mean, ca_std, cv_window, loco_total, ...
%     'VariableNames', {'WindowStart', 'TempMean', 'CaMean', 'CaStd', 'CV_window', 'LocoTotal'});
% T.AnimalID = repmat({animal_id}, height(T), 1);
% global_ca_mean = mean(ca_mean, 'omitnan');
% T.GlobalCaMean = repmat(global_ca_mean, height(T), 1);
% T.CV_norm = T.CaStd ./ abs(T.GlobalCaMean);
% 
% output_file = fullfile(output_dir, sprintf('%s_window_data.csv', animal_id));
% writetable(T, output_file);
% fprintf('窗口数据已保存至: %s\n', output_file);

%% 函数定义
function [F_F0, z_score] = calculate(curr_signal, fitted_reference)
    % 计算 F/F0
    F_F0 = curr_signal ./ fitted_reference;
    average_F_F0 = mean(F_F0);
    std_F_F0 = std(F_F0);
    z_score = (F_F0 - average_F_F0) ./ std_F_F0;
end

function [fit_coeffs, fitted_signal] = fitReferenceToSignal(curr_reference, fitting_signal, fitting_reference, order)
    fit_coeffs = polyfit(fitting_reference, fitting_signal, order);
    fitted_signal = polyval(fit_coeffs, curr_reference);
end