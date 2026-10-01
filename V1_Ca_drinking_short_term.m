%% 
% LZY script_20260714 - 优化 reference fitting
clear all; clc;

% ---- 读取数据 ----
filename = 'F:\PhD Thesis\Fiber photometry\Cre line_Ca2+ GCamP\Processing_raw_data\Continuous-2-hr\20260721_220200.csv'; 
data = csvread(filename,1,0);

%% ---- 通道选择 ----
channel = 3;   % 请根据实际修改
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

% ---- 截取前两小时 ----
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
    
% ---- 参考通道异常值修复（同原代码） ----
threshold = 0.2 * median(reference);
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

%% ---- 读取饮水事件 ----
drink_file = 'F:\PhD Thesis\Fiber photometry\Cre line_Ca2+ GCamP\Processing_raw_data\Continuous-2-hr\20260721_220200-Event.csv';
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

% ---- 读取数据 ----
try
    if has_header_drink
        T_drink = readtable(drink_file);
        drink_raw = table2array(T_drink(:, col_drink));
    else
        drink_raw = readmatrix(drink_file);
        if size(drink_raw,2) >= col_drink
            drink_raw = drink_raw(:, col_drink);
        else
            error('列数不足');
        end
    end
    drink_raw = drink_raw(:);
catch ME
    error('饮水文件读取失败: %s', ME.message);
end

% 时间轴（500 Hz）
fs_drink = 500;
drink_time = (0:length(drink_raw)-1)' / fs_drink;

% 截取到与钙信号相同长度
max_ca_time = max(time);
idx_drink = drink_time <= max_ca_time;
drink_time = drink_time(idx_drink);
drink_raw = drink_raw(idx_drink);
fprintf('饮水数据点数 = %d, 时间范围 %.2f 秒\n', length(drink_raw), drink_time(end));

% ---- 检测 lick 事件（0→1 跳变） ----
diff_drink = [0; diff(drink_raw)];
lick_onset_idx = find(diff_drink == 1);
lick_times = drink_time(lick_onset_idx);
fprintf('原始 lick 事件数: %d\n', length(lick_times));

if isempty(lick_times)
    fprintf('未检测到 lick 事件，跳过饮水分析\n');
    return;
end

%% ---- 分析 lick 间隔分布（帮助选择合并阈值） ----
if length(lick_times) > 1
    intervals = diff(lick_times);
    fprintf('\n--- Lick 间隔统计 ---\n');
    fprintf('间隔均值: %.3f 秒\n', mean(intervals));
    fprintf('间隔中位数: %.3f 秒\n', median(intervals));
    fprintf('间隔标准差: %.3f 秒\n', std(intervals));
    fprintf('间隔四分位数 (25%%-50%%-75%%): %.3f - %.3f - %.3f 秒\n', ...
        prctile(intervals,25), prctile(intervals,50), prctile(intervals,75));
    fprintf('最小间隔: %.3f 秒\n', min(intervals));
    fprintf('最大间隔: %.3f 秒\n', max(intervals));
    
    % 绘制间隔分布
    figure('Name', 'Lick interval distribution');
    subplot(2,1,1);
    histogram(intervals, 'BinWidth', 0.02, 'FaceColor', [0.6 0.6 0.8], 'EdgeColor', 'none');
    hold on;
    [f, xi] = ksdensity(intervals);
    plot(xi, f * length(intervals) * 0.02, 'k-', 'LineWidth', 2);
    xlabel('间隔 (秒)');
    ylabel('频数');
    title('Lick 间隔分布');
    xline(median(intervals), 'r--', sprintf('中位数=%.2fs', median(intervals)), 'LineWidth', 1.5);
    xline(mean(intervals), 'g--', sprintf('均值=%.2fs', mean(intervals)), 'LineWidth', 1.5);
    legend('频数', '核密度', '中位数', '均值');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    subplot(2,1,2);
    boxplot(intervals, 'Orientation', 'horizontal', 'Symbol', '+');
    xlabel('间隔 (秒)');
    title('间隔箱线图');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    
    fprintf('\n--- 阈值建议 ---\n');
    fprintf('建议阈值1 (0.5秒, 常用): 0.50 秒\n');
    fprintf('建议阈值2 (2倍中位数): %.2f 秒\n', 2*median(intervals));
    fprintf('建议阈值3 (75%%分位数): %.2f 秒\n', prctile(intervals,75));
    fprintf('请根据直方图谷底选择合适阈值。\n');
end

%% ---- 将 lick 合并为 bout（阈值 15 秒，可根据上述分析调整） ----
bout_gap = 15;   % 秒
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
    fprintf('合并为 %d 个 lick bout (间隔阈值 %.2f 秒)\n', n_bouts, bout_gap);
else
    n_bouts = 1;
    bout_onset = lick_times;
    bout_lick_count = 1;
    bout_duration = 0;
end

%% ---- 输出 bout 长度统计 ----
fprintf('\n--- 饮水 bout 长度统计 ---\n');
fprintf('Bout 总数: %d\n', n_bouts);
if n_bouts > 1
    fprintf('持续时间 (秒) 均值: %.3f, 中位数: %.3f, 标准差: %.3f\n', ...
        mean(bout_duration), median(bout_duration), std(bout_duration));
    fprintf('持续时间 四分位数 (25%%-50%%-75%%): %.3f - %.3f - %.3f\n', ...
        prctile(bout_duration,25), prctile(bout_duration,50), prctile(bout_duration,75));
    fprintf('持续时间 最小值: %.3f, 最大值: %.3f\n', min(bout_duration), max(bout_duration));
    fprintf('Lick 数 均值: %.2f, 中位数: %.2f, 标准差: %.2f\n', ...
        mean(bout_lick_count), median(bout_lick_count), std(bout_lick_count));
    fprintf('Lick 数 范围: %d - %d\n', min(bout_lick_count), max(bout_lick_count));
    
    % 绘制分布图
    figure('Name', 'Drinking bout duration and lick count');
    subplot(2,1,1);
    histogram(bout_duration, 'BinWidth', 0.1, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
    xlabel('持续时间 (秒)');
    ylabel('频数');
    title('饮水 bout 持续时间分布');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    
    subplot(2,1,2);
    histogram(bout_lick_count, 'BinWidth', 1, 'FaceColor', [0.8 0.4 0.6], 'EdgeColor', 'none');
    xlabel('Lick 数');
    ylabel('频数');
    title('饮水 bout Lick 数分布');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
else
    fprintf('仅有一个 bout，无法计算分布统计。\n');
end

% ---- 分析 bout 之间的间隔 ----
if n_bouts > 1
    bout_end = bout_onset + bout_duration;
    bout_intervals = bout_onset(2:end) - bout_end(1:end-1);
    fprintf('\n--- 饮水 bout 间隔统计 ---\n');
    fprintf('间隔数: %d\n', length(bout_intervals));
    fprintf('间隔 (秒) 均值: %.3f, 中位数: %.3f, 标准差: %.3f\n', ...
        mean(bout_intervals), median(bout_intervals), std(bout_intervals));
    fprintf('间隔 四分位数 (25%%-50%%-75%%): %.3f - %.3f - %.3f\n', ...
        prctile(bout_intervals,25), prctile(bout_intervals,50), prctile(bout_intervals,75));
    fprintf('间隔 最小值: %.3f, 最大值: %.3f\n', min(bout_intervals), max(bout_intervals));
    
    % 绘制间隔分布直方图
    figure('Name', 'Drinking bout interval distribution');
    histogram(bout_intervals, 'BinWidth', 1, 'FaceColor', [0.6 0.6 0.8], 'EdgeColor', 'none');
    xlabel('间隔 (秒)');
    %xlim([0,30]);
    ylabel('频数');
    title('饮水回合间隔分布');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
else
    fprintf('Bout 数不足2个，无法计算间隔。\n');
end
%% ---- 事件锁时参数（统一时间轴） ----
win_before = 30;   % 秒（前10秒用于拟合）
win_after  = 60;   % 秒（后10秒）
t_uniform = -win_before : 0.1 : win_after;   % 0.1秒分辨率

% ---- 在钙信号时间轴上定位事件 ----
event_indices = zeros(n_bouts, 1);
for i = 1:n_bouts
    [~, idx] = min(abs(time - bout_onset(i)));
    event_indices(i) = idx;
end

% ---- 筛选事件窗口完整且远离边界 ----
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
fprintf('有效 bout 数（窗口完整且远离边界）: %d\n', n_events);

if n_events == 0
    error('无有效 bout，请检查数据或调整窗口参数。');
end

%% ---- 对每个 bout 独立进行 reference fitting ----
event_dFF = zeros(n_events, length(t_uniform));

for i = 1:n_events
    t0 = time(event_indices(i));
    
    % 1. 拟合窗口：事件前 10 秒（不含 t0）
    fit_mask = (time >= t0 - win_before) & (time < t0);
    fit_idx = find(fit_mask);
    if length(fit_idx) < 5
        warning('Bout %d 拟合点数不足，跳过', i);
        continue;
    end
    sig_fit = signal(fit_idx);
    ref_fit = reference(fit_idx);
    
    % 2. 线性拟合（signal = a*ref + b）
    p = polyfit(ref_fit, sig_fit, 1);
    
    % 3. 整个事件窗口（前10～后10秒）
    evt_mask = (time >= t0 - win_before) & (time <= t0 + win_after);
    evt_idx = find(evt_mask);
    t_rel = time(evt_idx) - t0;
    sig_evt = signal(evt_idx);
    ref_evt = reference(evt_idx);
    
    % 4. 使用拟合参数预测 fitted_reference
    fitted_ref = polyval(p, ref_evt);
    
    % 5. 计算 ΔF/F (%) = (signal - fitted_ref) / fitted_ref * 100
    dFF = (sig_evt - fitted_ref) ./ fitted_ref * 100;
    
    % 6. 插值到统一时间轴
    dFF_uniform = interp1(t_rel, dFF, t_uniform, 'linear', 'extrap');
    if any(isnan(dFF_uniform))
        dFF_uniform = fillmissing(dFF_uniform, 'nearest');
    end
    
    event_dFF(i, :) = dFF_uniform;
end

% 剔除失败事件（全为NaN）
valid_events = ~isnan(event_dFF(:,1));
event_dFF = event_dFF(valid_events, :);
bout_onset = bout_onset(valid_events);
n_events = size(event_dFF, 1);
fprintf('成功处理 bout 数: %d\n', n_events);

if n_events == 0
    error('无成功处理的事件，请检查数据。');
end

%% ---- 检查前10个事件（每行4张）并导出数据为 CSV ----
time_axis = t_uniform;
n_check = max(3, n_events);
if n_check > 0
    n_rows = ceil(n_check / 4);
    figure('Name', 'Check drinking bout events with behavior (interpolated)');
    
    % 预分配用于导出的数据存储
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
        
        % ---- 收集数据用于导出 ----
        all_time = [all_time; time_axis(:)];
        all_dFF = [all_dFF; event_dFF(i, :)'];
        all_behavior = [all_behavior; beh_interp(:)];
        all_eventID = [all_eventID; repmat(i, length(time_axis), 1)];
    end
    
    % ---- 导出为 CSV ----
    export_table = table(all_eventID, all_time, all_dFF, all_behavior, ...
        'VariableNames', {'EventID', 'TimeFromOnset_sec', 'dFF_percent', 'Behavior_01'});
    output_dir = "F:\PhD Thesis\Fiber photometry\Cre line_Ca2+ GCamP\Processed_file-short-term\Water_intake_simul_photometry\Sample_traces";
    export_filename = fullfile(output_dir, 'J944_licking_events_export.csv');
    writetable(export_table, export_filename);
    fprintf('已导出 %d 个事件的数据到 CSV 文件: %s\n', n_check, export_filename);
end

%% ---- 平均响应 + SEM ----
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

%% ---- 热图（按峰值排序） ----
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

%% ---- 双Y轴图：全时程钙信号 + 饮水事件 ----
% 先计算全局 ΔF/F₀（如果未计算）
if ~exist('deltaF_F_global', 'var')
    fprintf('计算全局 ΔF/F₀ 用于总览图...\n');
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

% 插值饮水事件到钙信号时间轴（最近邻）
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