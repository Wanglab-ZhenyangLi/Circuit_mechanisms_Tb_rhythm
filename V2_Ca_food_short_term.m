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

%% ---- 读取进食事件 ----
eating_filename = 'F:\PhD Thesis\Fiber photometry\Cre line_Ca2+ GCamP\Processed_file-short-term\Food_intake_simul_photometry\food intake_recordings\corrected_J944-FED008_22-24_0721_extract_output.csv';
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
% 截取到与钙信号相同长度
max_ca_time = max(time);
idx_eat = eat_time <= max_ca_time;
eat_time = eat_time(idx_eat);
eat_raw = eat_raw(idx_eat);

% ---- 事件检测（0→1 跳变） ----
diff_eat = [0; diff(eat_raw)];
onset_idx = find(diff_eat == 1);
onset_time = eat_time(onset_idx);
fprintf('检测到 %d 次进食开始事件\n', length(onset_time));

%% ---- 合并连续事件为回合（bouts） ----
if length(onset_time) > 1
    bout_threshold = 60;   % 秒，间隔<=此值则合并为同一回合
    bout_start_times = [];
    bout_end_times = [];
    bout_event_counts = [];
    bout_idx = 1;
    current_bout_start = onset_time(1);
    current_bout_end = onset_time(1);
    current_bout_count = 1;
    
    for i = 2:length(onset_time)
        if onset_time(i) - onset_time(i-1) <= bout_threshold
            % 属于同一回合，更新结束时间和计数
            current_bout_end = onset_time(i);
            current_bout_count = current_bout_count + 1;
        else
            % 结束当前回合，记录
            bout_start_times(bout_idx) = current_bout_start;
            bout_end_times(bout_idx) = current_bout_end;
            bout_event_counts(bout_idx) = current_bout_count;
            bout_idx = bout_idx + 1;
            % 开始新回合
            current_bout_start = onset_time(i);
            current_bout_end = onset_time(i);
            current_bout_count = 1;
        end
    end
    % 记录最后一个回合
    bout_start_times(bout_idx) = current_bout_start;
    bout_end_times(bout_idx) = current_bout_end;
    bout_event_counts(bout_idx) = current_bout_count;
    
    % 计算每个回合的持续时间（秒）
    bout_durations = bout_end_times - bout_start_times;
    
    % 更新 onset_time 为回合开始时间
    onset_time = bout_start_times';
    
    fprintf('合并后回合数: %d (阈值=%.0f 秒)\n', length(onset_time), bout_threshold);
else
    % 只有一个事件，单独作为一个回合
    bout_start_times = onset_time;
    bout_end_times = onset_time;
    bout_event_counts = ones(size(onset_time));
    bout_durations = zeros(size(onset_time));
    fprintf('事件数不足2个，每个事件单独作为一个回合。\n');
end

%% ---- 输出回合长度统计 ----
fprintf('\n--- 进食回合长度统计 ---\n');
fprintf('回合总数: %d\n', length(onset_time));
if length(bout_durations) > 1
    fprintf('持续时间 (秒) 均值: %.2f, 中位数: %.2f, 标准差: %.2f\n', ...
        mean(bout_durations), median(bout_durations), std(bout_durations));
    fprintf('持续时间 四分位数 (25%%-50%%-75%%): %.2f - %.2f - %.2f\n', ...
        prctile(bout_durations,25), prctile(bout_durations,50), prctile(bout_durations,75));
    fprintf('持续时间 最小值: %.2f, 最大值: %.2f\n', min(bout_durations), max(bout_durations));
    fprintf('事件数 均值: %.2f, 中位数: %.2f, 标准差: %.2f\n', ...
        mean(bout_event_counts), median(bout_event_counts), std(bout_event_counts));
    fprintf('事件数 范围: %d - %d\n', min(bout_event_counts), max(bout_event_counts));
else
    fprintf('仅有一个回合，无法计算分布统计。\n');
end

% ---- 绘制回合长度分布（直方图） ----
if length(bout_durations) > 1
    figure('Name', 'Eating bout duration distribution');
    subplot(2,1,1);
    histogram(bout_durations, 'BinWidth', 10, 'FaceColor', [0.8 0.4 0.4], 'EdgeColor', 'none');
    xlabel('持续时间 (秒)');
    ylabel('频数');
    title('进食回合持续时间分布');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    
    subplot(2,1,2);
    histogram(bout_event_counts, 'BinWidth', 1, 'FaceColor', [0.4 0.6 0.8], 'EdgeColor', 'none');
    xlabel('事件数 (颗粒数)');
    ylabel('频数');
    title('进食回合事件数分布');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
end

% ---- 分析 bout 之间的间隔 ----
if length(bout_start_times) > 1
    bout_intervals = bout_start_times(2:end) - bout_end_times(1:end-1);
    fprintf('\n--- 进食 bout 间隔统计 ---\n');
    fprintf('间隔数: %d\n', length(bout_intervals));
    fprintf('间隔 (秒) 均值: %.2f, 中位数: %.2f, 标准差: %.2f\n', ...
        mean(bout_intervals), median(bout_intervals), std(bout_intervals));
    fprintf('间隔 四分位数 (25%%-50%%-75%%): %.2f - %.2f - %.2f\n', ...
        prctile(bout_intervals,25), prctile(bout_intervals,50), prctile(bout_intervals,75));
    fprintf('间隔 最小值: %.2f, 最大值: %.2f\n', min(bout_intervals), max(bout_intervals));
    
    % 绘制间隔分布直方图
    figure('Name', 'Eating bout interval distribution');
    histogram(bout_intervals, 'BinWidth', 30, 'FaceColor', [0.6 0.6 0.8], 'EdgeColor', 'none');
    xlabel('间隔 (秒)');
    ylabel('频数');
    title('进食回合间隔分布');
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
else
    fprintf('Bout 数不足2个，无法计算间隔。\n');
end


% ---- 重新计算事件索引（因为 onset_time 已更新为回合开始） ----
event_indices = zeros(length(onset_time),1);
for i = 1:length(onset_time)
    [~, idx] = min(abs(time - onset_time(i)));
    event_indices(i) = idx;
end

%% ---- 事件锁时参数 ----
win_before = 30;   % 秒
win_after  = 60;
t_uniform = -win_before : 0.1 : win_after;

% ---- 筛选事件窗口完整且远离边界 ----
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
fprintf('有效回合数（窗口完整且远离边界）: %d\n', n_events);

% ---- 对每个事件进行 reference fitting ----
event_dFF = zeros(n_events, length(t_uniform));

for i = 1:n_events
    t0 = time(event_indices(i));
    
    % 1. 拟合窗口：事件前 10 秒（不含 t0）
    fit_mask = (time >= t0 - win_before) & (time < t0);
    fit_idx = find(fit_mask);
    if length(fit_idx) < 5
        warning('事件 %d 拟合点数不足，跳过', i);
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

% 剔除失败事件（全为 NaN 的行）
valid_events = ~isnan(event_dFF(:,1));
event_dFF = event_dFF(valid_events, :);
onset_time = onset_time(valid_events);
n_events = size(event_dFF, 1);
fprintf('成功处理事件数: %d\n', n_events);

if n_events == 0
    error('无有效事件，请检查数据。');
end

%% ---- 检查前几个事件（绘图）并导出数据为 CSV ----
time_axis = t_uniform;
n_check = max(12, n_events);
if n_check > 0
    n_rows = ceil(n_check / 4);
    figure('Name', 'Check eating events with behavior (interpolated, reference fitted)');
    
    % 预分配用于导出的数据存储
    all_time = [];
    all_dFF = [];
    all_behavior = [];
    all_eventID = [];
    
    for i = 1:n_check
        % ---- 左轴：Ca 信号 ΔF/F ----
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
        
        % ---- 右轴：进食标记（0/1） ----
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
        
        % ---- 公共设置 ----
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
        
        % ---- 收集数据用于导出 ----
        all_time = [all_time; time_axis(:)];
        all_dFF = [all_dFF; event_dFF(i, :)'];
        all_behavior = [all_behavior; beh_interp(:)];
        all_eventID = [all_eventID; repmat(i, length(time_axis), 1)];
    end
    
    % ---- 导出为 CSV ----
    export_table = table(all_eventID, all_time, all_dFF, all_behavior, ...
        'VariableNames', {'EventID', 'TimeFromOnset_sec', 'dFF_percent', 'Behavior_01'});
    output_dir = "F:\PhD Thesis\Fiber photometry\Cre line_Ca2+ GCamP\Processed_file-short-term\Food_intake_simul_photometry\Sample traces";
    export_filename = fullfile(output_dir, 'J944_eating_events_export.csv');
    writetable(export_table, export_filename);
    fprintf('已导出 %d 个事件的数据到 CSV 文件: %s\n', n_check, export_filename);
end

%% ---- 平均响应与标准误 ----
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

% ---- 热图（按峰值排序） ----
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
    %caxis([-5, 10]);   % 注意：如果您的 MATLAB 版本不支持 clim，请使用 caxis
    grid off;
    ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    hold on;
    line([0, 0], [0.5, n_events+0.5], 'Color', 'w', 'LineWidth', 1.5, 'LineStyle', '--');
    hold off;
end

% %% ===== 窗口相关性分析：进食比例 vs Ca信号 =====
% % 使用已截取的 time, deltaF_F, eat_time, eat_raw
% win_sec = 180;      % 窗口长度（秒）
% step_sec = 60;     % 步长（秒）
% t_starts = 0:step_sec:(max(time)-win_sec);
% n_windows = length(t_starts);
% 
% eat_ratio = nan(n_windows, 1);
% ca_mean_win = nan(n_windows, 1);
% ca_std_win = nan(n_windows, 1);
% 
% for i = 1:n_windows
%     t0 = t_starts(i);8
%     t1 = t0 + win_sec;
%     
%     % 钙信号窗口
%     idx_ca = (time >= t0) & (time < t1);
%     if sum(idx_ca) < 10
%         continue;
%     end
%     ca_win = deltaF_F(idx_ca);
%     ca_mean_win(i) = mean(ca_win, 'omitnan');
%     ca_std_win(i) = std(ca_win, 'omitnan');
%     
%     % 进食数据窗口（使用 eat_time）
%     idx_eat = (eat_time >= t0) & (eat_time < t1);
%     if sum(idx_eat) < 1
%         continue;
%     end
%     eat_win = eat_raw(idx_eat);
%     eat_ratio(i) = sum(eat_win) / length(eat_win);   % 进食比例
% end
% 
% % 去除无效窗口
% valid = ~isnan(eat_ratio) & ~isnan(ca_mean_win) & ~isnan(ca_std_win);
% eat_ratio = eat_ratio(valid);
% ca_mean_win = ca_mean_win(valid);
% ca_std_win = ca_std_win(valid);
% fprintf('进食窗口分析有效窗口数 = %d\n', length(eat_ratio));
% 
% if length(eat_ratio) > 5
%     % 相关性分析
%     [rho_mean, p_mean] = corr(eat_ratio, ca_mean_win, 'Type', 'Spearman');
%     [rho_std, p_std] = corr(eat_ratio, ca_std_win, 'Type', 'Spearman');
%     fprintf('进食比例 vs Ca均值: rho = %.3f, p = %.4f\n', rho_mean, p_mean);
%     fprintf('进食比例 vs Ca波动: rho = %.3f, p = %.4f\n', rho_std, p_std);
%     
%     % 散点图：进食比例 vs Ca均值
%     figure('Name', 'Eating ratio vs Ca mean');
%     scatter(eat_ratio, ca_mean_win, 40, 'filled', 'MarkerFaceAlpha', 0.5);
%     xlabel('Eating ratio (window)');
%     ylabel('Ca mean (ΔF/F)');
%     title('Eating ratio vs Ca signal mean');
%     xlim([-0.05, 1.05]); ylim([-10, 20]);
%     grid off;
%     ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
%     coef = polyfit(eat_ratio, ca_mean_win, 1);
%     x_fit = linspace(0, 1, 100);
%     y_fit = polyval(coef, x_fit);
%     hold on; plot(x_fit, y_fit, 'k-', 'LineWidth', 2);
%     legend('Data', 'Linear trend', 'Location', 'best'); hold off;
%     
%     % 散点图：进食比例 vs Ca波动
%     figure('Name', 'Eating ratio vs Ca std');
%     scatter(eat_ratio, ca_std_win, 40, 'filled', 'MarkerFaceAlpha', 0.5);
%     xlabel('Eating ratio (window)');
%     ylabel('Ca std (ΔF/F)');
%     title('Eating ratio vs Ca fluctuation');
%     xlim([-0.05, 1.05]); ylim([-10, 20]);
%     grid off;
%     ax = gca; ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
%     coef2 = polyfit(eat_ratio, ca_std_win, 1);
%     y_fit2 = polyval(coef2, x_fit);
%     hold on; plot(x_fit, y_fit2, 'k-', 'LineWidth', 2);
%     legend('Data', 'Linear trend', 'Location', 'best'); hold off;
% else
%     fprintf('有效窗口太少，跳过窗口相关性分析\n');
% end


%% ---- 计算全局 ΔF/F₀（用于总览图） ----
percentile_fit = 50;   % 使用低于50%分位数的点作为基线
thresh_fit = prctile(signal, percentile_fit);
idx_low_fit = signal <= thresh_fit;
if sum(idx_low_fit) < 10
    warning('基线点不足，使用全部数据拟合');
    p_global = polyfit(reference, signal, 1);
else
    p_global = polyfit(reference(idx_low_fit), signal(idx_low_fit), 1);
end
fitted_ref_global = polyval(p_global, reference);
deltaF_F = (signal - fitted_ref_global) ./ fitted_ref_global * 100;
fprintf('全局 ΔF/F₀ 已计算。\n');

% ---- 绘制全时程双Y轴图 ----
% 该图展示整个两小时内的钙信号和进食事件
figure('Name', 'Calcium signal and eating events over time');

% 左轴：钙信号 ΔF/F
yyaxis left;
plot(time / 60, deltaF_F, 'g-', 'LineWidth', 0.8);
ylabel('ΔF/F (%)');
ylim([-10, 20]);  % 统一范围

% 右轴：进食标记（0/1）
yyaxis right;
% 将进食数据插值到钙信号的时间点（最近邻插值，保持0/1）
eat_interp = interp1(eat_time, eat_raw, time, 'nearest', 0);
% 用垂直线显示进食事件（在每个进食点画一条垂直线）
stem(time / 60, eat_interp, 'b-', 'LineWidth', 0.5, 'Marker', 'none');
ylabel('Eating (0/1)');
ylim([-0.1, 2.1]);   % 保持0/1范围

% 公共设置
xlabel('Time (minutes)');
title('Calcium signal and eating events over 2 hours');
xlim([0, 120]);
xticks(0:30:120);
grid off;
ax = gca;
ax.XAxis.TickDirection = 'out';
ax.YAxis(1).TickDirection = 'out';
ax.YAxis(2).TickDirection = 'out';
