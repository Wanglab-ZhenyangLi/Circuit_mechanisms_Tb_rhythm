% 清除环境
clear; close all; clc;

% 读取数据，第一行作为列名
data = readtable("F:\PhD Thesis\Fluoresecent images\Nikon\Prork2-iCre-Vglut2-Flpo_based\CV and correlation\4-Night_Combined_E800-E883-E884-E449-Time Trace.csv");

% 获取列名（细胞名称）
cell_names = data.Properties.VariableNames;

% 转换为数值矩阵
calcium_data = table2array(data);

% 只选取前150行数据
if size(calcium_data, 1) < 150
    warning('数据只有%d行，使用所有可用数据', size(calcium_data, 1));
    rows_to_use = size(calcium_data, 1);
else
    rows_to_use = 150;
end

calcium_data = calcium_data(1:rows_to_use, :);
%% 
%拟合校正bleach
% 获取数据基本信息
[num_frames, num_cells] = size(calcium_data);
fprintf('数据信息：%d个时间点，%d个细胞\n', num_frames, num_cells);
fprintf('时间范围：%.1f - %.1f 秒（假设帧率1Hz）\n\n', 0, num_frames-1);

%% 2. 数据预处理和可视化
% 创建时间轴（假设帧率为1Hz）
time_axis = (0:num_frames-1)';

% 绘制原始数据
figure('Position', [100, 100, 1200, 800]);

% 原始信号（子图1）
subplot(3, 3, [1, 2, 3]);
plot(time_axis, calcium_data, 'LineWidth', 1);
title('原始钙信号（所有细胞）', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('时间 (秒)', 'FontSize', 11);
ylabel('荧光强度 (a.u.)', 'FontSize', 11);
grid on;

% 前5个细胞的放大视图（子图2）
subplot(3, 3, [4, 5]);
num_to_show = min(5, num_cells);
colors = lines(num_to_show);
for i = 1:num_to_show
    plot(time_axis, calcium_data(:, i), 'Color', colors(i, :), 'LineWidth', 1.5);
    hold on;
end
title('前5个细胞的原始信号（放大）', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('时间 (秒)', 'FontSize', 10);
ylabel('荧光强度 (a.u.)', 'FontSize', 10);
legend(cell_names(1:num_to_show), 'Location', 'best', 'FontSize', 8);
grid on;

%% 3. 核心：为每个ROI进行光漂白校正
fprintf('开始光漂白校正...\n');

% 初始化结果存储
corrected_data = zeros(size(calcium_data));  % 除法校正后数据
baseline_traces = zeros(size(calcium_data)); % 拟合的漂白曲线
fitted_params = cell(num_cells, 1);          % 拟合参数
correction_factors = zeros(size(calcium_data)); % 校正因子

% 为每个细胞单独进行校正
for cell_idx = 1:num_cells
    fprintf('  校正细胞 %d/%d: %s\n', cell_idx, num_cells, cell_names{cell_idx});
    
    % 1. 获取信号 ------------------------------------------------
    raw_signal = calcium_data(:, cell_idx);
    
    % 处理缺失值
    if any(isnan(raw_signal))
        raw_signal = fillmissing(raw_signal, 'linear');
        fprintf('      -> 注: 信号中存在NaN，已线性插值填充\n');
    end
    
    fit_signal = raw_signal; % 使用原始信号，不做平滑
    
    % 2. 基线选择 ------------------------------------------------
    % 2.1 强制初始段 (前10%)
    force_init_len = max(10, round(0.1 * num_frames));
    force_init_idx = 1:force_init_len;
    
    % 2.2 智能选择最平稳段
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
    
    % 2.3 合并区域
    baseline_idx = union(force_init_idx, stable_idx);
    
    fprintf('      -> 基线合并选择:\n');
    fprintf('         强制初始段: 帧 1 到 %d\n', force_init_len);
    if ~isempty(stable_idx)
        fprintf('         智能平稳段: 帧 %d 到 %d (局部方差=%.4f)\n', ...
                stable_idx(1), stable_idx(end), min_var_mean);
    end
    fprintf('         合并后总基线帧数: %d\n', length(baseline_idx));
    
    % 指数衰减模型
    exp_model = @(p, t) p(1) * exp(-p(2) * t) + p(3);
    
    % 初始参数猜测
    baseline_signal = fit_signal(baseline_idx);
    baseline_time = time_axis(baseline_idx) - time_axis(baseline_idx(1));
    
    init_c = quantile(baseline_signal, 0.25);
    
    if mean(baseline_signal(1:round(end/3))) > mean(baseline_signal(round(2*end/3):end))
        init_a = mean(baseline_signal(1:round(end/3))) - init_c;
    else
        init_a = (max(baseline_signal) - min(baseline_signal)) * 0.5;
    end
    
    % 基于强制初始段估计衰减率
    init_segment = fit_signal(force_init_idx);
    init_time_segment = time_axis(force_init_idx) - time_axis(1);
    
    if length(init_segment) >= 5
        p = polyfit(init_time_segment, init_segment, 1);
        slope = p(1);
        
        if slope < 0 && init_segment(1) > 0
            init_b = abs(slope) / (init_segment(1) - init_c + eps);
            init_b = min(max(init_b, 0.001), 0.2);
            fprintf('      -> 基于初始段斜率估计衰减率: b=%.4f\n', init_b);
        else
            init_b = 0.01;
        end
    else
        init_b = 0.02;
    end
    
    fprintf('      -> 动态初始值: a=%.2f, b=%.4f, c=%.2f\n', init_a, init_b, init_c);
    
    % 参数边界
    lb = [0, 1e-5, 0.5*init_c];
    ub = [3*init_a, 0.5, max(fit_signal)*1.5];
    
    options = optimset('Display', 'off', 'TolFun', 1e-8, 'TolX', 1e-8, ...
                       'MaxIter', 2000, 'MaxFunEvals', 3000);
    
    % 第一层次：全段数据指数拟合 --------------------------------------
    fprintf('      -> 尝试全段数据指数拟合\n');
    try
        [params_full, resnorm_full, ~, exitflag_full] = lsqcurvefit(exp_model, ...
            [init_a, init_b, init_c], ...
            time_axis - time_axis(1), ...
            fit_signal, lb, ub, options);
        
        % 计算全段拟合质量指标
        bleach_curve_full = exp_model(params_full, time_axis - time_axis(1));
        R2_full = 1 - (resnorm_full / sum((fit_signal - mean(fit_signal)).^2));
        global_ratio_range_full = max(bleach_curve_full ./ fit_signal) - min(bleach_curve_full ./ fit_signal);
        
        % 判断条件1和3：R2<0.8或全局比率范围>1.0
        if R2_full >= 0.8 && global_ratio_range_full <= 1.0
            % 全段拟合成功
            fprintf('      -> 全段指数拟合成功 (R^2=%.3f, 比率范围=%.2f)\n', R2_full, global_ratio_range_full);
            bleach_curve = bleach_curve_full;
            fitted_params{cell_idx}.method = 'exponential_full';
            fitted_params{cell_idx}.params = params_full;
            fitted_params{cell_idx}.R2 = R2_full;
            fitted_params{cell_idx}.exitflag = exitflag_full;
            
            % 执行校正
            baseline_traces(:, cell_idx) = bleach_curve;
            correction_factors(:, cell_idx) = bleach_curve(1) ./ bleach_curve;
            corrected_data(:, cell_idx) = raw_signal .* correction_factors(:, cell_idx);
            
            continue; % 跳过后续步骤，处理下一个细胞
        else
            fprintf('      -> 全段指数拟合不佳 (R^2=%.2f, 比率范围=%.2f)\n', R2_full, global_ratio_range_full);
            % 继续尝试第二层次
        end
    catch ME
        fprintf('      -> 全段指数拟合失败: %s\n', ME.message);
        % 继续尝试第二层次
    end
    
   % 第二层次：初始段+10%分位数点指数拟合 -------------------------------
    fprintf('      -> 尝试初始段+10%%分位数点指数拟合\n');
    
    % 计算10%分位数点（下包络）
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
        lower_envelope(t) = quantile(window_data, 0.10); % 10%分位数
    end
    
    % 选择初始段和下包络线点作为拟合数据
    combined_time = [time_axis(force_init_idx); time_axis];
    combined_signal = [fit_signal(force_init_idx); lower_envelope];
    
    % 为第二层次重新估计初始参数
    init_segment_2 = fit_signal(force_init_idx);
    init_time_2 = time_axis(force_init_idx) - time_axis(1);
    
    % 重新估计参数
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
    
    fprintf('      -> 第二层次初始值: a=%.2f, b=%.4f, c=%.2f\n', init_a_2, init_b_2, init_c_2);
    
    % === 关键修复：清除可能的遗留标志 ===
    use_linear_fit = false; % 先明确设置为false
    
    try
        [params_combined, resnorm_combined, ~, exitflag_combined] = lsqcurvefit(exp_model, ...
            [init_a_2, init_b_2, init_c_2], ...
            combined_time - combined_time(1), ...
            combined_signal, lb, ub, options);
        
        % 计算拟合质量
        bleach_curve_combined = exp_model(params_combined, time_axis - time_axis(1));
        R2_combined = 1 - (resnorm_combined / sum((combined_signal - mean(combined_signal)).^2));
        global_ratio_range_combined = max(bleach_curve_combined ./ fit_signal) - min(bleach_curve_combined ./ fit_signal);
        
        % === 调试：打印实际判断值 ===
        fprintf('      [调试] 实际判断值: R2_combined=%.3f, 比率范围=%.3f\n', R2_combined, global_ratio_range_combined);
        
        % 判断条件1和3：R2<0.7或全局比率范围>1.0
        if R2_combined >= 0.7 && global_ratio_range_combined <= 1.0
            % 第二层次拟合成功
            fprintf('      -> 初始段+10%%分位数点指数拟合成功 (R^2=%.3f, 比率范围=%.2f)\n', R2_combined, global_ratio_range_combined);
            bleach_curve = bleach_curve_combined;
            fitted_params{cell_idx}.method = 'exponential_combined';
            fitted_params{cell_idx}.params = params_combined;
            fitted_params{cell_idx}.R2 = R2_combined;
            fitted_params{cell_idx}.exitflag = exitflag_combined;
            
            % === 关键修复：明确设置，防止进入第三层次 ===
            use_linear_fit = false;
            
            % 执行校正并跳过后续层次
            baseline_traces(:, cell_idx) = bleach_curve;
            correction_factors(:, cell_idx) = bleach_curve(1) ./ bleach_curve;
            corrected_data(:, cell_idx) = raw_signal .* correction_factors(:, cell_idx);
            
            continue; % === 关键：直接跳到下一个细胞，不再执行第三层次 ===
            
        else
            % 第二层次拟合不佳
            fprintf('      -> 初始段+10%%分位数点指数拟合不佳 (R^2=%.2f, 比率范围=%.2f)\n', R2_combined, global_ratio_range_combined);
            use_linear_fit = true; % 允许进入第三层次
        end
        
    catch ME
        % 指数拟合过程失败
        fprintf('      -> 初始段+10%%分位数点指数拟合失败: %s\n', ME.message);
        use_linear_fit = true; % 允许进入第三层次
    end
    % === 第二层次结束 ===
    
    % 第三层次：10%分位数线性拟合 ----------------------------------------
    if exist('use_linear_fit', 'var') && use_linear_fit
        fprintf('      -> 执行备用方案：基于10%%分位数点的线性拟合\n');
        % 计算5%分位数点（更低的分位点，更能排除钙瞬变）
        envelope_window = max(5, min(31, round(num_frames * 0.07)));
        if mod(envelope_window, 2) == 0
            envelope_window = envelope_window + 1;
        end
        
        lower_envelope = zeros(num_frames, 1);
        half_win = floor(envelope_window / 2);
        
        % 计算滚动5%分位数，构建更保守的"下包络"曲线
        for t = 1:num_frames
            start_idx = max(1, t - half_win);
            end_idx = min(num_frames, t + half_win);
            window_data = fit_signal(start_idx:end_idx);
            lower_envelope(t) = quantile(window_data, 0.01); % 改为5%分位数
        end
        % 准备拟合数据：使用下包络线上的点
        X = time_axis - time_axis(1);
        Y = lower_envelope;
        
        % 使用 polyfit 进行线性拟合，无条件接受其结果
        % p_simple(1) 是斜率，p_simple(2) 是截距
        [p_simple, S] = polyfit(X, Y, 1);
        
        % 计算拟合值
        [Y_fit_simple, delta] = polyval(p_simple, X, S);
        
        % 计算 R²
        residuals = Y - Y_fit_simple;
        ss_res_simple = sum(residuals .^ 2);
        ss_tot_simple = sum((Y - mean(Y)) .^ 2);
        
        if ss_tot_simple > eps
            R2_simple = 1 - (ss_res_simple / ss_tot_simple);
        else
            R2_simple = 0; % 数据无变化，R²设为0
        end
        
        final_slope = p_simple(1);
        final_intercept = p_simple(2);
        bleach_curve = Y_fit_simple;
        
        % 报告拟合结果
        fprintf('         线性拟合完成 (斜率=%.2e, 截距=%.2f, R²=%.3f)\n', ...
                final_slope, final_intercept, R2_simple);
        
        % 存储结果
        fitted_params{cell_idx}.method = 'linear_10percentile_polyfit';
        fitted_params{cell_idx}.slope = final_slope;
        fitted_params{cell_idx}.intercept = final_intercept;
        fitted_params{cell_idx}.R2_linear = R2_simple;
        
        % 可选：轻微平滑（根据您的"不做平滑"要求，此行已被注释）
        % bleach_curve = smoothdata(bleach_curve, 'movmean', 3);
        
    end
    % 结束第三层次
        
    % 执行除法校正 ----------------------------------------------------
    baseline_traces(:, cell_idx) = bleach_curve;
    correction_factors(:, cell_idx) = bleach_curve(1) ./ bleach_curve;
    corrected_data(:, cell_idx) = raw_signal .* correction_factors(:, cell_idx);
    
end
fprintf('漂白校正完成！\n');
%%
F_over_F0 = zeros(size(calcium_data));
for cell_idx = 1:num_cells
    raw = calcium_data(:, cell_idx);
    baseline = baseline_traces(:, cell_idx);
    
    % 计算F/F₀ = raw / baseline
    F_over_F0(:, cell_idx) = raw ./ baseline;
end
% 按平均F/F₀从大到小排序
mean_F_over_F0 = mean(F_over_F0, 1);  % 计算每列(细胞)的平均值
[~, sort_idx] = sort(mean_F_over_F0, 'descend');  % 获取排序索引
F_over_F0_sorted = F_over_F0(:, sort_idx);  % 按排序索引重排数据

%% 5. 钙信号校正因子热图 - Z-score标准化 (Science期刊标准)
% 创建新图窗，设置Science期刊推荐尺寸
fig_width_cm = 17.8; % Science双栏宽度
fig_height_cm = 10; % 适当高度
fig_width = fig_width_cm / 2.54 * 96; % 转换为像素
fig_height = fig_height_cm / 2.54 * 96;

fig_heatmap = figure('Position', [100, 100, fig_width, fig_height], ...
    'Color', 'white', ...
    'Units', 'inches', ...
    'PaperUnits', 'inches', ...
    'PaperSize', [fig_width_cm/2.54, fig_height_cm/2.54]);

% 准备热图数据 - 每行是一个细胞
heatmap_data = deltaF_over_F0';

% Z-score标准化：每行数据减去其均值，除以其标准差
heatmap_data_zscore = zeros(size(heatmap_data));
row_mean = zeros(size(heatmap_data, 1), 1);
row_std = zeros(size(heatmap_data, 1), 1);

for i = 1:size(heatmap_data, 1)
    row_data = heatmap_data(i, :);
    row_mean(i) = mean(row_data);
    row_std(i) = std(row_data);
    
    if row_std(i) > 1e-10  % 避免除以零
        heatmap_data_zscore(i, :) = (row_data - row_mean(i)) / row_std(i);
    else
        heatmap_data_zscore(i, :) = zeros(size(row_data)); % 如果标准差为零，则设为0
    end
end

% 按原始平均校正因子排序（从大到小）
mean_correction = mean(heatmap_data, 2); % 使用原始数据的平均值排序
[mean_correction_sorted, sort_idx] = sort(mean_correction, 'descend');

% 按排序结果重新排列数据
heatmap_data_sorted = heatmap_data_zscore(sort_idx, :);
row_mean_sorted = row_mean(sort_idx);
row_std_sorted = row_std(sort_idx);
cell_indices_sorted = 1:num_cells; % 原始细胞索引
cell_indices_sorted = cell_indices_sorted(sort_idx); % 排序后的索引

% 创建热图 - 使用适合钙信号的颜色映射
% 对于钙信号，我们通常使用暖色调（红、橙、黄）来突出钙瞬变
% 选项1: 'hot' - 从黑到红到黄到白，非常适合钙信号
% 选项2: 'parula' - MATLAB默认，对钙信号也适用
% 选项3: 'inferno' 或 'plasma' - 从黑到亮黄，非常突出瞬变
% 选项4: 自定义红色调色板

% 我推荐使用'inferno'或'plasma'，因为它们能很好地突出高值
if exist('inferno', 'file') || exist('inferno.m', 'file')
    cmap = inferno(256); % inferno从黑到亮黄，非常突出
elseif exist('plasma', 'file') || exist('plasma.m', 'file')
    cmap = plasma(256); % plasma从紫到黄，也非常好
else
    % 如果没有这些颜色映射，使用hot或创建自定义
    cmap = hot(256); % hot是MATLAB自带的，从黑到红到黄到白
end

% 或者创建自定义的红色调色板，专门为钙信号优化
% cmap = create_calcium_colormap();

imagesc(time_axis, 1:size(heatmap_data_sorted, 1), heatmap_data_sorted);
colormap(cmap);
cbar = colorbar;

% Science期刊要求的字体设置
set(gca, 'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'normal');
set(cbar, 'FontName', 'Arial', 'FontSize', 8);

% 设置坐标轴标签
xlabel('时间 (秒)', 'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'bold');
ylabel('细胞编号', 'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'bold');
cbar.Label.String = 'Z-score';
cbar.Label.FontName = 'Arial';
cbar.Label.FontSize = 9;
cbar.Label.FontWeight = 'bold';

% 设置颜色范围：钙信号通常关注正值（钙瞬变）
% 对于Z-score，负值表示低于基线，可能不是主要关注点
z_min = min(heatmap_data_sorted(:));
z_max = max(heatmap_data_sorted(:));

% 为了突出钙瞬变，我们可以将颜色范围重点放在正值
% 但保持一定的负值范围以显示基线波动
if z_max > 0
    % 如果数据有正值，将颜色范围设置为[-1, max(3, z_max)]
    color_min = max(-1, min(-0.5, z_min)); % 负值范围较小
    color_max = max(3, z_max); % 确保能看到高Z-score的钙瞬变
    caxis([color_min, color_max]);
    
    % 设置颜色条刻度
    if color_max <= 5
        cbar_ticks = [color_min, 0, 1, 2, 3, color_max];
        cbar_ticklabels = arrayfun(@(x) sprintf('%.1f', x), cbar_ticks, 'UniformOutput', false);
    else
        cbar_ticks = [color_min, 0, 1, 2, 3, 5, color_max];
        cbar_ticklabels = {sprintf('%.1f', color_min), '0', '1', '2', '3', '5', sprintf('%.1f', color_max)};
    end
else
    % 如果没有正值，使用全范围
    caxis([z_min, z_max]);
    cbar_ticks = linspace(z_min, z_max, 5);
    cbar_ticklabels = arrayfun(@(x) sprintf('%.1f', x), cbar_ticks, 'UniformOutput', false);
end

% 设置颜色条刻度
cbar.Ticks = cbar_ticks;
cbar.TickLabels = cbar_ticklabels;

% 设置X轴显示
if length(time_axis) <= 20
    xticks(time_axis);
else
    % 显示6个时间刻度
    num_ticks = min(6, length(time_axis));
    tick_indices = round(linspace(1, length(time_axis), num_ticks));
    xticks(time_axis(tick_indices));
    xticklabels(arrayfun(@(x) sprintf('%.0f', x), time_axis(tick_indices), 'UniformOutput', false));
end

% 设置Y轴显示
num_cells_total = size(heatmap_data_sorted, 1);
if num_cells_total <= 30
    % 细胞数≤30时显示所有细胞标签
    yticks(1:num_cells_total);
    yticklabels(cell_indices_sorted);
elseif num_cells_total <= 60
    % 每5个细胞显示一个标签
    ytick_interval = 5;
    yticks(1:ytick_interval:num_cells_total);
    yticklabels(cell_indices_sorted(1:ytick_interval:end));
else
    % 每10个细胞显示一个标签
    ytick_interval = 10;
    yticks(1:ytick_interval:num_cells_total);
    yticklabels(cell_indices_sorted(1:ytick_interval:end));
end

% 添加标题
title('钙信号校正因子Z-score热图 (按平均校正因子降序排列)', ...
    'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');

% 添加网格线
grid on;
set(gca, 'GridColor', [0.3, 0.3, 0.3], 'GridAlpha', 0.1, 'GridLineStyle', ':');

% 美化图形
set(gca, 'Box', 'on', 'LineWidth', 0.5);
set(gca, 'TickDir', 'out');
set(gca, 'TickLength', [0.008, 0.008]);

% 使用矢量图形渲染器
set(gcf, 'Renderer', 'painters');
%% 4. 结果可视化
figure('Position', [200, 100, 1400, 900]);

% 示例展示前4个细胞的详细分析
num_examples = min(12, num_cells);
for i = 1:num_examples
    subplot(4, 4, i);
    
    % 原始信号 vs 拟合基线
    plot(time_axis, calcium_data(:, i+24), 'b-', 'LineWidth', 1.5, 'DisplayName', '原始信号');
    hold on;
    plot(time_axis, baseline_traces(:, i+24), 'r--', 'LineWidth', 2, 'DisplayName', '拟合基线');
    
    title(sprintf('细胞 %d: %s', i, cell_names{i+24}), 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('时间 (秒)', 'FontSize', 9);
    ylabel('荧光强度', 'FontSize', 9);
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    
    % 显示拟合参数
    if length(fitted_params{i}) >= 3
        if fitted_params{i}(3) == 0  % 线性拟合
            text(0.05, 0.15, sprintf('线性校正\n斜率: %.3e', fitted_params{i}(1)), ...
                'Units', 'normalized', 'FontSize', 8, 'BackgroundColor', 'white');
        else  % 指数拟合
            text(0.05, 0.15, sprintf('指数校正\n衰减常数: %.4f', fitted_params{i}(2)), ...
                'Units', 'normalized', 'FontSize', 8, 'BackgroundColor', 'white');
        end
    end
    
%     % 校正后信号（ΔF/F0）
%     subplot(4, 4, i+4);
%     plot(time_axis, corrected_data(:, i), 'g-', 'LineWidth', 1.5);
%     title(sprintf('校正后 ΔF/F0 (细胞 %d)', i), 'FontSize', 10, 'FontWeight', 'bold');
%     xlabel('时间 (秒)', 'FontSize', 9);
%     ylabel('ΔF/F0', 'FontSize', 9);
%     grid on;
    
%     % 添加零线参考
%     hold on;
%     plot([time_axis(1), time_axis(end)], [0, 0], 'k--', 'LineWidth', 0.5);
end

% % 校正前后的平均信号对比
% subplot(4, 4, [11, 12, 15, 16]);
% % 计算平均信号
% mean_raw = mean(calcium_data, 2);
% mean_corrected = mean(corrected_data, 2);
% 
% yyaxis left;
% plot(time_axis, mean_raw, 'b-', 'LineWidth', 1.5);
% ylabel('平均原始荧光 (a.u.)', 'FontSize', 10);
% ylim([min(mean_raw)*0.9, max(mean_raw)*1.1]);
% 
% yyaxis right;
% plot(time_axis, mean_corrected, 'r-', 'LineWidth', 1.5);
% ylabel('平均 ΔF/F0', 'FontSize', 10);
% 
% title('校正前后平均信号对比', 'FontSize', 11, 'FontWeight', 'bold');
% xlabel('时间 (秒)', 'FontSize', 10);
% legend('原始平均', '校正后平均 ΔF/F0', 'Location', 'best', 'FontSize', 9);
% grid on;

% %% 5. 保存结果
% % 创建结果表格
% result_table = array2table([time_axis, corrected_dfof], ...
%     'VariableNames', ['Time_sec', cell_names]);
% 
% % 保存为CSV文件
% output_filename = 'calcium_bleach_corrected_results.csv';
% writetable(result_table, output_filename);
% fprintf('结果已保存至: %s\n', output_filename);
% 
% % 保存拟合参数
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
% fprintf('拟合参数已保存至: bleach_fitting_parameters.csv\n\n');

%% 
%计算统计指标
num_cells = size(calcium_data, 2);

% 初始化统计矩阵
stats_table = table();
stats_table.CellName = cell_names';
stats_table.CellIndex = (1:num_cells)';
stats_table.Mean = zeros(num_cells, 1);
stats_table.Std = zeros(num_cells, 1);
stats_table.DynamicRange = zeros(num_cells, 1); % 最大值-最小值
stats_table.CV = zeros(num_cells, 1); % 变异系数 = 标准差/均值
stats_table.Kurtosis = zeros(num_cells, 1);
stats_table.Skewness = zeros(num_cells, 1); % 额外添加偏度
stats_table.Min = zeros(num_cells, 1);
stats_table.Max = zeros(num_cells, 1);

% 计算每个细胞的统计量
for i = 1:num_cells
    %signal = corrected_data(:, i);
    signal = F_over_F0(:, i);
    % 移除NaN值
    signal = signal(~isnan(signal));
    
    if ~isempty(signal)
        stats_table.Mean(i) = mean(signal);
        stats_table.Std(i) = std(signal);
        stats_table.DynamicRange(i) = max(signal) - min(signal);
        
        % 变异系数 (避免除0)
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
% ========== 偏态CV分布的τ值提取 ==========
% 假设您已有细胞CV值数组: cell_cv_values

% 方法1：逆累积分布法（最稳健，推荐首选）
sorted_cv = sort(stats_table.CV);
n = length(sorted_cv);
cdf = (1:n) / n;

% 找到多个特征衰减点
tau_50 = interp1(cdf, sorted_cv, 0.50);  % 中位数
tau_63 = interp1(cdf, sorted_cv, 0.632); % 经典1-1/e衰减点
tau_80 = interp1(cdf, sorted_cv, 0.80);  % 80%衰减点
tau_90 = interp1(cdf, sorted_cv, 0.90);  % 90%衰减点

% 方法2：双指数混合模型拟合（适合偏态分布）
% 模型：p(x) = w1*exp(-x/τ1) + w2*exp(-x/τ2), w1+w2=1
double_exp_pdf = @(p, x) p(1)*exp(-x/p(2)) + (1-p(1))*exp(-x/p(3));

% 准备直方图数据
[counts, bin_centers] = hist(stats_table.CV, 50);
bin_width = bin_centers(2) - bin_centers(1);
pdf_values = counts / (sum(counts) * bin_width);

% 初始参数猜测
init_w1 = 0.7;  % 第一个成分的权重
init_tau1 = mean(stats_table.CV) * 0.5;  % 快衰减成分
init_tau2 = mean(stats_table.CV) * 2;    % 慢衰减成分

% 拟合双指数
try
    [params_double, resnorm] = lsqcurvefit(double_exp_pdf, ...
        [init_w1, init_tau1, init_tau2], ...
        bin_centers, pdf_values, ...
        [0.1, 0.001, 0.001], [0.9, Inf, Inf]); % 边界约束
    
    % 计算加权平均τ
    w1_fit = params_double(1);
    tau1_fit = params_double(2);
    tau2_fit = params_double(3);
    tau_weighted = w1_fit * tau1_fit + (1-w1_fit) * tau2_fit;
    
    fprintf('双指数拟合结果:\n');
    fprintf('  快成分τ1 = %.4f (权重=%.2f)\n', tau1_fit, w1_fit);
    fprintf('  慢成分τ2 = %.4f (权重=%.2f)\n', tau2_fit, 1-w1_fit);
    fprintf('  加权平均τ = %.4f\n', tau_weighted);
    
catch
    fprintf('双指数拟合失败，使用逆累积分布法\n');
    tau_weighted = tau_63; % 回退到方法1
end

% 方法3：Gamma分布拟合（专为偏态正数分布设计）
if exist('gamfit', 'file')
    [param_gam, ci_gam] = gamfit(stats_table.CV);
    shape_param = param_gam(1);  % 形状参数k
    scale_param = param_gam(2);  % 尺度参数θ
    
    % Gamma分布均值 = k*θ，可视为特征τ
    tau_gamma = shape_param * scale_param;
    
    fprintf('Gamma分布拟合:\n');
    fprintf('  形状参数k = %.4f, 尺度参数θ = %.4f\n', shape_param, scale_param);
    fprintf('  分布均值τ = %.4f\n', tau_gamma);
end

% 方法4：偏态修正τ（考虑分布偏度）
cv_skewness = skewness(stats_table.CV);
cv_mean = mean(stats_table.CV);
cv_median = median(stats_table.CV);

% 偏态修正公式：τ_skew = 中位数 * (1 + 偏度校正因子)
if cv_skewness > 0
    % 右偏分布：均值 > 中位数，使用修正
    skew_correction = min(0.5, cv_skewness/5); % 限制校正幅度
    tau_skew_adjusted = cv_median * (1 + skew_correction);
else
    % 左偏或对称，直接用中位数
    tau_skew_adjusted = cv_median;
end

% ========== 输出所有τ估计值 ==========
fprintf('\n========== CV分布τ值估计汇总 ==========\n');
fprintf('基本统计:\n');
fprintf('  均值 = %.4f, 中位数 = %.4f, 偏度 = %.4f\n', cv_mean, cv_median, cv_skewness);
fprintf('逆累积分布法:\n');
fprintf('  τ(50%%) = %.4f, τ(63%%) = %.4f\n', tau_50, tau_63);
fprintf('  τ(80%%) = %.4f, τ(90%%) = %.4f\n', tau_80, tau_90);
if exist('tau_weighted', 'var')
    fprintf('双指数模型: 加权τ = %.4f\n', tau_weighted);
end
if exist('tau_gamma', 'var')
    fprintf('Gamma分布: 均值τ = %.4f\n', tau_gamma);
end
fprintf('偏态修正: τ_skew = %.4f\n', tau_skew_adjusted);
%% 应用筛选标准：只考虑 变异系数 >= 0.03   %动态范围 >= 100 且
stats_table.IsActive = stats_table.CV >= 0.035;
stats_table.ActivityType = cell(num_cells, 1);
for i = 1:num_cells
    if stats_table.IsActive(i)
        stats_table.ActivityType{i} = '活跃';
    else
        stats_table.ActivityType{i} = '不活跃';
    end
end

% 获取活跃和不活跃细胞的索引
active_indices = find(stats_table.IsActive);
inactive_indices = find(~stats_table.IsActive);

fprintf('\n=== 细胞分类结果 ===\n');
fprintf('活跃细胞数量 (动态范围 >= 100 且 CV >= 0.02): %d\n', length(active_indices));
fprintf('不活跃细胞数量: %d\n', length(inactive_indices));

% 显示活跃细胞的详细信息
if ~isempty(active_indices)
    fprintf('\n--- 活跃细胞列表 ---\n');
    for i = 1:length(active_indices)
        idx = active_indices(i);
        fprintf('细胞 %d: %s (动态范围: %.2f, CV: %.4f)\n', ...
            idx, cell_names{idx}, stats_table.DynamicRange(idx), stats_table.CV(idx));
    end
end

% 显示不活跃细胞的详细信息
if ~isempty(inactive_indices)
    fprintf('\n--- 不活跃细胞列表 ---\n');
    for i = 1:min(length(inactive_indices), 10) % 只显示前10个不活跃细胞
        idx = inactive_indices(i);
        fprintf('细胞 %d: %s (动态范围: %.2f, CV: %.4f)\n', ...
            idx, cell_names{idx}, stats_table.DynamicRange(idx), stats_table.CV(idx));
    end
    if length(inactive_indices) > 10
        fprintf('... 还有 %d 个不活跃细胞\n', length(inactive_indices) - 10);
    end
end

%% 从活跃和不活跃细胞中各选择5个进行可视化
% 如果某类细胞不足5个，则使用所有可用细胞
num_to_select = 5;

if length(active_indices) >= num_to_select
    selected_active = active_indices(randperm(length(active_indices), num_to_select));
else
    selected_active = active_indices;
    fprintf('\n注意: 只有 %d 个活跃细胞，全部用于可视化\n', length(active_indices));
end

if length(inactive_indices) >= num_to_select
    selected_inactive = inactive_indices(randperm(length(inactive_indices), num_to_select));
else
    selected_inactive = inactive_indices;
    fprintf('注意: 只有 %d 个不活跃细胞，全部用于可视化\n', length(inactive_indices));
end

fprintf('\n=== 选择的细胞用于可视化 ===\n');
fprintf('活跃细胞: %s\n', mat2str(selected_active));
fprintf('不活跃细胞: %s\n', mat2str(selected_inactive));

%% 可视化活跃与不活跃细胞的对比
figure('Position', [100, 100, 1400, 1000]);

% 子图1: 活跃细胞的钙信号轨迹
subplot(2, 3, 1);
if ~isempty(selected_active)
    colors_active = parula(length(selected_active));
    for i = 1:length(selected_active)
        idx = selected_active(i);
        plot(corrected_data(:, idx), 'Color', colors_active(i, :), 'LineWidth', 2);
        hold on;
    end
    xlabel('时间点');
    ylabel('钙信号强度');
    title('活跃细胞钙信号轨迹');
    % 创建图例标签
    legend_labels_active = cell(length(selected_active), 1);
    for i = 1:length(selected_active)
        idx = selected_active(i);
        legend_labels_active{i} = sprintf('细胞 %d', idx);
    end
    legend(legend_labels_active, 'Location', 'eastoutside');
    grid on;
else
    text(0.5, 0.5, '无活跃细胞', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title('活跃细胞钙信号轨迹');
end

% 子图2: 不活跃细胞的钙信号轨迹
subplot(2, 3, 2);
if ~isempty(selected_inactive)
    colors_inactive = parula(length(selected_inactive)); 
    for i = 1:length(selected_inactive)
        idx = selected_inactive(i);
        plot(corrected_data(:, idx), 'Color', colors_inactive(i, :), 'LineWidth', 2);
        hold on;
    end
    xlabel('时间点');
    ylabel('钙信号强度');
    title('不活跃细胞钙信号轨迹');
    % 创建图例标签
    legend_labels_inactive = cell(length(selected_inactive), 1);
    for i = 1:length(selected_inactive)
        idx = selected_inactive(i);
        legend_labels_inactive{i} = sprintf('细胞 %d', idx);
    end
    legend(legend_labels_inactive, 'Location', 'eastoutside');
    grid on;
else
    text(0.5, 0.5, '无不活跃细胞', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title('不活跃细胞钙信号轨迹');
end

% 子图3: 动态范围分布比较
subplot(2, 3, 3);
if ~isempty(active_indices) && ~isempty(inactive_indices)
    % 创建分组数据
    group = [ones(length(active_indices), 1); 2*ones(length(inactive_indices), 1)];
    data_dr = [stats_table.DynamicRange(active_indices); stats_table.DynamicRange(inactive_indices)];
    
    boxplot(data_dr, group, 'Labels', {'活跃细胞', '不活跃细胞'});
    ylabel('动态范围');
    title('动态范围分布比较');
    grid on;
else
    text(0.5, 0.5, '数据不足', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title('动态范围分布比较');
end

% 子图4: 变异系数分布比较
subplot(2, 3, 4);
if ~isempty(active_indices) && ~isempty(inactive_indices)
    % 创建分组数据
    group = [ones(length(active_indices), 1); 2*ones(length(inactive_indices), 1)];
    data_cv = [stats_table.CV(active_indices); stats_table.CV(inactive_indices)];
    
    boxplot(data_cv, group, 'Labels', {'活跃细胞', '不活跃细胞'});
    ylabel('变异系数 (CV)');
    title('变异系数分布比较');
    grid on;
else
    text(0.5, 0.5, '数据不足', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title('变异系数分布比较');
end

% 子图5: 动态范围 vs 变异系数散点图
subplot(2, 3, 5);
if ~isempty(active_indices) && ~isempty(inactive_indices)
    scatter(stats_table.DynamicRange(active_indices), stats_table.CV(active_indices), ...
            50, 'g', 'filled', 'MarkerFaceAlpha', 0.7);
    hold on;
    scatter(stats_table.DynamicRange(inactive_indices), stats_table.CV(inactive_indices), ...
            50, 'r', 'filled', 'MarkerFaceAlpha', 0.7);
    
    % 添加分类线
    xline(100, '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    yline(0.02, '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
    
    xlabel('动态范围');
    ylabel('变异系数 (CV)');
    title('动态范围 vs 变异系数');
    legend('活跃细胞', '不活跃细胞', '分类阈值', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, '数据不足', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title('动态范围 vs 变异系数');
end

% 子图6: 两类细胞的统计特性对比
subplot(2, 3, 6);
if ~isempty(active_indices) && ~isempty(inactive_indices)
    % 定义要比较的统计指标
    stat_names = {'Mean', 'Std', 'DynamicRange', 'CV', 'Kurtosis', 'Skewness'};
    stat_labels = {'均值', '标准差', '动态范围', '变异系数', '峰度', '偏度'};
    
    % 计算两类细胞的平均统计值
    active_means = zeros(1, length(stat_names));
    inactive_means = zeros(1, length(stat_names));
    
    for i = 1:length(stat_names)
        active_means(i) = mean(stats_table.(stat_names{i})(active_indices), 'omitnan');
        inactive_means(i) = mean(stats_table.(stat_names{i})(inactive_indices), 'omitnan');
    end
    
    % 创建分组柱状图 - 确保X和Y长度相同
    X = 1:length(stat_names);
    bar_data = [active_means; inactive_means]';
    
    % 绘制柱状图
    h = bar(X, bar_data);
    
    % 设置颜色
    h(1).FaceColor = [0.2, 0.8, 0.2]; % 活跃细胞 - 绿色
    h(2).FaceColor = [0.8, 0.2, 0.2]; % 不活跃细胞 - 红色
    
    xlabel('统计指标');
    ylabel('平均值');
    title('两类细胞统计特性对比');
    set(gca, 'XTickLabel', stat_labels, 'XTickLabelRotation', 45);
    legend('活跃细胞', '不活跃细胞', 'Location', 'best');
    grid on;
    
    % 添加数值标签
    for i = 1:length(X)
        % 活跃细胞的数值标签
        if ~isnan(active_means(i))
            text(X(i)-0.18, active_means(i)+max(active_means)/50, sprintf('%.2f', active_means(i)), ...
                'FontSize', 7, 'HorizontalAlignment', 'center');
        end
        
        % 不活跃细胞的数值标签
        if ~isnan(inactive_means(i))
            text(X(i)+0.18, inactive_means(i)+max(inactive_means)/50, sprintf('%.2f', inactive_means(i)), ...
                'FontSize', 7, 'HorizontalAlignment', 'center');
        end
    end
else
    text(0.5, 0.5, '数据不足', 'HorizontalAlignment', 'center', 'Units', 'normalized');
    title('两类细胞统计特性对比');
end

sgtitle('活跃与不活跃细胞对比分析', 'FontSize', 14, 'FontWeight', 'bold');
%%
%%可视化
% 方式1作为示例，你可以修改这些索引
selected_indices = [4,12,31]; % 修改为你想要分析的细胞索引
fprintf('\n=== 指定细胞的详细统计参数 ===\n');
for i = 1:length(selected_indices)
    idx = selected_indices(i);
    fprintf('\n--- 细胞 %d: %s ---\n', idx, cell_names{idx});
    fprintf('均值: %.4f\n', stats_table.Mean(idx));
    fprintf('标准差: %.4f\n', stats_table.Std(idx));
    fprintf('最小值: %.4f\n', stats_table.Min(idx));
    fprintf('最大值: %.4f\n', stats_table.Max(idx));
    fprintf('动态范围: %.4f\n', stats_table.DynamicRange(idx));
    fprintf('变异系数(CV): %.4f\n', stats_table.CV(idx));
    fprintf('峰度: %.4f\n', stats_table.Kurtosis(idx));
    fprintf('偏度: %.4f\n', stats_table.Skewness(idx));
end
% 可视化指定细胞的信号和统计
figure('Position', [100, 100, 1400, 1000]);

% 子图1: 指定细胞的钙信号轨迹
subplot(2, 3, 1);
colors = hsv(length(selected_indices));
for i = 1:length(selected_indices)
    idx = selected_indices(i);
    plot(corrected_data(:, idx), 'Color', colors(i, :), 'LineWidth', 2);
    hold on;
end
xlabel('时间点');
ylabel('钙信号强度');
title('指定细胞的钙信号轨迹');
% 创建图例标签 - 使用索引和名称
legend_labels = cell(length(selected_indices), 1);
for i = 1:length(selected_indices)
    idx = selected_indices(i);
    legend_labels{i} = sprintf('细胞 %d: %s', idx, cell_names{idx});
end
legend(legend_labels, 'Location', 'eastoutside', 'Interpreter', 'none');
grid on;

% 子图2: 指定细胞的标准差比较
subplot(2, 3, 2);
bar(stats_table.Std(selected_indices), 'FaceColor', [0.2, 0.6, 0.8]);
% 使用索引作为x轴标签
x_labels = cell(length(selected_indices), 1);
for i = 1:length(selected_indices)
    x_labels{i} = num2str(selected_indices(i));
end
set(gca, 'XTickLabel', x_labels);
xlabel('细胞索引');
ylabel('标准差');
title('指定细胞的标准差比较');
grid on;

% 子图3: 指定细胞的动态范围比较
subplot(2, 3, 3);
bar(stats_table.DynamicRange(selected_indices), 'FaceColor', [0.8, 0.4, 0.2]);
set(gca, 'XTickLabel', x_labels);
xlabel('细胞索引');
ylabel('动态范围');
title('指定细胞的动态范围比较');
grid on;

% 子图4: 指定细胞的变异系数比较
subplot(2, 3, 4);
bar(stats_table.CV(selected_indices), 'FaceColor', [0.4, 0.8, 0.4]);
set(gca, 'XTickLabel', x_labels);
xlabel('细胞索引');
ylabel('变异系数 (CV)');
title('指定细胞的变异系数比较');
grid on;

% 子图5: 指定细胞的峰度比较
subplot(2, 3, 5);
bar(stats_table.Kurtosis(selected_indices), 'FaceColor', [0.8, 0.2, 0.8]);
hold on;
plot(xlim, [3, 3], 'r--', 'LineWidth', 2); % 正态分布参考线
set(gca, 'XTickLabel', x_labels);
xlabel('细胞索引');
ylabel('峰度');
title('指定细胞的峰度比较 (红线: 正态分布=3)');
grid on;

% 子图6: 指定细胞的偏度比较
subplot(2, 3, 6);
bar(stats_table.Skewness(selected_indices), 'FaceColor', [0.9, 0.7, 0.1]);
hold on;
plot(xlim, [0, 0], 'r--', 'LineWidth', 2); % 对称分布参考线
set(gca, 'XTickLabel', x_labels);
xlabel('细胞索引');
ylabel('偏度');
title('指定细胞的偏度比较 (红线: 对称分布=0)');
grid on;

sgtitle('指定细胞的钙信号统计分析', 'FontSize', 14, 'FontWeight', 'bold');

%% 可视化统计分布
figure('Position', [100, 100, 1400, 1000]);

% 子图1: 标准差分布
subplot(2, 3, 1);
histogram(stats_table.Std, 20, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'black');
xlabel('标准差');
ylabel('细胞数量');
title('标准差分布');
grid on;

% 子图2: 动态范围分布
subplot(2, 3, 2);
histogram(stats_table.DynamicRange, 20, 'FaceColor', [0.8, 0.4, 0.2], 'EdgeColor', 'black');
xlabel('动态范围');
ylabel('细胞数量');
title('动态范围分布');
grid on;

% 子图3: 变异系数分布
subplot(2, 3, 3);
histogram(stats_table.CV, 20, 'FaceColor', [0.4, 0.8, 0.4], 'EdgeColor', 'black');
xlabel('变异系数 (CV)');
ylabel('细胞数量');
title('变异系数分布');
grid on;

% 子图4: 峰度分布
subplot(2, 3, 4);
histogram(stats_table.Kurtosis, 20, 'FaceColor', [0.8, 0.2, 0.8], 'EdgeColor', 'black');
hold on;
% 标记正态分布参考线
y_limits = ylim;
plot([3, 3], y_limits, 'r--', 'LineWidth', 2);
xlabel('峰度');
ylabel('细胞数量');
title('峰度分布 (红线: 正态分布=3)');
legend('细胞分布', '正态分布', 'Location', 'best');
grid on;

% 子图5: 偏度分布
subplot(2, 3, 5);
histogram(stats_table.Skewness, 20, 'FaceColor', [0.9, 0.7, 0.1], 'EdgeColor', 'black');
hold on;
% 标记对称分布参考线
y_limits = ylim;
plot([0, 0], y_limits, 'r--', 'LineWidth', 2);
xlabel('偏度');
ylabel('细胞数量');
title('偏度分布 (红线: 对称分布=0)');
legend('细胞分布', '对称分布', 'Location', 'best');
grid on;

% 子图6: 所有信号的示例图
subplot(2, 3, 6);
% 随机选择几个代表性信号绘制
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
    xlabel('时间点');
    ylabel('钙信号强度');
    title('示例钙信号轨迹');
    legend(cell_names(cells_to_plot), 'Location', 'best', 'Interpreter', 'none');
    grid on;
end

sgtitle('钙信号统计特性分析', 'FontSize', 14, 'FontWeight', 'bold');

%%
%获取活跃细胞的原始钙信号数据
active_cell_indices = find(stats_table.IsActive);
active_cell_names = cell_names(active_cell_indices);

active_calcium_data = calcium_data(:, active_cell_indices);
%% 计算活跃细胞的Pearson相关系数矩阵
if ~isempty(active_calcium_data)
    fprintf('\n=== 计算活跃细胞的Pearson相关系数矩阵 ===\n');
    
    % 计算相关系数矩阵
    correlation_matrix = corr(active_calcium_data, 'Rows', 'pairwise');
    
    % 获取活跃细胞数量
    num_active_cells = size(active_calcium_data, 2);
    
    fprintf('活跃细胞数量: %d\n', num_active_cells);
    fprintf('相关系数矩阵大小: %d x %d\n', size(correlation_matrix));
    fprintf('相关系数范围: [%.4f, %.4f]\n', min(correlation_matrix(:)), max(correlation_matrix(:)));
    
    % 提取上三角部分（不包括对角线）用于分布分析
    triu_indices = triu(true(size(correlation_matrix)), 1);
    correlation_values = correlation_matrix(triu_indices);
    
    fprintf('相关系数统计:\n');
    fprintf('  平均值: %.4f\n', mean(correlation_values, 'omitnan'));
    fprintf('  中位数: %.4f\n', median(correlation_values, 'omitnan'));
    fprintf('  标准差: %.4f\n', std(correlation_values, 'omitnan'));
    
    % 计算显著相关的比例（例如 |r| > 0.5）
    strong_positive = sum(correlation_values > 0.5) / length(correlation_values) * 100;
    strong_negative = sum(correlation_values < -0.5) / length(correlation_values) * 100;
    moderate_positive = sum(correlation_values > 0.3 & correlation_values <= 0.5) / length(correlation_values) * 100;
    moderate_negative = sum(correlation_values < -0.3 & correlation_values >= -0.5) / length(correlation_values) * 100;
    
    fprintf('强正相关 (r > 0.5): %.2f%%\n', strong_positive);
    fprintf('中等正相关 (0.3 < r <= 0.5): %.2f%%\n', moderate_positive);
    fprintf('强负相关 (r < -0.5): %.2f%%\n', strong_negative);
    fprintf('中等负相关 (-0.5 <= r < -0.3): %.2f%%\n', moderate_negative);
    
    %% 可视化相关系数分布
    figure('Position', [100, 100, 1200, 500]);
    
    % 子图1: 相关系数分布直方图
    subplot(1, 2, 1);
    histogram(correlation_values, 50, 'FaceColor', [0.3, 0.6, 0.9], 'EdgeColor', 'black');
    hold on;

    % 添加参考线 - 先保存句柄用于图例
    y_limits = ylim;
    h_zero = plot([0, 0], y_limits, 'k--', 'LineWidth', 1.5);
    h_strong = plot([0.5, 0.5], y_limits, 'r--', 'LineWidth', 1);
    h_strong_neg = plot([-0.5, -0.5], y_limits, 'r--', 'LineWidth', 1);
    h_moderate = plot([0.3, 0.3], y_limits, 'g--', 'LineWidth', 1);
    h_moderate_neg = plot([-0.3, -0.3], y_limits, 'g--', 'LineWidth', 1);

    xlabel('Pearson 相关系数');
    ylabel('频率');
    title('活跃细胞相关系数分布');
    % 图例使用实际的线条句柄，确保颜色一致
    legend([h_zero, h_strong, h_moderate], ...
        '零相关', '强相关阈值 (±0.5)', '中等相关阈值 (±0.3)', 'Location', 'best');
    grid on;
    
    % 子图2: 相关系数箱线图
    subplot(1, 2, 2);
    boxplot(correlation_values, 'Orientation', 'horizontal');
    xlabel('Pearson 相关系数');
    title('相关系数箱线图');
    grid on;
    
    sgtitle('活跃细胞Pearson相关系数分析', 'FontSize', 14, 'FontWeight', 'bold');
    
    %% 进行层次聚类（获取cluster_order）
    fprintf('\n=== 进行层次聚类分析 ===\n');
    
    % 使用层次聚类对相关系数矩阵进行分群
    % 计算距离矩阵（1 - 相关系数）
    distance_matrix = 1 - correlation_matrix;
    
    % 初始化cluster_order
    cluster_order = 1:num_active_cells; % 默认顺序
    
    if num_active_cells > 1
        % 将距离矩阵转换为向量格式（pdist格式）
        distance_vector = squareform(distance_matrix, 'tovector');
        
        % 进行层次聚类
        linkage_tree = linkage(distance_vector, 'average');
        
        % 获取聚类顺序 - 使用dendrogram获取排序
        [~, ~, cluster_order] = dendrogram(linkage_tree, 0);
        cluster_order = cluster_order';
        fprintf('层次聚类完成，获得细胞排序顺序。\n');
    else
        % 如果只有一个细胞，直接使用1
        cluster_order = 1;
        linkage_tree = [];
        fprintf('只有一个活跃细胞，无需聚类。\n');
    end
    
    %% 创建带分群的热图
    figure('Position', [100, 100, max(800, num_active_cells*40), 800]);
    
    % 按照聚类顺序重新排列相关系数矩阵和细胞名称
    sorted_correlation_matrix = correlation_matrix(cluster_order, cluster_order);
    sorted_cell_names = active_cell_names(cluster_order);
    
    % 创建热图
    imagesc(sorted_correlation_matrix);
    colormap(jet); % 使用jet颜色图，蓝色表示负相关，红色表示正相关
    colorbar;
    caxis([-1, 1]); % 设置颜色轴范围为-1到1
    
    % 设置坐标轴标签 - 强制显示所有标签
    xticks(1:num_active_cells);
    yticks(1:num_active_cells);
    xticklabels(sorted_cell_names);
    yticklabels(sorted_cell_names);
    
    % 动态调整布局
    if num_active_cells > 20
        % 细胞很多时的布局
        xtickangle(90);
        ytickangle(0);
        font_size = max(6, 10 - num_active_cells/20);
        set(gca, 'FontSize', font_size);
        
        % 调整图形大小以适应所有标签
        fig_position = get(gcf, 'Position');
        set(gcf, 'Position', [100, 100, max(1000, num_active_cells*30), 800]);
    elseif num_active_cells > 10
        % 中等数量细胞的布局
        xtickangle(45);
        ytickangle(0);
        set(gca, 'FontSize', 9);
    else
        % 少量细胞的布局
        xtickangle(0);
        ytickangle(0);
        set(gca, 'FontSize', 10);
    end
    
    xlabel('细胞');
    ylabel('细胞');
    title('活跃细胞Pearson相关系数热图（带分群）');
    
    % 添加网格线以便更好地观察
    hold on;
    for i = 1:num_active_cells
        plot([0.5, num_active_cells+0.5], [i+0.5, i+0.5], 'k-', 'LineWidth', 0.5);
        plot([i+0.5, i+0.5], [0.5, num_active_cells+0.5], 'k-', 'LineWidth', 0.5);
    end
    
    % 在热图右侧单独显示树状图
    if num_active_cells > 1
        dendro_axes = axes('Position', [0.85, 0.1, 0.12, 0.8]);
        dendrogram(linkage_tree, 0, 'Orientation', 'right');
        set(dendro_axes, 'XTick', [], 'YTick', [], 'YAxisLocation', 'right');
        title(dendro_axes, '聚类树');
    end
    
    %% 创建网络可视化（如果有多个细胞）
    if num_active_cells > 1
        figure('Position', [100, 100, 1000, 800]);
        
        % 只显示较强的相关性（例如 |r| > 0.3）
        threshold = 0.3;
        strong_correlations = sorted_correlation_matrix;
        strong_correlations(abs(strong_correlations) < threshold) = 0;
        
        % 创建图对象
        G = graph(strong_correlations, sorted_cell_names, 'upper');
        
        % 计算节点大小（基于动态范围）
        node_sizes = zeros(length(sorted_cell_names), 1);
        for i = 1:length(sorted_cell_names)
            cell_idx = active_cell_indices(cluster_order(i));
            node_sizes(i) = 100 + stats_table.DynamicRange(cell_idx) / 10; % 缩放因子
        end
        
        % 计算节点颜色（基于平均钙信号）
        node_colors = zeros(length(sorted_cell_names), 1);
        for i = 1:length(sorted_cell_names)
            cell_idx = active_cell_indices(cluster_order(i));
            node_colors(i) = stats_table.Mean(cell_idx);
        end
        
        % 绘制网络图
        p = plot(G, 'Layout', 'force', 'UseGravity', true);
        
        % 设置节点属性
        p.MarkerSize = node_sizes / 20; % 进一步缩放
        p.NodeCData = node_colors;
        p.NodeLabel = sorted_cell_names;
        
        % 设置边属性
        edge_weights = G.Edges.Weight;
        positive_edges = edge_weights > 0;
        negative_edges = edge_weights < 0;
        
        % 正相关边用红色，负相关边用蓝色
        p.EdgeColor = zeros(length(edge_weights), 3);
        p.EdgeColor(positive_edges, :) = repmat([1, 0, 0], sum(positive_edges), 1);
        p.EdgeColor(negative_edges, :) = repmat([0, 0, 1], sum(negative_edges), 1);
        
        % 设置边宽度基于相关强度
        p.LineWidth = 2 * abs(edge_weights);
        
        colorbar;
        colormap(parula);
        title(sprintf('活跃细胞相关性网络 (|r| > %.1f)', threshold));
        xlabel('节点大小: 动态范围 | 节点颜色: 平均钙信号 | 边颜色: 红=正相关, 蓝=负相关');
    end
    
else
    fprintf('\n没有活跃细胞，无法计算相关系数矩阵。\n');
end

% function trend = robust_moving_average(signal, time_axis, num_frames)
%     % 对全段数据执行稳健的移动平均，用于估计漂白背景趋势
%     % 输入： signal - 原始信号， time_axis - 时间轴， num_frames - 总帧数
%     % 输出： trend - 估计的背景趋势曲线
%     
%     % 1. 动态计算窗宽：总时长的15%-25%，并确保为奇数且合理
%     total_time = time_axis(end) - time_axis(1);
%     avg_window_ratio = 0.2; % 默认窗宽占总时长的20%
%     window_in_seconds = total_time * avg_window_ratio;
%     
%     % 将时间窗宽转换为帧数（假设时间均匀采样）
%     avg_frame_duration = median(diff(time_axis));
%     window_in_frames = round(window_in_seconds / avg_frame_duration);
%     
%     % 确保窗宽是奇数（使窗口对称），并在合理范围内
%     window_in_frames = max(5, min(floor(num_frames/4), window_in_frames)); % 最小5帧，最大不超过总帧数1/4
%     if mod(window_in_frames, 2) == 0
%         window_in_frames = window_in_frames + 1; % 转为奇数
%     end
%     fprintf('         移动平均窗宽: %d 帧 (约%.1f秒)\n', window_in_frames, window_in_frames*avg_frame_duration);
%     
%     % 2. 使用 'movmean' 计算移动平均，处理边界 ('endpoints' 参数很关键)
%     % 'discard' 选项会使输出变短，我们使用 'fill' 用NaN填充，然后自己处理
%     raw_trend = movmean(signal, window_in_frames, 'Endpoints', 'fill');
%     
%     % 3. 处理边界点（开头和结尾 window_in_frames/2 的帧）
%     half_win = floor(window_in_frames / 2);
%     % 3.1 对于开头，用前 half_win 帧的平均值填充
%     raw_trend(1:half_win) = mean(signal(1:half_win));
%     % 3.2 对于结尾，用后 half_win 帧的平均值填充
%     raw_trend(end-half_win+1:end) = mean(signal(end-half_win+1:end));
%     
%     % 4. 二次平滑 (可选，使用更小的窗宽让曲线更平滑)
%     final_trend = smoothdata(raw_trend, 'movmean', 5);
%     
%     % 5. 最终合理性检查：确保趋势整体非上升 (漂白特性)
%     % 计算首尾差值，如果趋势是上升的，则强制赋予一个非常平缓的下降趋势
%     if final_trend(end) > final_trend(1) * 1.01 % 如果上升超过1%
%         fprintf('         警告: 移动平均趋势轻微上升，已强制平缓化\n');
%         % 赋予一个极平缓的线性下降，斜率基于信号初始值的0.5%
%         gentle_slope = -0.005 * final_trend(1) / (num_frames - 1);
%         final_trend = final_trend(1) + gentle_slope * (0:num_frames-1)';
%     end
%     
%     trend = final_trend;
% end