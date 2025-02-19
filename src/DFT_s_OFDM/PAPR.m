% DFT-s-OFDM 和 OFDMA 的 PAPR 比较（无工具箱）
clear; clc;

% 参数设置
N = 256;           % IFFT/FFT 大小
M = 64;            % 数据子载波数
num_symbols = 1e4; % 仿真符号数
mod_order = 16;    % 调制阶数（16QAM）

% 生成16QAM星座表（手动映射）
qam_table = (1/sqrt(10)) * [  % 归一化因子sqrt(10)
    -3-3i; -3-1i; -3+3i; -3+1i;   % 索引0-3
    -1-3i; -1-1i; -1+3i; -1+1i;   % 索引4-7
    3-3i;  3-1i;  3+3i;  3+1i;    % 索引8-11
    1-3i;  1-1i;  1+3i;  1+1i;    % 索引12-15
];

% 初始化 PAPR 存储
papr_dft_s_ofdm = zeros(num_symbols, 1);
papr_ofdma = zeros(num_symbols, 1);

% PAPR 计算函数
calculate_papr = @(x) 10*log10(max(abs(x).^2) / mean(abs(x).^2));

%% 仿真循环
for i = 1:num_symbols
    % 生成随机数据（0到15的索引）
    tx_data = randi([0 mod_order-1], M, 1);
    
    % 调制（16QAM）
    tx_symbols = qam_table(tx_data + 1); % MATLAB索引从1开始
    
    % === DFT-s-OFDM ===
    % DFT 预编码
    dft_precoded = fft(tx_symbols) / sqrt(M);
    
    % 子载波映射（中心对齐）
    ifft_input = zeros(N, 1);
    start_idx = floor((N-M)/2) + 1;
    ifft_input(start_idx:start_idx+M-1) = dft_precoded;
    
    % IFFT 变换
    tx_time_dft_s_ofdm = ifft(ifft_input) * sqrt(N);
    
    % 计算 PAPR
    papr_dft_s_ofdm(i) = calculate_papr(tx_time_dft_s_ofdm);
    
    % === OFDMA ===
    % 子载波映射（中心对齐）
    ifft_input = zeros(N, 1);
    ifft_input(start_idx:start_idx+M-1) = tx_symbols;
    
    % IFFT 变换
    tx_time_ofdma = ifft(ifft_input) * sqrt(N);
    
    % 计算 PAPR
    papr_ofdma(i) = calculate_papr(tx_time_ofdma);
end

%% 绘制 PAPR 分布
figure;
[count_dft, bin_dft] = hist(papr_dft_s_ofdm, 100);
[count_ofdma, bin_ofdma] = hist(papr_ofdma, 100);

% 归一化
count_dft = count_dft / num_symbols;
count_ofdma = count_ofdma / num_symbols;

% 绘制 PDF
plot(bin_dft, count_dft, 'b', 'LineWidth', 2);
hold on;
plot(bin_ofdma, count_ofdma, 'r', 'LineWidth', 2);
grid on;
xlabel('PAPR (dB)');
ylabel('Probability Density');
title('DFT-s-OFDM 和 OFDMA 的 PAPR 分布比较');
legend('DFT-s-OFDM', 'OFDMA');

% 计算并显示平均 PAPR
mean_papr_dft_s_ofdm = mean(papr_dft_s_ofdm);
mean_papr_ofdma = mean(papr_ofdma);
fprintf('DFT-s-OFDM 的平均 PAPR: %.2f dB\n', mean_papr_dft_s_ofdm);
fprintf('OFDMA 的平均 PAPR: %.2f dB\n', mean_papr_ofdma);