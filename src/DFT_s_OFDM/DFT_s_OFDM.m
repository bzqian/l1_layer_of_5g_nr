% DFT-s-OFDM 仿真（不使用通信工具箱）
clear; clc;

% 参数设置
M = 64;        % DFT大小（子载波数）
N = 256;       % IFFT大小
cp_len = 64;   % 循环前缀长度（通常取IFFT大小的1/4）
snr_dB = 30;   % 信噪比（dB）
num_symbols = 1;% 传输的DFT-s-OFDM符号数

% 生成16QAM星座表（手动映射）
qam_table = (1/sqrt(10)) * [  % 归一化因子sqrt(10)
    -3-3i; -3-1i; -3+3i; -3+1i;   % 索引0-3
    -1-3i; -1-1i; -1+3i; -1+1i;   % 索引4-7
    3-3i;  3-1i;  3+3i;  3+1i;    % 索引8-11
    1-3i;  1-1i;  1+3i;  1+1i;    % 索引12-15
];

%% 发射端
% 生成随机二进制数据
num_bits = M * 4 * num_symbols; % 每个符号4位
tx_bits = randi([0 1], num_bits, 1);

% 二进制转十进制索引（手动实现）
tx_groups = reshape(tx_bits, 4, [])';
tx_symbols_index = tx_groups * [8; 4; 2; 1];

% QAM调制
tx_symbols = qam_table(tx_symbols_index + 1); % MATLAB索引从1开始

% DFT预编码
dft_precoded = fft(tx_symbols)/sqrt(M);

% 子载波映射（中心对齐）
ifft_input = zeros(N, 1);
start_idx = floor((N-M)/2)+1;
ifft_input(start_idx:start_idx+M-1) = dft_precoded;

% IFFT变换
tx_time = ifft(ifft_input) * sqrt(N);

% 添加循环前缀
tx_signal = [tx_time(end-cp_len+1:end); tx_time];

%% 信道（AWGN）
% 计算信号功率
sig_power = mean(abs(tx_signal).^2);

% 计算噪声功率
noise_power = sig_power / (10^(snr_dB/10));

% 生成复高斯噪声
noise = sqrt(noise_power/2)*(randn(size(tx_signal)) + 1i*randn(size(tx_signal)));
rx_signal = tx_signal + noise;

%% 接收端
% 移除循环前缀
rx_time = rx_signal(cp_len+1 : cp_len+N);

% FFT变换
rx_freq = fft(rx_time)/sqrt(N);

% 子载波解映射
rx_dft = rx_freq(start_idx : start_idx+M-1);

% IDFT解码
rx_symbols = ifft(rx_dft) * sqrt(M);

% QAM解调（最近邻判决）
rx_symbols_index = zeros(size(rx_symbols));
for k = 1:length(rx_symbols)
    [~, idx] = min(abs(rx_symbols(k) - qam_table));
    rx_symbols_index(k) = idx-1; % MATLAB索引从1开始
end

% 十进制转二进制（手动实现）
rx_groups = zeros(length(rx_symbols_index),4);
for k = 1:length(rx_symbols_index)
    num = rx_symbols_index(k);
    rx_groups(k,1) = bitand(num,8)/8;   % 1000
    rx_groups(k,2) = bitand(num,4)/4;   % 0100
    rx_groups(k,3) = bitand(num,2)/2;   % 0010
    rx_groups(k,4) = bitand(num,1);     % 0001
end
rx_bits = reshape(rx_groups', [], 1);

% 计算误码率
ber = sum(tx_bits ~= rx_bits)/num_bits;
fprintf('误码率: %.4f\n', ber);