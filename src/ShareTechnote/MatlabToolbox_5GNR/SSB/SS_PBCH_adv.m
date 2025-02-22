% 定义gNB结构体
gNB = struct();

% 设置gNB参数
gNB.NDLRB = 100;               % 带宽，以资源块的数量表示
gNB.SubcarrierSpacing = 30;    % 子载波间隔 (kHz)
gNB.WaveformType = 'CP-OFDM';  % 波形类型: 'CP-OFDM', 'W-OFDM' 或 'F-OFDM'
gNB.CyclicPrefix = 'Normal';   % 循环前缀: 'Normal' 或 'Extended'
gNB.UseDCSubcarrier = 'Off';   % 是否使用DC子载波: 'On' 或 'Off'
gNB.NCellID = 1;               % 小区标识

% 设置SS Burst参数
gNB.SSBurst = struct();
gNB.SSBurst.BurstType = 'CaseC';   
gNB.SSBurst.SubcarrierOffset = 0;

% 位图指示在突发中传输的块
gNB.SSBurst.SSBTransmitted = [1 1 1 1 1 1 1 1];

% SS突发集的周期性 (ms)
gNB.SSBurst.SSBPeriodicity = 10;        

% 创建SS突发波形和信息结构
ssburst = gNB.SSBurst;
ssburst.DisplayBurst = true;

% 自定义函数生成SS突发波形和信息
[ssbWaveform, ssbInfo] = generateSSBurst(gNB, ssburst);

% 显示SS突发内容
if ssburst.DisplayBurst
    plotSSBurst(ssbWaveform, ssbInfo);
end


function [ssbWaveform, ssbInfo] = generateSSBurst(gNB, ssburst)
    % 参数提取
    NDLRB = gNB.NDLRB;               % 带宽，以资源块的数量表示
    SubcarrierSpacing = gNB.SubcarrierSpacing; % 子载波间隔 (kHz)
    NCellID = gNB.NCellID;           % 小区标识
    SSBTransmitted = ssburst.SSBTransmitted; % SSB传输位图
    SSBPeriodicity = ssburst.SSBPeriodicity; % SSB周期 (ms)
    SubcarrierOffset = ssburst.SubcarrierOffset; % 子载波偏移

    % 常量定义
    Nfft = 4096; % FFT点数
    NumSymbolsPerSSB = 4; % 每个SSB的OFDM符号数
    NumSubcarriersPerRB = 12; % 每个资源块的子载波数
    NumRBsPerSSB = 20; % 每个SSB占用的资源块数
    NumSubcarriersPerSSB = NumRBsPerSSB * NumSubcarriersPerRB; % 每个SSB的子载波数
    CPLengths = []; % 循环前缀长度

    % 根据子载波间隔确定OFDM参数
    switch SubcarrierSpacing
        case 15
            SampleRate = 30.72e6; % 采样率 (Hz)
            CPLengths = [160; 144]; % 循环前缀长度
        case 30
            SampleRate = 61.44e6; % 采样率 (Hz)
            CPLengths = [320; 288]; % 循环前缀长度
        case 60
            SampleRate = 122.88e6; % 采样率 (Hz)
            CPLengths = [640; 576]; % 循环前缀长度
        otherwise
            error('Unsupported subcarrier spacing.');
    end

    % 生成SSB的频域资源网格（增加DM-RS映射）
    ssbGrid = complex(zeros(NumSubcarriersPerSSB, NumSymbolsPerSSB));
    dmrsGrid = zeros(size(ssbGrid)); % 用于标记DM-RS位置

    % 生成PSS并记录位置
    pssSymbols = generatePSS(NCellID);
    ssbGrid(1:127, 1) = pssSymbols;

    % 生成SSS并记录位置
    sssSymbols = generateSSS(NCellID);
    ssbGrid(1:127, 3) = sssSymbols;

    % 生成PBCH及其DM-RS
    [pbchSymbols, dmrsIndices] = generatePBCH(NCellID, NumSubcarriersPerSSB);
    ssbGrid(:, [2,4]) = pbchSymbols;
    dmrsGrid(dmrsIndices, [2,4]) = 1; % 标记DM-RS位置

    % 将资源网格信息存入ssbInfo
    ssbInfo = struct();
    ssbInfo.ssbGrid = ssbGrid;
    ssbInfo.dmrsGrid = dmrsGrid;
    ssbInfo.NumSymbolsPerSSB = NumSymbolsPerSSB;
    ssbInfo.NumSubcarriersPerSSB = NumSubcarriersPerSSB;

    % OFDM调制
    for symbolIdx = 1:NumSymbolsPerSSB
        % 获取当前符号的循环前缀长度
        if symbolIdx == 1
            CP = CPLengths(1);
        else
            CP = CPLengths(2);
        end

        % 将频域信号映射到OFDM子载波
        ofdmSymbol = zeros(Nfft, 1);
        ofdmSymbol(SubcarrierOffset + (1:NumSubcarriersPerSSB)) = ssbGrid(:, symbolIdx);

        % IFFT变换
        timeDomainSymbol = ifft(fftshift(ofdmSymbol), Nfft);

        % 添加循环前缀
        timeDomainSymbol = [timeDomainSymbol(end-CP+1:end); timeDomainSymbol];

        % 将符号添加到SS突发波形中
        ssbWaveform((symbolIdx-1)*(Nfft+CP) + (1:Nfft+CP)) = timeDomainSymbol;
    end

    % 生成SSB信息结构
    ssbInfo.NDLRB = NDLRB;
    ssbInfo.SubcarrierSpacing = SubcarrierSpacing;
    ssbInfo.NCellID = NCellID;
    ssbInfo.SSBTransmitted = SSBTransmitted;
    ssbInfo.SSBPeriodicity = SSBPeriodicity;
    ssbInfo.SubcarrierOffset = SubcarrierOffset;
    ssbInfo.SampleRate = SampleRate;
end

% 生成PSS（主同步信号）
function pssSymbols = generatePSS(NCellID)
    % PSS序列生成（基于5G NR标准）
    n = 0:126; % PSS序列长度
    switch mod(NCellID, 3)
        case 0
            pssSymbols = exp(-1j * pi * n .* (n + 1) / 127);
        case 1
            pssSymbols = exp(-1j * pi * n .* (n + 1) / 127) .* exp(-1j * 2 * pi / 3);
        case 2
            pssSymbols = exp(-1j * pi * n .* (n + 1) / 127) .* exp(-1j * 4 * pi / 3);
    end
end

% 生成SSS（辅同步信号）
function sssSymbols = generateSSS(NCellID)
    % SSS序列生成（基于5G NR标准）
    n = 0:126; % SSS序列长度
    m0 = mod(NCellID, 31);
    m1 = floor(NCellID / 31);
    sssSymbols = exp(-1j * pi * (m0 * n .* (n + 1) / 127 + m1 * (n + 1) .* (n + 2) / 127));
end

% 修改后的PBCH生成函数（包含DM-RS）
function [pbchSymbols, dmrsIndices] = generatePBCH(NCellID, NumSubcarriers)
    % 生成PBCH QPSK符号
    pbchSymbols = (1 + 1j)/sqrt(2) * (2*randi([0 1], NumSubcarriers, 2) - 1);
    
    % 生成DM-RS位置（每4个子载波一个DM-RS）
    dmrsIndices = 1:4:NumSubcarriers;
    dmrsIndices = dmrsIndices(2:end-1); % 排除边缘
end

% 修改后的绘图函数
function plotSSBurst(ssbWaveform, ssbInfo)
    % 创建资源网格可视化
    figure;
    
    % 准备数据
    ssbGrid = abs(ssbInfo.ssbGrid);       % 主信号能量
    dmrsMask = ssbInfo.dmrsGrid;          % DM-RS位置
    [nSubcarriers, nSymbols] = size(ssbGrid);
    
    % 创建坐标网格
    [X,Y] = meshgrid(1:nSymbols, 1:nSubcarriers);
    
    % 绘制主信号
    imagesc(ssbGrid);
    colormap('parula');
    hold on;
    
    % 标记DM-RS位置（黄色）
    scatter(X(dmrsMask==1), Y(dmrsMask==1), 40, 'y', 'filled');
    
    % 标记PSS位置（红色）
    rectangle('Position',[0.5 0.5 1 127], 'EdgeColor','r', 'LineWidth',2)
    
    % 标记SSS位置（绿色）
    rectangle('Position',[2.5 0.5 1 127], 'EdgeColor','g', 'LineWidth',2)
    
    % 标注文字
    text(0.8, 60, 'PSS', 'Color','w', 'FontWeight','bold', 'Rotation',90)
    text(2.8, 60, 'SSS', 'Color','w', 'FontWeight','bold', 'Rotation',90)
    text(1.5, 180, 'PBCH', 'Color','w', 'FontWeight','bold', 'HorizontalAlignment','center')
    text(1.5, 200, 'DM-RS', 'Color','k', 'FontWeight','bold', 'HorizontalAlignment','center')
    
    % 坐标轴设置
    ylabel('Subcarriers');
    xlabel('OFDM Symbols');
    title(sprintf('SS Burst, SCS=%dkHz, NDLRB=%d',...
         ssbInfo.SubcarrierSpacing, ssbInfo.NDLRB));
    
    % 调整坐标轴范围
    xlim([0.5 nSymbols+0.5]);
    ylim([0.5 nSubcarriers+0.5]);
    
    % 添加色标
    colorbar;
    hold off;
end