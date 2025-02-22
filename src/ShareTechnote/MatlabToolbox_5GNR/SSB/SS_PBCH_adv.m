% 参数设置
gNB = struct();
gNB.NDLRB = 100;
gNB.SubcarrierSpacing = 30;
gNB.NCellID = 1;

ssburst = struct();
ssburst.SSBTransmitted = [1 1 1 1];
ssburst.SSBPeriodicity = 10;
ssburst.SubcarrierOffset = 0;
ssburst.DisplayBurst = true;

% 生成波形并绘图
[ssbWaveform, ssbInfo] = generateSSBurst(gNB, ssburst);
plotSSBurst(ssbWaveform, ssbInfo);

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

    % 初始化SS突发波形
    ssbWaveform = complex(zeros(NumSymbolsPerSSB * (Nfft + CPLengths(1)), 1));

    % 生成SSB的频域资源网格
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
    ssbInfo.ssbGrid = ssbGrid;
    ssbInfo.dmrsGrid = dmrsGrid;
    ssbInfo.NumSymbolsPerSSB = NumSymbolsPerSSB;
    ssbInfo.NumSubcarriersPerSSB = NumSubcarriersPerSSB;
    ssbInfo.SubcarrierSpacing = SubcarrierSpacing;
    ssbInfo.NDLRB = NDLRB;

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

% 生成PBCH（物理广播信道）
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
    colormap('jet'); % 使用jet颜色映射
    hold on;
    
    % 标记DM-RS位置（黄色）
    scatter(X(dmrsMask==1), Y(dmrsMask==1), 40, 'y', 'filled');
    
    % 标记PSS位置（红色）
    rectangle('Position',[0.5 0.5 1 127], 'EdgeColor','r', 'LineWidth',2)
    
    % 标记SSS位置（绿色）
    rectangle('Position',[2.5 0.5 1 127], 'EdgeColor','g', 'LineWidth',2)
    
    % 标注文字
    text(0.8, 60, 'PSS', 'Color','w', 'FontWeight','bold', 'Rotation',90, 'FontSize',12)
    text(2.8, 60, 'SSS', 'Color','w', 'FontWeight','bold', 'Rotation',90, 'FontSize',12)
    text(1.5, 180, 'PBCH', 'Color','w', 'FontWeight','bold', 'HorizontalAlignment','center', 'FontSize',12)
    text(1.5, 200, 'DM-RS', 'Color','k', 'FontWeight','bold', 'HorizontalAlignment','center', 'FontSize',12)
    
    % 坐标轴设置
    ylabel('Subcarriers', 'FontSize',12);
    xlabel('OFDM Symbols', 'FontSize',12);
    title(sprintf('SS Burst, SCS=%dkHz, NDLRB=%d',...
         ssbInfo.SubcarrierSpacing, ssbInfo.NDLRB), 'FontSize',14);
    
    % 调整坐标轴范围
    xlim([0.5 nSymbols+0.5]);
    ylim([0.5 nSubcarriers+0.5]);
    
    % 添加网格线
    grid on;
    set(gca, 'GridColor', 'k', 'GridAlpha', 0.5, 'LineWidth', 0.5);
    
    % 添加色标
    colorbar;
    hold off;
end