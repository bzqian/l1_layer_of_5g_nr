% 5G NR 资源网格可视化（不依赖 MATLAB Toolbox）
% 作者：通信助手
% 功能：生成包含 PDSCH/DM-RS/SSB 的资源网格图

%% 参数配置
clear; close all; clc;

% 全局参数
gNB.NDLRB = 20;                % 下行资源块数量 (6~275)
gNB.SubcarrierSpacing = 30;    % 子载波间隔 (kHz): 15,30,60,120,240
gNB.WaveformType = 'CP-OFDM';  % 波形类型: 'CP-OFDM'
gNB.CyclicPrefix = 'Normal';   % 循环前缀类型: 'Normal'
gNB.NCellID = 1;               % 物理小区 ID (0~1007)

% PDSCH 参数
pdsch.PRBSet = 0:19;           % PRB 分配 (0-based)
pdsch.SymbolSet = 0:13;        % 符号分配 (0-based)
pdsch.Modulation = '16QAM';    % 调制方式
pdsch.NLayers = 2;             % 传输层数
pdsch.PortSet = [0, 2];        % 天线端口
pdsch.DL_DMRS_typeA_pos = 2;   % DM-RS 起始符号位置 (2或3)
pdsch.DL_DMRS_config_type = 1; % DM-RS 配置类型 (1或2)
pdsch.DL_DMRS_max_len = 1;     % DM-RS 符号长度 (1或2)

% SSB 参数
ssb.BurstType = 'CaseB';       % SSB 突发类型 (30kHz)
ssb.SSBTransmitted = [1 0 0 0 0 0 0 0]; % SSB 传输位图 (仅第一个SSB激活)
ssb.SubcarrierOffset = 0;      % SSB 频域偏移

%% 手动实现 OFDM 参数计算
NSubcarriers = gNB.NDLRB * 12;    % 总子载波数
SymbolsPerSlot = 14;              % 每时隙符号数
FFTSize = 2^ceil(log2(NSubcarriers)); % FFT 大小
SampleRate = gNB.NDLRB * 12 * gNB.SubcarrierSpacing * 1e3; % 采样率

%% 生成资源网格
resourceGrid = complex(zeros(NSubcarriers, SymbolsPerSlot, pdsch.NLayers));
typeGrid = zeros(NSubcarriers, SymbolsPerSlot, pdsch.NLayers); % 资源类型标记

% 生成 PDSCH 资源索引
[pdschIndices, pdschInfo] = generatePDSCHIndices(gNB, pdsch);

% 生成 PDSCH 调制符号（示例用随机QAM符号）
pdschSymbols = (randn(numel(pdschIndices),1) + 1i*randn(numel(pdschIndices),1))/sqrt(2);

% 生成 DM-RS 资源索引和符号
dmrsIndices = generateDMRSIndices(gNB, pdsch);
dmrsSymbols = generateDMRSSymbols(gNB, pdsch, dmrsIndices);

% 生成 SSB 资源网格（复数符号）
ssbGrid = generateSSBGrid(gNB, ssb);

% 资源映射
resourceGrid(pdschIndices) = pdschSymbols;
typeGrid(pdschIndices) = 2; % PDSCH标记为2

resourceGrid(dmrsIndices) = dmrsSymbols;
typeGrid(dmrsIndices) = 3; % DM-RS标记为3

ssbIndices = find(abs(ssbGrid) > 0);
for layer = 1:pdsch.NLayers
    resourceGrid(ssbIndices + (layer-1)*numel(ssbGrid)) = ssbGrid(ssbIndices);
    typeGrid(ssbIndices + (layer-1)*numel(ssbGrid)) = 1; % SSB标记为1
end

%% 可视化资源映射
plotResourceGrid(typeGrid, gNB, pdsch, ssb);

%% 自定义函数实现
% ========================== 生成 PDSCH 索引 ==========================
function [indices, info] = generatePDSCHIndices(gNB, pdsch)
    NSubcarriers = gNB.NDLRB * 12;
    SymbolsPerSlot = 14;
    
    indices = [];
    for layer = 0:pdsch.NLayers-1
        for sym = pdsch.SymbolSet
            for prb = pdsch.PRBSet
                startSubcarrier = prb * 12;
                subcarriers = startSubcarrier + (0:11);
                for sc = subcarriers
                    linearIndices = sub2ind([NSubcarriers, SymbolsPerSlot, pdsch.NLayers],...
                                      sc + 1, sym + 1, layer + 1);
                    indices = [indices; linearIndices];
                end
            end
        end
    end
    info = struct('PRBSet', pdsch.PRBSet, 'SymbolSet', pdsch.SymbolSet);
end

% ========================== 生成 DM-RS 索引 ==========================
function indices = generateDMRSIndices(gNB, pdsch)
    NSubcarriers = gNB.NDLRB * 12;
    SymbolsPerSlot = 14;
    
    dmrsSymbols = pdsch.DL_DMRS_typeA_pos + (0:pdsch.DL_DMRS_max_len-1);
    dmrsSymbols = dmrsSymbols(dmrsSymbols < SymbolsPerSlot);
    
    indices = [];
    for layer = 0:pdsch.NLayers-1
        port = pdsch.PortSet(layer+1);
        for sym = dmrsSymbols
            for prb = pdsch.PRBSet
                k = mod(port, 2) + (0:2:11); % 类型2：每个PRB 6个RE
                subcarriers = prb*12 + k;
                for sc = subcarriers
                    linearIndices = sub2ind([NSubcarriers, SymbolsPerSlot, pdsch.NLayers],...
                                          sc + 1, sym + 1, layer + 1);
                    indices = [indices; linearIndices];
                end
            end
        end
    end
end

% ========================== 生成 DM-RS 符号 ==========================
function dmrsSymbols = generateDMRSSymbols(gNB, pdsch, dmrsIndices)
    % 根据 3GPP 38.211 7.4.1.1 生成 DM-RS 序列
    c_init = mod(gNB.NCellID * 2^16 + pdsch.PortSet(1) * 2^11 + 1, 2^31);
    r = nrPRBS(c_init, 2*numel(dmrsIndices)); % 生成伪随机序列
    r = 1/sqrt(2) * (1 - 2*r(1:2:end)) + 1i/sqrt(2) * (1 - 2*r(2:2:end)); % QPSK调制
    dmrsSymbols = r;
end

% ========================== 生成 SSB 资源网格 ==========================
function ssbGrid = generateSSBGrid(gNB, ssb)
    NSubcarriers = gNB.NDLRB * 12;
    SymbolsPerSlot = 14;
    ssbGrid = complex(zeros(NSubcarriers, SymbolsPerSlot));
    
    % SSB 位置参数（根据 3GPP 38.211 7.4.3.1）
    ssbStartSymbol = 2; % SSB 起始符号
    ssbSubcarriers = (0:239) + ssb.SubcarrierOffset + 1; % SSB 占240个子载波
    
    % 生成 PSS 序列（3GPP 38.211 7.4.2.2）
    pssSeq = generatePSS(gNB.NCellID);
    
    % 生成 SSS 序列（3GPP 38.211 7.4.2.3）
    sssSeq = generateSSS(gNB.NCellID);
    
    % 生成 PBCH DMRS 序列（3GPP 38.211 7.4.1.4）
    pbchDMRSSeq = generatePBCHDMRS(gNB.NCellID);
    
    % 映射 SSB 到资源网格
    validSubcarriers = ssbSubcarriers(ssbSubcarriers <= NSubcarriers);
    numValidSubcarriers = numel(validSubcarriers);
    
    % 如果 validSubcarriers 超过序列长度，截断或填充序列
    if numValidSubcarriers > 127
        pssSeq = [pssSeq, zeros(1, numValidSubcarriers - 127)]; % 填充 PSS（行向量）
        sssSeq = [sssSeq, zeros(1, numValidSubcarriers - 127)]; % 填充 SSS（行向量）
    end
    if numValidSubcarriers > 144
        pbchDMRSSeq = [pbchDMRSSeq; zeros(2*numValidSubcarriers - 144, 1)]; % 填充 PBCH DMRS
    end
    
    % 映射 PSS（符号2）
    ssbGrid(validSubcarriers, ssbStartSymbol + 0) = pssSeq(1:numValidSubcarriers);
    
    % 映射 SSS（符号3）
    ssbGrid(validSubcarriers, ssbStartSymbol + 1) = sssSeq(1:numValidSubcarriers);
    
    % 映射 PBCH DMRS（符号4和5）
    ssbGrid(validSubcarriers, ssbStartSymbol + 2) = pbchDMRSSeq(1:numValidSubcarriers);
    ssbGrid(validSubcarriers, ssbStartSymbol + 3) = pbchDMRSSeq(numValidSubcarriers+1:2*numValidSubcarriers);
end

% ========================== 生成 PSS 序列 ==========================
function pssSeq = generatePSS(NCellID)
    % 根据 3GPP 38.211 7.4.2.2 生成 PSS 序列
    nID2 = mod(NCellID, 3);
    pssSeq = exp(-1i * pi * nID2 * (0:126) .* ((0:126) + 1) / 127);
end

% ========================== 生成 SSS 序列 ==========================
function sssSeq = generateSSS(NCellID)
    % 根据 3GPP 38.211 7.4.2.3 生成 SSS 序列
    nID1 = floor(NCellID / 3);
    m0 = mod(nID1, 31);
    m1 = mod(m0 + floor(nID1 / 31) + 1, 31);
    
    % 初始化 x0 和 x1 序列
    x0 = zeros(1, 31);
    x1 = zeros(1, 31);
    x0(1) = 1; % x0(0) = 1
    x1(1) = 1; % x1(0) = 1
    
    % 生成 x0 和 x1 序列
    for n = 2:31
        x0(n) = mod(x0(n-1) + x0(max(n-2, 1)), 2); % 避免索引为0
        x1(n) = mod(x1(n-1) + x1(max(n-2, 1)), 2); % 避免索引为0
    end
    
    % 生成 SSS 序列
    sssSeq = (1 - 2 * x0(mod((0:126) + m0, 31) + 1)) .* (1 - 2 * x1(mod((0:126) + m1, 31) + 1));
end

% ========================== 生成 PBCH DMRS 序列 ==========================
function pbchDMRSSeq = generatePBCHDMRS(NCellID)
    % 根据 3GPP 38.211 7.4.1.4 生成 PBCH DMRS 序列
    c_init = mod(NCellID, 2^16);
    r = nrPRBS(c_init, 2*144); % 生成伪随机序列
    pbchDMRSSeq = 1/sqrt(2) * (1 - 2*r(1:2:end)) + 1i/sqrt(2) * (1 - 2*r(2:2:end));
end

% ========================== 可视化 ==========================
function plotResourceGrid(typeGrid, gNB, pdsch, ssb)
    figure('Name','5G NR 资源网格可视化','NumberTitle','off');
    colormap([1 1 1; 0.7 0.7 1; 1 0.7 0.7; 0.7 1 0.7]); % 白|蓝|红|绿
    
    % 确保 typeGrid 的维度正确
    if size(typeGrid, 3) < pdsch.NLayers
        error('typeGrid 的层数与 pdsch.NLayers 不匹配');
    end
    
    % 绘制每个层的资源网格
    for layer = 1:pdsch.NLayers
        % 全资源网格
        subplot(2, pdsch.NLayers, (layer-1)*2 + 1);
        imagesc(typeGrid(:,:,layer));
        axis xy;
        title(sprintf('Layer %d - 全资源网格', layer-1));
        xlabel('OFDM 符号');
        ylabel('子载波');
        set(gca,'XTick',1:14,'XTickLabel',0:13);
        
        % 缩放视图
        subplot(2, pdsch.NLayers, (layer-1)*2 + 2);
        imagesc(typeGrid(:,:,layer));
        axis xy;
        xlim([0.5 14.5]);
        ylim([pdsch.PRBSet(1)*12+0.5, (pdsch.PRBSet(end)+1)*12+0.5]);
        title(sprintf('Layer %d - 缩放视图', layer-1));
        xlabel('OFDM 符号');
        ylabel('子载波');
        set(gca,'XTick',1:14,'XTickLabel',0:13);
    end
    
    % 创建图例
    h = gobjects(3,1); % 创建图形对象数组
    h(1) = plot(NaN,NaN,'s','MarkerFaceColor',[0.7 0.7 1],'MarkerEdgeColor','k');
    hold on;
    h(2) = plot(NaN,NaN,'s','MarkerFaceColor',[1 0.7 0.7],'MarkerEdgeColor','k');
    h(3) = plot(NaN,NaN,'s','MarkerFaceColor',[0.7 1 0.7],'MarkerEdgeColor','k');
    hold off;
    
    % 添加图例
    legend(h, {'SSB','PDSCH','DM-RS'}, 'Location', 'southoutside','Orientation','horizontal');
    
    set(gcf,'Position',[100 100 800 600]);
    sgtitle(sprintf('5G NR 资源网格 (NDLRB=%d, SCS=%dkHz, Ports=%s)',...
            gNB.NDLRB, gNB.SubcarrierSpacing, num2str(pdsch.PortSet)));
end