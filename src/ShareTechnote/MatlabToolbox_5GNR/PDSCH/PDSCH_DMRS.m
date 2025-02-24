% 5G NR 资源网格可视化（不依赖 MATLAB Toolbox）
% 作者：通信助手
% 功能：生成包含 PDSCH/DM-RS/SSB 的资源网格图

%% 参数配置
clear; close all;

% 全局参数
gNB.NDLRB = 20;                % 下行资源块数量 (6~275)
gNB.SubcarrierSpacing = 30;    % 子载波间隔 (kHz): 15,30,60,120,240
gNB.WaveformType = 'CP-OFDM';  % 波形类型: 'CP-OFDM'
gNB.CyclicPrefix = 'Normal';   % 循环前缀类型: 'Normal'
gNB.NCellID = 1;               % 物理小区 ID (0~1007)

% PDSCH 参数
pdsch.PRBSet = 0:0;            % PRB 分配 (0-based)
pdsch.SymbolSet = 0:13;        % 符号分配 (0-based)
pdsch.Modulation = '16QAM';    % 调制方式
pdsch.NLayers = 2;             % 传输层数
pdsch.PortSet = [0, 2];        % 天线端口
pdsch.DL_DMRS_typeA_pos = 2;   % DM-RS 起始符号位置 (2或3)
pdsch.DL_DMRS_config_type = 2; % DM-RS 配置类型 (1或2)
pdsch.DL_DMRS_max_len = 1;     % DM-RS 符号长度 (1或2)

% SSB 参数
ssb.BurstType = 'CaseB';       % SSB 突发类型 (30kHz)
ssb.SSBTransmitted = [0 0 0 0 0 0 0 0]; % SSB 传输位图 (全关闭)
ssb.SubcarrierOffset = 0;      % SSB 频域偏移

%% 手动实现 OFDM 参数计算
NSubcarriers = gNB.NDLRB * 12;    % 总子载波数
SymbolsPerSlot = 14;              % 每时隙符号数
FFTSize = 2^ceil(log2(NSubcarriers)); % FFT 大小
SampleRate = gNB.NDLRB * 12 * gNB.SubcarrierSpacing * 1e3; % 采样率

%% 生成资源网格
grid = zeros(NSubcarriers, SymbolsPerSlot, pdsch.NLayers);

% 生成 PDSCH 资源索引
[pdschIndices, pdschInfo] = generatePDSCHIndices(gNB, pdsch);

% 生成 DM-RS 资源索引
dmrsIndices = generateDMRSIndices(gNB, pdsch);

% 生成 SSB 资源网格
ssbGrid = generateSSBGrid(gNB, ssb);

% 映射资源到网格
grid = mapResourcesToGrid(grid, pdschIndices, dmrsIndices, ssbGrid);

%% 可视化
plotResourceGrid(grid, gNB, pdsch, ssb);

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
                % 使用循环生成每个子载波的索引
                for sc = subcarriers
                    linearIndices = sub2ind([NSubcarriers, SymbolsPerSlot, pdsch.NLayers], ...
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
    
    % DM-RS 参数
    dmrsSymbols = pdsch.DL_DMRS_typeA_pos + (0:pdsch.DL_DMRS_max_len-1);
    dmrsSymbols = dmrsSymbols(dmrsSymbols < SymbolsPerSlot);
    
    indices = [];
    for layer = 0:pdsch.NLayers-1
        port = pdsch.PortSet(layer+1);
        for sym = dmrsSymbols
            for prb = pdsch.PRBSet
                % 根据 3GPP 38.211 7.4.1.1 生成 DM-RS 位置
                k = mod(port, 2) + (0:2:11); % 类型2：每个PRB 6个RE
                subcarriers = prb*12 + k;
                % 使用循环生成每个子载波的索引
                for sc = subcarriers
                    linearIndices = sub2ind([NSubcarriers, SymbolsPerSlot, pdsch.NLayers], ...
                                          sc + 1, sym + 1, layer + 1);
                    indices = [indices; linearIndices];
                end
            end
        end
    end
end

% ========================== 生成 SSB 资源网格 ==========================
function ssbGrid = generateSSBGrid(gNB, ssb)
    NSubcarriers = gNB.NDLRB * 12;
    SymbolsPerSlot = 14;
    ssbGrid = zeros(NSubcarriers, SymbolsPerSlot);
    
    % SSB 位置参数（简化实现）
    ssbStartSymbol = 2; % 假设 SSB 起始符号
    ssbSubcarriers = (0:239) + ssb.SubcarrierOffset + 1; % 假设 SSB 占240个子载波
    
    % 标记 SSB 位置（用随机数据模拟）
    validSubcarriers = ssbSubcarriers(ssbSubcarriers <= NSubcarriers);
    ssbGrid(validSubcarriers, ssbStartSymbol + (0:3)) = rand(numel(validSubcarriers),4);
end

% ========================== 资源映射 ==========================
function grid = mapResourcesToGrid(grid, pdschIndices, dmrsIndices, ssbGrid)
    % 映射 PDSCH（值=2）
    grid(pdschIndices) = 2;
    
    % 映射 DM-RS（值=3）
    grid(dmrsIndices) = 3;
    
    % 映射 SSB（值=1）
    ssbMask = abs(ssbGrid) > 0;
    for layer = 1:size(grid,3)
        grid(:,:,layer) = grid(:,:,layer) + ssbMask;
    end
end

% ========================== 可视化 ==========================
function plotResourceGrid(grid, gNB, pdsch, ssb)
    figure('Name','5G NR 资源网格可视化','NumberTitle','off');
    colormap([1 1 1; 0.7 0.7 1; 1 0.7 0.7; 0.7 1 0.7]); % 白|蓝|红|绿
    
    for layer = 1:pdsch.NLayers
        subplot(2,2,layer*2-1);
        imagesc(grid(:,:,layer));
        axis xy;
        title(sprintf('Layer %d - 全资源网格', layer-1));
        xlabel('OFDM 符号');
        ylabel('子载波');
        set(gca,'XTick',1:14,'XTickLabel',0:13);
        
        subplot(2,2,layer*2);
        imagesc(grid(:,:,layer));
        axis xy;
        xlim([0.5 14.5]);
        ylim([pdsch.PRBSet(1)*12+0.5, (pdsch.PRBSet(end)+1)*12+0.5]);
        title(sprintf('Layer %d - 缩放视图', layer-1));
        xlabel('OFDM 符号');
        ylabel('子载波');
        set(gca,'XTick',1:14,'XTickLabel',0:13);
    end
    
    % 创建图例
    % h = zeros(3,1);
    % h(1) = plot(NaN,NaN,'s','MarkerFaceColor',[0.7 0.7 1],'MarkerEdgeColor','k');
    % h(2) = plot(NaN,NaN,'s','MarkerFaceColor',[1 0.7 0.7],'MarkerEdgeColor','k');
    % h(3) = plot(NaN,NaN,'s','MarkerFaceColor',[0.7 1 0.7],'MarkerEdgeColor','k');
    % legend(h, {'SSB','PDSCH','DM-RS'}, 'Location', 'southoutside','Orientation','horizontal');
    
    set(gcf,'Position',[100 100 800 600]);
    sgtitle(sprintf('5G NR 资源网格 (NDLRB=%d, SCS=%dkHz, Ports=%s)',...
            gNB.NDLRB, gNB.SubcarrierSpacing, num2str(pdsch.PortSet)));
end