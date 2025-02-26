function dmrs_simulation()
    % 参数设置
    N_RB = 1;                % 资源块数量
    num_ports = 2;           % 天线端口数
    dmrs_symbols = [2, 11];  % DMRS符号位置（0-based）
    dmrs_additional_position = [0, 1]; % DMRS附加符号位置（0-based，取值范围：0, 1, 2, 3）
    comb_offset = [0, 1];    % 梳状偏移量（0或1）
    cyclic_shift = [0, 3];   % 循环移位量（子载波数）
    c_init = 123;            % 伪随机序列初始化值
    dmrs_type = 1;           % DMRS类型（1 或 2）

    % 检查 dmrs_additional_position 的合法性
    if any(dmrs_additional_position < 0 | dmrs_additional_position > 3)
        error('dmrs_additional_position must be in the range [0, 3]');
    end

    % 派生参数
    N = 6 * N_RB;            % 每个端口的DMRS子载波数
    num_subcarriers = 12 * N_RB; % 总子载波数
    num_symbols = 14;         % 时隙内的OFDM符号数
    
    % 生成伪随机序列
    c = generate_pr_seq(c_init, 2*N);
    
    % 生成基序列
    real_part = 1 - 2*c(1:2:2*N);
    imag_part = 1 - 2*c(2:2:2*N);
    base_seq = (real_part + 1i*imag_part)/sqrt(2);
    
    % 生成端口DMRS序列
    dmrs_sequences = zeros(num_ports, N);
    for port = 1:num_ports
        phase = exp(1i*2*pi*cyclic_shift(port)*(0:N-1)/N);
        dmrs_sequences(port,:) = base_seq .* phase;
    end

    % 应用 OCC
    if dmrs_type == 1
        % Type 1: OCC [1, 1] 和 [1, -1]
        occ = [1, 1; 1, -1];
    elseif dmrs_type == 2
        % Type 2: OCC [1, 1, 1, 1] 和 [1, -1, 1, -1]
        occ = [1, 1, 1, 1; 1, -1, 1, -1];
    else
        error('Unsupported DMRS type');
    end

    % 创建资源网格
    resource_grid = zeros(num_subcarriers, num_symbols);
    
    % 映射到资源网格
    dmrs_subcarriers = zeros(num_ports, N); % 用于记录 DMRS 子载波位置
    for port = 1:num_ports
        if dmrs_type == 1
            % Type 1: 梳状结构
            subcarriers = comb_offset(port) + 2*(0:N-1);
        elseif dmrs_type == 2
            % Type 2: 块状结构
            subcarriers = (port-1)*N + (1:N);
        end
        dmrs_subcarriers(port, :) = subcarriers; % 记录 DMRS 子载波位置
        for sym = dmrs_symbols
            % 应用 OCC
            if dmrs_type == 1
                % Type 1: OCC [1, 1] 和 [1, -1]
                occ = [1, 1; 1, -1];
                for port = 1:num_ports
                    % 将 dmrs_sequences 重塑为 2 x (N/2)
                    reshaped_seq = reshape(dmrs_sequences(port, :), [2, N/2]);
                    % 扩展 occ(port, :) 为 2 x (N/2)
                    occ_expanded = repmat(occ(port, :)', 1, N/2);
                    % 逐元素乘法
                    dmrs_sequences(port, :) = reshape(reshaped_seq .* occ_expanded, 1, N);
                end
            elseif dmrs_type == 2
                % Type 2: OCC [1, 1, 1, 1] 和 [1, -1, 1, -1]
                occ = [1, 1, 1, 1; 1, -1, 1, -1];
                for port = 1:num_ports
                    % 将 dmrs_sequences 重塑为 4 x (N/4)
                    reshaped_seq = reshape(dmrs_sequences(port, :), [4, N/4]);
                    % 扩展 occ(port, :) 为 4 x (N/4)
                    occ_expanded = repmat(occ(port, :)', 1, N/4);
                    % 逐元素乘法
                    dmrs_sequences(port, :) = reshape(reshaped_seq .* occ_expanded, 1, N);
                end
            else
                error('Unsupported DMRS type');
            end
            % 映射 DMRS 到资源网格
            for pos = dmrs_additional_position
                if sym + pos < num_symbols
                    resource_grid(subcarriers+1, sym+pos+1) = dmrs_sequences(port,:);
                end
            end
        end
    end

    % 绘制资源网格
    figure;
    imagesc(abs(resource_grid));
    xlabel('OFDM Symbols');
    ylabel('Subcarriers');
    title('DMRS Resource Grid');
    colorbar;

    % 添加子载波索引
    yticks(1:num_subcarriers);
    yticklabels(arrayfun(@num2str, num_subcarriers-1:-1:0, 'UniformOutput', false));

    % 添加 OFDM 符号索引
    xticks(1:num_symbols);
    xticklabels(arrayfun(@num2str, 0:num_symbols-1, 'UniformOutput', false));

    % 定义颜色
    colors = ['r', 'g', 'b', 'k']; % 为每个端口分配不同颜色

    % 标记 DMRS 符号和非 DMRS 符号
    hold on;
    for port = 1:num_ports
        color = colors(port); % 选择当前端口的颜色
        for sym = dmrs_symbols
            dmrs_sc = dmrs_subcarriers(port, :);
            % 标记 DMRS 符号
            for pos = dmrs_additional_position
                if sym + pos < num_symbols
                    plot(sym+pos+1, dmrs_sc+1, '.', 'Color', color, 'MarkerSize', 10); 
                end
            end
        end
    end
    % 标记非 DMRS 符号
    non_dmrs_symbols = setdiff(0:num_symbols-1, dmrs_symbols + dmrs_additional_position);
    for sym = non_dmrs_symbols
        plot(sym+1, 1:num_subcarriers, '.', 'Color', 'k', 'MarkerSize', 5); % 使用黑色点标记非 DMRS 符号
    end
    hold off;

    % 添加图例
    legend_entries = [arrayfun(@(p) ['Port ', num2str(p)], 1:num_ports, 'UniformOutput', false); repmat({'Non-DMRS Symbols'}, 1, num_ports)];   
    legend(legend_entries{:});
end

function seq = generate_pr_seq(c_init, len)
    % 生成伪随机序列（Gold序列）
    Nc = 1600;               % 初始垃圾比特数
    total_len = Nc + len;    % 总序列长度
    
    % 初始化x1寄存器
    x1 = zeros(1, total_len + Nc + 31);  % 扩展长度以避免索引越界
    x1(1) = 1;  % x1[0]=1, others=0
    
    % 初始化x2寄存器
    x2 = zeros(1, total_len + Nc + 31);
    for i = 0:30
        x2(i+1) = bitand(bitshift(c_init, -i), 1);
    end
    
    % 生成扩展序列
    for n = 31:(total_len + Nc + 31 - 1)
        x1(n+1) = mod(x1(n-31+1) + x1(n-3+1), 2);
        x2(n+1) = mod(x2(n-31+1) + x2(n-3+1) + x2(n-2+1) + x2(n-1+1), 2);
    end
    
    % 生成输出序列
    c = zeros(1, total_len);
    for n = 1:total_len
        c(n) = mod(x1(n+Nc) + x2(n+Nc), 2);
    end
    seq = c(Nc+1:Nc+len);
end