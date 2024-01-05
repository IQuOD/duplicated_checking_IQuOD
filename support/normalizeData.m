function normalized = normalizeData(data, newMin, newMax)
    % data: 输入数据
    % newMin: 新的最小值（例如 -1）
    % newMax: 新的最大值（例如 1）

    % 计算原始数据的最小值和最大值
    minData = nanmin(data);
    maxData = nanmax(data);
    
    for i=1:length(data(1,:))
        if maxData(i) ~= minData(i)
            % 归一化到 [0, 1]
            normalized(:,i) = (data(:,i) - minData(i)) ./ (maxData(i) - minData(i));

            % 缩放到 [newMin, newMax]
            normalized(:,i) = normalized(:,i) * (newMax - newMin) + newMin;
        else
            normalized(:,i)= minData(i);
        end
    end
end