function standardizedData = standardizeData(data,targetMean,targetStd)
    % data: 输入数据
    % targetMean: 目标均值
    % targetStd: 目标标准差

    % 计算数据的原始均值和标准差
    meanData = nanmean(data);
    stdData = nanstd(data);
    

    % 先进行原始数据的标准化（均值为0，标准差为1）
    standardizedData = (data - meanData) ./ stdData;

    % 将标准化后的数据调整到目标均值和标准差
    standardizedData = standardizedData * targetStd + targetMean;
end