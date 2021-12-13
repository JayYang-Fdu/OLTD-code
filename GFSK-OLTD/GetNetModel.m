function [Net,GMModel] = GetNetModel(y,label,numClasses)
yTrain = [real(y), imag(y)].';
v_fYcat = cell(size(yTrain, 1),1);
inputSize = size(yTrain, 1);

numHiddenUnits = 100;
LSTMLayer = lstmLayer(numHiddenUnits,'OutputMode','last'... 
    , 'RecurrentWeightsLearnRateFactor', 0 ...
    , 'RecurrentWeightsL2Factor', 0 ...
    );
LSTMLayer.RecurrentWeights = zeros(4*numHiddenUnits,numHiddenUnits);
layers = [ ...
    sequenceInputLayer(inputSize)
    LSTMLayer
    fullyConnectedLayer(floor(numHiddenUnits/2))
    reluLayer
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];
v_fXcat = categorical(label);
% v_fYcat = num2cell(yTrain');
for ii=1:length(yTrain)
    v_fYcat{ii} = yTrain(:,ii);
end
maxEpochs = 100;
miniBatchSize = 27;
learnRate = 0.01;
options = trainingOptions('adam', ... 
    'ExecutionEnvironment','cpu', ...
    'InitialLearnRate', learnRate, ...
    'LearnRateDropFactor',0.2, ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'Verbose',true);%,'Plots','training-progress'); % This can be unmasked to display training convergence
Net = trainNetwork(v_fYcat,v_fXcat,layers,options);

end