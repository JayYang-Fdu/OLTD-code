clear
snrdB = -8:-2;
theta = 0.1;0.1:0.2:1;
channelTap = 2;
h = 0.5;
alpha = 0.5;
tau = 0;
TrainBits = 32;
ValidBits = 5;
ModulationType = 8; 2;% 2:S=2;8:s=8;
if ModulationType == 2
    ModType = 2;
elseif ModulationType == 8
    ModType = 8;
end
SER = zeros(1, length(snrdB));
sequence = randperm(ModType*TrainBits);
sequence1 = randperm(ModType*ValidBits);
traindata = zeros((TrainBits)*ModType, 3);
validdata = zeros((ValidBits)*ModType, 3);
% s_fFrameSize = 50;
%%
for Hidx = 1:length(theta)
    channelH = exp(1j*theta(Hidx)*2*pi);
    Hidx
    channelH
     % generate the train data
    akTrain = round(rand(TrainBits,1));
    codedBits = FEC_Coding(akTrain,ModType);
    interleavedBits = codedBits(sequence);
    xkTrain = gfsk_modulation(interleavedBits,h,tau);  
    vkTrain = filter(channelH, 1, xkTrain);
    
    akValid = round(rand(ValidBits,1));
    codedBits = FEC_Coding(akValid,ModType);
    interleavedBits1 = codedBits(sequence1);
    xkValid = gfsk_modulation(interleavedBits1,h,tau);  
    vkValid = filter(channelH, 1, xkValid);
    
%     akTrainInf = round(rand(2*length(interleavedBits),1));
%     akTrainInf = randi([0 3],length(interleavedBits),1);
%     xkTrainInf = gfsk_modulation(akTrainInf,h,tau);
%     xkTrainInf = pammod(akTrainInf,4);
%     xkTrainInf = qammod(akTrainInf, 4, 'InputType', 'bit', 'UnitAveragePower', true);
    
%     akValidInf = round(rand(2*length(interleavedBits1),1));
%     akValidInf = randi([0 3],length(interleavedBits1),1);
%     xkValidInf = gfsk_modulation(akValidInf,h,tau);\
%     xkValidInf = qammod(akValidInf, 4, 'InputType', 'bit', 'UnitAveragePower', true);
%     xkValidInf = pammod(akValidInf,4);
    %% generate the train data
    for ii = 1:length(snrdB)
        fprintf([ '\n', 'SNR = %ddB ', datestr(now), '\n'], snrdB(ii));
%         interference = 2*xkInf;
        Gaussion_noise_train = (randn(size(vkTrain)) + 1j*randn(size(vkTrain)))/sqrt(2)*db2mag(-snrdB(ii)); % Gaussian nosie
        Gaussion_noise_valid = (randn(size(vkValid)) + 1j*randn(size(vkValid)))/sqrt(2)*db2mag(-snrdB(ii)); % Gaussian nosie

        %% CSI perfect
        yTrain = vkTrain + Gaussion_noise_train ;  
        yValid = vkValid + Gaussion_noise_valid ;    
%         yTrain = vkTrain ;  
%         yValid = vkValid ; 
        
        label = phaseCompute(2*interleavedBits-1);
        traindata(1:(TrainBits)*ModType, 1) = real(yTrain);
        traindata(1:(TrainBits)*ModType, 2) = imag(yTrain);
        traindata(1:(TrainBits)*ModType, 3) = label;
        
        label1 = phaseCompute(2*interleavedBits1-1);
        validdata(1:(ValidBits)*ModType, 1) = real(yValid);
        validdata(1:(ValidBits)*ModType, 2) = imag(yValid);
        validdata(1:(ValidBits)*ModType, 3) = label1;

        csvwrite(strcat('./train/traindata', num2str(theta(Hidx)), '-',num2str(snrdB(ii)),'dB.csv'),traindata);
        csvwrite(strcat('./train/validdata', num2str(theta(Hidx)), '-',num2str(snrdB(ii)),'dB.csv'),validdata);
    end
    
   
end


