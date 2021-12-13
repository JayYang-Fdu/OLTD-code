clear
% snrdB = [0.001,0.01,0.05,0.1:0.1:1];
snrdB = 16;
alpha = rand;
MonteCalo = 1;
M = 4; %modulation:4:qpsk;2:BPSK
channelTap = 2;
TrainBits = 400+log2(M);
SER = zeros(1, length(snrdB));
ISIh = 0.1:0.2:1;
numClasses = M^channelTap;
traindata = zeros((TrainBits-log2(M))/log2(M), 3);
% s_fFrameSize = 50;
%%
for Hidx = 1:length(ISIh)
    channelH = sqrt(exp(-ISIh(Hidx)*(0:1))/sum(exp(-ISIh(Hidx)*(0:1))));
    Hidx
    channelH
     % generate the train data
    akTrain = round(rand(TrainBits,1));
    xkTrain = qammod(akTrain, M, 'InputType', 'bit', 'UnitAveragePower', true);
%     vkTrain = filter(channelH + sqrt(0.05)*randn(size(channelH)), 1, xkTrain);
    vkTrain = filter(channelH, 1, xkTrain);
    
    akInf1 = randi([0 3],length(vkTrain),1);
    xkInf = pammod(akInf1,4);
    % generate the train data
    for ii = 1:length(snrdB)
        fprintf([ '\n', 'SNR = %ddB ', datestr(now), '\n'], snrdB(ii));
        interference = xkInf*exp(1j*alpha);
        Gaussion_noise = (randn(size(vkTrain)) + 1j*randn(size(vkTrain)))/sqrt(2)*db2mag(-snrdB(ii)); % Gaussian nosie
%         Laplace_nosie = (-sqrt(0.5)*sign(rand(size(vkTrain))-0.5).*log(1-2*abs(rand(size(vkTrain))-0.5)) ...
%             - 1j*sqrt(0.5)*sign(rand(size(vkTrain))-0.5).*log(1-2*abs(rand(size(vkTrain))-0.5)))...
%             /sqrt(2)*db2mag(-snrdB(ii)); % Laplace nosie -
%         pd = makedist("stable","alpha",0.5,"beta",0);    
%         Stable_noise = random(pd,size(vkTrain)); % stable-dist noise 
%         cauchy_noise = snrdB(ii)*tan(pi*(rand(size(vkTrain))-0.5)) + 1j*snrdB(ii)*tan(pi*(rand(size(vkTrain))-0.5)); %cauthy_noise
        %% CSI perfect
        yTrain = vkTrain + Gaussion_noise + interference;     
        yTrain = yTrain(2:end);
        label = GenLabel(akTrain, yTrain, M);
        traindata(1:(TrainBits-log2(M))/log2(M), 1) = real(yTrain);
        traindata(1:(TrainBits-log2(M))/log2(M), 2) = imag(yTrain);
        traindata(1:(TrainBits-log2(M))/log2(M), 3) = label;

%         csvwrite(strcat('./train/traindata', num2str(ISIh(Hidx)), '-',num2str(snrdB(ii)),'dB.csv'),traindata);
        csvwrite(strcat('./train/validdata', num2str(ISIh(Hidx)), '-',num2str(snrdB(ii)),'dB.csv'),traindata);
    end
    
   
end

