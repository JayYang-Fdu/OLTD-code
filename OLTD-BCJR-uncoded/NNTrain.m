clear
warning("off")
alpha = rand; %should be the same as the pilot
snrdB =8:2:20;
MonteCalo = 1000;
M = 4; %modulation qpsk
channelTap = 2;
% TrainBits = 1002;
TestBits = 1000+log2(M);
% SER = zeros(1, length(snrdB));
ISIh = 0.1:1;
v_nCurves = [...          % Curves
    1 ...                  
    1 ...                  
    ];

% m_fSERAvg = zeros(length(v_nCurves),length(snrdB));
algotype =  " BCJR";"Viterbi";%
s_nCurves = length(v_nCurves);
v_stProts = strvcat(  ...
    strcat('OLTD-based ' ,algotype), ...
    strcat('Model-based ',algotype));
%% 
BER = zeros(length(v_nCurves),length(snrdB));
for Hidx = 1:length(ISIh)
    Hidx
    channelH = sqrt(exp(-ISIh(Hidx)*(0:1))/sum(exp(-ISIh(Hidx)*(0:1))));
    % generate the test data
    akTest = round(rand(TestBits,1));
    xkTest = qammod(akTest, M, 'InputType', 'bit', 'UnitAveragePower', true);
    % generate co-channel interference
    akInf1 = randi([0 3],length(xkTest),1);
    xkInf = pammod(akInf1,4);
    interference = xkInf*exp(1j*alpha);
for ii = 1:length(snrdB)
    fprintf([ '\n', 'SNR = %ddB ', datestr(now), '\n'], snrdB(ii));
    if(v_nCurves(1)==1) 
        Net =  importKerasNetwork(strcat('./inteference-500/dnn_model', num2str(ISIh(Hidx)), '-',num2str(snrdB(ii)),'dB.h5'));
    end
    for mm = 1:MonteCalo
        vkTest = filter(channelH, 1, xkTest);
        Gaussion_noise=(randn(size(vkTest)) + 1j*randn(size(vkTest)))/sqrt(2)*db2mag(-snrdB(ii)); %Guassion noise
%         Laplace_nosie = (-sqrt(0.5)*sign(rand(size(vkTest))-0.5).*log(1-2*abs(rand(size(vkTest))-0.5)) ...
%             - 1j*sqrt(0.5)*sign(rand(size(vkTest))-0.5).*log(1-2*abs(rand(size(vkTest))-0.5)))...
%             /sqrt(2)*db2mag(-snrdB(ii)); % Laplace nosie
        
%         pd = makedist("stable","alpha",0.5,"beta",0);    
%         Stable_noise = random(pd,size(vkTest)); % stable-dist noise

%         cauthy_noise = snrdB(ii)*tan(pi*(rand(size(vkTest))-0.5)) + 1j*snrdB(ii)*tan(pi*(rand(size(vkTest))-0.5)); %cauthy_noise
        yTest = vkTest+Gaussion_noise+interference;
        yTest = yTest(2:end);
        if(v_nCurves(1)==1) 
            akHat1 = ApplyLANN(Net, yTest, algotype, M);
            BER(1,ii) = mean(akHat1 ~= akTest(log2(M)+1:end))+BER(1,ii);
        end
        if(v_nCurves(2)==1)
            akHat2 = BCJR(yTest, db2pow(-snrdB(ii)), channelH);
            BER(2,ii) = mean(akHat2 ~= akTest(log2(M)+1:end)) + BER(2,ii);
        end
    end

end
end
m_fSERAvg = BER/(length(ISIh)*MonteCalo);
Legend = [];
fig1 = figure;
set(fig1, 'WindowStyle', 'docked');
%% plot
v_stPlotType = strvcat( '-rs', '-bo', '-mx', '-mv');
for aa=1:s_nCurves
    if (v_nCurves(aa) ~= 0)
        Legend = strvcat(Legend,  v_stProts(aa,:));
        semilogy(snrdB, m_fSERAvg(aa,:), v_stPlotType(aa,:),'LineWidth',1,'MarkerSize',6);
        hold on;
    end
end
xlabel('SNR [dB]');
ylabel('BER');

grid on;
legend(Legend,'Location','SouthWest');
% legend(strcat('LANN-',algotype,' algorithm perfect CSI'), strcat(algotype,' algorithm'));