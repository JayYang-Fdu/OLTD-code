clear

warning("off")
snrdB = 16;
alpha = rand; % should be th same as the pilot
MonteCalo = 1000;
M = 4; %modulation qpsk
channelTap = 2;
% TrainBits = 1002;
TestBits = 1000+log2(M);
% SER = zeros(1, length(snrdB));
ISIh = 0.1:0.2:1;
v_nCurves = [...          % Curves
    1 ...                  
    0 ...
    0 ...
    ];
trainpilots = [50,100,200,500,1000,2000];
% m_fSERAvg = zeros(length(v_nCurves),length(snrdB));
% m_fSERAvg = zeros(length(v_nCurves),length(trainpilots));
algotype = "Viterbi";"BCJR"; %
s_nCurves = length(v_nCurves);
v_stProts = strvcat(  ...
    strcat('OLTD-',algotype, ' perfect CSI'), ...
    'ViterbiNet', ...
    strcat(algotype,' algorithm'));
%%
BER = zeros(length(v_nCurves),length(trainpilots));
for ss = 1:length(trainpilots)
    trainpilots(ss)
for Hidx = 1:length(ISIh)
    Hidx
    channelH = sqrt(exp(-ISIh(Hidx)*(0:1))/sum(exp(-ISIh(Hidx)*(0:1))));
    % generate the test data
    akTest = round(rand(TestBits,1));
    xkTest = qammod(akTest, M, 'InputType', 'bit', 'UnitAveragePower', true);
    akInf1 = randi([0 3],length(xkTest),1);
    xkInf = pammod(akInf1,4);
%     BER = zeros(length(v_nCurves),length(snrdB));
for ii = 1:length(snrdB)
    fprintf([ '\n', 'SNR = %ddB ', datestr(now), '\n'], snrdB(ii));
    if(v_nCurves(1)==1) 
        Net =  importKerasNetwork(strcat('./inteference-',num2str(trainpilots(ss)),'/dnn_model', num2str(ISIh(Hidx)), '-',num2str(snrdB(ii)),'dB.h5'));
    end
    interference = xkInf*exp(1j*alpha);
    for mm = 1:MonteCalo
        vkTest = filter(channelH, 1, xkTest);
        Gaussion_noise=(randn(size(vkTest)) + 1j*randn(size(vkTest)))/sqrt(2)*db2mag(-snrdB(ii)); %Guassion noise
%         Laplace_nosie = (-sqrt(0.5)*sign(rand(size(vkTest))-0.5).*log(1-2*abs(rand(size(vkTest))-0.5)) ...
%             - 1j*sqrt(0.5)*sign(rand(size(vkTest))-0.5).*log(1-2*abs(rand(size(vkTest))-0.5)))...
%             /sqrt(2)*db2mag(-snrdB(ii)); % Laplace nosie
        
%         pd = makedist("stable","alpha",0.5,"beta",0);    
%         Stable_noise = random(pd,size(vkTest)); % stable-dist noise

%         cauthy_noise = snrdB(ii)*tan(pi*(rand(size(vkTest))-0.5)) + 1j*snrdB(ii)*tan(pi*(rand(size(vkTest))-0.5)); %cauthy_noise
        yTest = vkTest+Gaussion_noise + interference;
        yTest = yTest(2:end);
        if(v_nCurves(1)==1) 
            akHat1 = ApplyLANN(Net, yTest, algotype, M);
            BER(1,ss) = mean(akHat1 ~= akTest(log2(M)+1:end))+BER(1,ss);
        end
        if(v_nCurves(2)==1) 
%             GMModelreal = fitgmdist(real(yTest),16,'RegularizationValue',0.01);
%             GMModelimag = fitgmdist(imag(yTest),16,'RegularizationValue',0.01);
            akHat2 = ViterbiNet(Net, yTest, algotype, M);
            BER(2,ss) = mean(akHat2 ~= akTest(log2(M)+1:end))+BER(2,ss);
        end
        if(v_nCurves(3)==1)
            akHat3 = Viterbi(yTest, db2pow(snrdB(ii)), channelH, M);
            BER(3,ss) = mean(akHat3 ~= akTest(log2(M)+1:end)) + BER(3,ss);
        end
    end

end
%     m_fSERAvg = m_fSERAvg + BER/MonteCalo;
end
end
m_fSERAvg = BER/(length(ISIh)*MonteCalo);
v_stLegend = [];
fig1 = figure;
set(fig1, 'WindowStyle', 'docked');
%% plot
v_stPlotType = strvcat( '--rs', '--ro', '--mx', 'mv', '-.b^');
for aa=1:s_nCurves
    if (v_nCurves(aa) ~= 0)
        v_stLegend = strvcat(v_stLegend,  v_stProts(aa,:));
        semilogy(trainpilots, m_fSERAvg(aa,:), v_stPlotType(aa,:),'LineWidth',1,'MarkerSize',10);
        hold on;
    end
end
xlabel('Pilots length');
ylabel('BER');
% title("500")
grid on;
legend(v_stLegend,'Location','SouthWest');
% legend(strcat('LANN-',algotype,' algorithm perfect CSI'), strcat(algotype,' algorithm'));