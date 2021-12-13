clear
warning("off")
snrdB =[0.001,0.01,0.05,0.1:0.1:1];
MonteCalo = 100;
M = 4; %modulation qpsk
channelTap = 2;
% TrainBits = 1002;
TestBits = 10000+log2(M);
SER = zeros(1, length(snrdB));
ISIh = 0.1:0.2:1;
v_nCurves = [...          % Curves
    0 ...
    1 ...
    ];

% m_fSERAvg = zeros(length(v_nCurves),length(snrdB));
algotype = "Viterbi"; "BCJR";%
s_nCurves = length(v_nCurves);
v_stProts = strvcat(  ...
    strcat('OLTD-',algotype, ' perfect CSI'), ...
    strcat('Model-based ',algotype));
%%
BER = zeros(length(v_nCurves),length(snrdB));
for Hidx = 1:length(ISIh)
    Hidx
    channelH = sqrt(exp(-ISIh(Hidx)*(0:1))/sum(exp(-ISIh(Hidx)*(0:1))));
    % generate the test data
    akTest = round(rand(TestBits,1));
    xkTest = qammod(akTest, M, 'InputType', 'bit', 'UnitAveragePower', true);
    
    
    for ii = 1:length(snrdB)
        fprintf([ '\n', 'SNR = %ddB ', datestr(now), '\n'], snrdB(ii));
        if(v_nCurves(1)==1)
            Net =  importKerasNetwork(strcat('./cauchy_model-500/dnn_model', num2str(ISIh(Hidx)), '-',num2str(snrdB(ii)),'dB.h5'));
        end
        for mm = 1:MonteCalo
            vkTest = filter(channelH, 1, xkTest);
            %         Gaussion_noise=(randn(size(vkTest)) + 1j*randn(size(vkTest)))/sqrt(2)*db2mag(-snrdB(ii)); %Gaussion noise
            %         Laplace_nosie = (-sqrt(0.5)*sign(rand(size(vkTest))-0.5).*log(1-2*abs(rand(size(vkTest))-0.5)) ...
            %             - 1j*sqrt(0.5)*sign(rand(size(vkTest))-0.5).*log(1-2*abs(rand(size(vkTest))-0.5)))...
            %             /sqrt(2)*db2mag(-snrdB(ii)); % Laplace nosie
            
            
            cauchy_noise = snrdB(ii)*tan(pi*(rand(size(vkTest))-0.5)) + 1j*snrdB(ii)*tan(pi*(rand(size(vkTest))-0.5)); %cauchy_noise
            yTest = vkTest+cauchy_noise;
            yTest = yTest(2:end);
            if(v_nCurves(1)==1)
                akHat1 = ApplyLANN(Net, yTest, algotype, M);
                BER(1,ii) = mean(akHat1 ~= akTest(log2(M)+1:end))+BER(1,ii);
            end
            if(v_nCurves(2)==1)
                if algotype=="Viterbi"
                    akHat2 = Viterbi(yTest, db2pow(snrdB(ii)), channelH, M, snrdB(ii));
                end
                if algotype=="BCJR"
                    akHat2 = BCJR(yTest, db2pow(-snrdB(ii)), channelH,snrdB(ii));
                end
                BER(2,ii) = mean(akHat2 ~= akTest(log2(M)+1:end)) + BER(2,ii);
            end
        end
        
    end
    %     m_fSERAvg = m_fSERAvg + BER/MonteCalo;
end
m_fSERAvg = BER/(MonteCalo*length(ISIh));
v_stLegend = [];
fig1 = figure;
set(fig1, 'WindowStyle', 'docked');
%% plot
v_stPlotType = strvcat( '-rs', '--ro', 'mx', 'mv', '-.b^');
for aa=1:s_nCurves
    if (v_nCurves(aa) ~= 0)
        v_stLegend = strvcat(v_stLegend,  v_stProts(aa,:));
        semilogy(-log10(snrdB), m_fSERAvg(aa,:), v_stPlotType(aa,:),'LineWidth',1,'MarkerSize',6);
        hold on;
    end
end
xlabel('SNR [dB]');
ylabel('BER');
grid on;
legend(v_stLegend,'Location','SouthWest');
% legend(strcat('LANN-',algotype,' algorithm perfect CSI'), strcat(algotype,' algorithm'));