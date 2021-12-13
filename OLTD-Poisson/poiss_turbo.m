%% =============================
% Turbo Equalization with MAP
% ==============================
profile on
clear
warning off
%% 1.set up the parameters
TestBits = 500;    % (bit)
M = 2; % the modulate order
ISIh = 0.1;0.1:1;
hLen = length(ISIh);
snrdB = 10:2:16;    % (dB)
MonteCalo = 10;
sequence = randperm(TestBits*2);
Pak = 1/2;
iternum = 3;

%% 2.compute the BER/snr & monte carlo for the BER
turErrBitNN = zeros(length(snrdB), 3);
turErrBitNN1 = zeros(length(snrdB), iternum);
fprintf([ '\n Simulation begins at ', datestr(now), '\n']);
for ss = 1:hLen
    ss
    channelH = sqrt(exp(-ISIh(ss)*(0:1))/sum(exp(-ISIh(ss)*(0:1))));
    akTest = [round(rand(TestBits-3,1));zeros(3,1)];
    codedBits = FEC_Coding(akTest,2);
    interleavedBits = codedBits(sequence);
    
    for ii = 1:length(snrdB)
        fprintf([ '\n', 'SNR = %ddB ', datestr(now), '\n'], snrdB(ii));
        
        Net =  importKerasNetwork(strcat('./model-500/dnn_model', num2str(ISIh(ss)), '-',num2str(snrdB(ii)),'dB.h5'));
        
        
        for mm = 1 : MonteCalo  % Monte Carlo
            vk = db2mag(snrdB(ii))*filter(channelH, 1, interleavedBits) + 1;  % signal go through the channel +randn(1,2)*sqrt(0.01)
            y = poissrnd(vk);
            turErrBitNN = ApplyNN(Net, y, sequence, iternum, turErrBitNN, akTest,ii, "OLTD"); % algo :algorithm type
            turErrBitNN1 = ApplyNN(Net, y, sequence, iternum, turErrBitNN1, akTest,ii, "BCJRNet"); 
            
        end
        
    end
end
BER = turErrBitNN/(MonteCalo*TestBits*hLen);
BER1= turErrBitNN1/(MonteCalo*TestBits*hLen);
%% 3.plot the curve of the BER-SNR
% Legend = [];
fig1 = figure;
set(fig1, 'WindowStyle', 'docked');
%% plot
PlotType = strvcat( '-rs', '-bo', '-mx', '-mv', '-b');
Legend1 = strvcat('OLTD-based turbo 0 iter', 'OLTD-based turbo 1 iter', 'OLTD-based turbo 2 iter');
for aa=1:3
%     Legend = strvcat(Legend,  NameProts(aa,:));
    semilogy(snrdB, BER(:,aa), PlotType(aa,:),'LineWidth',1,'MarkerSize',6);
    hold on;
end

PlotType1 = strvcat( '--rs', '--bo', '--mx', '--mv', '--b^');
Legend2 = strvcat('BCJRNet-based turbo 0iter', 'BCJRNet-based turbo 1iter', 'BCJRNet-based turbo 2iter');
for bb=1:3
%     Legend = strvcat(Legend,  NameProts1(bb,:));
    semilogy(snrdB, BER1(:,bb), PlotType1(bb,:),'LineWidth',1,'MarkerSize',6);
    hold on;
end
Legend = strvcat(Legend1(1:iternum,:),Legend2(1:iternum,:));
xlabel('SNR [dB]');
ylabel('SER');

grid on;
legend(Legend,'Location','SouthWest');