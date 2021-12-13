clear
warning("off")
snrdB =-8:-2;
MonteCalo = 100;
alpha = 0.5;
ModulationType = 2; 2;% 2:S=2;8:s=8;
if ModulationType == 2
    ModType = 2;
elseif ModulationType == 8
    ModType = 8;
end
h = 0.5;
tau = 0;
TestBits = 512;
ISIh = 0.1:1;
v_nCurves = [...          % Curves
    1 ...                  
    1 ...                  
    ];
iternum = 3;
% m_fSERAvg = zeros(length(snrdB), iternum);
sequence = randperm(ModType*TestBits);
% sequence = interleaving(1:TestBits*2,16);
algotype = "Turbo";"BCJR";"Viterbi"; %

%%
turErrBitNN = zeros(length(snrdB), iternum);
turErrBit = zeros(length(snrdB), iternum);
for Hidx = 1:length(ISIh)
    ISIh(Hidx)
    channelH = exp(1j*ISIh(Hidx)*2*pi);
    % generate the test data
    akTest = [round(rand(TestBits-3,1));zeros(3,1)];
    
    codedBits = FEC_Coding(akTest,ModulationType);
    interleavedBits = codedBits(sequence);
    
    
    xkTest = gfsk_modulation(interleavedBits,h,tau);

for ii = 1:length(snrdB)
    fprintf([ '\n', 'SNR = %ddB ', datestr(now), '\n'], snrdB(ii));
    if (v_nCurves(1)==1) 
        Net =  importKerasNetwork(strcat('./S2_Model/dnn_model', num2str(ISIh(Hidx)), '-',num2str(snrdB(ii)),'dB.h5'));


    end
    vkTest = filter(channelH, 1, xkTest);
    for mm = 1:MonteCalo
        
        Gaussion_noise=(randn(size(vkTest)) + 1j*randn(size(vkTest)))/(sqrt(2)/db2mag(-snrdB(ii))); %Guassion noise 
        yTest = vkTest + Gaussion_noise ;
        if (v_nCurves(1)==1)
            turErrBitNN = ApplyOLTD(Net, yTest, algotype, sequence, iternum, turErrBitNN, akTest, ii);
        end
        if (v_nCurves(2)==1)
            turErrBit = Turbo_iter(yTest, sequence,channelH, iternum, turErrBit, akTest, db2pow(-snrdB(ii)),ii);
        end
    end

end
end
m_fSERAvg = turErrBitNN / (MonteCalo*TestBits*length(ISIh));
m_fSERAvg1 = turErrBit / (MonteCalo*TestBits*length(ISIh));

%% plot
fig1 = figure;
v_stLegend = [];
set(fig1, 'WindowStyle', 'docked');
v_stPlotType = strvcat( '-rs', '-bo', '-mx', '-mx', '-bo','-rs');
if (v_nCurves(1) == 1)
    v_stLegend1 = strvcat( "OLTD-based turbo 0iter","OLTD-based turbo 1iter","OLTD-based turbo 2iter");
    for aa = 1:iternum 
        semilogy(snrdB, m_fSERAvg(:,aa), v_stPlotType(aa,:),'LineWidth',1,'MarkerSize',6);
        hold on;
    end
    v_stLegend = strvcat(v_stLegend,v_stLegend1);
end
% legend(v_stLegend1,'Location','SouthWest');
v_stPlotType1 = strvcat( '--rs', '--bo', '--mx', '--mx', '--bo','--rs');
if (v_nCurves(2) == 1)
    v_stLegend2 = strvcat( "Model-based turbo 0iter","Model-based turbo 1iter","Model-based turbo 2iter");
    for aa = 1:iternum
        semilogy(snrdB, m_fSERAvg1(:,aa), v_stPlotType1(aa,:),'LineWidth',1,'MarkerSize',6);
        hold on;
    end
    v_stLegend = strvcat(v_stLegend,v_stLegend2);
end
legend(v_stLegend,'Location','SouthWest');
xlabel('SNR(dB)')
ylabel('BER')
grid on
