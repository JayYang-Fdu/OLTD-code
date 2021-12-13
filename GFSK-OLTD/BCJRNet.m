clear
snrdB = -2:2:20;
MonteCalo = 10000;
M = 4; %modulation qpsk
% state = 12;
numBits = 1000;
SER = zeros(1, length(snrdB));
ISIh = [0.866 0.5];
for ii = 1:length(snrdB)
    fprintf([ '\n', 'SNR = %ddB ', datestr(now), '\n'], snrdB(ii));
    for mm = 1:MonteCalo
        ak = round(rand(numBits,1));
        xk = qammod(ak, M, 'InputType', 'bit', 'UnitAveragePower', true);
        vk = filter(ISIh, 1, xk);
        noise = (randn(size(vk)) + 1j*randn(size(vk)))/sqrt(2)*db2mag(-snrdB(ii));
        yk = vk + noise;
        Lext = zeros(1,2*length(yk));
        hHat = channelEstISI(snrdB(ii),ISIh); 
%         akHat = BCJR(yk, Lext, db2pow(-snrdB(ii)), hHat);
        akHat = Viterbi(yk, db2pow(-snrdB(ii)), hHat);
        SER(ii)= SER(ii) + sum(akHat ~= ak);
%         xkHat = qammod(LLR, M, 'InputType', 'bit', 'UnitAveragePower', true);
%         SER(ii)= SER(ii) + sum(xkHat ~= xk);
    end 
    SER(ii) = SER(ii) / (MonteCalo*(numBits));
    disp(SER(ii))
end
%% plot the curve of the BER-SNR

figure
semilogy(snrdB,SER,'-*')
% semilogy(snrdB,Errbits_vitb,'-*')
% legend('0 iter','1 iter','2 iter')
% legend('viterb-dec')
title(strcat('ÐÅµÀ¹À¼Æ100samples'));
xlabel('SNR(dB)')
ylabel('BER')
grid on