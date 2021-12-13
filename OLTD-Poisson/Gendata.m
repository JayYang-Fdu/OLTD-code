profile on
%% 1.set up the parameters
numBits = 51;    % (bit)
M = 2; % the modulate order
% h = [1 0]; % 
ISIh = 0.1:0.2:1.5;
hLen = length(ISIh);
snrdB = 30:2:38;    % (dB)
runtime = 1;
Pak = 1/2;
SER = zeros(length(snrdB), 1);
traindata = zeros(numBits-1,2);
%% 2.compute the BER/snr & monte carlo for the BER
fprintf([ '\nSimulation begins at ', datestr(now), '\n']);
for ss = 1:hLen
    ss
    channelH = sqrt(exp(-ISIh(ss)*(0:1))/sum(exp(-ISIh(ss)*(0:1))));
for ii =  1 : length(snrdB)    % (dB)
    fprintf([ '\n', 'SNR = %ddB ', datestr(now), '\n'], snrdB(ii));
    for ll = 1 : runtime  % Monte Carlo
        ak = randi([0, 1],numBits, 1);
%         vk = db2mag(snrdB(ii))*filter(ISIh+[randn*sqrt(0.02) randn*sqrt(0.02)], 1, ak) + 1;  % signal go through the channel without noise
        vk = db2mag(snrdB(ii))*filter(channelH, 1, ak) + 1;  % signal go through the channel without noise
        label = gen_label(ak);
        yk = poissrnd(vk(2:end));
        traindata((ll-1)*(numBits-1)+1:(ll)*(numBits-1), 1) = yk;
        traindata((ll-1)*(numBits-1)+1:(ll)*(numBits-1), 2) = label; % the training label
        
%         informationBits((ll-1)*(numBits-1)+1:ll*(numBits-1),ii) = ak(2:end);
%         codedBits((ll-1)*(numBits-1)+1:(ll)*(numBits-1),1) = real(yk);

        
    end 
    csvwrite(strcat('./train/validdata', num2str(ISIh(ss)), '-',num2str(snrdB(ii)),'dB.csv'),traindata);
%     csvwrite(strcat('./test/testdata_CSI',num2str(snrdB(ii)),'dB.csv'),codedBits);
    
end
end
% fid = strcat('./test/informationBits_CSI-',num2str(numBits-1),'bits.mat');
% save(fid,'informationBits');
