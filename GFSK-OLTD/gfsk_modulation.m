function x = gfsk_modulation(bits,h,tau)
    numBits = length(bits);
%     upSampledBits = zeros(1,length(bits)*upSampRate);
    modulBits = bits*2-1;
%     gt = pulse_shape(t,B,T); %这个函数里面采样，来表示这个函数
%     H = [0.053, 0.447, 0.5*ones(1,numBits-2)]; % for upsample t=T/2,3T/2
%     H = [0, 0.25, 0.5*ones(1,numBits-2)]; % for upsample t=0,T
%     H = [0.053, 0.447, 0.5*ones(1,numBits-2)]; % for upsample t=3T/4,7T/4
%     H = [0.25, 0.5, 0.5*ones(1,numBits-2)];%for upsample t=T,2T
    T = 1;
    B = 0.5;
    %(n+tau)T,(n+1+tau)T
    sample1 = gfsk_qt(1+tau,T,B);
    sample2 = gfsk_qt(2+tau,T,B);
    pulseShaper = [sample1,sample2,0.5*ones(1,numBits-2)];
    faiT = 2*pi*h*filter(pulseShaper,1,modulBits);
    % faiT = 2*pi*h*filter(gt,1,upSampledBits);
    % theta = cumsum(faiT/upSampRate);  % 相当于把积分符号放到外面来求和了
    x = exp(1j*faiT);
end