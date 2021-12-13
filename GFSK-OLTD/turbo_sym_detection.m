%% ===========================================
% Function: Equalization forward/backward algorithm 
% ================================================= 
% this algorithm is designed for sequence detection 
%  Input : x (length 2^N),lecky<=>prior probability
% Output : result（soft decision）
% priority_p    the prior LLR of x
% noiseVar    the noise variance 噪声方差
%%
function Lcky = turbo_sym_detection(x, Leckp, noiseVar, channelH)
%% upsample = 1;
    len = length(x);
    Lcky = zeros(len, 1);
    trellisOutput = channelH.*[0 0.7071+0.7071i 0 0.7071-0.7071i;
                     0.7071+0.7071i 0 -0.7071+0.7071i 0;
                     0 -0.7071+0.7071i 0 -0.7071-0.7071i;
                     0.7071-0.7071i 0  -0.7071-0.7071i 0]; %t=T、2T处采样

    Apos = [0 1 0 0;
            0 0 1 0;
            0 0 0 1;
            1 0 0 0]; %  A(+1)


    Aneg = [0 0 0 1;
            1 0 0 0;
            0 1 0 0;
            0 0 1 0]; %  A(-1)

   matrices = [0 1 0 1;
               1 0 1 0;
               0 1 0 1;
               1 0 1 0];    % for the Pk

    noiseVar2 = noiseVar+db2pow(-1); %干扰功率加噪声功率
%     noiseVar2 = noiseVar;
%     sqrtnoiseVar = pi*noiseVar2;
    gammakpost = zeros(4,4,len);
    PkPost = zeros(4,4,len);
    bk1 =  ones(4,1,len+1);  %bk


    for ii = len:-1:1
        gammakpost(:,:,ii) = exp(Apos.*Leckp(ii))/(1 + exp(Leckp(ii))).*matrices; % P(xk=xi,j)
        PkPost(:,:,ii) = gammakpost(:,:,ii).* (exp(-(abs(x(ii)-trellisOutput)).^2/noiseVar2)); % Pk,gamma构成的P矩阵 
        
        bk1(:,:,ii) = PkPost(:,:,ii)*bk1(:,:,ii+1);
        bk1(:,:,ii) = bk1(:,:,ii)/sum(bk1(:,:,ii));
       
    end
    fk = ones(4,1); %f0
    Pkinit = rand(4,4).*matrices;
    fk = Pkinit'*fk; %f1
    Lcky(1) = log(fk'*(Apos.*Pkinit)*bk1(:,:,2)/(fk'*(Aneg.*Pkinit)*bk1(:,:,2)));
    for kk = 2 : len
        
        gammakprior = exp(Apos.*Leckp(kk-1))/(1 + exp(Leckp(kk-1))).*matrices; % P(xk=xi,j)
        PkPrior = gammakprior.*(exp(-(abs(x(kk-1)-trellisOutput)).^2/noiseVar2));  % P,gamma构成的P矩阵
        fk = PkPrior'*fk;   % 前向更新 fk
        fk = fk / sum(fk);
        gammak = exp(Apos.*Leckp(kk))/(1 + exp(Leckp(kk))).*matrices;
        Pk=gammak.*exp(-(abs(x(kk)-trellisOutput)).^2/noiseVar2);
        Lcky(kk) = log(fk'*(Apos.*Pk)*bk1(:,:,kk+1)/(fk'*(Aneg.*Pk)*bk1(:,:,kk+1)));
        
    end
%     Lcky = Lcky(2:end-1);
   
end

