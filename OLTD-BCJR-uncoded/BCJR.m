function LLR = BCJR(yk,  noiseVar, ISIh)
%======INPUT==================
%yk:the channel output
%Lext: initial LLR 0
%noiseVar: the variance of noise
%ISIh: channel gain
%======OUTPUT==================
%LLR:the LLR of every bit
% m = 2; % for qpsk
%     gamma = mm;
    len = length(yk);
    Lext = zeros(len*2,1);
    Lcky = zeros(len*2, 1);
    tmp = [1-1i, 1+1i, -1+1i, -1-1i]/sqrt(2);
    trellisOutput =zeros(4,4); % probability transfer matrix
    trellisOutput(:,1) = ISIh(1)*tmp(1) + ISIh(2)*tmp; 
    trellisOutput(:,2) = ISIh(1)*tmp(2) + ISIh(2)*tmp;
    trellisOutput(:,3) = ISIh(1)*tmp(3) + ISIh(2)*tmp;
    trellisOutput(:,4) = ISIh(1)*tmp(4) + ISIh(2)*tmp;

    Afirst1 = [1 1 0 0;
               1 1 0 0;
               1 1 0 0;
               1 1 0 0]; %  A(+1)


    Afirst0 = [0 0 1 1;
               0 0 1 1;
               0 0 1 1;
               0 0 1 1]; %  A(-1)

   Asecond1 = [1 0 0 1;
               1 0 0 1;
               1 0 0 1;
               1 0 0 1];    % for the Pk
   Asecond0 = [0 1 1 0;
               0 1 1 0;
               0 1 1 0;
               0 1 1 0];

%     noiseVar2 = noiseVar+1; %干扰功率加噪声功率
    noiseVar2 = noiseVar+5;
%     sqrtnoiseVar = sqrt(noiseVar);
    gammakpost = ones(4,4,len);
    PkPost = zeros(4,4,len);
    bk1 =  ones(4,1,(len+1));  %bk


    for ii = len:-1:1
        
        gammakpost(:,:,ii) = exp(Afirst1.*Lext(2*ii-1))/(1 + exp(Lext(2*ii-1))).*...
        (exp(Asecond1.*Lext(2*ii))/(1 + exp(Lext(2*ii))));
%         tmp = (yk(ii)-trellisOutput);
        
        PkPost(:,:,ii) = gammakpost(:,:,ii).* (exp(-(abs(yk(ii)-trellisOutput)).^2/noiseVar2)); % Pk,gamma构成的P矩阵
%         PkPost(:,:,ii) = gammakpost(:,:,ii).*(exp(-(abs(yk(ii)-trellisOutput))/(sqrtnoiseVar/2))); %Laplace\
%         PkPost(:,:,ii) = gammakpost(:,:,ii).*(gamma^2/(pi^2)./((((real(tmp))^2+gamma^2)).*(((imag(tmp))^2+gamma^2)))); %cauthy
        bk1(:,:,ii) = PkPost(:,:,ii)*bk1(:,:,ii+1);
        bk1(:,:,ii) = bk1(:,:,ii)/sum(bk1(:,:,ii));
       
    end
    fk = ones(4,1); %f0
    Pkinit = rand(4,4);
    fk = Pkinit'*fk; %f1
    
    Lcky(1) = log((fk'*(Afirst1.*Pkinit)*bk1(:,:,2))/(fk'*(Afirst0.*Pkinit)*bk1(:,:,2)));
    Lcky(2) = log((fk'*(Asecond1.*Pkinit)*bk1(:,:,2))/(fk'*(Asecond0.*Pkinit)*bk1(:,:,2)));
    for kk = 2 : len
        
        gammakprior = exp(Afirst1.*Lext(2*(kk-1)-1))/(1 + exp(Lext(2*(kk-1)-1))).*... 
            exp(Asecond1.*Lext(2*(kk-1)))/(1 + exp(Lext(2*(kk-1)))); % P(xk=xi,j)
        PkPrior = gammakprior.*(exp(-(abs(yk(kk-1)-trellisOutput)).^2/noiseVar2));  % P,gamma构成的P矩阵
%         tmp1 = (yk(kk-1)-trellisOutput);
%         PkPrior = gammakprior.*(gamma^2/(pi^2)./((((real(tmp1))^2+gamma^2)).*(((imag(tmp1))^2+gamma^2)))); %cauthy
%         PkPrior = gammakprior.*(exp(-(abs(yk(kk-1)-trellisOutput))/(sqrtnoiseVar/2)));
        fk = PkPrior'*fk;   % 前向更新 fk
        fk = fk / sum(fk);
        gammak = exp(Afirst1.*Lext(2*kk-1))/(1 + exp(Lext(2*kk-1))).*... 
            exp(Asecond1.*Lext(2*kk))/(1 + exp(Lext(2*kk)));
        Pk=gammak.*exp(-(abs(yk(kk)-trellisOutput)).^2/noiseVar2);
%         Pk=gammak.*(exp(-(abs(yk(kk)-trellisOutput))/(sqrtnoiseVar/2)));
%         tmp2 = (yk(kk)-trellisOutput);
%         Pk=gammak.*(gamma^2/(pi^2)./((((real(tmp2))^2+gamma^2)).*(((imag(tmp2))^2+gamma^2)))); %cauthy
        Lcky(2*kk-1) = log(fk'*(Afirst1.*Pk)*bk1(:,:,kk+1)/(fk'*(Afirst0.*Pk)*bk1(:,:,kk+1)));
        Lcky(2*kk) = log(fk'*(Asecond1.*Pk)*bk1(:,:,kk+1)/(fk'*(Asecond0.*Pk)*bk1(:,:,kk+1)));
    end
%     Lcky = Lcky(2:end-1);
   LLR = Lcky;
   LLR(LLR>0)=1;
   LLR(LLR<0)=0;
end
