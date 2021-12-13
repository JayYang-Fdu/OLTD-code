function LLR = Apply_LANNBCJR(yk, yprob)
    len = length(yk);
    Lcky = zeros(len, 1);
    Lext = zeros(len*2,1);
%     tmp = [1-1i, 1+1i, -1+1i, -1-1i]/sqrt(2);
%     trellisOutput =zeros(4,4); % probability transfer matrix
%     trellisOutput(:,1) = ISIh(1)*tmp(1) + ISIh(2)*tmp; 
%     trellisOutput(:,2) = ISIh(1)*tmp(2) + ISIh(2)*tmp;
%     trellisOutput(:,3) = ISIh(1)*tmp(3) + ISIh(2)*tmp;
%     trellisOutput(:,4) = ISIh(1)*tmp(4) + ISIh(2)*tmp;
    gamma = zeros(4,4,len);
    for mm = 1:len
        gamma(:,:,mm) =  reshape(yprob(mm,:),4,4);
    end
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
%     noiseVar2 = noiseVar;
%     sqrtnoiseVar = pi*noiseVar;

    f0 = ones(4,1); %f0
    Pkinit = rand(4,4);
    fk = Pkinit'*f0; %f1
    

%% 
    gammakpost = ones(4,4,len);
    PkPost = zeros(4,4,len);
    bk1 =  ones(4,1,(len+1));  %bk
    for ii = len:-1:1
        
        gammakpost(:,:,ii) = exp(Afirst1.*Lext(2*ii-1))/(1 + exp(Lext(2*ii-1))).*...
        (exp(Asecond1.*Lext(2*ii))/(1 + exp(Lext(2*ii))));
        
        
        PkPost(:,:,ii) = gammakpost(:,:,ii).* gamma(:,:,ii); % Pk,gamma构成的P矩阵
%         PkPost(:,:,ii) =  gamma(:,:,ii);
        
        bk1(:,:,ii) = PkPost(:,:,ii)*bk1(:,:,ii+1);
        bk1(:,:,ii) = bk1(:,:,ii)/sum(bk1(:,:,ii));
       
    end
    Lcky(1) = log((fk'*(Afirst1.*Pkinit)*bk1(:,:,2))/(fk'*(Afirst0.*Pkinit)*bk1(:,:,2)));
    Lcky(2) = log((fk'*(Asecond1.*Pkinit)*bk1(:,:,2))/(fk'*(Asecond0.*Pkinit)*bk1(:,:,2)));
    for kk = 2 : len
        
        gammakprior = exp(Afirst1.*Lext(2*(kk-1)-1))/(1 + exp(Lext(2*(kk-1)-1))).*... 
            exp(Asecond1.*Lext(2*(kk-1)))/(1 + exp(Lext(2*(kk-1)))); % P(xk=xi,j)
        PkPrior = gammakprior.*gamma(:,:,kk-1);  % P,gamma构成的P矩阵
%         PkPrior = gamma(:,:,kk-1); 
        fk = PkPrior'*fk;   % 前向更新 fk
        fk = fk / sum(fk);
        gammak = exp(Afirst1.*Lext(2*kk-1))/(1 + exp(Lext(2*kk-1))).*... 
            exp(Asecond1.*Lext(2*kk))/(1 + exp(Lext(2*kk)));
        Pk=gammak.*gamma(:,:,kk);
        Lcky(2*kk-1) = log(fk'*(Afirst1.*Pk)*bk1(:,:,kk+1)/(fk'*(Afirst0.*Pk)*bk1(:,:,kk+1)));
        Lcky(2*kk) = log(fk'*(Asecond1.*Pk)*bk1(:,:,kk+1)/(fk'*(Asecond0.*Pk)*bk1(:,:,kk+1)));
    end
%     Lcky = Lcky(2:end-1);
   LLR = Lcky;
   LLR(LLR>0)=1;
   LLR(LLR<0)=0;

%%

% bn = ones(4,1);
% pf = repmat(eye(4),[1,1,len]);
% pb = repmat(eye(4),[1,1,len+1]);
% for ii=1:len
% 
%     if ii == 1
%         pf(:,:,ii)= Pkinit*gamma(:,:,ii);
%         pb(:,:,len) = gamma(:,:,len);
%     else
%         pf(:,:,ii) = pf(:,:,ii-1)*gamma(:,:,ii);
%         pb(:,:,len+1-ii) = gamma(:,:,len+1-ii)*pb(:,:,len+2-ii);
%     end
% end
% Lcky(1) = log((fk'*(Afirst1.*Pkinit)*pb(:,:,2)*bn)/(fk'*(Afirst0.*Pkinit)*pb(:,:,2)*bn));
% Lcky(2) = log((fk'*(Asecond1.*Pkinit)*pb(:,:,2)*bn)/(fk'*(Asecond0.*Pkinit)*pb(:,:,2)*bn));
% for ii = 2:len
%     Lcky(2*ii-1) = log((fk'*(pf(:,:,ii).*Afirst1)*pb(:,:,ii+1)*bn)/(fk'*(pf(:,:,ii).*Afirst0)*pb(:,:,ii+1)*bn));
%     Lcky(2*ii) = log((fk'*(pf(:,:,ii).*Asecond1)*pb(:,:,ii+1)*bn)/(fk'*(pf(:,:,ii).*Asecond0)*pb(:,:,ii+1)*bn));
% end
% 
%    LLR = Lcky;
%    LLR(LLR>0)=1;
%    LLR(LLR<0)=0;



end