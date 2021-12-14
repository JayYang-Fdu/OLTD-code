function  LLR = ApplyBCJRNet(yk, LeBkP,yprob)
%  %LeCkP = LLR_Old;LeBkP=LLR_New
%==================================
len = length(yk);
Lcky = zeros(len, 1);
gamma = zeros(2,2,len);


Psk = 1/2*ones(2,1,len);
for ss = 2:len
    Psk(1,1,ss) = exp(LeBkP(ss-1))/(1 + exp(LeBkP(ss-1)));
    Psk(2,1,ss) = 1/(1 + exp(LeBkP(ss-1)));
end
for mm = 1:len
    if mm == 1
        gamma(:,:,mm) =  reshape(yprob(mm,:),2,2);
    else
        gammatmp=  reshape(yprob(mm,:),2,2);
        gamma(1,:,mm) = gammatmp(1,:)/Psk(1,1,mm-1);
        gamma(2,:,mm) = gammatmp(2,:)/Psk(2,1,mm-1);
    end
    
end



Apos = [1 0 ;
    1 0 ]; %  A(+1)
Aneg = [0 1 ;
    0 1 ]; %  A(-1)

% gammakpost = ones(2,2,len);
PkPost = zeros(2,2,len);
bk1 =  ones(2,1,(len+1));  %bk

% for mm = 1:len
%     gammakpost(:,:,mm) = exp(Apos.*LeBkP(mm))/(1 + exp(LeBkP(mm)));
% end
for ii = len:-1:1
    
    % PkPost(:,:,ii) = exp(Apos.*LeBkP(ii))/(1 + exp(LeBkP(ii))).* gamma(:,:,ii); % Pk,gamma构成的P矩阵
    PkPost(:,:,ii) = gamma(:,:,ii); % Pk,gamma构成的P矩阵
    bk1(:,:,ii) = PkPost(:,:,ii)*bk1(:,:,ii+1);
    bk1(:,:,ii) = bk1(:,:,ii)/sum(bk1(:,:,ii));
    
end
fk = ones(2,1); %f0
Pkinit = rand(2,2);
fk = Pkinit'*fk; %f1

Lcky(1) = log((fk'*(Apos.*Pkinit)*bk1(:,:,2))/(fk'*(Aneg.*Pkinit)*bk1(:,:,2)));

for kk = 2 : len
    
    %gammakprior = exp(Apos.*LeBkP(kk-1))/(1 + exp(LeBkP(kk-1))); % P(xk=xi,j)
    PkPrior = gamma(:,:,kk-1);  % P,gamma构成的P矩阵
    fk = PkPrior'*fk;   % 前向更新 fk
    fk = fk / sum(fk);
    %gammak = exp(Apos.*LeBkP(kk))/(1 + exp(LeBkP(kk)));
    Pk=gamma(:,:,kk);
    
    Lcky(kk) = log(fk'*(Apos.*Pk)*bk1(:,:,kk+1)/(fk'*(Aneg.*Pk)*bk1(:,:,kk+1)));
end
% Lcky = Lcky(2:end-1);
LLR = Lcky;
% LLR(LLR>0)=1;
% LLR(LLR<0)=0;
end

