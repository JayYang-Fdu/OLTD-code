function LLR = ApplyOLTD(yk, Lext, yprob)
len = length(yk);
Lcky = zeros(len, 1);
gamma = zeros(2,2,len);
for mm = 1:len
    gamma(:,:,mm) =  reshape(yprob(mm,:),2,2);
end
Apos = [1 0 ;
    1 0 ]; %  A(+1)
Aneg = [0 1 ;
    0 1 ]; %  A(-1)

PkPost = zeros(2,2,len);
bk1 =  ones(2,1,(len+1));  %bk


for ii = len:-1:1
    
    PkPost(:,:,ii) = exp(Apos.*Lext(ii))/(1 + exp(Lext(ii))).* gamma(:,:,ii); % Pk,gamma构成的P矩阵
    
    bk1(:,:,ii) = PkPost(:,:,ii)*bk1(:,:,ii+1);
    bk1(:,:,ii) = bk1(:,:,ii)/sum(bk1(:,:,ii));
    
end
fk = ones(2,1); %f0
Pkinit = rand(2,2);
fk = Pkinit'*fk; %f1

Lcky(1) = log((fk'*(Apos.*Pkinit)*bk1(:,:,2))/(fk'*(Aneg.*Pkinit)*bk1(:,:,2)));

for kk = 2 : len
    
    gammakprior = exp(Apos.*Lext(kk-1))/(1 + exp(Lext(kk-1))); % P(xk=xi,j)
    PkPrior = gammakprior.*gamma(:,:,kk-1);  % P,gamma构成的P矩阵
    fk = PkPrior'*fk;   % 前向更新 fk
    fk = fk / sum(fk);
    gammak = exp(Apos.*Lext(kk))/(1 + exp(Lext(kk)));
    Pk=gammak.*gamma(:,:,kk);
    
    Lcky(kk) = log(fk'*(Apos.*Pk)*bk1(:,:,kk+1)/(fk'*(Aneg.*Pk)*bk1(:,:,kk+1)));
end
%    Lcky = Lcky(2:end-1);
LLR = Lcky;
%    LLR(LLR>0)=1;
%    LLR(LLR<0)=0;
end