function Lcky = BCJRmodulator(y,LCk_n,yprob)
    len = length(y);
    Lcky = zeros(len, 1);
    gammakpost = zeros(4,4,len);
    M = [0 2 0 8;
         5 0 3 0;
         0 6 0 4;
         1 0 7 0];
    R = (1:1:8);
    [~,pos]=ismember(R,M);
    A = zeros(size(M));
    for mm = 1:len
        A(pos)=yprob(mm,:);
        gammakpost(:,:,mm) = A;
    end

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

    noZeroPos = find(matrices~=0);
    PkPost = zeros(4,4,len);
    bk1 =  ones(4,1,len+1);  %bk
    
    gammak = zeros(4,4,len);
    for mm =1:len
%         a = exp(Apos.*LCk(mm))/(1 + exp(LCk(mm))).*matrices;
        b = exp(Apos.*LCk_n(mm))/(1 + exp(LCk_n(mm))).*matrices;
        c = b(noZeroPos)./1;
        c = c/sum(c(1:2));
        d = gammak(:,:,mm);
        d(noZeroPos) = c;
        gammak(:,:,mm) = d;
    end
    
    
    
    for ii = len:-1:1
        PkPost(:,:,ii) = gammakpost(:,:,ii).*gammak(:,:,ii); % Pk,gamma构成的P矩阵 
        bk1(:,:,ii) = PkPost(:,:,ii)*bk1(:,:,ii+1);
        bk1(:,:,ii) = bk1(:,:,ii)/sum(bk1(:,:,ii));
       
    end
    fk = ones(4,1); %f0
    Pkinit = rand(4,4).*matrices;
    fk = Pkinit'*fk; %f1
    Lcky(1) = log(fk'*(Apos.*Pkinit)*bk1(:,:,2)/(fk'*(Aneg.*Pkinit)*bk1(:,:,2)));
    for kk = 2 : len
        
        PkPrior = gammakpost(:,:,(kk-1)).*gammak(:,:,kk-1);  % P,gamma构成的P矩阵
        fk = PkPrior'*fk;   % 前向更新 fk
        fk = fk / sum(fk);
        Pk=gammakpost(:,:,kk).*gammak(:,:,kk);
        Lcky(kk) = log(fk'*(Apos.*Pk)*bk1(:,:,kk+1)/(fk'*(Aneg.*Pk)*bk1(:,:,kk+1)));
        
    end
   
end