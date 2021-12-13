function DecodeBit = Viterbi(yk, noiseVar, ISIh, M)
pd = makedist("stable","alpha",0.5,"beta",0.75);    
sqrtnoiseVar = sqrt(noiseVar);
LenISIh = length(ISIh);
if M==4
    len = length(yk);
    Path = zeros(M^(LenISIh-1),1);
    tmp = [1-1i, 1+1i, -1+1i, -1-1i]/sqrt(log2(M));
    trellisOutput =zeros(M^(LenISIh-1),M^(LenISIh-1)); % probability transfer matrix
    trellisOutput(:,1) = ISIh(1)*tmp(1) + ISIh(2)*tmp; 
    trellisOutput(:,2) = ISIh(1)*tmp(2) + ISIh(2)*tmp;
    trellisOutput(:,3) = ISIh(1)*tmp(3) + ISIh(2)*tmp;
    trellisOutput(:,4) = ISIh(1)*tmp(4) + ISIh(2)*tmp;
    
    PathMetric = [0;0;0;0];
    NewPathMetric = zeros(4,1);
    for ii = 1:len
    %     yprob = (exp(-(abs(yk(ii)-trellisOutput)).^2/noiseVar))/(pi*noiseVar);%Guassion
        yprob = (exp(-(abs(yk(ii)-trellisOutput))/(sqrtnoiseVar/2)))/noiseVar; % laplace
        Branchmetric = reshape(-log(yprob),4,4);
        [Onemin,Oneclo] = min(Branchmetric(:,1)+PathMetric);
        NewPathMetric(1) = Onemin;
        NewP1 = [Path(Oneclo,:) ,  Oneclo];

        [Twomin,Twoclo] = min(Branchmetric(:,2)+PathMetric);
        NewPathMetric(2) = Twomin;
        NewP2 = [Path(Twoclo,:) ,  Twoclo];

        [Threemin,Threeclo] = min(Branchmetric(:,3)+PathMetric);
        NewPathMetric(3) = Threemin;
        NewP3 = [Path(Threeclo,:) ,  Threeclo];

        [Fourmin,Fourclo] = min(Branchmetric(:,4)+PathMetric);
        NewPathMetric(4) = Fourmin;
        NewP4 = [Path(Fourclo,:) ,  Fourclo];

        Path = [Path, zeros(4,1)];
        Path(1,:) = NewP1;
        Path(2,:) = NewP2;
        Path(3,:) = NewP3;
        Path(4,:) = NewP4;
        PathMetric = NewPathMetric;
    end
        [~,index] = min(PathMetric);
        FinalPath = [Path(index,3:end), index];
        DecodeBit = zeros(len*2,1);
        for jj =  1: len
            if FinalPath(jj) == 1
                DecodeBit(2*jj-1:2*jj) = [1;1];
            elseif FinalPath(jj) == 2
                DecodeBit(2*jj-1:2*jj) = [1;0];
            elseif FinalPath(jj) == 3
                DecodeBit(2*jj-1:2*jj) = [0;0];
            else
                DecodeBit(2*jj-1:2*jj) = [0;1];
            end
        end
elseif M==2
        len = length(yk);
        Path = zeros(M^(LenISIh-1),1);
        tmp = [1 -1]/sqrt(log2(M));
        trellisOutput =zeros(M^(LenISIh-1),M^(LenISIh-1)); % probability transfer matrix
        trellisOutput(:,1) = ISIh(1)*tmp(1) + ISIh(2)*tmp; 
        trellisOutput(:,2) = ISIh(1)*tmp(2) + ISIh(2)*tmp;

        % noiseVar = noiseVar;
        PathMetric = [0;0];
        NewPathMetric = zeros(2,1);
        for ii = 1:len
%             yprob = (exp(-(abs(yk(ii)-trellisOutput)).^2/(2*noiseVar)))/(sqrt(2*pi*noiseVar));%Gaussian
            yprob = (exp(-(abs(yk(ii)-trellisOutput))/(sqrtnoiseVar/2)))/(noiseVar);%Laplace
%             yprob = pdf(pd, yk(ii)-trellisOutput*sqrt(noiseVar));
%             yprob = pdf('stable', yk(ii)-trellisOutput*sqrt(noiseVar),0.5, 0, 1, 0);
            Branchmetric = reshape(-log(yprob),2,2);
            [Onemin,Oneclo] = min(Branchmetric(:,1)+PathMetric);
            NewPathMetric(1) = Onemin;
            NewP1 = [Path(Oneclo,:) ,  Oneclo];

            [Twomin,Twoclo] = min(Branchmetric(:,2)+PathMetric);
            NewPathMetric(2) = Twomin;
            NewP2 = [Path(Twoclo,:) ,  Twoclo];



            Path = [Path, zeros(2,1)];
            Path(1,:) = NewP1;
            Path(2,:) = NewP2;

            PathMetric = NewPathMetric;
        end
        [~,index] = min(PathMetric);
        FinalPath = [Path(index,3:end), index];
        DecodeBit = zeros(len,1);
        for jj =  1: len
            if FinalPath(jj) == 1
                DecodeBit(jj) = 1;

            else
                DecodeBit(jj) = 0;
            end
        end
end

end