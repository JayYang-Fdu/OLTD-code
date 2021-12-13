function DecodeBit = Apply_LANNViterbi(yk,yprob, M)

if M == 4
    len = length(yk);
    Path = zeros(4,1);


    PathMetric = [0;0;0;0];
    NewPathMetric = zeros(4,1);
    for ii = 1:len
        Branchmetric = reshape(-log(yprob(ii,:)),4,4);
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
elseif M ==2
    len = length(yk);
    Path = zeros(2,1);
    PathMetric = [0;0];
    NewPathMetric = zeros(2,1);
    for ii = 1:len
        Branchmetric = reshape(-log(yprob(ii,:)),2,2);
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