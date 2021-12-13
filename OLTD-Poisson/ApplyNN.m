function turErrBit = ApplyNN(Net, y, sequence, iternum,turErrBit, akTest, mm, algo)


Likelihood = predict(Net, y);
LeBkP = zeros(length(y), 1);

if algo == "OLTD"
    for nn = 1 : iternum %均衡器循环
        LeCkP = LeBkP;
        LCkY = ApplyOLTD(y, LeCkP, Likelihood);
        LeCkY = LCkY - LeCkP;  % Lext(ck|y)
        LeCkY(LeCkY>12) = 12;
        LeCkY(LeCkY<-12) = -12;% the treshold have an impact on results
        LeBkY = deinterleaving(LeCkY, sequence);    % deinterleaving, Lext(bk|y),bkSoft
        [LAkP, LBkP] = turbo_decode(LeBkY, 1/2);% [ akSoft , L(bk|p) ]
        LeBkP = LBkP - LeBkY;  %Lext(bk|p)
        LeBkP = LeBkP(sequence);  %Lext(ck|p) LLR_new
        
        if nn == 1
            errbit0 = sum(LAkP~= akTest(1:end-3));
            %                 errbit0 = length(find(LAkP~= ak));
            turErrBit(mm,nn) = turErrBit(mm,nn)+errbit0;
        end
        if nn ==2
            %                 errbit1 = length(find(LAkP~= ak));
            errbit1 = sum(LAkP~= akTest(1:end-3));
            turErrBit(mm,nn) = turErrBit(mm,nn)+errbit1;
        end
        if nn ==iternum
            %                 errbit1 = length(find(LAkP~= ak));
            errbit2 = sum(LAkP~= akTest(1:end-3));
            turErrBit(mm,3) = turErrBit(mm,3)+errbit2;
        end
    end
elseif algo == "BCJRNet"
    for nn = 1 : iternum %均衡器循环
        LeCkP = LeBkP;
        LCkY = ApplyBCJRNet(y, LeCkP, Likelihood);
        LeCkY = LCkY - LeCkP;  % Lext(ck|y)
        LeCkY(LeCkY>12) = 12;
        LeCkY(LeCkY<-12) = -12;% the treshold have an impact on results
        LeBkY = deinterleaving(LeCkY, sequence);    % deinterleaving, Lext(bk|y),bkSoft
        [LAkP, LBkP] = turbo_decode(LeBkY, 1/2);% [ akSoft , L(bk|p) ]
        LeBkP = LBkP - LeBkY;  %Lext(bk|p)
        LeBkP = LeBkP(sequence);  %Lext(ck|p) LLR_new
        %     LCkY = BCJRdemodulator([1;y;1], LeCkP, [0;LeBkP;0],[[1,0,0,0];Likelihood;[1,0,0,0]],algo); %LeCkP = LLR_Old;LeBkP=LLR_New
        
        if nn == 1
            errbit0 = sum(LAkP~= akTest(1:end-3));
            %                 errbit0 = length(find(LAkP~= ak));
            turErrBit(mm,nn) = turErrBit(mm,nn)+errbit0;
        end
        if nn ==2
            %                 errbit1 = length(find(LAkP~= ak));
            errbit1 = sum(LAkP~= akTest(1:end-3));
            turErrBit(mm,nn) = turErrBit(mm,nn)+errbit1;
        end
        if nn ==iternum
            %                 errbit1 = length(find(LAkP~= ak));
            errbit2 = sum(LAkP~= akTest(1:end-3));
            turErrBit(mm,3) = turErrBit(mm,3)+errbit2;
        end
    end
end





end