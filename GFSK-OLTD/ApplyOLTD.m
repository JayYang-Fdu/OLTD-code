function turErrBit = ApplyOLTD(Net, y, algotype, sequence, iternum,turErrBit, akTest, mm)

yTest = [real(y), imag(y)];
Likelihood = predict(Net, yTest); 


if algotype == "Turbo"
   
    LeBkP = zeros(length(y), 1);
    for nn = 1 : iternum %¾ùºâÆ÷Ñ­»·

        LeCkP = LeBkP;  % Lext(ck|p)
        LCkY = BCJRmodulator(y, LeCkP, Likelihood); 
        LeCkY = LCkY - LeCkP;  % Lext(ck|y)
        LeCkY(LeCkY>50) = 50;
        LeCkY(LeCkY<-50) = -50;
        LeBkY = deinterleaving(LeCkY, sequence);    % deinterleaving, Lext(bk|y),bkSoft
        LeBY = defunc(LeBkY);
        [LAkP, LBkP] = turbo_decode(LeBY, 1/2);% [ akSoft , L(bk|p) ]
        LBkP = func(LBkP);
        LeBkP = LBkP - LeBkY;  %Lext(bk|p)
        LeBkP = LeBkP(sequence);  %Lext(ck|p)
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
        if nn == 3
    %                 errbit2 = length(find(LAkP~= ak));
            errbit2 = sum(LAkP~= akTest(1:end-3));
            turErrBit(mm,nn) = turErrBit(mm,nn)+errbit2; 
        end

    end

end




end