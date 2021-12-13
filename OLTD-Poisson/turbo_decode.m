%% ===========================================
% Function: MAP decoding method
% Last Mofified: 2020/2/15
% ================================================= 
% this algorithm is designed for  
% x    The LLR : x = ln((P(s = +1))/(P(s = -1)))
% y    the output with soft decision
% priority_p    the priority probability of x
% dec_type    the decision type: hard or soft    
function [y2, y1]  = turbo_decode(Lebkp, akPrior)
    len = length(Lebkp);
    y2 = zeros(len/2, 1);
    y1 = zeros(len, 1);
    trellisInput = [1 0 0 0 0 0 0 1;
                    1 0 0 0 0 0 0 1;
                    0 1 0 0 0 0 1 0;     % ����ͼȫ������
                    0 1 0 0 0 0 1 0;
                    0 0 1 0 0 1 0 0;
                    0 0 1 0 0 1 0 0;
                    0 0 0 1 1 0 0 0;
                    0 0 0 1 1 0 0 0];    % describe the trellis connection: 1 connected and 0 unconnected 
    trellisOutOdd = [0 0 0 0 0 0 0 1;
                     1 0 0 0 0 0 0 0;
                     0 0 0 0 0 0 1 0;
                     0 1 0 0 0 0 0 0;
                     0 0 0 0 0 1 0 0;
                     0 0 1 0 0 0 0 0;
                     0 0 0 0 1 0 0 0;
                     0 0 0 1 0 0 0 0];    % output first code ��֧·���
    trellisOutEven = [0 0 0 0 0 0 0 1;
                      1 0 0 0 0 0 0 0;
                      0 0 0 0 0 0 1 0;
                      0 1 0 0 0 0 0 0;
                      0 0 1 0 0 0 0 0;
                      0 0 0 0 0 1 0 0;
                      0 0 0 1 0 0 0 0;
                      0 0 0 0 1 0 0 0];   % output second code
%% the connection relation between the restrict bit, the following is coded bit  
% trellisOut ���Ϊb�������Ӧ��Լ������Ӧ���b
    AbOddZeros = [1 0 0 0 0 0 0 0;
                  0 0 0 0 0 0 0 1;
                  0 1 0 0 0 0 0 0;   % trellisOutOdd�����Ϊ0��λ�ã���1����ʾfirst code���Ϊ0������
                  0 0 0 0 0 0 1 0;
                  0 0 1 0 0 0 0 0;
                  0 0 0 0 0 1 0 0;
                  0 0 0 1 0 0 0 0;
                  0 0 0 0 1 0 0 0];  
    AbOddOnes =  [0 0 0 0 0 0 0 1;
                  1 0 0 0 0 0 0 0;
                  0 0 0 0 0 0 1 0;
                  0 1 0 0 0 0 0 0;
                  0 0 0 0 0 1 0 0;
                  0 0 1 0 0 0 0 0;
                  0 0 0 0 1 0 0 0;
                  0 0 0 1 0 0 0 0];    %trellisOutOdd�����Ϊ 1 ��λ�ã���1��  ... 1 ...
    AbEvenZeros = [1 0 0 0 0 0 0 0;
                   0 0 0 0 0 0 0 1;
                   0 1 0 0 0 0 0 0;
                   0 0 0 0 0 0 1 0;
                   0 0 0 0 0 1 0 0;
                   0 0 1 0 0 0 0 0;
                   0 0 0 0 1 0 0 0;
                   0 0 0 1 0 0 0 0]; % trellisOutEven... 0 ��λ�ã���1����ʾsecond code���Ϊ 0 ������
    AbEvenOnes = [0 0 0 0 0 0 0 1;
                  1 0 0 0 0 0 0 0;
                  0 0 0 0 0 0 1 0;
                  0 1 0 0 0 0 0 0;
                  0 0 1 0 0 0 0 0;
                  0 0 0 0 0 1 0 0;
                  0 0 0 1 0 0 0 0;
                  0 0 0 0 1 0 0 0];% trellisOutEven... 1 ��λ�ã���1��  ... 1 ...
%% the connection relation between the restrict bit, the following is uncoded bit     
% trellisOut  ������a�������Ӧ��Լ������Ӧ����a
    AaOnes = [0 0 0 0 0 0 0 1;
              0 0 0 0 0 0 0 1;
              0 0 0 0 0 0 1 0;
              0 0 0 0 0 0 1 0;
              0 0 0 0 0 1 0 0;
              0 0 0 0 0 1 0 0;
              0 0 0 0 1 0 0 0;
              0 0 0 0 1 0 0 0];  % ����a=1 ��λ�ã���1
    AaZeros = [1 0 0 0 0 0 0 0;
               1 0 0 0 0 0 0 0;
               0 1 0 0 0 0 0 0;
               0 1 0 0 0 0 0 0;
               0 0 1 0 0 0 0 0;
               0 0 1 0 0 0 0 0;
               0 0 0 1 0 0 0 0;
               0 0 0 1 0 0 0 0];  % ����a=0 ��λ�ã���1
%%

    bii = 1;
    for kk = 1 : 2 : len
        if kk == 1
            fk = [1; zeros(7, 1)];
        else
            PbOddpriortmp = exp(trellisOutOdd*Lebkp(kk-2)) / (1 + exp(Lebkp(kk-2)));  % P(b2k-1 = b1,i,j | y)
            PbOddprior = PbOddpriortmp.*trellisInput;
            PbEventmp = exp(trellisOutEven*Lebkp(kk-1)) / (1 + exp(Lebkp(kk-1))); % P(b2k =b2,i,j | y)
            PbEvenprior = PbEventmp.*trellisInput;
            PkPrior = akPrior*PbOddprior.*PbEvenprior; % gamma,Pk
            fk = PkPrior'*fk;
            fk = fk / sum(fk);                
        end % if kk == 1 end
        count = len - 1;    % the reverse time slot of the trellis
        bk = [1; zeros(7, 1)];
        while count > kk
            PbOddposttmp = exp(trellisOutOdd*Lebkp(count)) / (1 + exp(Lebkp(count)));
            PbOddpost = PbOddposttmp.*trellisInput;
            PbEvenposttmp = exp(trellisOutEven*Lebkp(count+1)) / (1 + exp(Lebkp(count+1)));
            PbEvenpost = PbEvenposttmp.*trellisInput;
            PkPost = akPrior*PbOddpost.*PbEvenpost;
            bk = PkPost * bk;
            bk = bk / sum(bk);
            count = count - 2;
        end % while end
        PbOddtmp = exp(trellisOutOdd*Lebkp(kk)) / (1 + exp(Lebkp(kk)));
        PbOdd = PbOddtmp.*trellisInput;
        PbEventmp = exp(trellisOutEven*Lebkp(kk+1)) / (1 + exp(Lebkp(kk+1)));
        PbEven = PbEventmp.*trellisInput;
        Pk = akPrior*PbOdd .* PbEven;
        y2((kk+1)/2) = log(fk'*(AaOnes.*Pk)*bk/(fk'*(AaZeros.*Pk)*bk));
        y1(bii) = log(fk'*(AbOddOnes.*Pk)*bk/(fk'*(AbOddZeros.*Pk)*bk));
        y1(bii+1) = log(fk'*(AbEvenOnes.*Pk)*bk/(fk'*(AbEvenZeros.*Pk)*bk));
        bii = bii + 2;
    end % for kk end
    y2((y2)>0) = 1;
    y2((y2)<0) = 0;
    y2 = y2(1:end-3);
end % function end