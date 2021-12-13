function qx = gfsk_qt(x,T,B)
%   x = kT;
%    gfsk: BT = 0.5;
%     T = 1;
    N = 2;%gfsk
    delt = T/1e5;
    if x <2
        t = -N*T:delt:(x-1)*T;
        qx = sum( 1/(2*T)*(qfunc(2*pi*B*(t-T/2)/sqrt(log(2)))...
             -qfunc(2*pi*B*(t+T/2)/sqrt(log(2))))*delt);
    else
        qx = 0.5;
    end
end
% function qx = gfsk_qt(T,B)
% %   x = kT;
% %    gfsk: BT = 0.5;
% %     T = 1;
%     N = 2;%gfsk
%     delt = T/1e6;
%     a = [];
%     for ii = 0:0.01:3
%         if ii <2
%             t = -N*T:delt:(ii-1)*T;
%             qx = sum( 1/(2*T)*(qfunc(2*pi*B*(t-T/2)/sqrt(log(2)))...
%                  -qfunc(2*pi*B*(t+T/2)/sqrt(log(2))))*delt);
% %         else
% %             qx = 0.5;
%         end
%        a = [a qx];
%     end
% %     plot(0:0.01:3,a)
% end

