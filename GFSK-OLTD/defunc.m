function [y] = defunc(x)
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    len = length(x);
    y = zeros(len/4,1);
    
    for mm = 1:len/4
        y(mm) = (x(4*(mm-1)+1)+x(4*(mm-1)+2)-x(4*(mm-1)+3)-x(4*(mm-1)+4));
    end
end

