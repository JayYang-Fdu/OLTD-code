function [y] = func(x)
%UNTITLED4 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
len = length(x);
y = zeros(4*len, 1);
for mm = 1:len
    y(4*(mm-1)+1:4*(mm-1)+2) = x(mm);
    y(4*(mm-1)+3:4*(mm-1)+4) = -x(mm);
end
end

