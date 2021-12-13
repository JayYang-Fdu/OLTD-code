function label = gen_label(ak)
% input:  ak: the input bit
%    
% output: label: the ground truth
len = length(ak);
label = zeros(len-1,1);
for ii= 1:len-1
    if ak(ii:ii+1) == [1;1]
        label(ii)=0;
    elseif ak(ii:ii+1) == [0;1]
        label(ii)=1;
    elseif ak(ii:ii+1) == [1;0]
        label(ii)=2;
    elseif ak(ii:ii+1) == [0;0]
        label(ii)=3;
    end
end

end

