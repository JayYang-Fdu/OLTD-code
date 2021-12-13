function [encoded] = FEC_Coding(data,S)
dataT = reshape(data,[],1);
% trellis = poly2trellis(4,{'1+x+x^2+x^3','1+x+x^3'});
trellis = poly2trellis(4,[17, 13]);
codedDataT = convenc(dataT,trellis);
codedData = reshape(codedDataT,[],1);
switch S
    case 8
    for i =1:length(codedData)
        if codedData(i) == 0
            encoded((4*i-3):(4*i),1) = [0 0 1 1];
        end
        if codedData(i) == 1
            encoded((4*i-3):(4*i),1) = [1 1 0 0];
        end
    end
    
    case 2
        encoded = codedDataT;
        
    otherwise
        disp('false FEC mode');
        
end
end