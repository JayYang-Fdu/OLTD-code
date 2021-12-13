function y = phaseCompute(x)
    len = length(x);
    y = zeros(len,1);
    summ = 0;
    for ii = 1:len
        
        if summ==0 && x(ii)==1
            y(ii) = 1;
        elseif summ==0 && x(ii)==-1
            y(ii) = 7;
           
        elseif summ==1 && x(ii)==1
            y(ii) = 2;    
        elseif summ==1 && x(ii)==-1
            y(ii) = 4;
        elseif summ==2 && x(ii)==1
            y(ii) = 3; 
        elseif summ==2 && x(ii)==-1
            y(ii) = 5;
        elseif summ==3 && x(ii)==1
            y(ii) = 0; 
        elseif summ==3 && x(ii)==-1
            y(ii) = 6; 
        end
        summ = mod(summ + x(ii),4);
    end
    
    
    
    
end