function y = GenLabel(ak, yk, M)
len = length(yk);
y = zeros(len, 1);
if M==4
    for ii = 1:len
        if ak(2*ii-1:2*ii+2)==[1; 1; 1; 1]
            y(ii) = 0;
        elseif ak(2*ii-1:2*ii+2)==[1; 0; 1; 1]
             y(ii) = 1;   
        elseif ak(2*ii-1:2*ii+2)==[0; 0; 1; 1]
             y(ii) = 2; 
        elseif ak(2*ii-1:2*ii+2)==[0; 1; 1; 1]
             y(ii) = 3; 


        elseif ak(2*ii-1:2*ii+2)==[1; 1; 1; 0]
             y(ii) = 4;   
        elseif ak(2*ii-1:2*ii+2)==[1; 0; 1; 0]
             y(ii) = 5; 
        elseif ak(2*ii-1:2*ii+2)==[0; 0; 1; 0]
             y(ii) = 6;
        elseif ak(2*ii-1:2*ii+2)==[0; 1; 1; 0]
             y(ii) = 7;   


        elseif ak(2*ii-1:2*ii+2)==[1; 1; 0; 0]
             y(ii) = 8; 
        elseif ak(2*ii-1:2*ii+2)==[1; 0; 0; 0]
             y(ii) = 9; 
        elseif ak(2*ii-1:2*ii+2)==[0; 0; 0; 0]
             y(ii) = 10;   
        elseif ak(2*ii-1:2*ii+2)==[0; 1; 0; 0]
             y(ii) = 11; 


        elseif ak(2*ii-1:2*ii+2)==[1; 1; 0; 1]
             y(ii) = 12;
        elseif ak(2*ii-1:2*ii+2)==[1; 0; 0; 1]
             y(ii) = 13;   
        elseif ak(2*ii-1:2*ii+2)==[0; 0; 0; 1]
             y(ii) = 14; 
        elseif ak(2*ii-1:2*ii+2) ==[0; 1; 0; 1]
             y(ii) = 15;
        end
    end
elseif M==2
    for ii = 1:len
        if ak(ii:ii+1) == [1;1]
            y(ii) = 0;
        elseif ak(ii:ii+1) == [0;1]
             y(ii) = 1;   
        elseif ak(ii:ii+1) == [1;0]
             y(ii) = 2; 
        elseif ak(ii:ii+1) == [0;0]
             y(ii) = 3; 
        end
    end
    
end

end