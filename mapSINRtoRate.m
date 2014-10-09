function rt = mapSINRtoRate(sinr)

    
    if sinr > 24
        rt = 65e6;
 
    elseif sinr > 22
        rt = 58e6;
    elseif sinr > 20
        rt = 54e6;       
    elseif sinr > 17
        rt = 36e6;
        
    elseif sinr > 13 
        rt = 24e6;
    elseif sinr > 9 
        rt = 12e6;
    elseif sinr > 7
        rt = 6e6;
    else 
        %rt = 6e6; %cannot be received, but TX will try to send?
        rt = 0;
    end
