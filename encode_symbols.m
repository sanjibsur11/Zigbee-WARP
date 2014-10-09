% Encode input 4 bit symbol to a 32 bit chip sequence
% 
% Author: Sanjib Sur
% Institution: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 01/14/2014
% 
% Comments: 


function [output] = encode_symbols(symbol)

global CHIP_ENCODE_MAPPING;

for cs = 1:length(symbol)
    assert(symbol(cs) < 16);
    output(cs) = CHIP_ENCODE_MAPPING(symbol(cs)+1);
end

end