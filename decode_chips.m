% Decode input 32 bit chip sequence to a 4 bit symbol
% 
% Author: Sanjib Sur
% Institution: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 01/14/2014
% 
% Comments:

function [output, error_count] = decode_chips(chips)

global DEFAULT_ERROR_THRESHOLD CHIP_DECODE_MAPPING;

% Best match when all bits are correct
best_match = hex2dec('FF');

% Ignore the first and last chip bit since it depends on the last chip
modified_chip = bitand(chips, hex2dec('7FFFFFFE'));
modified_decode_mapping = bitand(CHIP_DECODE_MAPPING, hex2dec('7FFFFFFE'));

% Replicate input chips to compare with the chip decoder vector
err_vec = de2bi(bitxor(repmat(modified_chip, 16, 1), ...
                                        modified_decode_mapping));

err_count = sum(err_vec.');

% Find the match with minimum error count
[threshold_value output] = min(err_count);

if threshold_value > DEFAULT_ERROR_THRESHOLD
    output = -1;
else
    % Convert array index into symbol value
    output = output - 1;    
end % end if threshold_value

% Symbols recovered with number of errors (if any)
error_count = threshold_value;  

end