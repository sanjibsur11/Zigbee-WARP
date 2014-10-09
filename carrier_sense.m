% Carrier sensing function on received samples. This is required to
% remove extraneous samples at the beginning which might cause trouble for
% the FM demodulator, DC offset corrector and/or M&M circuit
% 
% Author: Sanjib Sur
% Institution: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 01/14/2014
% 
% Comments: 

function [idx] = carrier_sense(rxData)

global padding_size CARRIER_SENSE_WIN_SIZE RSSthresh;

% Find the current window energy starting from the beginning
prev_win_energy = 0;
ns = 1;

% Start scanning through all the samples for energy going beyond the
% previous energy
while ns < length(rxData) - CARRIER_SENSE_WIN_SIZE
    curr_win_energy = sum(abs(rxData(ns:ns+CARRIER_SENSE_WIN_SIZE)).^2);
    
    % Check if the current windows energy is significantly higher than the 
    % previous windows energy level. 
    if curr_win_energy - prev_win_energy > RSSthresh
        % We found a match now
        break;
    else
        ns = ns + CARRIER_SENSE_WIN_SIZE;
        prev_win_energy = curr_win_energy;
    end
end

% We did not found any signal
if ns >= length(rxData)
    idx = 0;
else
    % We found some significant signal
    if ns > CARRIER_SENSE_WIN_SIZE
        idx = ns - CARRIER_SENSE_WIN_SIZE;
    else
        idx = ns;
    end
end

end