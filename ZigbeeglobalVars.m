% This file contains global variables which are used in difference function
% of Zigbee PHY code
% 
% Author: Sanjib Sur
% Institution: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 01/17/2014
% 
% Comments: 


function f = ZigbeeglobalVars()

global USESIM LTF_shift LTF_t LTFscale STF_t LTF_up2 STF_up2 osamp...
    SM_CHAN repeat_preamble;


% -------- Generate 802.11 STF ---------
STF = [0 0 0 0 0 0 0 0 1+1i 0 0 0 -1+1i 0 0 0 -1-1i 0 0 0 1-1i 0 0 0 ...
    -1-1i 0 0 0 1-1i 0 0 0 0 0 0 0 1-1i 0 0 0 -1-1i 0 0 0 1-1i 0 0 0 ...
    -1-1i 0 0 0 -1+1i 0 0 0 1+1i 0 0 0 0 0 0 0].';
% Before ifft, ensure the zero-frequency on index 1; minus-frequency begins
% from half (so the samples near the half-freq should be a sequence of 0s)
STF_t = ifft(fftshift(STF)).';
STF_t_short = STF_t(1:16);
STF_10 = repmat(STF_t_short,1,10);
STF_I = real(STF_10);
STF_Q = imag(STF_10);
% Upsample by osamp so the standard preamble occupies a bandwidth of 
% +-40/osamp MHz (computed for a sampling frequency of 40 MHz)
[STF_I_up2] = interp(STF_I, osamp);
[STF_Q_up2] = interp(STF_Q, osamp);
% Scale to span -1,1 range of DAC
scale_ShortSyms = max([ max(abs(STF_I_up2)), max(abs(STF_Q_up2)) ]);
[STF_I_up2] = (1/scale_ShortSyms)*STF_I_up2;
[STF_Q_up2] = (1/scale_ShortSyms)*STF_Q_up2;
STF_up2 = (STF_I_up2 + sqrt(-1)*STF_Q_up2);


LTF_freq_bot = [0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 ...
    1 -1 1 -1 1 1 1 1]';
LTF_freq_top = [1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 ...
    -1 1 1 1 1 0 0 0 0 0]';
LTF_freq = [LTF_freq_bot ; 0 ; LTF_freq_top];


LTF_shift = fftshift(LTF_freq);
LTF_t = ifft(LTF_shift).';

LTF_time_up2 = interp(LTF_t, osamp); % Upsample 
LTFscale = 1/max([max(abs(real(LTF_time_up2))), max(abs(imag(LTF_time_up2))) ]);
LTF_time_up2 = LTFscale * LTF_time_up2; % Scale to span -1,1 range of DAC 

% Concatenate two long training symbols and add cyclic prefix
LTF_up2 = [LTF_time_up2(32*osamp+1:64*osamp) repmat(LTF_time_up2,1,2)];


% Adding channel estimation preamble
repeat_preamble = 4;
SM_CHAN = [];
preamble_symbols = [];

for ps = 1:repeat_preamble
    for k = 1:4
        preamble_symbols = [preamble_symbols qpsk_half_sine_mod(k-1)];
    end
end

SM_CHAN = preamble_symbols;

end