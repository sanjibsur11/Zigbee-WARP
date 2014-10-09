function f = MUMIMOglobalVars()

global VLTF_shift LTF_shift LTF_t LTFscale STF_t numTxAntenna LTF_up2 STF_up2 VLTF_t_up VLTF_freq USESIM osamp numOFDMsymb otherTo48 numRxNode;


osamp = 2; % Oversampling rate or Number of samples per symbol, 40Mhz/osamp wide band



numOFDMsymb = 10;

if USESIM 
    osamp = 1;
end

% -------- Generate 802.11 STF ---------
STF = [0 0 0 0 0 0 0 0 1+1i 0 0 0 -1+1i 0 0 0 -1-1i 0 0 0 1-1i 0 0 0 -1-1i 0 0 0 1-1i 0 0 0 0 0 0 0 1-1i 0 0 0 -1-1i 0 0 0 1-1i 0 0 0 -1-1i 0 0 0 -1+1i 0 0 0 1+1i 0 0 0 0 0 0 0].';
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


LTF_freq_bot = [0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1]';
LTF_freq_top = [1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0]';
LTF_freq = [LTF_freq_bot ; 0 ; LTF_freq_top];



LTF_shift = fftshift(LTF_freq);
LTF_t = ifft(LTF_shift).';

LTF_time_up2 = interp(LTF_t, osamp); % Upsample 
LTFscale = 1/max([max(abs(real(LTF_time_up2))), max(abs(imag(LTF_time_up2))) ]);
LTF_time_up2 = LTFscale * LTF_time_up2; % Scale to span -1,1 range of DAC 

% Concatenate two long training symbols and add cyclic prefix
LTF_up2 = [LTF_time_up2(32*osamp+1:64*osamp) repmat(LTF_time_up2,1,2)];

% ------ VHT-LTF -------
VLTF_freq_bot(1,:) = [0 0 0 0 0 0 1 -1 -1 1 1 1 -1 -1 -1 1 -1 1 1 -1 1 -1 1 1 1 -1 1 -1 -1 1 1 1]';
VLTF_freq_top(1,:) = [-1 -1 -1 1 -1 1 1 -1 1 -1 -1 1 -1 -1 1 1 -1 -1 1 -1 -1 -1 1 1 -1 1 0 0 0 0 0]';

VLTF_freq_bot(2,:) = [0 0 0 0 0 0 -1 -1 1 -1 1 1 -1 1 -1 -1 -1 1 -1 1 1 -1 1 1 1 -1 1 -1 -1 -1 1 1]';
VLTF_freq_top(2,:) = [1 -1 1 -1 -1 -1 1 1 -1 1 -1 -1 1 -1 -1 1 1 -1 1 1 -1 1 1 1 -1 -1 0 0 0 0 0]';

% --- for rx 3 ---
VLTF_freq_bot(3,:) = [0 0 0 0 0 0 1 -1 -1 -1 -1 1 -1 1 1 -1 1 1 1 -1 -1 1 -1 -1 1 -1 -1 1 1 1 -1 -1]';
VLTF_freq_top(3,:) = [-1 1 1 1 -1 1 1 -1 -1 1 -1 1 1 -1 1 -1 1 -1 1 1 1 -1 1 -1 -1 1 0 0 0 0 0]';

% --- for rx 4 ---
VLTF_freq_bot(4,:) = [0 0 0 0 0 0 -1 1 -1 1 1 -1 -1 -1 -1 1 1 -1 1 1 1 -1 1 1 -1 -1 1 1 -1 1 1 -1]';
VLTF_freq_top(4,:) = [1 -1 -1 1 1 -1 -1 1 1 1 1 -1 -1 -1 1 -1 -1 1 -1 -1 -1 1 1 -1 1 -1 0 0 0 0 0]';

for rxn = 1:4
	VLTF_freq(rxn,:) = [VLTF_freq_bot(rxn,:) 0 VLTF_freq_top(rxn,:)];
end
VLTF_freq(1,:) = [VLTF_freq_bot(4,:) 0 VLTF_freq_top(4,:)];
VLTF_freq(4,:) = [VLTF_freq_bot(1,:) 0 VLTF_freq_top(1,:)];

for rxn = 1:numRxNode
	VLTF_freq(rxn,:) = LTF_freq;
	VLTF_shift(:,rxn) = fftshift(VLTF_freq(rxn,:));
end

%012013: init
for txa = 1:numTxAntenna
    VLTF_t(txa, :) = zeros(1, 64+16); 
    VLTF_t_up(txa, :) = interp(VLTF_t(txa, :), osamp);
end

%map subcarrier format for precoding to 48 data subc ID
otherTo48 =[ 0     0     0     0     0     0    25    26    27    28    29     0    30    31    32    33    34    35    36    37    38 ...
    39    40    41    42     0    43    44    45    46    47    48     0     1     2     3     4     5     6     0     7     8 ...
     9    10    11    12    13    14    15    16    17    18    19     0    20    21    22    23    24     0     0     0     0     0];


