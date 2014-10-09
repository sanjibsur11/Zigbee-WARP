% Frequency compensation using LTF style preamble sequence
% 
% Author: Sanjib Sur
% Institution: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 01/17/2014
% 
% Comments: 


function [rxData_comp] = freq_compensation(rxData)

global USESIM LTF_t;
sigin = rxData;

ns = 0;
% Get the exact position of LTF via crosscorr 
k = 1;
for t = ns+10:ns+30
    crosscorrLTF(k) = abs(sum(sigin(t+1:t+16).*conj(LTF_t(1:16).')));
    k = k + 1;
end   
    
[~, peakPos] = max(crosscorrLTF);

LTFpeakPos = ns + 10 - 1 + peakPos;

% Estimate the frequency offset
co = sum(sigin((LTFpeakPos+1):(LTFpeakPos+64)).*...
    conj(sigin((LTFpeakPos+64+1):(LTFpeakPos+64+64))));

iFreq = angle(co)/64;

% Compensate for the frequency offset
ns = LTFpeakPos;
rxData_comp = zeros(1, length(rxData) - LTFpeakPos);
k = 1;
while ns < length(sigin)
    rxData_comp(k) = sigin(ns+1) * exp(1i*iFreq*(k-1));
    ns = ns + 1;
    k = k + 1;
end

rxData_comp = rxData_comp.';

end