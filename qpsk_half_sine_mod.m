% Implementation of QPSK half-sine pulse shape modulator
% Converts 2 bits symbol to 4 samples half sine pulse
% 
% Author: Sanjib Sur
% Institution: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 01/14/2014
% 
% Comments:


function [output] = qpsk_half_sine_mod(input)

% 2 bits to QPSK constellation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                              %
%       1     |     3          %
%        *    |    *           %
%             |                %
%     -----------------        %
%             |                %
%        *    |    *           %
%       0     |     2          %
%                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
constellation = [-1-1j -1+1j 1-1j 1+1j];
% QPSK modulation
qpsk_samples = constellation(input+1);

% Half-sine pulse shape filter
half_sine_samples(1) = 0;
half_sine_samples(2) = qpsk_samples / sqrt(2);
half_sine_samples(3) = qpsk_samples;
half_sine_samples(4) = qpsk_samples / sqrt(2);

output = half_sine_samples;

end