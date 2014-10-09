% Frequency offset estimation and removal of DC offset using single pole
% IIR Filter
% 
% Author: Sanjib Sur
% Institution: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 01/14/2014
% 
% Comments: 


function [raw_samples] = remove_dco(alpha, raw_samples_wDCO)

assert(alpha < 1);

% Estimate DC offset per sample, based on previous samples
freq_offset = single_pole_iir_filter(alpha, raw_samples_wDCO);
% Remove DC offset from the samples
raw_samples = raw_samples_wDCO - freq_offset;

end