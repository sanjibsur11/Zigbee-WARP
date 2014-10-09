% MMSE Finite Impulse Response interpolator for the input FM demodulated
% samples. Derived from gnuradio block mmse_fir_interpolator_ff
% 
% Author: Sanjib Sur
% Institution: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 01/14/2014
% 
% Comments:

function [filtered_out] = mmse_fir_interpolator(raw_samples, mu)

global MM_FILTER_NSTEPS MM_FILTER_NTAPS FILTER_TAPS;

% Make sure the number of raw samples are equal to the number of taps in
% the MMSE filter
assert(length(raw_samples) == MM_FILTER_NTAPS);

% Select correct tap array from mu value
imu = round(mu * MM_FILTER_NSTEPS) + 1;
% Dot product of input samples with filter taps
filtered_out = FILTER_TAPS(imu, :) * raw_samples;

end