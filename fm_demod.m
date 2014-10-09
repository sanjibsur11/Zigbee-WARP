% FM Demodulation function using gnuradio block quadrature_demod_cf
% 
% Author: Sanjib Sur
% Institution: University of Wisconsin - Madison
% Version: 0.0.1
% Last modified: 01/14/2014
% 
% Comments: 

function [delTheta] = fm_demod(quad_samples, gain)
    
    for qs = 2:length(quad_samples)
        complex_product = quad_samples(qs) * conj(quad_samples(qs-1));
        if real(complex_product) == 0
            delTheta(qs-1) = 0;
        else
            delTheta(qs-1) = gain * atan(imag(complex_product)/...
                                         real(complex_product));
        end
    end

    delTheta = delTheta.';
end