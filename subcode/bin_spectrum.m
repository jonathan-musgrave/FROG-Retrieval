function [var_out] = bin_spectrum(var_in, f1, f2)

% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019

var_out = interp1(f1, var_in.', f2,'pchip');

var_out = var_out.';

end
   
 