function out = super_gaussian_1d(in , coeff, power)

% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019

    in = in(:).';
    [~, n] = size(in);
    h = -n/2: n/2-1;
    wn = n/2 * coeff;
 	
    sgauss = exp(-(h.^2/wn^2 ).^power);
    sgauss = quickscale(sgauss);
    out = in .* sgauss;

end