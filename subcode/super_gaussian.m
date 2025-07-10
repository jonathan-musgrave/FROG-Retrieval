function out = super_gaussian(in , coeff_n, coeff_m, power)
  
% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019


    [n, m] = size(in);
    hn = -n/2: n/2-1;
    hm = -m/2: m/2-1;
    [M , N] = meshgrid(hm,hn);

    wn = n/2 * coeff_n;
    wm = m/2 * coeff_m;
    sgauss = exp(-(N.^2/wn^2 + M.^2/wm^2).^power);
    sgauss = sgauss +1e-6; 
    sgauss = quickscale(sgauss);

    out = in .* sgauss;

end