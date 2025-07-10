function [Elam,lam_eq]=Ew2l(Ew_in,f_in)

% Author: Rana Jafari
% Georgia Tech
% email: rjafari7@gatech.edu
% June 2019 

% w, frequency axis in units of s^(-1)
    c = 3e8;        

  
    f = f_in(f_in>0);
    Ew = Ew_in(f_in>0);
    n = length(f);
    f = f(:);
    Ew = Ew(:);

    [f, order] = sort(f);
    Ew = Ew(order);

    Sw = abs(Ew).^2;
    Pw = unwrap(angle(Ew));

    Sw = Sw.*f.^2;


    lam = (c./f);


    lam_min = c/max(f);
    lam_max = c/min(f);
    lam_eq = linspace(lam_min,lam_max,n);

    Slam = interp1(lam, Sw, lam_eq, 'pchip');
    Slam(~isfinite(Slam)) = 0;

    Plam = interp1(lam, Pw, lam_eq, 'pchip');
    
    Elam = abs(sqrt(Slam)).*exp(1i*Plam);
    
    Elam = quickscale(Elam);

end