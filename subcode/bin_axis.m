function var_out = bin_axis(coeff, dim, dt, option)

% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019

switch option
    
    case 1 % time
        h = -dim/2:dim/2-1;
        var_out = (coeff*dt)*h;
    
    case 2 % freq
        h = -dim/2:dim/2-1;
        df = (1/dim/dt/coeff);
        var_out = df*h;        
end

end