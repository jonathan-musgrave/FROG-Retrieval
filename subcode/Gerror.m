function g = Gerror(N)

% Authors: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019

e = 1e-2 ; % G error for 32 x 32 trace
g = e *(32/N)^(1/2);

end