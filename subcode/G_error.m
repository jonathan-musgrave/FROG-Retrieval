function G =  G_error(Asig, Et,npad)

% Author: Rana Jafari
% Georgia Tech
% email: rjafari7@gatech.edu
% Feb 2020

Esig = CalcEsig(Et,Et);
    
Esig = fft_FROG(Esig,npad);
    
[G, a] =  MinGerr(Esig, Asig);


end