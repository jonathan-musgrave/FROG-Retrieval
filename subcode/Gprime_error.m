function Gp =  Gprime_error(Asig, Et,npad)

% Author: Rana Jafari
% Georgia Tech
% email: rjafari7@gatech.edu
% Feb 2020

Esig = CalcEsig(Et,Et);
    
Esig = fft_FROG(Esig,npad);
    
[G, a] =  MinGerr(Esig, Asig);

       
Gp = sqrt(trapz(trapz(abs(Asig.^2 - a *abs(Esig).^2).^2))/trapz(trapz(Asig.^4)));
 


end