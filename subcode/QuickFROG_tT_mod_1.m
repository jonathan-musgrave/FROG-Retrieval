function [EtB, Et, EtBp, G_f ] = QuickFROG_tT_mod_1(Asig, Et, t, w, options, Gpcutoff)
 
 
%   Default options.
defaultopt = struct('Display',      'all', ...
    'DisplayStep',  25, ...
    'DisplayLimit', 100, ...
    'Domain',       'time', ...
    'CorrectG',     'on', ...
    'KeepBest',     'g', ...
    'MaxIter',      500, ...
    'MinimumG',     0.0, ...
    'MinimumZ',     0.0,...
    'npad',         []);
 
% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(Asig,'defaults')
    EtB = defaultopt;
    return
end
 
if nargin < 5; options = []; end
 
n = FROG_optimget(options, 'MaxIter', defaultopt, 'fast') - 1;
Gcutoff = FROG_optimget(options, 'MinimumG', defaultopt, 'fast');
Zcutoff = FROG_optimget(options, 'MinimumZ', defaultopt, 'fast');
stepsz = FROG_optimget(options, 'DisplayStep', defaultopt, 'fast');
Klimit = FROG_optimget(options, 'DisplayLimit', defaultopt, 'fast');
CorrectG = strcmp(FROG_optimget(options, 'CorrectG', defaultopt, 'fast'), 'on');
npad =  FROG_optimget(options, 'npad', defaultopt, 'fast');

switch FROG_optimget(options, 'Display', defaultopt, 'fast')
case 'all'
    dsp = 1;
    prnt = 1;
case 'print'
    prnt = 1;
    dsp = 0;
case 'graph'
    dsp = 1;
    prnt = 0;
otherwise
    dsp = 0;
    prnt = 0;
end
 
switch FROG_optimget(options, 'KeepBest', defaultopt, 'fast')
case 'g'
    KeepBest = 1;
otherwise
    KeepBest = 0;
end
 
 
k = 0;
G(1,1) = Inf;   Gprime(1,1) = Inf;
Z(1)   = Inf;
Gbest  = Inf;
Zbest  = Inf;
Gpbest = Inf; 
 
% x = [-128:127]; Et = exp(-(x./10).^2);
Esig = CalcEsig(Et,Et); % Caluclate Esig of FRG
 
Esig = fft_FROG(Esig,npad);
 
Esig = MagRepl(Esig,Asig);
 
Esig = ifft_FROG(Esig,npad);
 
error_cond  = G(1,1) > Gcutoff;
 

while error_cond &&  Z(end) > Zcutoff && k < n
 
    k = k+1;
    
    dZ = -dZdE_shg(Esig,Et);
 
    
    [Et, Z(k)] = MinZerr_shg(Esig, Et, dZ);
 
    if ~mod(k-1,stepsz)
        Et = center(Et,'moment'); %remove centering for nsFROG simulation,SR
    end
    
    Esig = CalcEsig(Et,Et);
    
    Esig = fft_FROG(Esig,npad);
    
    [G(k,1),a ] = MinGerr(Esig,Asig);

       
    Gprime(k,1) = sqrt(trapz(trapz(abs(Asig.^2 - a *abs(Esig).^2).^2))/trapz(trapz(Asig.^4)));
 
    
    if CorrectG
        Et = Et * sqrt(sqrt(a));
    end 
 
    
	if KeepBest && G(k,1) <= Gbest
		Gbest = G(k,1);
		Zbest = Z(k);
		EtB = Et;
	elseif Z(k) <= Zbest
		Gbest = G(k,1);
		Zbest = Z(k);
		EtB = Et;
	end
    
    if KeepBest && Gprime(k,1) <= Gpbest
        Gpbest = Gprime(k,1);
        EtBp = Et;
    end
    
    if dsp && ~mod(k-1,stepsz)
        Ew = fftc(Et,npad);
        DisplayFROG(Asig, Esig, Et, Ew, t, w, k, Klimit, G, Z, Gbest, Zbest, Gcutoff, Zcutoff)
    end
    
    if prnt && ~mod(k-1,stepsz)  
        PrintFROG(k, G, Z, Gcutoff, Zcutoff, toc)
    end
    
    Esig = MagRepl(Esig,Asig);
    
    Esig = ifft_FROG(Esig,npad);
        
 
    error_cond  = min(Gprime) > Gpcutoff &  min(G) > Gcutoff;

    if k > 50
        stag = mean(abs(diff(G(k-10:k-1))));
        if (stag < 1e-6 && stag < Gcutoff/100); break; end
    end


end
 
k = k+1;
 
dZ = -dZdE_shg(Esig,Et);
 
[Et, Z(k)] = MinZerr_shg(Esig, Et, dZ);
 
Esig = CalcEsig(Et,Et);
 
Esig = fft_FROG(Esig,npad);
 
[G(k,1),a] = MinGerr(Esig,Asig);
 
Et = Et * sqrt(sqrt(a));
 
if dsp
    Ew = fftc(Et,npad);
    DisplayFROG(Asig, Esig, Et, Ew, t, w, k, Klimit, G, Z, Gbest, Zbest, Gcutoff, Zcutoff)
end
 
if prnt
    PrintFROG(k, G, Z, Gcutoff, Zcutoff, toc)
end
 
G_f = Gbest;

% Gp_f = min(Gprime);
 
beep
 

