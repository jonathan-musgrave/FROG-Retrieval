function [EtB, Et, Gbest] = QuickFROG_tT_SD(Asig, Et, t, w, w2, options)

Asig = quickscale(Asig);

%	Default options.
defaultopt = struct('Display',		'all', ...
	'DisplayStep',	25, ...
	'DisplayLimit',	100, ...
	'Domain',		'time', ...
	'CorrectG',		'on', ...
	'KeepBest',		'g', ...
	'MaxIter',		1200, ...
	'MinimumG',		1e-5, ...
	'MinimumZ',		0.0 );

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
G(1,1) = Inf;
Z(1)   = Inf;
Gbest = Inf;
Zbest = Inf;

tic;

Esig = CalcEsig(conj(Et), Et.*Et);

Esig = fft_FROG(Esig);

Esig = MagRepl(Esig,Asig);

Esig = ifft_FROG(Esig);

while G(end,1) > Gcutoff & Z(end) > Zcutoff & k < n
	k = k+1;
	
	dZ = -dZdE_sd(Esig,Et);
	
	[Et, Z(k)] = MinZerr_sd(Esig, Et, dZ);
	
	if ~mod(k-1,stepsz)
		Et = center(Et,'max');
    end
	
    Esig = CalcEsig(conj(Et), Et.*Et);
	
	Esig = fft_FROG(Esig);
	
	[G(k,1),a(k)] = MinGerr(Esig,Asig);
	
	if CorrectG
		Et = Et * (a(end))^(1/6);
	end
	
	if KeepBest & G(k,1) <= Gbest
		Gbest = G(k,1);
		Zbest = Z(k);
		EtB = Et;
	elseif Z(k) <= Zbest
		Gbest = G(k,1);
		Zbest = Z(k);
		EtB = Et;
	end
	
	if dsp & ~mod(k-1,stepsz)
		Ew = fftc(Et);
		DisplayFROG(Asig, Esig, Et, Ew, t, w, k, Klimit, G, Z, Gbest, Zbest, Gcutoff, Zcutoff)
	end
	
	if prnt & ~mod(k-1,stepsz)	
		PrintFROG(k, G, Z, Gcutoff, Zcutoff, toc)
	end
	
	Esig = MagRepl(Esig,Asig);
	
	Esig = ifft_FROG(Esig);
end

k = k+1;

dZ = -dZdE_sd(Esig,Et);

[Et, Z(k)] = MinZerr_sd(Esig, Et, dZ);

Esig = CalcEsig(conj(Et), Et.*Et);

Esig = fft_FROG(Esig);

[G(k,1),a(k)] = MinGerr(Esig,Asig);

Et = Et * a(end)^(1/6);

if dsp
	Ew = fftc(Et);
	DisplayFROG(Asig, Esig, Et, Ew, t, w, k, Klimit, G, Z, Gbest, Zbest, Gcutoff, Zcutoff)
end

if prnt
	PrintFROG(k, G, Z, Gcutoff, Zcutoff, toc)
end

beep

