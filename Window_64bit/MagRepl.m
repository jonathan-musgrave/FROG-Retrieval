function Esig = MagRepl(Esig, Asig)
% MagRepl.m
% MAGREPL Replace the magnitude of a FROG trace.
%	MAGREPL(ESIG,ASIG) Returns a matrix with the magnitude of ASIG and the
%	phase of ESIG.  Where ASIG < 0 this function returns ESIG.

%	$Id: MagRepl.m,v 1.1 2006-11-11 00:15:29 pablo Exp $
% MAGREPL Replace the magnitude of a FROG trace.
%	MAGREPL(ESIG,ASIG) Returns a matrix with the magnitude of ASIG and the
%	phase of ESIG.  Where ASIG < 0 this function returns ESIG.

% The MEX version of this file found some issues in the compile process.
% The problem is caused by the standard library included in gcc conflict
% with some definition in the header file.
% This m-file do the same operation.

%	$Id: MagRepl.m,v 1.1 2006-11-11 00:15:29 pablo Exp $
%
% By Jeff Wong (GaTech) - 2012-06-02, 0129

% By Rana Jafari (GaTech) - 2019-03-17, 


temp = abs(Esig);
Asig(real(Asig)<0) = abs(Esig(real(Asig)<0));

Esig(temp>0) = Asig(temp>0) .* Esig(temp>0) ./ temp(temp>0);
Esig(temp==0) = Asig(temp==0);


end