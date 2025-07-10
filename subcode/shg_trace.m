function Isig_wT = shg_trace(et,varargin)
   

% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019

    if isempty(varargin)
        npad = [];
    else
        npad = varargin{1};
    end
    Esig_tT = CalcEsig(et,et);
    Esig_wT = fft_FROG(Esig_tT,npad);
    Isig_wT = abs(Esig_wT).^2;
    Isig_wT = quickscale(Isig_wT);
    
end

