function Isig_wT = shg_trace_rec(et)
   

% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019


Esig_tT = CalcEsig(et,et);
Esig_wT = fft_FROG(Esig_tT);
Isig_wT = abs(Esig_wT).^2;
    


end

