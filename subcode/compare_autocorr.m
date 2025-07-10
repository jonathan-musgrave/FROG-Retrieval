function rms_err = compare_autocorr(ew_int, M_tau_trace)

% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019

    field = ifftc(ew_int);       
    
    Isig_tT = CalcEsig(abs(field).^2, abs(field).^2);   
    
    M_tau = sum(Isig_tT);       M_tau = abs(quickscale(M_tau));

    rms_err = sqrt(mean(abs(M_tau - M_tau_trace).^2));      
    
end