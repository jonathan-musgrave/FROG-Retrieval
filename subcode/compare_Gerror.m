function rms_err = compare_Gerror(ew_int, Asig,varargin)
   
% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019
    if isempty(varargin)
        npad = [];
    else
        npad = varargin{1};
    end
    


    field = ifftc(ew_int,npad);  
    
    Asig_wT = quickscale(abs(fft_FROG(CalcEsig(field, field),npad)));   
    
    [rms_err, ~] = MinGerr(Asig_wT, Asig);

end