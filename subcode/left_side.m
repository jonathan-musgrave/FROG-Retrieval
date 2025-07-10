function [Swc, rms_err1, rms_err2] = left_side(Mw, St, pos_neg,varargin)

   
    if isempty(varargin)
        npad = [];
    else
        npad = varargin{1};
    end

    
% Authors: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019


    N = length(pos_neg); 
    
    signs = pos_neg;    
    signs(N/2+2:end) = fliplr(pos_neg(2:N/2));

    Stc = St.* signs;
    Swc = fftc(Stc,npad);
    Swc = super_gaussian_1d(Swc, .65,5);
    
    A = abs(Swc);  
    autoS = quickscale(abs(conv(A, A,'same'))); 
    rms_err1 = mean(abs(autoS- Mw).^2); 
    
    [~ , indmax]  = max(abs(real(Swc)));
    Swc = sign(real(Swc(indmax))) * Swc; 
    
    negative_values = zeros(size(Swc));
    negative_values(real(Swc) < 0) = Swc(real(Swc) < 0);
    rms_err2 = trapz(abs(negative_values));
    
    Swc = quickscale(abs(Swc));
 
end