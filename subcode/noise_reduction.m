function [x, x_S] = noise_reduction(x,varargin)
 
% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019
  
   
    if isempty(varargin)
        npad = [];
    else
        npad = varargin{1};
    end
    
    x = super_gaussian(x,2,2,2);
 
    fx = ifft_FROG(x,npad);
    fx = super_gaussian(fx,.7,100,10);
    x = fft_FROG(fx,npad);
     
     
    fx = fftc(x,[],2);      
    fx = super_gaussian(fx,100,.7,10);
    x = abs(ifftc(fx,[],2));      
   
    perm = 3;
    x = frogbacksub_mean( x,  [perm perm] , [perm perm]);
    x = remove_neg(x);
 
 
    x_S = quickscale(x);  
    
    
    fx = ifft_FROG(x,npad);
    fx = super_gaussian(fx,.7,100,10);
    x = fft_FROG(fx,npad);
 
 
    fx = fftc(x,[],2);      
    fx = super_gaussian(fx,100,.7,10);
    x = abs(ifftc(fx,[],2));      
 
    perm = 4;
    
    x = frogbacksub( x,  [perm perm] , [perm perm]);
 
    x = remove_neg(x);
    
    x = frogbacksub_mean( x,  [perm perm] , [perm perm]);
 
    x = remove_neg(x);
 
 
    x = quickscale(x);
 
end


  