function EsigOut = PadFROG_fft(Esig,npad,dim);
    % Center padds array on the dim 1 and returns the padded array for an
    % fftc 
    % dim 1 adds npad/2 0's on the rows thus increasing the row size too
    % npad+nt;
    M = size(Esig);
    
    nt = max(M);
    if dim == 1;
       EsigOut = padarray(Esig,[floor((npad-nt)/2),0],0,'both');
    else
       EsigOut = padarray(Esig,[0,floor((npad-nt)/2)],0,'both');
    end
end