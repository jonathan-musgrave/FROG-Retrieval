function [Isig_wT, t ,f , Swr] = multi_grid(Isig_full, t_full, f_full, Swr_full,varargin)


% Authors: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019

% % GENERATING BINNED TRACES, TIME AXIS, FREQUENCY AXIS, AND SPECTRUM FOR
% % INITIAL GUESS

    Isig_wT{1} = Isig_full; t{1} = t_full;  f{1} = f_full;   Swr{1} = Swr_full;

    [N,M] = size(Isig_full);
    
    if isempty(varargin{:})
        npad = N;
    else
        npad = varargin{1};
    end

    a = 2;          % binnig coefficent : sqrt(a)
    k = 0;
    dt = mean(diff(t_full));
    df = mean(diff(f_full));

    while k < 2
        
        k = k+1;
        
        coeff = a^(k/2);

        t{k+1} = bin_axis(coeff,N/2^k, dt,1);
        
        f{k+1} = bin_axis(coeff*N/npad,N/2^k, dt,2);
                
        Isig_wT{k+1} = bin_trace(Isig_wT{1}, t{1}, f{1}, t{k+1}, f{k+1});
        
        Swr{k+1} = bin_spectrum(Swr{1}, f{1}, f{k+1});

    end

end