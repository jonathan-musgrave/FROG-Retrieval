function [iter_array, tot_initial] = iteration_number(N, snr_lev)

% [N x N, N/2 x N/2, N/4 x N/4]


if N <= 64
    if snr_lev > 51
        iter_array = [600, 10, 20]; 
        tot_initial = [4, 8, 16];
    
    elseif snr_lev > 21
        iter_array = [600, 25, 25]; 
        tot_initial = [4, 12, 16];  
        
    elseif snr_lev > 16
        iter_array = [600, 25, 25]; 
        tot_initial = [4, 12, 20];
    else
    	iter_array = [600, 30, 35]; 
        tot_initial = [4, 20, 32];
    end

    
elseif N == 128
    
	if snr_lev > 51
        iter_array = [600, 15, 15]; 
        tot_initial = [4, 12, 20];
        
	elseif snr_lev > 21
        iter_array = [600, 25, 25];
        tot_initial = [4, 12, 24];
       
    
	elseif snr_lev > 16
    	iter_array = [600, 30, 35]; 
        tot_initial = [4, 16, 32];
        
	else
        iter_array = [600, 30, 35]; 
        tot_initial = [4, 20, 40];
	end
 
   
elseif N == 256
    
	if snr_lev > 51
        iter_array = [600, 15, 15]; 
        tot_initial = [4, 12, 24];
        
	elseif snr_lev > 21
        iter_array = [600, 30, 30]; 
        tot_initial = [4, 16, 32];
    
	elseif snr_lev > 16
    	iter_array = [600, 30, 35]; 
        tot_initial = [4, 20, 36];
        
	else 
        iter_array = [600, 30, 35]; 
        tot_initial = [4, 24, 48];
	end
     
    
elseif N == 512
    
	if snr_lev > 51
        iter_array = [600, 25, 30]; 
        tot_initial = [4, 20, 32];
        
    elseif snr_lev > 21
        iter_array = [600, 30, 35];
        tot_initial = [4, 20, 36];
        
    elseif snr_lev > 16
        iter_array = [600, 30, 35];
        tot_initial = [4, 28, 44];
        
    else 
        iter_array = [600, 30, 35];
        tot_initial = [4, 32, 52];
        
	end
    
elseif N == 1024
    
    if snr_lev > 51
        iter_array = [600, 25, 30]; 
        tot_initial = [4, 20, 32];
     
    elseif snr_lev > 21
        iter_array = [600, 30, 30]; 
        tot_initial = [4, 24, 36];
    
    elseif snr_lev > 16
    	iter_array = [600, 30, 35]; 
        tot_initial = [4, 24, 40];
        
    else
        iter_array = [600, 30, 35]; 
        tot_initial = [8, 32, 52];
    end
    

elseif N == 2048  
    
    if snr_lev > 51
        iter_array = [600, 30, 30];
        tot_initial = [4, 28, 42];
    
     
    elseif snr_lev > 21
        iter_array = [600, 30, 30]; 
        tot_initial = [4, 28, 48];
    
    elseif snr_lev > 16
    	iter_array = [800, 30, 35]; 
        tot_initial = [4, 28, 52];
        
    else
        iter_array = [800, 35, 35]; 
        tot_initial = [8, 36, 60];
    end
    
    
elseif N == 4096
    
    iter_array = [800, 35, 40];
    tot_initial = [8, 36, 60];
    
else
    
    iter_array = [1000, 35, 40];
    tot_initial = [8, 40, 80];
    
end
    
    