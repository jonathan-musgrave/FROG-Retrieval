function [Et_final, Gout, Gpout, Sout] = RANA(I_FROG, time, G_cut, Gprime_cut) 


% Author: Rana Jafari
% Georgia Tech
% email: rjafari7@gatech.edu
% Feb 2020
%
% Revised by Shu-Wei Huang for ECEN4006/5006
%


% inputs: FROG traces, I_FROG more noise is removed compared to I_FROG_Mw
% time: delay vector, freq: freq vector, 
% G_cut: cutoff value for G error, Gprime_cut: cutoff value G prime error 

I_FROG_noisy = I_FROG;

[I_FROG, I_FROG_Mw] = noise_reduction(I_FROG);
% I_FROG_Mw = I_FROG;

tic
I_FROG = quickscale(I_FROG);
[N,~] = size(I_FROG);
N_org = N;    


SNR_level = snr(I_FROG, I_FROG_noisy - I_FROG);


[iter_array, tot_initial] = iteration_number(N, SNR_level);    % total number of iterations and number of initial guesses based on the trace size

Isig_wT{1} = I_FROG;

time = sort(time);      
dt = mean(diff(time));
h = -N_org/2:N_org/2-1;
freq = 1/N_org/dt*h; % SHG Frog 

time = time(:).';       t{1} = time;              
freq = freq(:).';       f{1} = freq;


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% RETRIVING SPECTRUM FROM THE FREQUENCY MARGINAL

Mw = trapz(quickscale(I_FROG_Mw),2);

Swr{1} = spec_ret_marginal(Mw, 1, 0); % Mw: frequency marginal, 0/1 : noise, 0/1: 4 spectra/ 1 spectrum 

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% GENERATING BINNED ARRAYS FOR EACH GRID

[Isig_wT, t, f, Swr] =   multi_grid(Isig_wT{1}, t{1}, f{1}, Swr{1}) ;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

k = max(size(t));      % number of grids
N = length(t{k});

% INITIAL GUESSES FOR SMALLEST GRID
num_spec = size(Swr{1},1);

% initial guesses:

for jj =  tot_initial(k):-num_spec:1
   if num_spec == 1 
        et03(jj,:) = ifftc(sqrt(abs(Swr{k})).'.*exp(1i*2*pi*rand(1,N))); 
   else
        et03(jj-num_spec+1:jj,:) = ifftc(sqrt(abs(Swr{k})).*exp(1i*2*pi*rand(num_spec ,N)),[],2);  
   end
end

        
% MULTI-GRID

for i_ret =  k :-1:1

	if i_ret > 3
    	disp('Too many grids')

	elseif i_ret == 3   % N/4 by N/4 grid
                iter = iter_array(i_ret);
               
                G = 1e-6;
                options = struct('Display',		'none', ...
                    'DisplayStep',	10, ...
                    'DisplayLimit',	100, ...
                    'Domain',		'time', ...
                    'CorrectG',		'on', ...
                    'KeepBest',		'g', ...
                    'MaxIter',		iter, ...
                    'MinimumG',		G, ...
                    'MinimumZ',     0.0 );
                
                shg_amp = abs(sqrt(Isig_wT{i_ret}));
                
                parfor i_initial = 1 :tot_initial(i_ret)
                    [E_out(:,i_initial), ~, ~, Gout] =  ...
                        QuickFROG_tT_mod_1(shg_amp, et03(i_initial,:), [], [], options, Gprime_cut);  % Phase-retrieval algorithm 
                    
                    measure(i_initial) = Gout;
                end       


	elseif i_ret == 2   % N/2 by N/2 grid
                % % % % % % % % % 
                iter = iter_array(i_ret); 
                
                G = 1e-6;
                options = struct('Display',		'none', ...
                    'DisplayStep',	10, ...
                    'DisplayLimit',	100, ...
                    'Domain',		'time', ...
                    'CorrectG',		'on', ...
                    'KeepBest',		'g', ...
                    'MaxIter',		iter, ...
                    'MinimumG',		G, ...
                    'MinimumZ',     0.0 );
                
                Gout = [];      measure = [];       E_out = [];
                shg_amp = abs(sqrt(Isig_wT{i_ret}));
                et02 =  E_t0{i_ret};
               
                parfor i_initial = 1 :tot_initial(i_ret)
                    [E_out(:,i_initial), ~, ~, Gout] =  ...
                            QuickFROG_tT_mod_1(shg_amp, et02(i_initial,:), [], [], options, Gprime_cut);	% Phase-retrieval algorithm 
                        
                     measure(i_initial) = Gout;
                    
                end   
          

	elseif i_ret == 1 % N by N grid
                % iter = 100;
                G = G_cut;

                options = struct('Display',		'print', ...
                    'DisplayStep',	100, ...
                    'DisplayLimit',	400, ...
                    'Domain',		'time', ...
                    'CorrectG',		'on', ...
                    'KeepBest',		'g', ...
                    'MaxIter',		Inf, ...
                    'MinimumG',		G, ...
                    'MinimumZ',     0.0 );

                Gout = [];      measure = [];       E_out = [];

                et01 =  E_t0{i_ret};        shg_amp = abs(sqrt(Isig_wT{i_ret}));
                
                for i_initial = 1: tot_initial(1)
                     if condition == 0
                        [Et_retr_b(:,i_initial), Et_retr(:,i_initial), Etp_retr(:,i_initial), G_f(i_initial), Gp_f(i_initial), k1(i_initial)] = ...
                                QuickFROG_tT_mod_2( shg_amp, et01(i_initial,:), t{i_ret}, f{i_ret}, options, Gprime_cut);	
                                % Phase-retrieval algorithm 
                        
                        if G_f(i_initial) <= G || Gp_f(i_initial) <= Gprime_cut
                            condition = 1;
                            break
                        else
                            condition = 0;
                        end
                     end
                end
                [val, ind1] = sort(G_f);
                G_out_dum = val(1);
                kout = k1(ind1(1));

                Gpout = Gprime_error(shg_amp, Et_retr(:,ind1(1)));
                
                [valp, indp1] = sort(Gp_f);
                Gp_sure = valp(1);
        
                E_in = Et_retr_b;     Et_retr_dum = Et_retr; Et_retr = [];
                
                if G_f(i_initial) <= G
                    Et_final = E_in(:,ind1(1));
                     num_guess = ind1(1);
                    
                elseif  Gp_sure <= Gprime_cut
                    Et_final = Etp_retr(:,indp1(1));
                    Gpout = Gp_sure;
                    G_out_dum = G_error(shg_amp, Et_final);

                    
                else
                    Et_final = Et_retr_b(:,ind1(1));

                end

                
                iter = iter_array(1)-40;
                iter = Inf;
                G_f = [];   
                options = struct('Display',		'none', ...
                    'DisplayStep',	400, ...
                    'DisplayLimit',	100, ...
                    'Domain',		'time', ...
                    'CorrectG',		'on', ...
                    'KeepBest',		'g', ...
                    'MaxIter',		iter, ...
                    'MinimumG',		G, ...
                    'MinimumZ',     0.0 );
               
                if condition ~= 1
                	parfor i_initial = 1: tot_initial(1)
                    	[Et_retr(:,i_initial),~, Etp_retr(:,i_initial), G_f(i_initial), Gpbest(i_initial), k2(i_initial)] = ...
                                        QuickFROG_tT_mod_2( shg_amp, Et_retr_dum(:,i_initial), t{i_ret}, f{i_ret}, options, Gprime_cut);   
                                        % Phase-retrieval algorithm 
                    end
                    
                    
                    for i_initial = 1: tot_initial(1)
                        
                        if G_f(i_initial) < G_out_dum
                            Et_final = Et_retr(:,i_initial);
                            G_out_dum = G_f(i_initial);
                            Gpout = Gprime_error(shg_amp, Et_final);

                        end
                         
                        if G_f(i_initial) <= G 
                            G_out_dum = G_f(i_initial);
                            Et_final = Et_retr(:,i_initial);
                            Gpout = Gprime_error(shg_amp, Et_final);
                            kout = k2(i_initial);
                            condition = 1;

                        elseif Gpbest(i_initial) <= Gprime_cut
                             
                            Gpout =  Gpbest(i_initial);
                            Et_final = Etp_retr(:,i_initial);
                            G_out_dum = G_error(shg_amp, Et_final);
                            condition = 1;
                            
                        end
                    end
                end

                
	end
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %                
    %                   % NEXT INITIAL GUESS %
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
                if i_ret > 1
                    et0 = [];
                    [~, ind_mes] = sort(measure);
    
                    
                    for i_case = 1: tot_initial(i_ret-1)

                        E_t0 = [];
              
                        E_proc = E_out(:,ind_mes(i_case));       E_proc = center(E_proc,'moment');
                        
                        E_proc_w = fftc(E_proc);

                        real_dum = interp1( f{i_ret}, real(E_proc_w), f{i_ret-1}, 'linear');
                        imag_dum = interp1( f{i_ret}, imag(E_proc_w), f{i_ret-1}, 'linear');
                        
                        real_dum(~isfinite(real_dum)) = 0;
                        imag_dum(~isfinite(imag_dum)) = 0;


                        phase = angle(real_dum+1i*imag_dum);
                        S_dum = Swr{i_ret-1};
                        
                        rms_a = [];         ew_int_a = [];
                        
                   
                        ew_int_b = real_dum+1i*imag_dum;
                        rms_b = compare_Gerror(ew_int_b,sqrt(Isig_wT{i_ret-1}));

                        for ii = 1: num_spec
                                ew_int_a(ii,:) = abs(sqrt(S_dum(ii,:))).* exp(1i*phase); 
                                rms_a(ii) = compare_Gerror(ew_int_a(ii,:),abs(sqrt(Isig_wT{i_ret-1}))); 
                                % determining whether to reapply the retrieved spectrum to the initial guess
                        end
                        
                              
                        [~, ind_rms] = min(rms_a);
                        
                        if rms_a(ind_rms) < rms_b
                            
                                ew_int = ew_int_a(ind_rms,:);  
                        else 
                                ew_int = ew_int_b;    
                        end
                        et0(i_case,:) = ifftc(ew_int);

                    end
                    
                    E_t0{i_ret-1} = et0;                
                    N = 2*N;
                    condition = 0;
                end
end

Sout = Swr{1};

Gout = G_error(abs(shg_amp), Et_final);
Gpout = Gprime_error(abs(shg_amp), Et_final);


end
     



   

