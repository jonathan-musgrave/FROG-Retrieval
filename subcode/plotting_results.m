function [] = plotting_results(Ishg, et_ret, time, lam0)

% Author: Rana Jafari
% Georgia Tech
% email: rjafari7@gatech.edu
% June 2019

Ishg = quickscale(Ishg);

et_ret = center(et_ret,'moment');       % Retrieved temporal field

ew_ret = fftc(et_ret);

N = length(et_ret);

Ishg_ret = shg_trace_rec(et_ret);       % Reconstructed trace


fs = 1e-15;
THz = 1e12;
nm = 1e-9;

time=sort(time);
dt = mean(diff(time));
dt = dt*fs;

freq = ((-N/2:N/2-1)/N/dt + 3e8/(lam0*1e-9))/THz;
freq2 =((-N/2:N/2-1)/N/dt + 2*3e8/(lam0*1e-9))/THz;

[Elam, lam] = Ew2l(ew_ret, freq*THz);

lam = lam/nm;

fig_num = 101;
hFig = figure(fig_num);
set(hFig, 'Position', [100 100 1300 550]); 


subplot(2,3,1);	plot_trace(Ishg, time, freq2);                  title('Measured trace');    

subplot(2,3,2);	plot_trace(Ishg_ret, time, freq2);              title('Retrieved trace');  

subplot(2,3,3);	plot_trace(Ishg-Ishg_ret, time, freq2);         title('Difference');             
 
subplot(2,3,4);	plot_Et(et_ret, time);

subplot(2,3,5);	plot_Ew(ew_ret, freq);

subplot(2,3,6);	plot_Elam(Elam, lam);

