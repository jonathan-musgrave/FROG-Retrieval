%
%  Simulated pulses from file
%
%  Author: Rana Jafari
%  Georgia Tech
%  email: rjafari7@gatech.edu
%  June 2019
%
%  Revised by Shu-Wei Huang for ECEN4006/5006

clear all;
addpath('Field\');
addpath('subcode\');
addpath('Window_64bit\');
addpath('FROG_v1.2.1\lib\frog\');
addpath('FROG_v1.2.1\lib\frog\Error\');
addpath('FROG_v1.2.1\lib\Misc\');
%%

% fs = 1e-15; THz = 1e12;
% c = 299792458;
% load('FROG_for_Rana.mat');
% Isig = FROG_Rana;
% time = delay;
% N_org = length(Isig);


%Isig = [Isig(:,1:N_org/2+1),flip(Isig(:,1:N_org/2),2)];
%%
% figure(1);clf;
% ax1 = subplot(1,2,1);
% imagesc(delay,freq,FROG); 
% subplot(1,2,2);
% imagesc(time,f,Isig); 
% xlim(get(ax1,'XLIM'))
% ylim(get(ax1,'YLIM'))


load('t_5.mat')
load('f_5.mat')
load('Et_5.mat')
Et(1,:) = exp(-(((1:128)-64.5)/10).^2);
i_pulse = 21; 

    nt = 2^7;
    npad = 2.5*nt;
    twin = 800e-15;
    dt = 2*twin./nt;
    df = 1/(2*dt*npad);
    Time = [-nt/2:nt/2-1]*dt;
    f = df*[-nt/2:nt/2-1];
    nn = 1:size(Et,1);
    [TT,NN] = meshgrid(Time,nn);
    Et = interp2(time*4,nn,Et, TT,NN);
    Et(isnan(Et)) =0;
    time = Time;
    
et = Et(i_pulse,:);

Isig = shg_trace(et,npad);
% Isig = shg_trace(et);
% npad = [];
%% new pulse
nt = 2^7;
npad = 8*nt;
twin = 800e-15;
dt = 2*twin./nt;
df = 1/(2*dt*npad);
time = [-nt/2:nt/2-1]*dt;
f = df*[-nt/2:nt/2-1];

tau_fw = 200e-15;
Et = sech(time./tau_fw);
Isig = shg_trace(Et,npad);
% 



G_cutoff = 1e-6;     Gp_cutoff = 1e-7;

% Et_final: retrieved field, Gout: Gerror, Gpout: G' error 
[Et_final, Gout, Gpout, Sout] = RANA_Experimental(abs(Isig), time, G_cutoff, Gp_cutoff,npad);
% [Et_final, Gout, Gpout, Sout] = RANA(abs(Isig), time, G_cutoff, Gp_cutoff);

Sret = Sout; % retrieved spectra from frequency marginal

%%
et_in = center(zeros(length(Et_final),1),'max');
et_ret = center(Et_final,'max');
et_ret_rev = center(conj(flipud(et_ret)),'max');


% Isig = Isig0;
freq = f;
N = length(time);
I_in = quickscale(abs(et_in).^2);       
phase_in = angle(et_in);        phase_in(I_in<max(I_in)*5e-3) = NaN;         
phase_in = unwrap(phase_in);  



Ishg_ret = shg_trace(et_ret,npad);
I_ret = quickscale(abs(et_ret).^2);     
phase_ret = angle(et_ret);     phase_ret(I_in<max(I_in)*5e-3) = NaN;      
phase_ret = unwrap(phase_ret)-phase_ret(length(phase_ret)/2+1);  
phase_ret = phase_ret-phase_ret(length(phase_ret)/2+1);

I_ret_rev = quickscale(abs(et_ret_rev).^2);     
phase_ret_rev = angle(et_ret_rev);     phase_ret_rev(I_in<max(I_in)*5e-3) = NaN;      
phase_ret_rev = unwrap(phase_ret_rev);  
phase_ret_rev = phase_ret_rev-phase_ret_rev(length(phase_ret_rev)/2+1);

fs = 1e15;
THz = 1e-12;

fig_num = 101;
trange = 600*2;
frange = 1.5*4;
phi_range = 2;
hFig = figure(fig_num);clf;
set(hFig, 'Position', [100 100 950 550]) 

subplot(2,2,1);     plot_trace(Isig, time*fs, (freq+freq(N/2+1)*0)*THz);            title('Simulated trace'); 
xlim([-1,1]*trange);
ylim([-1,1]*frange);
subplot(2,2,2);     plot_trace(Ishg_ret, time*fs, (freq+freq(N/2+1)*0)*THz);        title('Retrieved trace'); 
xlim([-1,1]*trange);
ylim([-1,1]*frange);
subplot(2,2,3);    % plot_Et_sim(et_in, et_ret ,time*fs);          
yyaxis('left');
plot(time*fs,abs(et_ret).^2,'r');
set(gca,'YCOLOR','r')
yyaxis('right'); plot(time*fs,phase_ret,'-b');
set(gca,'YCOLOR','b')
title('Retrieved field');
xlim([-1,1]*trange*2);
ylim([-1,1]*phi_range);
subplot(2,2,4);    % plot_Et_sim(et_in, et_ret ,time*fs);           
yyaxis('left')
plot(time*fs,abs(et_ret_rev).^2,'-r');
set(gca,'YCOLOR','r')
yyaxis('right'); plot(time*fs,phase_ret_rev,'-b');
set(gca,'YCOLOR','b')
title('Retrieved field');
xlim([-1,1]*trange*2);
ylim([-1,1]*phi_range);
