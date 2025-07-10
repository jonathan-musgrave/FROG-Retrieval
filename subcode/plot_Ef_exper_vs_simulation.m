function [] = plot_Ef_exper_vs_simulation(Et_ret,Et_sim,freq_exper,freq_sim)
% Author: Jonathan Musgrave
% University of Colorado Boulder
% email: jomu3154@colorado.edu
% Oct 2024


p = 4;
N = length(Et_ret);
Ef_ret = fftshift(fft(fftshift(Et_ret)))./N;
Ef_sim = fftshift(fft(fftshift(Et_sim)))./N;

I_ret = abs(Ef_ret).^2;             I_ret = quickscale(I_ret);
phase_f_exper = -unwrap(angle(Ef_ret));
phase_f_exper(I_ret<max(I_ret)*5e-3) = NaN;
phase_f_exper = phase_f_exper - (phase_f_exper(N/2+1));

I_sim = abs(Ef_sim).^2;             I_sim = quickscale(I_sim);
phase_f_sim = -unwrap(angle(Ef_sim));
phase_f_sim(I_sim<max(I_sim)*5e-3) = NaN;
phase_f_sim = phase_f_sim - (phase_f_sim(N/2+1));



yyaxis left
hold all;
ax = gca;
alph = 0.5;
plot(freq_exper, I_ret,'Color',[1,0,0,alph],'LineWidth', 2,'Displayname','I(t)')
plot(freq_sim, I_sim,':','Color',[0,0,1,alph],'LineWidth', 2,'Displayname','I^{sim}(t)')

ylabel('Intensity (a.u.)','Color',[1,0,0]);
xlabel('Frequency (THz)'); 
ax.YColor = [1,0,0];

axis tight;
ax.XLim = [min(freq_sim) max(freq_sim)];
ax.YTick = [0,1];
ax.YLim = [0,1];
yyaxis right
hold all;
plot(freq_exper, phase_f_exper,'Color',[1, 0.5, 0],'LineWidth', 2,'handlevisibility','off')
plot(freq_sim, phase_f_sim,':','Color',[0 0.5 1],'LineWidth', 2,'handlevisibility','off')
ax = gca;
ylabel({'Phase (\pi rad)'},'Color',[0,0,1]); 

axis tight
ax.YColor = [0,0,1];
ax.LineWidth = 1.5;
ax.FontSize = 14; %font size
ax.YLim = [min([phase_f_sim';phase_f_exper]) max([phase_f_sim';phase_f_exper])];
ax.YTick = min([phase_f_sim';phase_f_exper]):4*pi:max([phase_f_sim';phase_f_exper]);
ax.YTickLabel=[-p:4:p];
box on;

end