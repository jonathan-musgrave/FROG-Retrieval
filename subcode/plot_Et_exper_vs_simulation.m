
function [] = plot_Et_exper_vs_simulation(Et_ret,Et_sim,time_exper,time_sim)
% Author: Jonathan Musgrave
% University of Colorado Boulder
% email: jomu3154@colorado.edu
% Oct 2024

p = 4;
N = length(Et_ret);


I_ret = abs(Et_ret).^2;             I_ret = quickscale(I_ret);
phase_t_ret = -unwrap(angle(Et_ret));
phase_t_ret(I_ret<max(I_ret)*5e-3) = NaN;
phase_t_ret = phase_t_ret - phase_t_ret(N/2+1);

I_sim = abs(Et_sim).^2;             I_sim = quickscale(I_sim);
phase_t_sim = -unwrap(angle(Et_sim));
phase_t_sim(I_sim<max(I_sim)*5e-3) = NaN;
phase_t_sim = phase_t_sim - phase_t_sim(N/2+1);



yyaxis left
hold all;
ax = gca;
alph = 0.75;
plot(time_exper, I_ret,'Color',[1,0,0,alph],'LineWidth', 2,'Displayname','I(t)')
plot(time_sim, I_sim,':','Color',[0,0,1,alph],'LineWidth', 2,'Displayname','I^{sim.}(t)')

ylabel('Intensity (a.u.)','Color',[1,0,0]);
xlabel('Time (fs)'); 
ax.YColor = [1,0,0];

axis tight;
ax.XLim = [min(time_sim) max(time_sim)];
ax.YTick = [0,1];
ax.YLim = [0,1];
yyaxis right
hold all;
plot(time_exper, phase_t_ret,'-','Color',[1, 0.5, 0],'LineWidth', 2,'handlevisibility','off')
plot(time_sim, phase_t_sim,':','Color',[0 0.5 1],'LineWidth', 2,'handlevisibility','off')
ax = gca;
ylabel({'Phase (\pi rad)'},'Color',[0,0,1]); 

axis tight
ax.YColor = [0,0,1];
ax.LineWidth = 1.5;
ax.FontSize = 14; %font size
ax.YLim = [-p*pi p*pi];
ax.YTick = -p*pi:4*pi:p*pi;
ax.YTickLabel=[-p:4:p];
box on;

end