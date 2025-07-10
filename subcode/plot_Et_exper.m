
function [] = plot_Et_exper(Et_ret,Et_ret_rev,time)
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

I_ret_rev = abs(Et_ret_rev).^2;             I_ret_rev = quickscale(I_ret_rev);
phase_t_ret_rev = -unwrap(angle(Et_ret_rev));
phase_t_ret_rev(I_ret_rev<max(I_ret_rev)*5e-3) = NaN;
phase_t_ret_rev = phase_t_ret_rev - phase_t_ret_rev(N/2+1);



yyaxis left
hold all;
ax = gca;
alph = 0.75;
plot(time, I_ret,'Color',[1,0,0,alph],'LineWidth', 2,'Displayname','I(t)')
plot(time, I_ret_rev,'-','Color',[0,0,1,alph],'LineWidth', 2,'Displayname','I^{rev.}(t)')

ylabel('Intensity (a.u.)','Color',[1,0,0]);
xlabel('Time (fs)'); 
ax.YColor = [1,0,0];

axis tight;
ax.XLim = [min(time) max(time)];
ax.YTick = [0,1];
ax.YLim = [0,1];
yyaxis right
hold all;
plot(time, phase_t_ret,':','Color',[1, 0.5, 0],'LineWidth', 2)
plot(time, phase_t_ret_rev,':','Color',[0 0.5 1],'LineWidth', 2)
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