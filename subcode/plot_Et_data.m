
function [] = plot_Et_sim(Et_sim, Et_ret ,time)

% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019


p = 4;
N = length(Et_sim);

I_sim = abs(Et_sim).^2;             I_sim = quickscale(I_sim);

phase_t_sim = -unwrap(angle(Et_sim));
phase_t_sim(I_sim<max(I_sim)*5e-3) = NaN;
phase_t_sim = phase_t_sim - phase_t_sim(N/2+1);



I_ret = abs(Et_ret).^2;             I_ret = quickscale(I_ret);

phase_t_ret = -unwrap(angle(Et_ret));
phase_t_ret(I_ret<max(I_ret)*5e-3) = NaN;
phase_t_ret = phase_t_ret - phase_t_ret(N/2+1);


yyaxis left
hold all;
ax = gca;

plot(time, I_sim,'Color',[1,0.5,0],'LineWidth', 2)
plot(time, I_ret,':','Color',[1,0,0],'LineWidth', 2)

ylabel('Intensity (a.u.)','Color',[1,0,0]);
xlabel('Time (fs)'); 
ax.YColor = [1,0,0];

axis tight;
ax.XLim = [min(time) max(time)];
ax.YTick = [0,1];
ax.YLim = [0,1];
yyaxis right
hold all;
plot(time, phase_t_sim,'Color',[0 1 1],'LineWidth', 2)
plot(time, phase_t_ret,':','Color',[0 0 1],'LineWidth', 2)
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


