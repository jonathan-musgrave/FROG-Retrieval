function [] = plot_Et(Et, time)

% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019


Et = Et(:); 

p = 4;

N =length(Et);

I = conj(Et).* Et;      I = quickscale(I);

phase_t = -unwrap(angle(Et));

phase_t(I<max(I)*5e-3) = NaN;	% phase blanking

phase_shift = (max(phase_t)+min(phase_t))/2;

phase_t = phase_t - phase_shift;


yyaxis left
ax = gca;

plot(time, I,'Color',[1,0,0],'LineWidth', 2)
fprintf(['FWHM: ',num2str(round(fwhm(I,time),0)),' fs\n'])
ylabel('Intensity (a.u.)','Color',[1,0,0]);
xlabel('Time (fs)'); 
ax.YColor = [1,0,0];
axis tight;
ax.XLim = [min(time) max(time)];
ax.YTick = [0,1];
ax.YLim = [0,1];

yyaxis right
plot(time, phase_t,'Color',[0 0 1],'LineWidth', 2)
ax = gca;
ylabel({'Phase (\pi rad)'},'Color',[0,0,1]); 

axis tight
ax.YColor = [0,0,1];
ax.LineWidth = 1.5;
ax.FontSize = 14; %font size
ax.YLim = [-p*pi p*pi];
ax.YTick = -p*pi:4*pi:p*pi;
ax.YTickLabel=[-p:4:p];


