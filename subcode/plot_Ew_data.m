function [] = plot_Ew(Ew, freq)

% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019


Ew = Ew(:); 

p = 4;

S = conj(Ew).* Ew;      S = quickscale(S);

phase_w = -unwrap(angle(Ew));

phase_w(S<max(S)*5e-3) = NaN;   % phase blanking

phase_shift = (max(phase_w)+min(phase_w))/2;

phase_w = phase_w -phase_shift ;


yyaxis left
ax = gca;

plot(freq, S,'Color',[0 0.5 0],'LineWidth', 2)

ylabel('Intensity (a.u.)','Color',[0 0.5 0]);
xlabel('Frequency (THz)'); 
ax.YColor = [0 0.5 0];
axis tight;
ax.XLim = [min(freq) max(freq)];
ax.YTick = [0,1];
ax.YLim = [0,1];

yyaxis right
plot(freq, phase_w,'Color',[0.6 0 0.8],'LineWidth', 2)
ax = gca;
ylabel({'Phase (\pi rad)'},'Color',[0.6 0 0.8]); 

axis tight
ax.YColor = [0.6 0 0.8];
ax.LineWidth = 1.5;
ax.FontSize = 14; %font size
ax.YLim = [-p*pi p*pi];
ax.YTick = -p*pi:4*pi:p*pi;
ax.YTickLabel=[-p:4:p];


