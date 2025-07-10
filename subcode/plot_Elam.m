function [] = plot_Elam(Elam, lam)

% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019


Elam = Elam(:); 

p = 4;

S = conj(Elam).* Elam;      S = quickscale(S);

phase_lam = -unwrap(angle(Elam));

phase_lam(S<max(S)*5e-3) = NaN;   % phase blanking

phase_shift = (max(phase_lam)+min(phase_lam))/2;

phase_lam = phase_lam -phase_shift ;


yyaxis left
ax = gca;

plot(lam, S,'Color',[0 0.5 0],'LineWidth', 2)

ylabel('Intensity (a.u.)','Color',[0 0.5 0]);
xlabel('Wavelength (nm)'); 
ax.YColor = [0 0.5 0];
axis tight;
ax.XLim = [min(lam) max(lam)];
ax.YTick = [0,1];
ax.YLim = [0,1];

yyaxis right
plot(lam, phase_lam,'Color',[0.6 0 0.8],'LineWidth', 2)
ax = gca;
ylabel({'Phase (\pi rad)'},'Color',[0.6 0 0.8]); 

axis tight
ax.YColor = [0.6 0 0.8];
ax.LineWidth = 1.5;
ax.FontSize = 14; %font size
ax.YLim = [-p*pi p*pi];
ax.YTick = -p*pi:4*pi:p*pi;
ax.YTickLabel=[-p:4:p];


