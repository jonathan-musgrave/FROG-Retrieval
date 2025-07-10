function [] = plot_trace(Ifrog, time, freq)

% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019

imagesc(time, freq, Ifrog);

ax = gca;

colormap(cmap_frog());
c = colorbar;
axpos = ax.Position;
c.Position(3) = 0.6*c.Position(3);
ax.Position = axpos;

axis tight;
axis square;

fig=gcf;
set(findall(fig,'-property','FontSize'),'FontName', 'Helvetica','FontSize',14);
xlabel('Delay (fs)');
ylabel('Frequency (THz)');
ax.LineWidth = 1.5;

end


