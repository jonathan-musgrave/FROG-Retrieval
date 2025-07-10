function plot_ComparisonExperimentFROGs_Input_Output(Frog_data_output,Frog_data_sim,freq,time,uuLz,uuLz_rev,varargin)
for i = 1:length(varargin)
    switch varargin
        case strcmpi(varargin,'fignum')
            fignum = varargin{i+1}
    end
end
if~exist('fignum','var');fignum = 3;end
fs = 1e15;
THz = 1e-12;

hFig = figure(fignum);clf;
set(hFig, 'Position', [100 700 950 550]) 

subplot(2,2,1);     plot_trace(Frog_data_output{1}, time{1}*fs, (freq{1}{1}+freq{1}{2})*THz);            title('Measured Output Trace'); 

subplot(2,2,2);     plot_trace(Frog_data_output{2}, time{1}*fs, (freq{1}{1}+freq{1}{2})*THz);        title('Retrieved Output trace'); 

subplot(2,2,3);     plot_trace(Frog_data_sim{1}, time{2}*fs, (freq{2}{1}+freq{2}{2})*THz);          title('Simulated Output Trace');

subplot(2,2,4);     plot_trace(Frog_data_sim{2}, time{2}*fs, (freq{2}{1}+freq{2}{2})*THz);       title('Simulated Output Trace (Input rev.)');


hFig = figure(fignum+1);clf;
set(hFig, 'Position', [100 700 950 550]) 

subplot(2,2,1);     plot_Et_exper_vs_simulation(uuLz_rev{1},uuLz{2},time{1}*fs,time{2}*fs);            title('Output Field Exp. rev. Sim forward'); 

subplot(2,2,2);     plot_Et_exper_vs_simulation(uuLz{1},uuLz_rev{2},time{1}*fs,time{2}*fs);       title('Output Field Exp. forward. Sim reversed'); 
legend
subplot(2,2,3);     plot_Et_exper_vs_simulation(uuLz{1},uuLz{2},time{1}*fs,time{2}*fs);       title('Output Field Exp. forward. Sim Forward'); 

subplot(2,2,4);     plot_Et_exper_vs_simulation(uuLz_rev{1},uuLz_rev{2},time{1}*fs,time{2}*fs);   title('Output Field Exp. reversed. Sim reversed'); 

% subplot(2,2,4);     plot_Ef_exper_vs_simulation(uuLz_rev{1},uuLz_rev{2},freq{1}{1},freq{2}{1});      % title('Output Field rev');

end
