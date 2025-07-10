function plot_FROGs_Input_Output(Frog_data_input,Frog_data_output,freq,time,uu,uu_rev,varargin)
for i = 1:length(varargin)
    switch varargin
        case strcmpi(varargin,'fignum')
            fignum = varargin{i+1}
    end
end
if~exist('fignum','var');fignum = 1;end
fs = 1e15;
THz = 1e-12;
hFig = figure(fignum);clf;
set(hFig, 'Position', [100 100 950 550]) 

subplot(2,2,1);     plot_trace(Frog_data_input{1}, time{1}*fs, (freq{1}{1}+freq{1}{2})*THz);            title('Measured Input Trace'); 

subplot(2,2,2);     plot_trace(Frog_data_input{2}, time{1}*fs, (freq{1}{1}+freq{1}{2})*THz);        title('Retrieved Input trace'); 

subplot(2,2,3);     plot_Et_exper(uu{1}, uu_rev{1},time{1}*fs);           title('Retrieved Input Field I(t)');

subplot(2,2,4);     plot_Ef_exper(uu{1}, uu_rev{1} ,freq{1}{1}*THz);       title('Retrieved Input Field I(f)');

hFig = figure(fignum+1);
set(hFig, 'Position', [100+1000 100 950 550]) 

subplot(2,2,1);     plot_trace(Frog_data_output{1}, time{2}*fs, (freq{2}{1}+freq{2}{2})*THz);            title('Measured Output Trace'); 

subplot(2,2,2);     plot_trace(Frog_data_output{2}, time{2}*fs, (freq{2}{1}+freq{2}{2})*THz);        title('Retrieved Output trace'); 

subplot(2,2,3);     plot_Et_exper(uu{2}, uu_rev{2},time{2}*fs);           title('Retrieved Input Field I(t)');

subplot(2,2,4);     plot_Ef_exper(uu{2}, uu_rev{2} ,freq{2}{1}*THz);       title('Retrieved Output Field I(f)');
end
