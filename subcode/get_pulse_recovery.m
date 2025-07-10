
%  Author: Jonathan Musgrave

function [FROG_data, uu0_rec, uuLz_rec, tau, omega, T0,time_in,time_out,freq_in,freq_out] = get_pulse_recovery(Input_FROG_folder,Output_FROG_folder, options)

folderDir = fullfile(cd,'FROG_data');

% Do recovery on the input trace
options.return_fc = 'True';
[FROG_data_input, FROG_rec_input, uu0, uu0_rev, freq_in, time_in] =...
    RANA_recovery_from_data(fullfile(folderDir,Input_FROG_folder),options);

% Do recovery on the output trace
[FROG_data_output, FROG_rec_output, uuLz, uuLz_rev, freq_out, time_out] =...
    RANA_recovery_from_data(fullfile(folderDir,Output_FROG_folder),options);

uu0 = SubtractRelativePhase(uu0);


uu0_rev = center(conj(flipud(uu0)),'max');

%         FROG_data_output = FROG_data_input;
%         FROG_rec_output = FROG_rec_input;
%         uuLz = uu0;
%         uuLz_rev = uu0_rev;
%         freq_out = freq_in;
%         time_out = time_in;


plot_FROGs_Input_Output({FROG_data_input,FROG_rec_input}...
    ,{FROG_data_output,FROG_rec_output}...
    ,{freq_in,freq_out}...
    ,{time_in,time_out}...
    ,{uu0,uuLz}...
    ,{uu0_rev,uuLz_rev});

N = options.n;  
dt = mean(diff(time_in));
h = -N/2:N/2-1;
freq_FF = 1/(dt*N)*(h); % Fundamental frequencies
freq = 1/N/dt*h; % SHG Frog 

%% Save relevant Params for SSFM
    FWHM = fwhm(abs(uu0).^2,time_in);
    T0 = FWHM./(2*log(1+sqrt(2))); % T0 to T_FWHM for sech Pulse
UU0 = fftshift(fft(fftshift(uu0)));
uuFTL = fftshift(ifft(fftshift(abs(UU0))));
T0FTL = fwhm(abs(uuFTL).^2,time_in)./(2*sqrt(log(2)));
% T0 = fwhm(abs(uu0).^2,time_in);
tau   = time_in./T0;
omega = fftshift(2*pi*freq_FF.*T0FTL);
T0 = (T0./1e-15);
T0FTL = T0FTL./(1e-15);
T0 = {T0,T0FTL};
P0 = max((abs(uuFTL).^2));

uu0  = uu0./sqrt(P0);
uuLz = uuLz./sqrt(P0);

%% Organize data
FROG_data{1} = FROG_data_input;
FROG_data{2} = FROG_rec_input;
FROG_data{3} = FROG_data_output;
FROG_data{4} = FROG_rec_output;
uuLz_rec{1} = uuLz;
uuLz_rec{2} = uuLz_rev;
uu0_rec{1} = uu0;
uu0_rec{2} = uu0_rev;



end


