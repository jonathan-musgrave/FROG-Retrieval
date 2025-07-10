function [freq,time,Isig] = load_FROGdata(folderpath,options);
cdOrg = cd;
cd(folderpath);
c =  299792458;      % speed of light
if ~isfield(options,'wavCROP'); options.wavCROP = [500,600]*1e-9;end
if ~isfield(options,'TimeThreshold'); options.TimeThreshold = 0.01;end
if ~isfield(options,'FROG_Threshold'); options.FROG_Threshold = 0.01;end;
if ~isfield(options,'FROG_Threshold'); options.FROG_Threshold = 0.01;end;
if ~isfield(options,'return_fc'); options.return_fc = 'False';end;
if ~isfield(options,'nPad'); options.nPad =0; end
if ~isfield(options,'n'); options.n = 128;end % Number of frequency and delay sample points that  data is interpolated on
%%%%%%%%%%%%%%%%%%%%%% Frog spectrogram raw data
wavelength = load('SHGWavelength_NIR.lvm');
wavelength = wavelength *1E-9;
wavelength = wavelength(1,:);
wavCROP = sort(options.wavCROP);
[~,inx1] = min(abs(wavCROP(1)-wavelength));
[~,inx2] = min(abs(wavCROP(2)-wavelength));
wavCROP = inx1:inx2;
wavelength = wavelength(wavCROP)';
Freq_NL = c./wavelength;
% Load AC

AC_taken = isfile('AC.txt')
if AC_taken
    load('AC.txt')
    tau = AC(:,1);
    % Get rid of the background (approximate)
    DC = mean(AC(1:100,2));
    AC = AC-DC;
    AC = AC(:,2)./max(AC(:,2));

    % Get the inx for croppping time based on threshold values
    INX = AC>options.TimeThreshold;
        inx1 = find(INX == 1,1);
        inx2 = find((INX)==1,1,'last');
    
    % Fit AC to gaussian:
%         fitfunc = fittype('a.*exp(-((x-b)/c).^2)+d');
%         [fitted_curve,gof] = fit(tau,AC,fitfunc,'StartPoint',[.9,tau(floor(length(tau)/2)),max(tau)/8,0]);
%         gauss = @(a,b,c,d,tau) a.*exp(-((tau-b)./c).^2)+d;
%         
%         [yfit,gof] = fit(tau',AC','gauss1')
%         figure(1);clf;
%         plot(yfit,tau,AC)
    % Set Max delay point (MaxD) as COM of AC
        [~,MaxD] = max(AC);
        FWHM = fwhm(AC,tau);
        n = 60;
        [~,inx2] = min(abs(0.5.*n.*FWHM-tau));
        [~,inx1] = min(abs(0.5.*n.*FWHM+tau));
        % [~,MaxD] = min(abs(tau-sum(AC.*tau)./(sum(AC))));
    
    % Recenter data for max(AC) to be at nt/2+1
    delayTau1 = max(abs(MaxD-inx1),abs(MaxD-inx2));
        if delayTau1>length(tau)/2
            delayTau1 = min(abs(MaxD-inx1),abs(MaxD-inx2));
        end
    timCROP = MaxD-delayTau1:MaxD+delayTau1-1;
    tau = tau(timCROP).*1e-15 - tau(MaxD).*1e-15;
    dtau = mean(diff(tau));
    AC = AC(timCROP);
    
    figure(1);clf;
    plot(tau,AC);
    drawnow;
    spectrogram = load('FROGSpectrum_NIR.lvm');
else
end

%% Process the Spectrogram
dark = mean(mean(spectrogram));%load('spectrometer noise stepsize50 steprate500 movement10 integration50ms.lvm');
spect_ORG = spectrogram;
    spectrogram = spect_ORG(:,:);
%dark = ones(a,1)*dark(1,:);
spectrogram = (spectrogram - dark)';
spectrogram = spectrogram(wavCROP, timCROP);
% Crop tau
spectrogram = spectrogram.*(spectrogram>(options.FROG_Threshold*max(max(spectrogram))));
spectrogram = spectrogram/max(max(spectrogram));
spectrogram = padarray(spectrogram',options.nPad,0,'both')';
spectrogram = padarray(spectrogram,options.nPad,0,'both');


% Recenter data on the carrier frequency
for ind = 1:length(tau)
    Fspe(:,ind) = spectrogram(:,ind).*(wavelength.^2);  % f(lamda)dlambda = f(freq)dfreq ; dlamda=c*lamda^2*dfreq
end
Fspe = Fspe/max(max(Fspe));
freqSH = (c./wavelength);
Cspe = sum(Fspe,2);
freqSH_center = sum(Cspe.*freqSH)/sum(Cspe);
freqSH = freqSH - freqSH_center;

% NOT FOURIER AXES OF EACH OTHER
freq = linspace(c./options.wavCROP(2),c./options.wavCROP(1),options.n);
time = linspace(min(tau),max(tau),options.n);
freq = freq-freqSH_center;

% IS FOURIER AXES OF EACH OTHER
nt = options.n;
tWin = abs(min(tau));
dt = 2*tWin./nt;
dw = pi/tWin;
df = dw./(2*pi);
N = (-nt/2:nt/2-1);
time = dt.*N;
freq = df.*N;


[T,F] = meshgrid(time,freq);
FROG = interp2(tau,freqSH,spectrogram,T,F);
FROG(isnan(FROG)) = 0;
Isig = FROG./(max(max(FROG)));
cd(cdOrg);

if strcmpi(options.return_fc,'True')
    freq_temp = freq;
    freq = cell(2,1);
    freq{1} = freq_temp;
    freq{2} = freqSH_center;
end


end