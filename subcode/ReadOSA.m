% Read OSA code for  data processing

% note the axes or in Hz not in 2*pi*Hz
function [sA0,w0,lam] = ReadOSA(fname, fc, Ntrace, freq_center, threshold)
if ~exist('threshold','var'), threshold = 5e-5;end % Default threshold 
if ~exist('fc','var'), JustPlot = 1;else, JustPlot = 0;end % Withoutt Interpolation axes it just plots eve
if ~exist('freq_center','var'), CalcFreqC = 1;
else CalcFreqC = 0; end % Calculate the center frequenct of EACH trace. (interpolation will be different
if ~exist('Ntrace','var'), Ntrace = inf; end;
c =  299792458;      % speed of light
%%%%%%%%%%%%% spectrum raw data
Data = readmatrix(fname,'NumHeaderLines',319); % Get all the data (OSA 1.5-3um)
% Data = readmatrix(fname,'NumHeaderLines',278); % Get all the data (OSA 0.6-2um)
NSweeps = length(Data(~isnan(Data(1,:))))/2;
TraceName = {'A','B','C','D','E','F','G'};
sA0 = cell(length(TraceName),1);
lam = cell(length(TraceName),1);
w0 = cell(length(TraceName),1);
% w0 = zeros(NSweeps,1);
if JustPlot == 1
    figure();
    for i = 1:min(NSweeps,Ntrace); % Interpolate data onto the centered frequency value
        frr = Data(:,2*i);
        SpectRaw = frr;
        frr = 10.^(frr/10);
        lambda = Data(:,2*i-1);
        
        % freq = c*1e9./lambda; % NL frequency axis
        % Get the spectrum
        % frr = frr.*(lambda.^2);
        frr = ((frr)./max(frr));
        frr = frr.*(frr>threshold); % Thresholding
        frr = ((frr)./max(frr));
        %frr(isnan(frr)) = -100;
        lambda0 = sum(frr(~isnan(frr)).*lambda(~isnan(lambda)))/sum(frr(~isnan(frr)));
       
        subplot(1,2,1);
        plot(lambda,SpectRaw,'DisplayName',num2str(round(lambda0)));hold on;
        YLmax(i) = max(SpectRaw);
        minLam(i) = min(lambda);
        maxLam(i) = max(lambda);
        subplot(1,2,2); 
        plot(lambda,frr,'DisplayName',TraceName{i});hold on;
        sA0{i} = frr(~isnan(frr));
        lam{i} = lambda(~isnan(frr));
        w0{i} = 2*pi*c./(lambda0);
        
    end

    subplot(1,2,1); grid on;
    xlabel('Wavelength (nm)')
    ylabel('dBm')
    ylim([-90,max(YLmax)+10]);
    xlim([min(minLam),max(maxLam)])
    l =legend;
    l.Title.String = '\lambda_c'
    
    subplot(1,2,2); grid on;
    xlabel('Wavelength (nm)')
    ylabel('Linear units')
    title('Normalized Linear Spectrum')
    xlim([min(minLam),max(maxLam)])
    l =legend;
    ylim([0,1.2])
    l.Title.String = 'Trace';
    
else
    figure();
    sA0 = zeros(NSweeps,length(fc));
    for i = 1:min(NSweeps,Ntrace); % Interpolate data onto the centered frequency value
        frr = Data(:,2*i);
        if max((~isnan(frr))) == 0
        else
            frr = frr(~isnan(frr));
            frr = 10.^(frr/10);
            lambda = Data(:,2*i-1);
            lambda = lambda(~isnan(lambda));
            freq = c*1e9./lambda; % NL frequency axis
            % Get the spectrum
            frr = frr.*(lambda.^2);
            frr = ((frr)./max(frr));
            frr = frr.*(frr>threshold); % Thresholding
            frr = ((frr)./max(frr));

            if CalcFreqC

                freq_center = sum(frr.*freq)/sum(frr);
            end
            freq = freq - freq_center;

            sA0(i,:) = interp1(freq,frr,fc);
            sA0(i,isnan(sA0(i,:))) = 0;
            w0(i) = 2*pi*freq_center;

            subplot(1,2,1);
            plot(lambda,10*log10(frr));hold on;
            subplot(1,2,2); 
            plot(2*pi*fc./(1e12),sA0(i,:),'DisplayName',TraceName{i});hold on;
        end
        i
    end
        subplot(1,2,1)
        title('Normalized Raw Data')
        subplot(1,2,2)
        title('Interpolated Data on \omega (2\piTHz)')
    end

end