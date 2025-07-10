function [u_3d_2, v_3d_2, uu_3d_2, spect_3d_2, z, LD,LNL,N2]=Run_SSFM_from_data(exper_notes, uu0_rec, tau, omega, T0cell, options)

% Calculate Experimental Values
Exper = readmatrix(exper_notes);

T0 = T0cell{1}./(2*log((1+sqrt(2))));
T0 = T0cell{1};
Input_Power = Exper(1,2);
Input_repRate = Exper(2,2); % in MHz
Lz = Exper(3,2); % Fiber length under test in m
E0 = Input_Power./Input_repRate; % Pulse Energy in nJ
P0 = E0./T0; % Peak Power in MW
% P0 = E0./(fwhm(abs(uu0_rec{1}).^2,tau.*T0));

I0 = abs(uu0_rec{1}).^2./sum(abs(uu0_rec{1}).^2.*mean(diff(tau.*T0))).*E0;

% P0 = max(I0);
% P0 = E0./T0;
s = options.lam0/(2*pi*3e2*T0)*options.self_steepening;  % self-steepening

LD = T0.^2/(abs(options.beta2)./((1e-15).^2)); % in m
LNL = 1./(options.gamma.*P0*1e6); % in m
distance = Lz/LD;   % normalized to dispersion length
N = sqrt(LD./LNL);          % soliton order
N2 = N.^2;
N = 2;
distance = pi/2;
% N = 1;
sbeta2 = options.sbeta;    % sign of beta_2 (positive is normal negative is anomolous)
delta3= 0.00;   % beta3/16*T0*abs(beta2))
delta4 = 0.0;   % beta4/(24*T0^2*abs(beta2))

%---set silumation parameters
dtau = mean(diff(tau));
nt = length(tau);
Tmax = dtau*nt/2;       % FFT points and window size
% step_num = 99500*2*0 + 100000;   %round(40*distance*N^2); 
zstep = options.nz;      % no of z steps
step_num = zstep*100;
ng = step_num/zstep;        % store data for 3D.
deltaz = distance/step_num; % step size in z
z = 0:deltaz*ng:distance;


if options.RamOn
    fR = 0.245;  fb=0.21;
    tau1=12.2/T0; tau2=32/T0; tau_b=96/T0;
    h = (1-fb)*(tau1^2 + tau2^2)/(tau1*tau2^2)*exp(-tau/tau2).*sin(tau/tau1)...
            +fb*((2*tau_b-tau)/tau_b^2).*exp(-tau/tau_b);
    h(1:nt/2) = 0;      % causality

else
   h = zeros(1,nt); fR=0; 
end
Ram_res = h;

%----store dispersive phase shifts to speedup code
dispersion = exp(1i*deltaz*(1/2*sbeta2*omega.^2 +...
    delta3*omega.^3 + delta4*omega.^4));
hhz = 1i*N^2*deltaz;

for i = 1:2
   UU0_rec{i} = ifft(fftshift(uu0_rec{i})); 
end
for i = 1:2
    uu0 = uu0_rec{i}';
    C = 0;
    uu0 = uu0./sqrt(max(abs(uu0).^2));
    % uu0 = sech(tau).*exp(-1j*C.*tau.^2/2);
    % uu0 = real(uu0);
    % uu0 = sech(tau);
    uu_3d(1,:)=uu0.*conj(uu0); % Intensity in time domain
    u_3d(1,:) = real(uu0);
    v_3d(1,:) = imag(uu0);
    temp = (ifft(fftshift(uu0)));
    spect_3d(1,:)=abs(temp).^2; % spectrum au unit
    
    uu = uu0;
    %*********[ Beginning of MAIN Loop]***********
    % scheme: 1/2N -> D -> 1/2N; first half step nonlinear
    P=uu.*conj(uu);
    %   convl=(nt*dtau)*fft(ifft(ifftshift(h)).*ifft(ifftshift(P))); % Convolution of raman response with Power (Intensity) normalized to windowing 
    %   convl=fftshift(convl); 
    %   sst3=deriv1(convl,dtau);% Take derivative of convolution wrt to time 
    %   sstnew=1i*s*fR*(sst3 + (conj(uu).*convl)./(P+eps).*sst2); % Adjust sst with this value?
    %   temp = uu.*exp(((1-fR)*(abs(uu).^2)+ fR*convl).*hhz/2); % Apply full step on nonlinearity
    temp = uu.*exp(((abs(uu).^2)).*hhz/2); % Apply full step on nonlinearity

    for n=1:step_num
        f_temp = ifft(temp).*dispersion; % 1 full step hhz of dispersion    
        uu = fft(f_temp); % uu is in t domain
        %   P=uu.*conj(uu);
        %   convl=(nt*dtau)*fft(ifft(ifftshift(h)).*ifft(ifftshift(P))); % Convolution of raman response with Power (Intensity) normalized to windowing 
        %     convl=fftshift(convl);
        %   % Convolution terms based on JM math
        %   duu_dt = deriv1(uu,dtau);
        %   duuc_dt = deriv1(conj(uu),dtau);
        %   dconvl_dt = deriv1(convl,dtau);
        %   SST = (2*(1-fR)*1j*s*(duu_dt.*P + 1/2*(duuc_dt).*uu.*uu) + ...
        %       1j*s*fR*(duu_dt.*convl + dconvl_dt.*uu))./(uu+eps);

        %temp = uu.*exp(((1-fR)*(abs(uu).^2)+ fR*convl + SST).*hhz); % Apply full step on nonlinearity
        temp = uu.*exp(((abs(uu).^2)).*hhz); % Apply full step on nonlinearity
        spect=(ifft(fftshift(temp)));
        if mod(n,ng)==0
            u_3d(n/ng+1,:) = real(temp);
            v_3d(n/ng+1,:) = imag(temp);
            uu_3d((n/ng+1),:)=temp.*conj(temp); % temporal power
            spect_3d((n/ng+1),:)=spect.*conj(spect); % spectral power
        end
    end
    uuLz = uu.*exp(((abs(uu).^2)).*-hhz/2); % Apply - step/2 of nonlinearity
    u_3d(end,:) = real(uuLz);
    v_3d(end,:) = imag(uuLz);
    uu_3d(end,:) = abs(uuLz).^2;
    spect_3d(end,:) = abs(ifft(fftshift(uuLz))).^2;
    
    u_3d_2{i} = u_3d;
    v_3d_2{i} = v_3d;
    uu_3d_2{i} = uu_3d;
    spect_3d_2{i} = spect_3d;
end
    
end