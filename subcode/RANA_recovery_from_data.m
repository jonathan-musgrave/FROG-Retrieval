function [Isig, Ishg_ret, et_ret, et_ret_rev, freq,time] = RANA_recovery_from_data(folder,options)

[freq,time,Isig] = load_FROGdata(folder, options);

figure(1);clf;
plot(freq{1},sum(Isig,2))
xlim([-7,7]*1e12)

% adjust G/G' error based on the amount of noise and area of trace, 
% in the presence of additive noise G error is almost independent of the trace dimension:
G_cutoff = options.G_cutoff;     Gp_cutoff = options.Gp_cutoff;

[~,D] = fileparts(folder);
fprintf(strcat('\nRunning Recovery on data from: %s')); fprintf(D);fprintf('\n\n');
[Et_final, Gout, Gpout, Sout] = RANA(abs(Isig), time, G_cutoff, Gp_cutoff);

%% 
et_ret = center(Et_final,'max');

et_ret_rev = center(conj(flipud(et_ret)),'max');

Ishg_ret = shg_trace(et_ret);


end