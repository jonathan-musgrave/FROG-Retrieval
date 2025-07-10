function [A2] = bin_trace(A, t1, f1, t2, f2)

% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019


[T1, F1] = meshgrid(t1,f1);
[T2, F2] = meshgrid(t2,f2);

A2 = interp2(T1, F1, A, T2, F2, 'linear');

A2(~isfinite(A2)) = 1e-50;

end
   



