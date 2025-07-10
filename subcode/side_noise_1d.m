function X = side_noise_1d(X, n)


% Author: Rana Jafari
% Georgia Institute of Technology
% email: rjafari7@gatech.edu
% June 2019


    X = X(:).';

    a = mean([X(1:n) X(end-n+1:end)]);

    X = X - a;

    X(X<0) = 1e-5;

end