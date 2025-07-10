# FROG-Retrieval
FROG retrieval based on generalized projections and using the RANA algorithm for initial guess

% Author: Rana Jafari
% Georgia Tech
% email: rjafari7@gatech.edu
% Feb 2020
%
% Revised by Jon Musgrave to use any frequency axis
% I_FROG should still be  a [nt,nt] array which is cropped such that the
% ifftc was taken at npad. This code will do the cropping down to the
% original [-nt/2:nt/2-1] size.
% In this case it is important that the user considers the cropped
% frequency window which will be df*nt where df = 1/npad/dt
