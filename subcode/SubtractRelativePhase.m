function Eout = SubtractRelativePhase(Ein)

phi = unwrap(angle(Ein));
Eout = abs(Ein).*exp(1j*(phi - phi(length(phi)/2+1)));


end