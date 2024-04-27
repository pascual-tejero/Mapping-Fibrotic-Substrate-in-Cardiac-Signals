function RE=RelativeError(phi_meas,phi_calc)
% function RE=RelativeError(phi_meas,phi_calc)
% Eq. 34
nom=sqrt(sum((phi_meas-phi_calc).^2));
den=sqrt(sum(phi_meas.^2));
RE=nom/den;