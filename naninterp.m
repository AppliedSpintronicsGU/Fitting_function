function Xi = naninterp(X)
% Interpolate over NaNs
% See INTERP1 for more info
dd_x=find(~isnan(X));
dd_y=X(~isnan(X));
Xi=interp1(dd_x,dd_y,1:length(X))';
return