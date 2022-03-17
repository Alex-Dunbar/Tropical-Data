function [val] = trop_polyval(x,coeffs)
% Evaluates the single variable tropical polynomial 
% 
% q(t) = max(t + coeffs(2), 2*t + coeffs(3), ...
%                                 d*t + coeffs(d+1), coeffs(1));
%
% Inputs: x      - vector of evaluation points
%         coeffs - vector of coefficients 
%
% Outputs: val   - vector of size(x) with entries q(x)
%
% Example usage:
%
% val  = trop_polyval([0;1],[0;0]); returns val = [0;1];


d = numel(coeffs) - 1; N = numel(x);

%evaluate without forming Vandermonde
val = coeffs(1)*ones(N,1);
for t = 1:d
    val = max(val,t*x + coeffs(t+1));
end


end