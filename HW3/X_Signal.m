function [x,v] = X_Signal(N,b,q)

% Construct the real discrete signal x[k], k = 1, ..., N, derived
% as the output of a MA-q process with coefficients b, driven by
% white non-Gaussian noise v[k], from exponential distribution
% with mean value of 1.

v = exprnd(1,[1,N]);  % input: White non-Gaussian noise
v = v - mean(v);      % remove mean, better NRMSE

% Output of MA-q process
x = zeros(1,N);
 for k=1:N
   for j=0:q
     if k>j
        x(k) = x(k) + b(j+1)*v(k-j);
     end
   end
 end
end