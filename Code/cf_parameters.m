% Caratheodory-Fejer method for constructing a uniform rational approximation 
% of the function f of type (n,n) on the negative real axis. 
% This is a modified version of a code published in 
% Trefethen, Weideman, Schmelzer
% Talbot Quadratures and Rational Approximations
% BIT, 2006
% zi      -- poles of the rational function
% ci      -- residues associated with those poles
% r_inf   -- r at infinity
% decay   -- SVD of the Hankel matrix associated with the problem
% s       -- error introduced by the (n,n) approximation.
function [zi,ci,r_inf,s,decay] = cf_parameters(fct,n,varargin)
  K = 75;                                 % no of Cheb coeffs
  nf = 1024;                              % no of pts for FFT
  w = exp(2i*pi*(0:nf-1)/nf);             % roots of unity
  t = real(w);                            % Cheb pts (twice over)
  scl = 9;                                % scale factor for stability
  F = zeros(size(t));
  g = (t~=-1);
  F(g) = feval(fct,scl*(t(g)-1)./(t(g)+1),varargin{:});
  %phiS(scl*(t(g)-1)./(t(g)+1),l);  % evaluate Phi_l on neg. axis
  c = real(fft(real(F)))/nf;              % Cheb coeffs of F
  f = polyval(c(K+1:-1:1),w);             % analytic part f of F
  [U,S,V] = svd(hankel(c(2:K+1)));        % SVD of Hankel matrix
  s = S(n+1,n+1);                         % singular value
  decay = diag(S);                        % decay of the error
  u = U(K:-1:1,n+1)'; v = V(:,n+1)';      % singular vectors
  zz = zeros(1,nf-K);                     % zeros for padding
  b = fft([u zz])./fft([v zz]);           % finite Blaschke product
  rt = f-s*w.^K.*b;                       % extended function r-tilde
  zr = roots(v); qj = zr(abs(zr)>1);      % poles
  qc = poly(qj);                          % coeffs of denominator
  pt = rt.*polyval(qc,w);                 % numerator
  ptc = real(fft(pt)/nf);                 % coeffs of numerator
  ptc = ptc(n+1:-1:1); ci = 0*qj;
  for k = 1:n                             % calculate residues
    q = qj(k); q2 = poly(qj(qj~=q));
    ci(k) = polyval(ptc,q)/polyval(q2,q);
  end
  zi = scl*(qj-1).^2./(qj+1).^2;          % poles in z-plane
  ci = 4*ci.*zi./(qj.^2-1);               % residues in z-plane
  r_inf = 0.5*(feval(fct,0,varargin{:})+sum(ci./zi));    % r at infinity
  [m, order]=sort(imag(zi));              % sort poles
  zi=zi(order);                           % reorder poles
  ci=ci(order);                           % reorder residues
end
