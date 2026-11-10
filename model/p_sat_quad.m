function [y, dydx] = p_sat_quad(x, varargin)
%P_SAT_QUAD  Quadratic saturation: linear on an inner range, parabolic tails into [lo, hi].
%
%   y = p_sat_quad(x)                         slope 1, eps = 1, limits [0 100]
%   y = p_sat_quad(x, Name, Value, ...)       set the options listed below
%   [y, dydx] = p_sat_quad(...)               also return the slope s'(x)
%
%   This is the 'quad' kind of P_SAT on its own. With knees a = lo + eps
%   and b = hi - eps,
%       s(x) = x                        a <= x <= b
%       s(x) = b + eps * G((x - b)/eps) x > b
%       s(x) = a - eps * G((a - x)/eps) x < a
%   where
%       G(u) = u - u^2/4   for 0 <= u < 2,   G(u) = 1 for u >= 2
%       g(u) = G'(u) = 1 - u/2   for 0 <= u <= 2,   g = 0 beyond
%   s is C1, reaches lo and hi exactly at u = 2 (x = lo - eps and
%   x = hi + eps when k = 1), and is flat beyond.
%
%   The slope gain k gives s_k(x) = s(c + k*(x - c)), with c = (lo + hi)/2.
%   s_k is linear with slope k on the inner range, keeps the same tails and
%   limits, and satisfies 0 <= s_k' <= k. k = 1 is the identity on [a, b].
%
%   Options (Name, Value):
%     'eps'   tail width / headroom, scalar or [epsLo epsHi]   (default 1)
%     'k'     slope of the linear part, k > 0                  (default 1)
%     'lims'  output limits [lo hi]                            (default [0 100])
%
%   Example:
%     x = linspace(-30, 130, 1601);
%     [y, d] = p_sat_quad(x, 'eps', 10, 'k', 1.5);
%     plot(x, y, x, d)
%
%   See also P_SAT, P_SAT_SMOOTH.

opts = struct('eps', 1, 'k', 1, 'lims', [0 100]);
if mod(numel(varargin), 2) ~= 0
    error('p_sat_quad:args', 'Options must be given as Name, Value pairs.');
end
for i = 1:2:numel(varargin)
    name = lower(varargin{i});
    if ~isfield(opts, name)
        error('p_sat_quad:args', 'Unknown option ''%s''.', varargin{i});
    end
    opts.(name) = varargin{i+1};
end

epsv = opts.eps;
if isscalar(epsv), epsv = [epsv epsv]; end
k = opts.k;  lims = opts.lims;

assert(numel(lims) == 2 && lims(1) < lims(2), 'p_sat_quad:lims', 'lims must be [lo hi] with lo < hi.');
assert(numel(epsv) == 2 && all(epsv > 0), 'p_sat_quad:eps', 'eps must be positive (scalar or [epsLo epsHi]).');
assert(sum(epsv) <= lims(2) - lims(1), 'p_sat_quad:eps', 'eps is too large: need epsLo + epsHi <= hi - lo.');
assert(isscalar(k) && k > 0, 'p_sat_quad:k', 'Slope k must be a positive scalar.');

lo = lims(1);  hi = lims(2);
c = (lo + hi)/2;              % fixed point of the gain map
a = lo + epsv(1);             % lower knee
b = hi - epsv(2);             % upper knee

z = c + k * (double(x) - c);  % apply the slope gain
y = z;
g = ones(size(z));
up = z > b;  dn = z < a;

[G, g(up)] = tail_profile((z(up) - b) / epsv(2));
y(up) = b + epsv(2) * G;

[G, g(dn)] = tail_profile((a - z(dn)) / epsv(1));
y(dn) = a - epsv(1) * G;

dydx = k * g;
end


function [G, g] = tail_profile(u)
% Normalised tail profile G(u) and its slope g(u) = G'(u), for u >= 0.
uc = min(u, 2);               % G = 1 and g = 0 for u >= 2, including Inf
G = uc - uc.^2/4;
g = 1 - uc/2;
end
