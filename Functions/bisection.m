function [X] = bisection(Fun, Xlow, Xhigh, tol, ItMax)
%BISECTION   Solve FUN(X) = 0 using bisection
%
%   [X] = bisection(FUN, XLOW, XHIGH, TOL, ITMAX) Solves FUN(X) = 0 using
%   bisection. The search is computed using the initial range
%   XLOW < X < XHIGH and is completed when abs(FUN(X)) < TOL. If FUN(XLOW) 
%   and FUN(XHIGH) have the same sign then the initial boundaries are 
%   iteratively moved, maintaining the same range, to try and find
%   boundaries which span FUN(X) = 0. 
%
%   ITMAX is an optional input with default 100.
%
%   BISECTION returns a warning if the maximum number of iterations is
%   reached and a satisfactory solution of X has not been found.
%
%   Examples
%     FUN can be specified using @:
%       x = bisection(@myfun, 0, 10, 0.001)
%
%     where myfun is a MATLAB function such as:
%       function F = myfun(x)
%           F = x^2 - 4;
%       end
%
%     FUN can also be an anonymous function:
%       x = bisection(@(x) x^2 - 4, 0, 10, 0.001)
%
%     If FUN is parameterized, you can use anonymous functions to capture 
%     the problem-dependent parameters:
%       function F = myfun(x,c)
%           F = 2*x^2 - x - c;
%       end
%
%       x = bisection(@(x) myfun(x,15), 0, 10, 0.001)
%   
%   See also: FSOLVE.

% set default ItMax if required
if ~exist('ItMax', 'var')
    ItMax = 100;
end

% flip bounds if required
if Fun(Xhigh) < Fun(Xlow)
    tempHigh = Xhigh;
    Xhigh = Xlow;
    Xlow = tempHigh;
end

Iteration = 0;

% check range
Range = Xhigh - Xlow;
while Fun(Xhigh) < 0  && Iteration < ItMax;
    Xhigh = Xhigh + Range;
    Xlow = Xlow + Range;
    Iteration = Iteration + 1;
end
while Fun(Xlow)>0  && Iteration < ItMax;
    Xlow = Xlow - Range;
    Xhigh = Xhigh - Range;
    Iteration = Iteration + 1;
end

% find solution
X = (Xhigh+Xlow)/2;
err = Fun(X);
while abs(err) > tol && Iteration < ItMax;
    Iteration = Iteration + 1;
    if err > 0
        Xhigh = X;
    else
        Xlow = X;
    end
    X = (Xhigh + Xlow) / 2;
    err = Fun(X);
end

% check solution
if abs(err) > tol
    warning(['max number of iterations reached in routine bisection. ', ...
             'Error (%g) exceeds tolerance (%g).'], err, tol)
end

end

