function [X] = bisection(f, Xhigh, Xlow, tol, ItMax)
%BISECTION   Solve f(X) = 0 using bisection
%
%   [X] = bisection(f, Xhigh, Xlow, tol, ItMax) Solves f(X) = 0 using the
%   bisection method. The search is computed using the initial range
%   Xlow < X < Xhigh and is completed when abs(f(X)) < tol. If f(Xlow)and 
%   f(Xhigh) have the same sign then the initial boundaries are 
%   iteratively moved, maintaining the same range, to try and find
%   boundaries which span f(X) = 0.
%
%   BISECTION returns a warning if the maximum number of iterations is
%   reaches and a satisfactory solution of X has not been found.

% flip bounds if required
if f(Xhigh) < f(Xlow)
    tempHigh = Xhigh;
    Xhigh = Xlow;
    Xlow = tempHigh;
end

Iteration = 0;

% check range
Range = Xhigh - Xlow;
while f(Xhigh) < 0  && Iteration < ItMax;
    Xhigh = Xhigh + Range;
    Xlow = Xlow + Range;
    Iteration = Iteration + 1;
end
while f(Xlow)>0  && Iteration < ItMax;
    Xlow = Xlow - Range;
    Xhigh = Xhigh - Range;
    Iteration = Iteration + 1;
end

% find solution
X = (Xhigh+Xlow)/2;
err = f(X);
while abs(err) > tol && Iteration < ItMax;
    Iteration = Iteration + 1;
    if err > 0
        Xhigh = X;
    else
        Xlow = X;
    end
    X = (Xhigh + Xlow) / 2;
    err = f(X);
end

% check solution
if abs(err) > tol
    warning(['max number of iterations reached in routine bisection. ', ...
             'Error (%g) exceeds tolerance (%g).'], err, tol)
end

end

