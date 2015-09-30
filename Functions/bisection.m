function [solution] = bisection(f,high,low,tol)
% simple bisection solver
% f = function where we desire to find f(solution) = 0
% high = initial high guess i.e. f(high>0)
% low = initial low guess i.e. f(low<0)
% tol = tolerance
% Can handle solution outside initial range provided initial range is reasonably
% close...

% check range
while f(high)<0
    high = high + (high-low);
end
while f(low)>0
    low = low - (high-low);
end

% find solution
solution = (high+low)/2;
err = f(solution);
while abs(err)>tol;
    if err>0
        high = solution;
    else
        low = solution;
    end
    solution = (high+low)/2;
    err = f(solution);
end

end

