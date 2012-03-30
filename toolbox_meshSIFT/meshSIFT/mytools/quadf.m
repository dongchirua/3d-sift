function f = quadf(c)
% this function returns the equation of the implicit quadratic function as
% a string
% the equation: z = ax^2+by^2+cxy+dx+ey+f with
% c = [a b c d e f]
% Author: Chris Maes
% 2009/09

f = sprintf('%d *x^2+ %d *y^2+ %d *x*y+ %d *x+ %d *y+ %d - z',c);
end
