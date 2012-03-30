function pp = quadinfl(coef)
% c are the coefficients of 
% the equation: z = ax^2+by^2+cxy+dx+ey+f with
% c = [a b c d e f]
x = (-2*coef(2)*coef(4)-coef(5)*coef(3))/(coef(3)^2+4*coef(1)*coef(2));
y = (2*coef(1)*x+coef(4))/coef(3);
z = coef(1)*x^2+coef(2)*y^2+coef(3)*x*y+coef(4)*x+coef(5)*y+coef(6);
pp = [x y z];
end