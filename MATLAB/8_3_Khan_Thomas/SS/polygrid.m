function x = polygrid(m,a,b)

kvec = [1:m]';
z = -cos(pi/2/m*(2*kvec-1));
x = (z+1).*(b-a)./2 + a.*ones(m,1);