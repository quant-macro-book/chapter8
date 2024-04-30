function [pnew Thetanew Kvecnew Yagg Iagg Cagg Nagg wnew] = bisectp(pL,pH,cmat,knotsk,knotsm,mp,znow,Theta,Kvec)

global critbp;

% bisection
diff = 1e+4;
iter = 0;

while (diff>critbp)

    p0 = (pL+pH)/2;
    [pnew Thetanew Kvecnew Yagg Iagg Cagg Nagg wnew] = pricemapmy(p0,cmat,knotsk,knotsm,mp,znow,Theta,Kvec);
    B0 = p0-pnew; % g(w) = w-f(w)

    if B0<0; % B0*BL>0
        pL = p0;
%         BL = B0;
    else
        pH = p0;
    end

    diff = pH-pL;
    iter = iter + 1;
%         s = sprintf('  bisection %2d,  pH-pL = %6.10f',iter,diff);
%         disp(s);

end