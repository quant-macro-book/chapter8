function k0 = newton(fname,k0,kL,kH,varargin)

[FL DFL] = feval(fname,kL,varargin{:});
[FH DFH] = feval(fname,kH,varargin{:});
        
% testvalues = [kL:0.01:kH]';
% [f0 df0 d2f0] = vfuncsp(testvalues,cfz,yterm,Kgrid);
% m = figure;
% subplot(211);
% plot(testvalues,df0);
% subplot(212);
% plot(testvalues,d2f0);
% pause;
% close(m); 

if (DFL<0)
    k0 = kL;
elseif (DFH>0)
    k0 = kH;
else

    klow = kL;
    khigh = kH;

    diff = 1e+4;
    crit = 1e-5;
    s1 = 0;
    bt = 'N';
     
    % Newton-Rhapson
    while (diff>crit)
        
        s1 = s1 + 1;
        [f0 df0 d2f0] = feval(fname,k0,varargin{:});
        NewtonStep = df0/d2f0;
        k1 = k0 - NewtonStep;
        
        if df0>0
            klow = k0;
        else
            khigh = k0;
        end
        
        % backtrack
        if k1>khigh
            k1 = (klow+khigh)/2;
            bt = 'H';
        elseif k1<klow
            k1 = (klow+khigh)/2;
            bt = 'L';
        end
        
%         % apply bisection if Newton Step doe not reduce FOC
%         [f1 df1] = vfuncsp(k1,cfz,yterm,Kgrid);
%         if (abs(df1-df0)>2*abs(df0))
%             k1 = (klow+khigh)/2;
%         end

        diff = abs(k1-k0);
        k0 = k1;

    end

%     s = sprintf( ' kprime = %6.4f  Euler: %4.2f  BT: %s  Iter: %2d', ...
%       k0, df0, bt, s1);

end

% disp(s)