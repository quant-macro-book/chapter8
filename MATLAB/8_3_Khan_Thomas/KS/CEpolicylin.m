function [BetaK, Betap] = CEpolicylin(zvec, ivec, kvec, kprime, p, znum, simLength)

% Exogenous Laws of Motion for Capital and prices in the Discretised Baseline Model:
% Quadratic Form Regressions%
% [BetaK, Betap] = CEpolicylin(zvec, ivec, kvec, kprime, p, znum, simLength)
%
% simLength = length(zvec);
% zvec is 1 x simLength history of productivity shock levels
% ivec is the associated index level of productivity
% kvec is associated history of current capital 
% kprime is associated history of next period's stock of capital
% p is associated history of current relative price of output
% znum = length(Zgrid) the number of discrete shock levels
% BetaK is the estimated relation logk' = Betak*[1 logk logk^2 logk^3]
% Betap is the estimated relation logp  = Betap*[1 logk logk^2 logk^3]
% where rows of Betak and Betap are associated with distinct levels of
% current productivity z.

% Julia Thomas, Spring 2000
%disp(' ');
%disp ( ' CEpolicylin results with linear form ' )     

for i = 1:1:znum 
   zloc = find(ivec(1:simLength) == i);
  
   zobs(i) = length(zloc);
   kz = kvec(zloc);
   kzprime = kprime(zloc);
   pz = p(zloc);
  
   klog = log(kz');
   plog = log(pz);
   kprimelog = log(kzprime);
   unit = ones(zobs(i),1);
   time = linspace(1,zobs(i), zobs(i));
      
   % independent and dependent variables for second order multidimensional linear regression in logs
   Y1 = [kprimelog' plog'];   
   X2 = [unit klog];
      
   % Regression: Note that these are multi-dimensional multivariate linear regressions.  More than 
   % one log linear relationship is being estimated.  
   % Beta3(:,:,i) = inv(X3'*X3)*(X3'*Y3);
   Beta2(:,:,i) = (X2'*X2)\(X2'*Y1);
   
   % Estimates
   estimates2 = (X2*Beta2(:,:,i))';
   
   N = length(klog) - 2;
   
   % Residuals: sum of squared residuals, sum of squared deviations and R2 = 1 - ssr/ssd
   % for each dependent variable
   residkt = estimates2(1,:) - kprimelog;					% second order capital
   ssrkt = (residkt*residkt')/N;
   ssdkt = kprimelog - mean(kprimelog); 
   ssdkt = (ssdkt*ssdkt')/N;
   R2kt = 1 - ssrkt/ssdkt;									% This is also the unadjusted R-squared
   ssrkt = sqrt(ssrkt);										% This yields standard error for logk
   
   residpt = estimates2(2,:) - plog;						% second order price
   ssrpt = (residpt*residpt')/N;
   ssdpt = plog - mean(plog); 
   ssdpt = (ssdpt*ssdpt')/N;
   R2pt = 1 - ssrpt/ssdpt;
   ssrpt = sqrt(ssrpt);										% This yields standard error for logp
   
   % minimum, maximum and sum of squared residual errors for residuals of first and second order regressions
   error2(1,:,i) = [min(abs(residkt)) max(abs(residkt)) ssrkt R2kt];
   error2(2,:,i) = [min(abs(residpt)) max(abs(residpt)) ssrpt R2pt];

   % Tables of Results   

   disp ( ' ' );    
   disp ( ' -------------------------------------------------------------------------------- ' )
   disp ([ ' znum=',num2str(i),':  Regression Coefficients based on log(x) = Beta0 + Beta1*log(k) ' ])
   s1 = sprintf ( ' kf    %+10.4f  %+10.4f  %+10.4f  ', Beta2(1,1,i), Beta2(2,1,i));
   s2 = sprintf ( '  p    %+10.4f  %+10.4f  %+10.4f  ', Beta2(1,2,i), Beta2(2,2,i));   
   disp ( ' 	       Beta0 	   Beta1 	   ' )
   disp (s1)
   disp (s2)     
   disp ( ' ' )
      
   disp ([ ' znum=',num2str(i),':  Regression Statistics based on ', num2str(zobs(i)), ' observations '])
   s1 = sprintf( ' kf     %8.4e   %8.4e   %8.4e    %12.10f ',error2(1,1,i), error2(1,2,i), error2(1,3,i), error2(1,4,i));
   s2 = sprintf( '  p     %8.4e   %8.4e   %8.4e    %12.10f ',error2(2,1,i), error2(2,2,i), error2(2,3,i), error2(2,4,i));
      
   disp ( '         minimum        maximum        S. E.  	     R-SQUARED ' )
   disp (s1)
   disp (s2)
   disp ( ' -------------------------------------------------------------------------------- ' )
    
   % separate and reorganize the laws of motion   
   BetaK(i,:) = Beta2(:,1,i)';
   Betap(i,:) = Beta2(:,2,i)';

end
