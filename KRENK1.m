
function [phi1,phim1] = KRENK1(phi,cas)
%KRENK interpolation.
%   [phi1,phim1] = KRENK1(phi,cas) gives the KRENK interpolated values
%   'phi1' and 'phim1' at the points 1 and -1 for the given function
%   'phi', obtained as a solution of an SIE.
%
%   The variable 'cas' denotes the case for the weight function behaviour
%   w(s)=(1-s)^alpha*(1+s)^beta:
%
%        Case | Behaviour at -1 | Behaviour at +1
%        ----------------------------------------
%          I  |     Singular    |     Singular
%         II  |     Singular    |     Bounded
%        III  |     Bounded     |     Singular
%         IV  |     Bounded     |     Bounded
%       
%
%   University of Oxford 
%   Department of Engineering Science
%   Jhonatan Da Ponte Lopes, PhD 
%   May, 2017; Last revision: 2017-05-30

%-------------------------------------------------------------------
%                           CASE DEFINITION
%-------------------------------------------------------------------

N=length(phi); % Number of quadrature points

auxi=1:N;
switch cas
    case 1
        aux1=sin((2.*auxi-1)./(4*N)*pi*(2*N-1))./sin((2.*auxi-1)./(4*N)*pi);
        auxm1=aux1;
        W1=1/N;
        Wm1=W1;
        
    case 2
        aux1=sin((auxi*pi)./(2*N+1)*(2*N-1))./sin((auxi*pi)./(2*N+1));
        auxm1=cot((2*auxi-1)./(2*N+1).*pi/2).*sin((2*auxi-1)/(2*N+1)*N*pi);
        W1=1;
        Wm1=2./(2*N+1);
        
    case 3
        aux1=cot((2.*auxi-1)./(2.*N+1).*pi/2).*sin((2.*auxi-1)/(2.*N+1).*N.*pi);
        auxm1=sin((auxi.*pi)./(2.*N+1).*(2.*N-1))./(sin((auxi.*pi)./(2.*N+1)));
        W1=2./(2.*N+1);
        Wm1=1;
        
    case 4
        aux1=cot(auxi./(N+1).*pi/2).*sin(auxi./(N+1)*N*pi);
        auxm1=aux1;
        W1=1;
        Wm1=W1;
end
clear auxi

%-------------------------------------------------------------------
%                           INTERPOLATION
%-------------------------------------------------------------------

sum1=0;
summ1=0;
for i=1:N
  sum1=sum1+aux1(i)*phi(i);
  summ1=summ1+auxm1(i)*phi(N+1-i);
end

phi1=W1.*sum1;
phim1=Wm1.*summ1;

end