
function phik = KRENK2(phi,sint,cas)
%KRENK interpolation.
%   phik = KRENK1(phi,sint,cas) gives the KRENK interpolated value
%   'phik' at the desired point 'sint' for the given function
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
%   June, 2017; Last revision: 2017-06-30

%-------------------------------------------------------------------
%                           CASE DEFINITION
%-------------------------------------------------------------------

N=length(phi); % Number of quadrature points

auxi=1:N;
switch cas
    case 1
        gam=acos(sint);
        aux=0;
        for j=1:N-1
            aux=aux+cos((2.*auxi-1)./(2.*N).*j.*pi).*cos(j.*gam);
        end
        aux=1/2+aux;
        W=2/N;
        
    case 2
        gam=acos(sqrt((1+sint)/2));
        aux=0;
        for j=0:N-1
            aux=aux+sin(auxi.*pi./(2.*N+1)).*...
                sin(auxi.*pi./(2*N+1).*(2.*j+1)).*sin(2.*(j+1).*gam)./...
                (sin(gam));
        end
        W=4/(2.*N+1);
        
    case 3
        gam=acos(sqrt((1+sint)/2));
        aux=0;
        for j=0:N-1
            aux=aux+cos((2.*auxi-1)./(2.*N+1).*pi/2).*...
                cos((2.*auxi-1)./(2.*N+1).*pi./2.*(2.*j+1)).*...
                sqrt(2./(1+sint)).*cos((2.*j+1).*gam);
        end
        W=4/(2.*N+1);
        
    case 4
        gam=acos(sint);
        aux=0;
        for j=0:N-1
            aux=aux+sin(auxi.*pi./(N+1)).*...
                sin(auxi.*pi./(N+1).*(j+1)).*sin((j+1).*gam)./...
                (sin(gam));
        end
        W=2/(N+1);
        
end
clear auxi

%-------------------------------------------------------------------
%                           INTERPOLATION
%-------------------------------------------------------------------

sum=0;
for i=1:N
  sum=sum+aux(i).*phi(i);
end

phik=W.*sum;

end