function [s,t,W] = GCHEB_POINTS(N,cas)
%GCHEB_POINTS Gauss Chebyshev collocation and integration points.
%   [s,t,W] = GCHEB_POINTS(N,cas) gives the GAUSS CHEBYSHEV collocation s,
%   integration t points and weights W for solving an SIE.
%
%   The variable 'cas' denotes the case for the weight function behaviour
%   w(s)=(1-s)^alpha*(1+s)^beta:
%
%        Case  | Behaviour at -1 | Behaviour at +1
%        ----------------------------------------
%       1 (I)  |     Singular    |     Singular
%       2 (II) |     Singular    |     Bounded
%       3 (III)|     Bounded     |     Singular
%       4 (IV) |     Bounded     |     Bounded
%       
%
%   University of Oxford 
%   Department of Engineering Science
%   Jhonatan Da Ponte Lopes, PhD 
%   Feb, 2019; Last revision: 2019-02-12

%-------------------------------------------------------------------
%                           CASE DEFINITION
%-------------------------------------------------------------------

auxi=1:N;
switch cas
    case 1 % Singular/Singular
        s=cos(pi.*(2.*auxi-1)./(2.*N));
        auxk=1:N-1;
        t=cos(pi.*auxk./N);
        W=1./N.*ones(1,N);
        
    case 2 % Singular/Bounded
        s=cos(pi.*(2.*auxi)./(2.*N+1));
        auxk=1:N;
        t=cos(pi.*(2.*auxk-1)./(2.*N+1));
        W=2.*(1-s)./(2.*N+1); 
        
    case 3 % Bounded/Singular
        s=cos(pi.*(2.*auxi-1)./(2.*N+1));
        auxk=1:N;
        t=cos(pi.*(2.*auxk)./(2.*N+1));
        W=2.*(1+s)./(2.*N+1);  
        
    case 4 % Bounded/Bounded
        s=cos(pi.*(auxi)./(N+1));
        auxk=1:N+1;
        t=cos(pi.*(2.*auxk-1)./(2.*N+2));
        W=(1-s.^2)./(N+1);  
end

end