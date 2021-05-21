function [sd,tk,Wd] = GLAG_POINTS_CALC(N,cas)
%GLAG_POINTS_CALC Gauss-Laguerre collocation and integration points.
%   [sd,Wd,tk,lk] = GLAG_POINTS_CALC(N,cas) calculates the GAUSS LAGUERRE
%   integration sd, collocation tk points and weights Wd for solving an SIE
%   and write the values to a file to be used later.
%
%   The variable 'cas' denotes the case for the weight function behaviour
%   w(s)=s^alpha*e^(-s):
%
%        Case  | Behaviour at s=0 | Value of alpha
%        ----------------------------------------
%       1 (I)  |     Bounded      |       +1/2
%       2 (II) |     Singular     |       -1/2
%
%   Ioakimidis, N.I. 1980. Application of the Gauss- and Radau-Laguerre
%       quadrature rules to the numerical solution of cauchy type singular
%       integral equations.
%
%   University of Oxford 
%   Department of Engineering Science
%   Jhonatan Da Ponte Lopes, PhD 
%   April, 2019; Last revision: 2019-04-10

% Tested with alpha=-0.3, N=10; (Compared with Ioakimidis)

%-----------------------------------------------------------------------
%                            CASE DEFINITION
%-----------------------------------------------------------------------

if cas==1
    a=sym(1/2);
elseif cas==2
    a=sym(-1/2);
else
    error('Case not recognized!');
end

%-----------------------------------------------------------------------
%                     INTEGRATION POINTS AND WEIGHTS
%-----------------------------------------------------------------------

% Define symbolic variables
syms s 'real'
Ns=sym(N);

% Laguerre polynomial
Ln=laguerreL(Ns,a,s);
p=coeffs(Ln,s);
p=flip(p);

% Integration points
si=roots(p);
sd=double(si);

% Weights
W=gamma(Ns+a+1).*si./(factorial(Ns).*((Ns+1).*(laguerreL(Ns+1,a,si))).^2);
Wd=double(W);

%-----------------------------------------------------------------------
%                           COLLOCATION POINTS
%-----------------------------------------------------------------------

% Initialize vectors
tk=zeros(N+1,1);
lk=tk;

% Find roots
for i=1:N+1
    % Display iteration
    fprintf('\n\nk = %i \n', i);
    if i==N+1
        [tk(i), lk(i)]=newtonraphson_glag(N,cas,1.05*sd(i-1));
        %[tk(i), lk(i), flag(i), ~]=fzero(lambda,1.05*sd(i-1),Options);
    elseif i==1
        [tk(i), lk(i)]=newtonraphson_glag(N,cas,0.5*sd(i));
        %[tk(i), lk(i), flag(i), ~]=fzero(lambda,sd(i),Options);
    else
        [tk(i), lk(i)]=newtonraphson_glag(N,cas,(sd(i)+sd(i-1))/2);
        %[tk(i), lk(i), flag(i), ~]=fzero(lambda,[sd(i) sd(i-1)],Options);
    end
end

%-----------------------------------------------------------------------
%                                 CHECKS
%-----------------------------------------------------------------------

uniq_s=numel(uniquetol(sd));
uniq_W=numel(uniquetol(Wd,min(Wd)));
uniq_t=numel(uniquetol(tk));

fprintf('%i unique elements in s\n',uniq_s);
fprintf('%i unique elements in W\n',uniq_W);
fprintf('%i unique elements in t\n',uniq_t);

errvar=0;

% Check number of s
if uniq_s~=N
    warning('Did not find all the integration points s!');
    errvar=errvar+1;
end

% Check number of W
if uniq_W~=N
    warning('Did not find all the weights W!');
    errvar=errvar+1;
end

% Check number of t
if uniq_t~=N+1
    warning('Did not find all the collocation points tk!');
    errvar=errvar+1;
end

% Check t's are solutions
l_nonconv=numel(find(lk>1e-10));
if l_nonconv>0
    warning('Some tk did not converge!');
    errvar=errvar+1;
end

%-----------------------------------------------------------------------
%                              PRINT VALUES
%-----------------------------------------------------------------------

% Print to screen
T=table((1:N+1)',[sd;NaN],[Wd;NaN],tk,...
    'VariableNames',{'i','si','Wi','tk'});
disp(T);

% Print if there are no errors!
if errvar==0
    % Matrix with values of calculated s,W,t
    A=zeros(N+1,4);
    A(:,1)=1:N+1;
    A(1:end-1,2)=sd;
    A(end,2)=NaN;
    A(1:end-1,3)=Wd;
    A(end,3)=NaN;
    A(:,4)=tk;
    % Save to file
    fname=horzcat(cd,'\DATA\GLAG_POINTS_CASE',num2str(cas,'%u'),...
        '_N',num2str(N,'%03.0f'),'.dat');
    fid=fopen(fname,'w+'); % Open file, discards any existing content
    for i=1:N+1
        if i==N+1
            fprintf(fid,'%i \t %16.16e \t\t %16.16e \t\t %16.16e \r\n',A(i,:));
        else
            fprintf(fid,'%i \t %16.16e \t %16.16e \t %16.16e \r\n',A(i,:));
        end
    end
    fclose(fid);
end



end

%-----------------------------------------------------------------------
%                           NEWTON-RAPHSON SOLVER
%-----------------------------------------------------------------------

function [tzero, fzero]=newtonraphson_glag(N,cas,t0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial values

% Get a
if cas==1
    a=sym(1/2);
elseif cas==2
    a=sym(-1/2);
else
    error('Case not recognized!');
end

% Symbolic N for computations:
Nsym=sym(N);

% We are solving for z, z=-t!

% tk's are zeros of lambda function

f= @(z) hypergeom(Nsym+1,1-a,vpa(z));                 % Lambda Function
fp= @(z) 2*(Nsym+1)*hypergeom(Nsym+2,1-a+1,vpa(z));   % Lambda Derivative

zold=-t0;   % Initial guess
k=0;        % Number of iterations
err=1;      % Initial error
maxk=100;    % Maximum number of iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize headers

disp(' ');
fprintf('_______________________________________________________________');
fprintf('\n\nGauss-Laguerre Quadrature\n');
fprintf('\n\nCollocation Point Search\n');
fprintf('\nInitial point: t = %6.4e \n', t0);
fprintf('\n...........\n');

fprintf('\t k \t\t x \t\t\t\t f(x) \t\t\t err\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton-Raphson solver

while (err>=eps && k<=maxk)
    fold=f(zold);           % Function at current step
    fpold=fp(zold);         % Derivative at current step
    znow=zold-fold/fpold;   % Next iteration of z
    fnow=f(znow);           % Update function
    err=abs(znow-zold);     % Error (difference between z's in two steps)
    zold=znow;              % Update z
    k=k+1;                  % Update k
    fprintf('\t %i \t %6.4e \t %+6.4e \t %6.4e\n',k,-znow,fnow,err);
end

tzero=-znow;
fzero=fnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print solved value

fprintf('\nConvergence after %i iterations\n',k);
fprintf('\nThe root is: \n');
fprintf('\n\t tk = %6.4e\n',tzero);
fprintf('_______________________________________________________________');



end