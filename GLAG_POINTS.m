function [s,t,W] = GLAG_POINTS(N,cas)
%GLAG_POINTS_CALC Gauss-Laguerre collocation and integration points.
%   [s,t,W] = GLAG_POINTS(N,cas) reads the GAUSS LAGUERRE
%   integration s, collocation t points and weights W for solving an SIE.
%   The calculation of the points is done by a different programme, 
%   GLAG_POINTS_CALC.
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


%-----------------------------------------------------------------------
%                             READ FILE
%-----------------------------------------------------------------------

% file name
p=mfilename('fullpath'); % Get the folder where function is installed
[fpath,~,~]=fileparts(p); % Directory without file name
fname=horzcat(fpath,'\DATA_GLAG\GLAG_POINTS_CASE',num2str(cas,'%u'),...
        '_N',num2str(N,'%03.0f'),'.dat');
    
% Matrix with data
try
    A=dlmread(fname);
catch
    errmsg=horzcat('File not found for this value of N or case! ',...
        'Please use GLAG_POINTS_CALC to calculate quadrature points!');
    error(errmsg);
end

% Vectors
s=A(1:end-1,2)';
W=A(1:end-1,3)';
t=A(:,4)';
  
end
