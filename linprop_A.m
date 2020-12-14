function output = ex2_A(info_type,varargin)

% 0 --c1--> S. (unobservable)
% S --c2--> 0. (unobservable)
% S --c3--> S + A. (observable) 
% X = # S
% Y = # A

if strcmp(info_type,'nu')
 nu1 = [1; 0]; 
 nu2 = [-1; 0]; 
 nu3 = [0; 1];
 
 nu = [nu1 nu2 nu3];
 output = nu;
 
elseif strcmp(info_type,'x0')
   output = [5; 0];
elseif strcmp(info_type,'T')  
   output = 40;
elseif strcmp(info_type,'prop')
   x = varargin{1};
   c = varargin{2};
   output = prop(x,c);
else output=[];   
end
   
function a = prop(x,c)

a = [c(1); c(2)*x(1); c(3)*x(1)];
