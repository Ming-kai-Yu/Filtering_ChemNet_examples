function output = ex2(info_type,varargin)
% S ---> S + A. (unobservable) 
% 0 --> S, S ---> 0. (observable)
% X = # A
% Y = # S

if strcmp(info_type,'nu')
 nu1 = [1; 0]; 
 nu2 = [0; 1]; 
 nu3 = [0; -1];
 
 nu = [nu1 nu2 nu3];
 output = nu;
 
elseif strcmp(info_type,'x0')
   output = [0; 5];
elseif strcmp(info_type,'T')  
   output = 20;
elseif strcmp(info_type,'prop')
   x = varargin{1};
   c = varargin{2};
   output = prop(x,c);
else output=[];   
end
   
function a = prop(x,c)

a = [c(1)*x(2); c(2); c(3)*x(2)];
