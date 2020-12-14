function output = circuit(info_type,varargin)

% DA+A --c1---> DA'.
% DA' ---c2---> D + A.
% DA  ---c3---> DA + A.
% DA' ---c4---> DA' + A.
% A   ---c5---> 0.

% Z = (DA, DA', A).
% X = (DA, DA').
% Y = A.

if strcmp(info_type,'nu')
 nu1 = [-1; 1; -1];
 nu2 = [1; -1; 1];
 nu3 = [0; 0; 1];
 nu4 = [0; 0; 1];
 nu5 = [0; 0; -1];
 
 nu = [nu1, nu2, nu3, nu4, nu5];
 output = nu;
 
elseif strcmp(info_type,'x0')
   output = [3; 0; 15];
elseif strcmp(info_type,'T')  
   output = 100;
elseif strcmp(info_type,'prop')
   x = varargin{1};
   c = varargin{2};
   output = prop(x,c);
else output=[];   
end
end

function a = prop(x,c)
a = [c(1)*x(1)*x(3); c(2)*x(2); c(3)*x(1); c(4)*x(2); c(5)*x(3)];
end