function output = toggle(info_type,varargin)

% 0 --a1--> U. (unobservable)
% U --a2--> 0. (unobservable)
% 0 --a3--> V. (observable)
% V --a4--> 0. (observable)
% Y = (U, V) 

if strcmp(info_type,'nu')
 nu1 = [1; 0];
 nu2 = [-1; 0];
 nu3 = [0; 1];
 nu4 = [0; -1];
 nu = [nu1 nu2 nu3 nu4];
 output = nu;
 
elseif strcmp(info_type,'x0')
   output = [0; 0];
elseif strcmp(info_type,'T')  
   output = 600;
elseif strcmp(info_type,'prop')
   x = varargin{1};
   c = varargin{2};
   output = prop(x,c);
else output=[];   
end
end

function a = prop(x,c)
% c = [alpha1, alpha2, beta, gamma]
% e.g. c = [50, 16, 2.5, 1]
a = [c(1)/(1+x(2)^c(3)); x(1); c(2)/(1+x(1)^c(4)); x(2)];
end