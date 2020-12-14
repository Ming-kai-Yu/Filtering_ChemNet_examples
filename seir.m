function output = seir(info_type,varargin)

% S + I --c1--> E + I. (unobservable)
% E --c2--> I. (observable)
% I --c3--> R. (observable) 
% X = (S, E, R)
% Y = (I)

if strcmp(info_type,'nu')
 nu1 = [-1; 1; 0; 0];
 nu2 = [0; -1; 0; 1];
 nu3 = [0; 0; 1; -1];
 
 nu = [nu1 nu2 nu3];
 output = nu;
 
elseif strcmp(info_type,'x0')
   output = [500; 20; 0; 5];
elseif strcmp(info_type,'T')  
   output = 40;
elseif strcmp(info_type,'prop')
   x = varargin{1};
   c = varargin{2};
   output = prop(x,c);
else output=[];   
end
end

function a = prop(x,c)
a = [c(1)*x(1)*x(4)/sum(x); c(2)*x(2); c(3)*x(4)];
end