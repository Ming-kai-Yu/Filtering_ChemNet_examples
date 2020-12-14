function [V,w]=CTMC_filter_jump(V_,w_,y_,y,sys,para,n_unobs)
% 


if w_==0 
    w=0; V=V_;
else
nu = feval(sys,'nu'); [n m] = size(nu);
lambda = feval(sys,'prop',[V_;y_],para);


J = (nu(n_unobs+1:n,:) == y-y_);
J = find(prod(J,1));  % J contains the indices of the observable reactions
                      % cosistent with y - y_.

                      
if length(J)==1
    w = w_*lambda(J);
    if lambda(J)>0
        V = V_ + nu(1:n_unobs,J);
    else
        V = V_;
    end
else
    r = rand*length(J);
    q = 1:length(J);
    ind=find(q>r,1);
    w = w_*lambda(J(ind));
    if lambda(J(ind))>0
        V = V_ + nu(1:n_unobs,J(ind));
    else
        V = V_;
    end
end
end
    
    
    
    