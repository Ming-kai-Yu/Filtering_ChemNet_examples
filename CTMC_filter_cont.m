function [tarr,Varr,warr]=CTMC_filter_cont(V0,w0,t0,tf,y,sys,para,n_unobs,m_unobs)
% continuous evolution 

nu = feval(sys,'nu');


V=V0; w=w0;Varr=[V]; warr=[w];
t=t0;  tarr=[t];


while(t<tf)
    
    lambda = feval(sys,'prop',[V;y],para);
    lambda0 = sum(lambda);
    lambda_unobs = sum(lambda(1:m_unobs));
    lambda_obs = lambda0 - lambda_unobs;
    
    %tau = exprnd(1/lambda_unobs);
    tau = log(1/rand)/lambda_unobs;
    if (t + tau <= tf)
      r = rand*lambda_unobs;  
      q = cumsum(lambda(1:m_unobs));
      i=1;
      while (i<m_unobs && r > q(i))
          i = i+1;
      end
      
      V=V+nu(1:n_unobs,i);
      w=w*exp(-tau*lambda_obs);
      
      if sum(V<0) 
          V
      end
      
      t = t+tau;
      tarr = [tarr t];
      Varr = [Varr V];
      warr = [warr w];
    else 
      w = w*exp(-(tf-t)*lambda_obs);
      t = tf;
      tarr = [tarr t];
      Varr = [Varr V];
      warr = [warr w];
    end        
end %while
end %function