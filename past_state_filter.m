function [V, w] = past_state_filter...
(n_unobs, n_obs, m_unobs, x0, Ns, Nk, tobs, T_interm, T, yobs, sys, c)

tobs_before = tobs(tobs<=T_interm);
yobs_before = yobs(yobs<=T_interm);
tobs_after = tobs(tobs>T_interm);
yobs_after = yobs(tobs>T_interm);

N1 = length(tobs_before);
N2 = length(tobs_after);

V = zeros(n_unobs*2,Ns); 
Vnew = zeros(n_unobs*2,Ns); 
w = zeros(1,Ns);

% Initialize.
t = 0; y =x0(n_unobs+1:end);
for i=1:Ns
    V(1:n_unobs,i)=x0(1:n_unobs); w(i) = 1;
end


% From t=0 to the most recent jump before T_interm
for j=1:N1
    for i=1:Ns
        [tarr_s,Varr_s,warr_s]=CTMC_filter_cont(V(1:n_unobs,i),w(i),t,tobs(j),y,sys,c,n_unobs,m_unobs);
        V_ = Varr_s(:,end); w_ = warr_s(end); y_ = y;
        [V(1:n_unobs,i), w(i)]=CTMC_filter_jump(V_,w_,y_,yobs(:,j),sys,c,n_unobs);
    end
    
    wsum = sum(w);
    w = w/wsum;
    % Implement the branching algorithm.
    offs = offsprings(w);
    i=1; ind=1;
    for i=1:Ns
        w(i)=1;
        for l=1:offs(i)
            Vnew(:,ind) = V(:,i);
            ind = ind+1;
        end
    end
    V = Vnew;
    
    t = tobs(j);
    y = yobs(:,j);
end

% From tobs_before(end) to T_interm
for i=1:Ns
    [tarr_s,Varr_s,warr_s]=CTMC_filter_cont(V(1:n_unobs,i),w(i),t,T_interm,y,sys,c,n_unobs,m_unobs);
    V(1:n_unobs,i) = Varr_s(:,end); w(i) = warr_s(end);
    t = T_interm;
end
V(n_unobs+1:n_unobs*2,:) = V(1:n_unobs,:);

% From T_interm to tobs_after(end)
for j=1:N2
    for i=1:Ns
        [tarr_s,Varr_s,warr_s]=CTMC_filter_cont(V(1:n_unobs,i),w(i),t,tobs_after(j),y,sys,c,n_unobs,m_unobs);
        V_ = Varr_s(1:n_unobs,end); w_ = warr_s(end); y_ = y;
        [V(1:n_unobs,i), w(i)]=CTMC_filter_jump(V_,w_,y_,yobs_after(:,j),sys,c,n_unobs);
    end
    
    wsum = sum(w);
    w = w/wsum;
    % Implement the branching algorithm.
    offs = offsprings(w);
    i=1; ind=1;
    for i=1:Ns
        w(i)=1;
        for l=1:offs(i)
            Vnew(:,ind) = V(:,i);
            ind = ind+1;
        end
    end
    V = Vnew;
    
    t = tobs_after(j);
    y = yobs_after(:,j);
end

% From tobs_after(end) to T
for i=1:Ns
    [tarr_s,Varr_s,warr_s]=CTMC_filter_cont(V(1:n_unobs,i),w(i),t,T,y,sys,c,n_unobs,m_unobs);
    V(1:n_unobs,i) = Varr_s(:,end); w(i) = warr_s(end);
    t = T;
end
end