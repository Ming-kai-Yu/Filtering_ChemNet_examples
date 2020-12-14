%% Past state estimation, Algorithm 6 in manuscript
%% This script 'past_state_scheme.m' computes P[X(T_interm)|Y(T)]
%% where 0<= T_interm <=T



sys = @seir; 
%Z = (#S, #E, #R, #I)
c = [0.05; 0.2; 0.05];
n_unobs = 3; m_unobs = 1;
n_theta = 0;

x_ind = [2]; 
v_ind = 5; % interested in V(v_ind,:)

% initial condition 
x0 = feval(sys,'x0');
Ns = 10000;



nu = feval(sys,'nu'); [n, m] = size(nu);
n_obs = n-n_unobs;
x0 = feval(sys,'x0');
T = feval(sys,'T');

%% load observation
fid = fopen('obs_data.bin','r');
Nk = fread(fid,1, 'double');
tobs = fread(fid,Nk,'double');
yobs = fread(fid,Nk*n_obs,'double');
yobs = reshape(yobs,[n_obs  Nk]);
fclose(fid);


fid = fopen('full_data.bin', 'r');
N = fread(fid, 1, 'double');
tarr = fread(fid, N,'double');
xarr = fread(fid,N*n,'double');
xarr = reshape(xarr, [n, N]);
fclose(fid);

%% -- Overwrite T to get conditional pmf at time T.
T_interm = 0;
T = 0; 
tobs_trunc = tobs(tobs<=T); %edit for different truncated time
tobs = tobs_trunc;

tobs_before = tobs(tobs<=T_interm);
yobs_before = yobs(yobs<=T_interm);
tobs_after = tobs(tobs>T_interm);
yobs_after = yobs(tobs>T_interm);


N1 = length(tobs_before);
N2 = length(tobs_after);

resampling_status = 1;
% resampling_status = 0, never resample
% resampling_status = 1, resample after each jump
% resampling_status = 2, adaptive resampling
num_zero = 10;
w_ratio = 1000;
r_count = 0;



%%
tic;


% Initialize.
t = 0; y =x0(n_unobs+1:end);

% for seir example, x0 follows a binomial prior
x2 = binornd(x0(1)+x0(2), 0.04, 1, Ns); 
x1 = x0(1) + x0(2) - x2;
V = zeros(n_unobs*2,Ns); 
Vnew = zeros(n_unobs*2,Ns); 
w = ones(1,Ns);
V(1,:) = x1;
V(2,:) = x2;
V(3,:) = ones(1,Ns)*x0(3);

%{
% for 
for i=1:Ns
    V(:,i)=x0(1:n_unobs); w(i) = 1;
end
%}



% From t=0 to the most recent jump before T_interm
for j=1:N1
    for i=1:Ns
        [tarr_s,Varr_s,warr_s]=CTMC_filter_cont(V(1:n_unobs,i),w(i),t,tobs(j),y,sys,c,n_unobs,m_unobs);
        V_ = Varr_s(:,end); w_ = warr_s(end); y_ = y;
        [V(1:n_unobs,i), w(i)]=CTMC_filter_jump(V_,w_,y_,yobs(:,j),sys,c,n_unobs);
    end
    
    if resampling_status == 1 % reample after each jump
        [V, w] = resampling(V, w);
    elseif resampling_status == 0 % no resampling, average weight renormalized to be one
        w = w/sum(w)*Ns;
    elseif resampling_status == 2 % adaptive resampling
        wnz = w(w ~= 0); % nonzero elements of w
        if sum(w == 0) > num_zero || max(wnz)/min(wnz)> w_ratio
            [V, w] = resampling(V, w);
            r_count = r_count + 1;
        else
            w = w/sum(w)*Ns;
        end
    end
    
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
    
   if resampling_status == 1 % reample after each jump
        [V, w] = resampling(V, w);
    elseif resampling_status == 0 % no resampling, average weight renormalized to be one
        w = w/sum(w)*Ns;
    elseif resampling_status == 2 % adaptive resampling
        wnz = w(w ~= 0); % nonzero elements of w
        if sum(w == 0) > num_zero || max(wnz)/min(wnz)> w_ratio
            [V, w] = resampling(V, w);
            r_count = r_count + 1;
        else
            w = w/sum(w)*Ns;
        end
    end
    
    t = tobs_after(j);
    y = yobs_after(:,j);
end

% From tobs_after(end) to T
for i=1:Ns
    [tarr_s,Varr_s,warr_s]=CTMC_filter_cont(V(1:n_unobs,i),w(i),t,T,y,sys,c,n_unobs,m_unobs);
    V(1:n_unobs,i) = Varr_s(:,end); w(i) = warr_s(end);
    t = T;
end
toc;  


%% conditional pmf for P[X(T_interm)|Y(T)]
x_min = 0;
x_max = 40;
x_states = x_min:x_max;
x_count = zeros(1, x_max-x_min+1);
for i = 1:x_max-x_min+1
    vx_ind = (V(v_ind,:) == x_states(i));
    x_count(i) = vx_ind * w';
end

pmf = x_count/sum(w);





%% Ploting
%{
if T == 0
    pmf_0 = pmf;
elseif T == 10
    pmf_10 = pmf;
elseif T == 20
    pmf_20 = pmf;
elseif T == 40
    pmf_40 = pmf;
elseif T == 50
    pmf_50 = pmf;
elseif T == 100
    pmf_100 = pmf;
elseif T == 200
    pmf_200 = pmf;
end


xtrunc = xarr(x_ind,tarr<=T_interm); %edit for different truncated time
x_true = xtrunc(end);


figure
hold on
plot(x_states, pmf_0, '-x', 'LineWidth', 2)
plot(x_states, pmf_20, '-+', 'LineWidth', 2)
plot(x_states, pmf_40, '-*', 'LineWidth', 2)
line([x_true, x_true],ylim, 'Color','r')
%plot(xarr(1,end), 0, 'r*')
%plot(conf_intv, zeros(1,20),'-g')
%plot(m1, 0, 'g*')
%plot(m1-var_x, 0, '<g')
%plot(m1+var_x, 0, '>g')
%xlabel('the number of susceptible')
%xlabel('the number of exposed')
xlabel('P[X(0) = x | Y(T)]')
ylabel('conditional probability mass function')
%legend('conditional pmf from filtering', 'true exposed at time 40', ...
%    'confidence interval')
lgd = legend('P[X(0)| Y(0)]', 'P[X(0)| Y(20)]','P[X(0)| Y(40)]');
lgd.FontSize = 14;
hold off
saveas(gcf, 'seir_x0_cpmf', 'epsc')
%}

