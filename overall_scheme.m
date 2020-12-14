%% Overall algorithm, Algorithm 1 in manuscript
%% This script 'overall_scheme.m' computes
%% the conditional probability mass function pi(T,x)


% Specify the system
ex = 1;

if ex == 1
    sys = @linprop;
    c = [1, 5, 1];
    n_unobs = 1; m_unobs = 1;
elseif ex == 2
    sys = @linprop_A;
    c = [5, 1, 1];
    n_unobs = 1; m_unobs = 2;
elseif ex == 3
    sys = @circuit;
    c = [0.3, 3, 0.5, 0.2, 0.06];
    n_unobs = 2; m_unobs = 0;
elseif ex == 4
    sys = @toggle;
    c =  [50, 16, 2.5, 1];
    %c = [20, 9, 2.5, 1];
    n_unobs = 1; m_unobs = 2;
elseif ex == 5
    sys = @seir;
    c = [0.05; 0.2; 0.05];
    n_unobs = 3; m_unobs = 1;
end

nu = feval(sys,'nu'); [n, m] = size(nu);
n_obs = n-n_unobs;
x0 = feval(sys,'x0');
T = feval(sys,'T');

Ns = 10000;

resampling_status = 1;
% resampling_status = 0, never resample
% resampling_status = 1, resample after each jump
% resampling_status = 2, adaptive resampling
num_zero = 10;
w_ratio = 1000;
r_count = 0;

%%
fid = fopen('obs_data.bin','r');
Nk = fread(fid,1, 'double');
tobs = fread(fid,Nk,'double');
yobs = fread(fid,Nk*n_obs,'double');
yobs = reshape(yobs,[n_obs  Nk]);

fid = fopen('full_data', 'r');
N = fread(fid, 1, 'double');
tarr = fread(fid, N,'double');
xarr = fread(fid,N*n,'double');
xarr = reshape(xarr, [n, N]);
fclose(fid);


%% -- Overwrite T --
%T = 20;  %edit for different truncated time
tobs_trunc = tobs(tobs<=T); 
tobs = tobs_trunc;
Nk = length(tobs);

%%

tic;

V = zeros(n_unobs,Ns); Vnew = zeros(n_unobs,Ns); w = zeros(1,Ns);
% Initialize.
t = 0; y =x0(n_unobs+1:end);
for i=1:Ns
    V(:,i)=x0(1:n_unobs); w(i) = 1;
end

for j=1:Nk
    for i=1:Ns
        [tarr_s,Varr_s,warr_s]=CTMC_filter_cont(V(:,i),w(i),t,tobs(j),y,sys,c,n_unobs,m_unobs);
        V_ = Varr_s(:,end); w_ = warr_s(end); y_ = y;
        [V(:,i),w(i)]=CTMC_filter_jump(V_,w_,y_,yobs(:,j),sys,c,n_unobs);
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

for i=1:Ns
    [tarr_s,Varr_s,warr_s]=CTMC_filter_cont(V(:,i),w(i),t,T,y,sys,c,n_unobs,m_unobs);
    V(:,i) = Varr_s(:,end); w(i) = warr_s(end);
end

toc;

%% conditional pmf

%interested in the ith component of X(t)
ind = 1;

x_min = min(V(ind,:));
x_max = max(V(ind,:));
x_states = x_min:x_max;
x_count = zeros(1, x_max-x_min+1);
for i = 1:x_max-x_min+1
    v_ind = (V(ind,:) == x_states(i));
    x_count(i) = v_ind * w';
end

pmf = x_count/sum(w);
m1 = V(ind,:) * w'/sum(w);
m2 = V(ind,:).^2 * w'/sum(w);
var_x = sqrt(m2 - m1^2);
xtrunc = xarr(ind,tarr<=T); 
x_true = xtrunc(end);
conf_intv = linspace(m1-var_x,m1+var_x,20);


figure
plot(x_states, pmf, '-x')
hold on
line([x_true, x_true],ylim, 'Color','r')
%plot(xarr(1,end), 0, 'r*')
plot(conf_intv, zeros(1,20),'-g')
plot(m1, 0, 'g*')
plot(m1-var_x, 0, '<g')
plot(m1+var_x, 0, '>g')
xlabel('the copy number of species S1')
ylabel('conditional probability mass function')
legend('conditional pmf from filtering', 'true exposed at time 40', ...
    'confidence interval')
%saveas(gcf, 'state_pmf.png')
