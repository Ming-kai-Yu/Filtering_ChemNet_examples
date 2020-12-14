%% The script 'Bayesian_para_scheme.m' performs
%% the overall scheme for Bayesian inference of parameter one at a time

%% Specify the system
ex = 4;

if ex == 1
    sys = @linprop;
    c = [1, 5, 1];
    n_unobs = 1; m_unobs = 1;
    c_ind = [3];
    n_theta = length(c_ind);
    theta_min = 0;
    theta_max = 2;
elseif ex == 2
    sys = @elinprop_A;
    c = [5, 1, 1];
    n_unobs = 1; m_unobs = 2;
    c_ind = [1];
    n_theta = length(c_ind);
    theta_min = 3;
    theta_max = 8;
elseif ex == 3
    sys = @circuit;
    c = [0.3, 3, 0.5, 0.2, 0.06];
    n_unobs = 2; m_unobs = 0;
    
    n_theta = 1;
    c_ind = 2;
    theta_min = 0;
    theta_max = 6;
    c_true = c(c_ind);
elseif ex == 4
    sys = @toggle;
    %c =  [50, 16, 2.5, 1];
    c = [20, 9, 2.5, 1];
    n_unobs = 1; m_unobs = 2;
   
    n_theta = 1;
    c_ind = 1;
    theta_min = 10;
    theta_max = 30;
    c_true = c(c_ind);
elseif ex == 5
    sys = @seir;
    c = [0.05; 0.2; 0.05];
    n_unobs = 3; m_unobs = 1;
    
    n_theta = 1;
    c_ind = 2;
    theta_min = 0.1;
    theta_max = 0.3;
    c_true = c(c_ind);
end

nu = feval(sys,'nu'); [n m] = size(nu);
n_obs = n-n_unobs;
x0 = feval(sys,'x0');
T = feval(sys,'T');

Ns = 1000;
resampling_status = 1;
% resampling_status = 0, never resample
% resampling_status = 1, resample after each jump
% resampling_status = 2, adaptive resampling
num_zero = 10;
w_ratio = 1000;
r_count = 0;


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



%% -- Overwrite T --
T = 20; % edit for different truncated time
tobs_trunc = tobs(tobs<=T); 
tobs = tobs_trunc;
Nk = length(tobs);

%%
v_ind = 2; % interested in V(v_ind,:)

V = zeros(n_unobs + n_theta,Ns); 
Vnew = zeros(n_unobs + n_theta,Ns); 
w = zeros(1,Ns);


tic;

% Initialize.
t = 0; y =x0(n_unobs+1:end);

for i=1:Ns
    V(1:n_unobs,i)=x0(1:n_unobs); w(i) = 1;
end

if n_theta ~= 0
    prior = theta_min + (theta_max - theta_min)*rand(Ns, 1);
    V(n_unobs+1,:) = prior;
end



for j=1:Nk
    for i=1:Ns
        if n_theta ~= 0
          c(c_ind) = V(n_unobs+1, i);
        end
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
        end       
    end
    
    t = tobs(j);
    y = yobs(:,j);
    
end

for i=1:Ns
    if n_theta ~= 0
          c(c_ind) = V(n_unobs+1:n_unobs+n_theta, i);
    end
    [tarr_s,Varr_s,warr_s]=CTMC_filter_cont(V(1:n_unobs,i),w(i),t,T,y,sys,c,n_unobs,m_unobs);
    V(1:n_unobs,i) = Varr_s(:,end); w(i) = warr_s(end);
   
end
tobs = [tobs; T];

toc;  

fprintf('resampling happened %d times,', r_count);
fprintf('observable reaction happened %d times.\n', length(tobs));


%% posterior pdf

nbins = 40;
params_intv = linspace(theta_min, theta_max, nbins+1);
freq_count = zeros(1, nbins);
bin_width = params_intv(2)-params_intv(1);

for i = 1:nbins
     ind = (params_intv(i) <= V(v_ind,:) & V(v_ind,:)< params_intv(i+1));
     freq_count(i) = ind * w';
 end
pmf = freq_count/(sum(w)*bin_width);



m1 = V(v_ind,:) * w'/sum(w);
m2 = V(v_ind,:).^2 * w'/sum(w);
std_x = sqrt(m2 - m1^2);
conf_intv = linspace(m1-std_x,m1+std_x,20);
fprintf('T = %g\n. ', T);
fprintf('95 percent confidence interval is [%g, %g].\n', m1-2*std_x, m1+2*std_x);

if T == 0
    pmf_0 = pmf;
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
elseif T == 600
    pmf_600 = pmf;
elseif T == 1000
    pmf_1000 = pmf;
elseif T == 2000
    pmf_2000 = pmf;
end

%%

figure
hold on
params_mesh = params_intv(1:end-1) + 0.5*(params_intv(2)-params_intv(1));

plot(params_mesh, pmf_0, ':', 'LineWidth', 4);
plot(params_mesh, pmf_20, '-.', 'LineWidth', 4);
%plot(params_mesh, pmf_40, '-+', 'LineWidth', 2);
%plot(params_mesh, pmf_50, '-','LineWidth', 2);
plot(params_mesh, pmf_100, '-', 'LineWidth', 3);
plot(params_mesh, pmf_600, '-+', 'LineWidth', 2);
%plot(params_mesh, pmf_2000, '-', 'LineWidth', 2);
hold on
line([c_true, c_true],ylim, 'Color','r', 'LineWidth', 1.5)


%xlabel('the number of susceptible')
xlabel('parameter $\alpha_1$', 'FontSize', 14, 'interpreter', 'latex')
ylabel('conditional distribution', 'FontSize', 14)
lgd = legend('T=0', 'T=20', 'T=100', 'T=600', 'true parameter');
%lgd = legend('T = 0', 'T=20', 'T=50', 'T=100', 'true parameter');
%lgd = legend('T = 0', 'T=20', 'T=40', 'true parameter');
lgd.FontSize = 14;
lgd.Location = 'northwest';
hold off
saveas(gcf,'toggle-alpha1-cpdf-20-9','epsc')

