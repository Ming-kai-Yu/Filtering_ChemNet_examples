%% Specify the system
ex = 3;


if ex == 1
    sys = @linprop;
    c = [1, 5, 1];
    n_unobs = 1; m_unobs = 1;
elseif ex == 2
    sys = @ex2_A;
    c = [5, 1, 1];
    n_unobs = 1; m_unobs = 2;
elseif ex == 3
    sys = @circuit;
    c = [0.3, 3, 0.5, 0.2, 0.06];
    n_unobs = 2; m_unobs = 0;
elseif ex == 4
    sys = @toggle;
    %c =  [50, 16, 2.5, 1];
    c = [20, 9, 2.5, 1];
    n_unobs = 1; m_unobs = 2;
elseif ex == 5
    sys = @seir;
    c = [0.05; 0.2; 0.05];
    n_unobs = 3; m_unobs = 1;
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


fid = fopen('full_data.bin', 'r');
N = fread(fid, 1, 'double');
tarr = fread(fid, N,'double');
xarr = fread(fid,N*n,'double');
xarr = reshape(xarr, [n, N]);
fclose(fid);

T = 3;
tobs = tobs(tobs<=T);
yobs = yobs(tobs<=T);
xarr = xarr(:, tobs<=T);
tarr = tarr(tobs<=T);
Nk = length(tobs);

%%
v_ind = 1; % interested in V(v_ind,:)
V_est = zeros(Nk+1, 1);
V_err = zeros(Nk+1, 1);

V = zeros(n_unobs,Ns); 
Vnew = zeros(n_unobs,Ns); 
w = zeros(1,Ns);


tic;

% Initialize.
t = 0; y =x0(n_unobs+1:end);

for i=1:Ns
    V(1:n_unobs,i)=x0(1:n_unobs); w(i) = 1;
end



for j=1:Nk
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
    
    m1 = V(v_ind,:)*w'/sum(w);
    V_est(j) = m1;
    m2 = V(v_ind,:).^2 *w'/sum(w);
    V_err(j) = sqrt(m2 - m1^2);
end

for i=1:Ns
    if n_theta ~= 0
          c(c_ind) = V(n_unobs+1:n_unobs+n_theta, i);
    end
    [tarr_s,Varr_s,warr_s]=CTMC_filter_cont(V(1:n_unobs,i),w(i),t,T,y,sys,c,n_unobs,m_unobs);
    V(1:n_unobs,i) = Varr_s(:,end); w(i) = warr_s(end);
    
    m1 = V(v_ind,:)*w'/sum(w);
    V_est(Nk+1) = m1;
    m2 = V(v_ind,:).^2 *w'/sum(w);
    V_err(Nk+1) = sqrt(m2 - m1^2);
end
tobs = [tobs; T];

toc;  


%%
figure
hold on
errorbar(tobs, V_est, V_err,'*')
%plot(tobs, V_est)
%plot(tarr, xarr(v_ind,:))
%line([0 T], [c_true, c_true], 'Color', 'r' )
%line([0 T], [20, 20], 'Color', 'r' )
%plot(tobs, yobs)
stairs(tarr, xarr(1,:), 'LineWidth', 2)
xlabel('time')
%ylabel('population')
ylabel('parameter')
%legend('filtered susceptible', 'actual susceptible')
%legend('filtered exposed', 'actual exposed')
%legend('filtered S1', 'actual S1', 'S2')
%legend('Filtered alpha2', 'Actual alpha2')
%legend('Filtered Exposed at initial time', 'Actual Initial Exposed')
%legend('actual', 'filtered', 'sample')
%legend('filtered kappa','actual kappa')
legend('filtered D_A', 'actual D_A')
hold off
saveas(gcf,'seir-kappa-traj','epsc')

