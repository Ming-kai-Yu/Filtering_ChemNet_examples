% Specify the system
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
m_obs = m-m_unobs;
x0 = feval(sys,'x0');
T = feval(sys,'T');


t=0; x=x0;
tarr = t; xarr = x0;
tobs = []; yobs = [];

while(t<T)
    lambda = feval(sys,'prop',x,c);
    lambda0 = sum(lambda);
    tau = exprnd(1/lambda0);
    if (t + tau <= T)
      r = rand*lambda0;  
      q = cumsum(lambda);
      i=1;
      while (r > q(i))
          i = i+1;
      end
      x = x + nu(:,i);
      t = t + tau;
      if i > m_unobs
          tobs = [tobs t]; yobs = [yobs x(n_unobs+1:n)];
      end
    else
        t = T;
        x = x;
    end
    tarr = [tarr t];
    xarr = [xarr x];
end

%% save trajectory
fido = fopen('obs_data.bin','w');
fidf = fopen('full_data.bin','w');

Nobs = length(tobs);
fwrite(fido,Nobs, 'double');
fwrite(fido,tobs,'double');
fwrite(fido,yobs,'double');
fclose(fido);
N = length(tarr);
fwrite(fidf,N, 'double');
fwrite(fidf,tarr,'double');
fwrite(fidf,xarr,'double');
fclose(fidf);

%%

figure
plot(tarr,xarr(1,:))
hold on
plot(tarr, xarr(2,:))
%plot(tarr, xarr(3,:))
%plot(tarr, xarr(4,:))
xlabel('time')
legend('S1','S2')
%legend('Susceptible', 'Exposed', 'Recovered', 'Infectious')
%saveas(gcf,'full_sys', 'epsc')

%% zoom in the trajectory on [0, Ts], Ts<=T
%{
Ts = 200;
tarrT = tarr(tarr<=Ts);
xarrT = xarr(:, tarr<=Ts);
figure
hold on
plot(tarrT, xarrT(1,:);
plot(tarrT, xarrT(2,:);
%}
