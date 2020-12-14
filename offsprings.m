function offs = offsprings(w)
% The branching process algorithm 9.2.1 from Bain and Crisan.
% offs is an array consisting of the number of offsprings 
% offs(i) = the number of offsprings of the ith trajectory whose 
% weight is w(i). 
% The algorithm produces offs(i) randomly so that 
% E{offs(i)} = w(i) for all i 
% and for any k between 1 and Ns (Ns = length of w) 
% Var(\sum_i=^k offs(i)) = 

Ns = length(w);
u = rand([1 Ns-1]);
g = Ns; h = Ns;

for j=1:Ns-1
    Nw = Ns*w(j); Nwfloor = floor(Nw);
    Nwfrac = Nw - Nwfloor;
    gminusNwfrac = (g-Nw) - floor(g-Nw);
    gfloor = floor(g);
    gfrac = g-gfloor;
    
    if Nwfrac + gminusNwfrac < 1 
        if u(j) < 1 - Nwfrac/gfrac 
            offs(j) = Nwfloor;
        else
            offs(j) = Nwfloor + (h-gfloor);
        end
    else
        if u(j) < 1 - (1 - Nwfrac)/(1-gfrac)
            offs(j) = Nwfloor + 1;
        else
            offs(j) = Nwfloor + (h - gfloor);
        end
    end
    g = g - Nw;
    h = h - offs(j);
end
offs(Ns) = h;

            