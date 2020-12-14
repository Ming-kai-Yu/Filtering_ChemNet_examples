function [Vnew, wnew] = resampling(V, w)
    wsum = sum(w);
    w = w/wsum;
    [n_unobs, Ns] = size(V);
    Vnew = zeros(n_unobs, Ns);
    wnew = zeros(1, Ns);
    % Implement the branching algorithm.
    offs = offsprings(w);
    i=1; ind=1;
    for i=1:Ns
        wnew(i)=1;
        for l=1:offs(i)
            Vnew(:,ind) = V(:,i);
            ind = ind+1;
        end
    end
end