function h = ScaleFreeBuilder(n,gamma,edges)
flag = true;
if edges == (n^2)/2 %If all edges are known, the function is very expensive, simply set all ones
    h = ones(n,n);
    flag = false;
end
if flag
    BAG = zeros(n);
    
    kmin = 5;
    c = (gamma-1)*(kmin^(gamma-1));
    alpha = 1/(1-gamma);
    eta = zeros(1,n);
    for i = 1:104
        eta(1,i) = c/(i^alpha);
    end
    etamean = mean(eta);
    links = 0;
    flag = 0;
    for ii = 1:n
        for jj = 1:n
            if flag == 0
                p_ij = (eta(1,ii)*eta(1,jj))/(etamean*n);
                r = rand();
                if r < p_ij && BAG(ii,jj) == 0
                    BAG(ii,jj) = 1;
                    BAG(jj,ii) = 1;
                    links = links + 1;
                    if links == edges
                        flag = 1;
                    end
                end
            else
                break;
            end
        end
    end
    h = zeros(n);
    I = randperm(n);
    for ii = 1:n
        for jj = 1:n
            if BAG(I(ii), I(jj)) == 1
                h(ii,jj) = 1;
            end
        end
    end
end
end