function p = getpressure(U)

% U is size 3 x Nt
nphase = getnphase(U);
N = size(U,2);
p = zeros(N,1);

for k = 1:N
    if (nphase(k) == 1)
        p(k) = CO2.p_rhoT(U(1,k),U(3,k));
    elseif (nphase(k) == 2)
        p(k) = CO2.pVap(U(3,k));
    end
end


