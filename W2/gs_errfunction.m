function [ Y_initial ] = gs_errfunction( P0, Xobs )
%GS_ERRFUNCTION

N = length(Xobs)/2;

Y_initial = (Xobs(N+1:end) - P0 * Xobs(1:N))^2 + (P0 * Xobs(N+1:end) + Xobs(1:N))^2;

end

