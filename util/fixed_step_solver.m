%% Model configurations
% Pflag=0 and Vtrfl=0 and Freqflag=0
function sol = fixed_step_solver(Fun, x0, t)
    N = length(t);
    dt = t(2) - t(1); 
    nx = length(x0);

    X = zeros(nx, N);
    X(:,1) = x0;

    for k = 2:N
        dx = Fun(t(k-1), X(:,k-1), k-1);  % k-1 because MATLAB is 1-based
        X(:,k) = X(:,k-1) + dt * dx;
    end

    sol.t = t;
    sol.X = X;
end