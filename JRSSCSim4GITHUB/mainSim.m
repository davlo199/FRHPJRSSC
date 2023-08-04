function SIM = mainSim(NSIM, alpha, beta, gamma, mu, M0, MMAX, lambda, theta,c)
    Out = zeros(2, NSIM);
    t = 0;
    n = 0;
    epsilon = 1e-10;
    SimPoints = [];
    MAG = [];

    % Simulate Points
    while length(SimPoints) < 1
        M = mu;
        E = exprnd(1 / M, 1, 1);
        t = t + E;
        U = rand;

        if (U < (mu / M))
            n = n + 1;
            SimPoints = [SimPoints, t];
            magexp = exprnd(theta, 1, 1) + M0; % %W-a random variable
            magPar = gprnd(1 / lambda, M0/lambda, M0, 1, 1);
            mag = min(MMAX, min(magexp, magPar));
            MAG = [MAG, mag];
        end

    end

    while length(SimPoints) < NSIM
        M = mu + alpha * ((exp(gamma * (MAG - M0))) ...
            * (MLapp(t+epsilon-SimPoints,c,beta,beta).*(t+epsilon-SimPoints).^(beta-1)).');
        E = exprnd(1 / M, 1, 1);
        t = t + E;
        U = rand;

        if (U < (mu + alpha * ((exp(gamma * (MAG - M0))) * (MLapp(t-SimPoints,c,beta,beta) ...
                .* (t - SimPoints).^(beta - 1)).')) / M)
            n = n + 1;
            SimPoints = [SimPoints, t];
            magexp = exprnd(theta, 1, 1) + M0; % %W-a random variable
            magPar = gprnd(1 / lambda, M0/lambda, M0, 1, 1);
            mag = min(MMAX, min(magexp, magPar));
            MAG = [MAG, mag];
        end

    end

    Out(1, :) = SimPoints; %Times of events
    Out(2, :) = MAG; %Magnitude of events
    SIM = Out;
end
