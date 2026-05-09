%% Dual Projected Gradient Descent with Advanced Mitigation
% Updated with:
% 1. Increased Restarts (numRestarts)
% 2. Early Exit on stagnation (to save time)
% 3. Annealed Weighting for binary potential (Soft Start)
% 4. Mild Jittering

clear; close all; clc;
rng(0);
algoName = "Dual PGD + Momentum + Advanced Mitigation";

%% Configuration
M = 61; N = 61;
[X,Y] = meshgrid(linspace(-1,1,N), linspace(-1,1,M));

% Target: 61x61 ICIP logo
aperture = zeros(M, N);
pad = [
  1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 0 0;
  1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 0;
  0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 0 0 0 1 1 1;
  0 0 0 0 1 1 0 0 0 0 0 1 1 1 0 0 0 0 1 1 1;
  0 0 0 0 1 1 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0;
  0 0 0 0 1 1 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0;
  0 0 0 0 1 1 0 0 0 0 0 1 1 1 0 0 0 0 1 1 1;
  0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 0 0 0 1 1 1;
  1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 0;
  1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 0 0;
  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
  1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 0 0;
  1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 0;
  0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 1 1 1;
  0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 1 1;
  0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 1 1 1;
  0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0;
  0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0;
  0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0;
  1 1 1 1 1 1 1 1 1 1 0 0 1 1 0 0 0 0 0 0 0;
  1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 0 0 0 0 0;
];
aperture(21:41, 21:41) = pad;

%% Target Fourier magnitude
diffraction = fftshift(fft2(aperture) / sqrt(M * N));
diffPattern = abs(diffraction);

%% Operators
g  = @(input) gradient_fn(input);
PA = @(y) proja_fn(y, M, N, diffPattern);
error_fn = @(input) norm(abs(fftshift(fft2(input) / sqrt(M * N))) - diffPattern, 'fro');
score_fn = @(input) error_fn(double(input > 0.5));

%% Optimization settings
stepsize0      = 0.05;
stepDecay      = 2e-4;
maxIter        = 10000;
beta           = 0.90;      
numTrials      = 100;       % Focus on success rate quality
numRestarts    = 15;        % [STRATEGY] Increased restarts to explore more basins
jitter0        = 5e-4;      % [STRATEGY] Mild noise to kick out of shallow traps
errorThreshold = 1e-10;
stagnationTol  = 1e-7;      % [STRATEGY] If score change is smaller than this, give up early
successCount   = 0;

allErrors = zeros(numTrials, 1);
tic;

for trial = 1:numTrials
    bestTrialScore = inf;
    bestTrialU = [];

    for rr = 1:numRestarts
        U = make_init(rr, diffPattern, M, N, PA);
        v = zeros(size(U));
        bestLocalScore = inf;
        
        lastScore = inf;
        
        for i = 1:maxIter
            % [STRATEGY] Annealed binary weight: Slowly increase the importance of binary constraint
            % This helps find the global phase first before forcing binary values.
            lambda = min(1, i / 2000); 
            
            mu = max(1e-3, stepsize0 / (1 + stepDecay * (i - 1)));
            grad = lambda * g(U); 
            
            v = beta * v + (1 - beta) * grad;
            U = U - mu * v;

            if jitter0 > 0
                U = U + (jitter0 / sqrt(i)) * randn(size(U));
            end

            U = PA(U);

            if mod(i, 100) == 0
                curScore = score_fn(U);
                
                % [STRATEGY] Early Exit / Stagnation Detection
                % If we are not improving, don't waste time on this restart.
                if abs(lastScore - curScore) < stagnationTol
                    break; 
                end
                lastScore = curScore;

                if curScore < bestLocalScore
                    bestLocalScore = curScore;
                    bestLocalU = U;
                end
                if curScore < errorThreshold, break; end
            end
        end

        if bestLocalScore < bestTrialScore
            bestTrialScore = bestLocalScore;
            bestTrialU = bestLocalU;
        end
        if bestTrialScore < errorThreshold, break; end
    end

    U_final = double(bestTrialU > 0.5);
    finalError = error_fn(U_final);
    allErrors(trial) = finalError;

    if finalError < errorThreshold
        successCount = successCount + 1;
    end
    
    if mod(trial, 10) == 0
        fprintf('Trial %d: Current Success Rate = %.2f%%\n', trial, (successCount/trial)*100);
    end
end

totalTime = toc;

%% Results Output
fprintf('\n=== RESULTS ===\n');
fprintf('Success rate: %.2f%%\n', (successCount / numTrials) * 100);
fprintf('Mean error: %.2e\n', mean(allErrors));
fprintf('Total Time: %.2f seconds\n', totalTime);

%% ===== Functions =====
function g = gradient_fn(input)
    U = input;
    g = 2 * (2 * U.^3 - 3 * U.^2 + U);
end

function y = proja_fn(y, M, N, diffPattern)
    x = fftshift(fft2(y) / sqrt(M * N));
    x = diffPattern .* exp(1i * angle(x));
    y = real(ifft2(ifftshift(x)) * sqrt(M * N));
end

function U0 = make_init(restartID, diffPattern, M, N, PA)
    if restartID == 1
        U0 = PA(rand(M, N));
        return;
    end
    % Random-phase initialization
    phase0 = 2 * pi * rand(M, N);
    z0 = diffPattern .* exp(1i * phase0);
    U0 = real(ifft2(ifftshift(z0)) * sqrt(M * N));
    U0 = min(max(U0, 0), 1);
    if restartID >= 5
        U0 = U0 + 0.1 * randn(M, N); % More diversity for late restarts
        U0 = min(max(U0, 0), 1);
    end
    U0 = PA(U0);
end