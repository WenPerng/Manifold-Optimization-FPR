% Projected Gradient Descent
clear; close all; clc;

%% Initialization configurations
% Test Image
M = 201; % size of images
N = 201;

[X,Y] = meshgrid(linspace(-1,1,N), linspace(-1,1,M));

% aperture = rand(M, N);

% aperture = (abs(abs(X) - 0.1) < 0.05) .* (abs(Y) < 0.2);

% aperture = (abs(abs(X) - 0.05) < 0.01) .* (abs(abs(Y) - 0.05) < 0.01);

vx = [-0.43, -0.05, -0.10,  0.13,  0.05,  0.43,  0.18,  0.23, -0.22, -0.15];
vy = [ 0.00,  0.40,  0.05,  0.13, -0.38,  0.00,  0.43, -0.08, -0.3, -0.43];
in = inpolygon(X, Y, vx, vy);
aperture = double(in);

% aperture = (X .^ 2 + Y .^ 2) <= 0.2 .^ 2;

diffraction  = fftshift(fft2(aperture) / sqrt(M * N));
diffPattern = abs(diffraction);

% Gradient Descent --------------------------------------------------------
g = @(input) gradient(input);
PA = @(y) proja(y, M, N, diffPattern);

% analysis
error_fn = @(input) norm(abs(fftshift(fft2(input) / sqrt(M * N))) - diffPattern, 'fro');

%% Gradient Descent
stepsize = 0.05;
maxIter  = 2000;
beta     = 0.0;
eps0 = 1e-12;

numTrials = 1;

0;
successCount = 0;

allErrors = zeros(numTrials, 1);
errorThreshold = 1e-10;

tic;
for trial = 1:numTrials
    fprintf('Trial %d/%d\n', trial, numTrials);
    % Random initialization
    U = rand(M , N);
    grad   = zeros(size(U));
    Err = zeros(maxIter+1, 1);

    for i = 1 : maxIter
        grad = beta * grad + (1 - beta) * g(U);
        U = U - stepsize * grad;
        U = PA(U);
        Err(i+1) = error_fn(U);
    end

    U = U .* (U > 0);
    U = ones(size(U)) .* (U > 0.5);

    % Check success
    finalError = error_fn(U);
    allErrors(trial) = finalError;
    if finalError < errorThreshold
        successCount = successCount + 1;
    end
end

totalTime = toc;

%% Results
fprintf('\n=== RESULTS ===\n');
fprintf('Success: %d / %d\n', successCount, numTrials);
fprintf('Success rate: %.2f%%\n', (successCount/numTrials)*100);
fprintf('Min error: %.2e\n', min(allErrors));
fprintf('Max error: %.2e\n', max(allErrors));
fprintf('Mean error: %.2e\n', mean(allErrors));
fprintf('Time: %.2f seconds\n', totalTime);

%% Function
function g = gradient(input)
    U = input;
    g = 2 * (2 * U .^ 3 - 3 * U .^ 2 + U);
end

function y = proja(y, M, N, diffPattern)
    x = fftshift(fft2(y) / sqrt(M * N));
    x = diffPattern .* exp(1i * angle(x));
    y = real(ifft2(ifftshift(x) * sqrt(M * N)));
end