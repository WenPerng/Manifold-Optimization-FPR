% Projected Gradient Descent
% Run 2000 trials and compute success rate
% Success: final error < 1e-10
clear; close all; clc;
algoName = "Projected Gradient Descent (Fourier magnitude loss + 0/1 prior projection)";

%% Configuration
M = 101;
N = 101;

[X,Y] = meshgrid(linspace(-1,1,N), linspace(-1,1,M));

% Double slit
% aperture = (abs(abs(X) - 0.1) < 0.05) .* (abs(Y) < 0.2);

% Polygon
vx = [-0.43, -0.05, -0.10,  0.13,  0.05,  0.43,  0.18,  0.23, -0.22, -0.15];
vy = [ 0.00,  0.40,  0.05,  0.13, -0.38,  0.00,  0.43, -0.08, -0.3, -0.43];
in = inpolygon(X, Y, vx, vy);
aperture = double(in);

% 61x61 NTUEE
% aperture = zeros(M, N);
% pad = [
%   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%   1 1 1 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1;
%   0 1 1 0 0 0 1 1 1 1 1 1 1 1 1 1 0 0 1 1 0;
%   0 1 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0;
%   0 1 1 1 1 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0;
%   0 1 1 1 1 1 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0;
%   0 1 1 0 1 1 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0;
%   0 1 1 0 0 1 1 1 0 0 1 1 0 0 1 1 1 1 1 1 0;
%   1 1 1 0 0 0 1 1 0 0 1 1 0 0 0 1 1 1 1 0 0;
%   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%   0 0 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 1 0 0;
%   0 0 0 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 0 0;
%   0 0 0 1 1 0 0 0 1 0 0 0 0 1 1 0 0 0 1 0 0;
%   0 0 0 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0;
%   0 0 0 1 1 1 1 0 0 0 0 0 0 1 1 1 1 0 0 0 0;
%   0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0;
%   0 0 0 1 1 0 0 0 1 0 0 0 0 1 1 0 0 0 1 0 0;
%   0 0 0 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 0 0;
%   0 0 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 1 0 0;
%   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
% ];
% aperture(21 : 41, 21 : 41)=pad;

% Blobs
% r = 0.08;
% cx = [-0.3, -0.25, -0.3, -0.35, 0.25, 0.3, 0.25, -0.05, 0, 0.05, 0];           
% cy = [0.25, 0.3, 0.2, 0.25, 0.1, 0.15, 0.05,-0.25, -0.3, -0.25, -0.2];
% D2min = min( (X - reshape(cx,1,1,[])).^2 + (Y - reshape(cy,1,1,[])).^2 , [], 3);
% aperture = double(D2min <= r^2);

% Target
diffraction  = fftshift(fft2(aperture) / sqrt(M * N));
diffPattern = abs(diffraction);

% Functions
g = @(input) gradient_fn(diffPattern, input);
error_fn = @(input) norm(abs(fftshift(fft2(input) / sqrt(M * N))) - diffPattern, 'fro');

%% Run Trials
stepsize = 0.05;
maxIter  = 10000;
beta     = 0.0;
numTrials = 2000;
errorThreshold = 1e-10;
successCount = 0;

allErrors = zeros(numTrials, 1);

tic;

for trial = 1:numTrials
    if mod(trial, 50) == 0
        fprintf('Trial %d/%d\n', trial, numTrials);
    end
    
    % Random initialization
    U = rand(size(aperture));
    grad = zeros(size(U));
    
    % Gradient Descent
    for i = 1:maxIter
        grad = beta * grad + (1 - beta) * g(U);
        U = U - stepsize * grad;
        % U = U .* (U > 0);
        U = min(max(U,0),1);

        % Early stopping
        if mod(i, 100) == 0 && i > 500
            U_temp = U .* (U > 0); U_temp = ones(size(U_temp)) .* (U_temp > 0.5);
            if error_fn(U_temp) < errorThreshold
                break; 
            end
        end
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
fprintf('Algorithm: %s\n', algoName);
fprintf('Success: %d / %d\n', successCount, numTrials);
fprintf('Success rate: %.2f%%\n', (successCount/numTrials)*100);
fprintf('Min error: %.2e\n', min(allErrors));
fprintf('Max error: %.2e\n', max(allErrors));
fprintf('Mean error: %.2e\n', mean(allErrors));
fprintf('Time: %.2f seconds\n', totalTime);

%% Function
function g = gradient_fn(target, input)
    M = size(target, 1);
    N = size(target, 2);
    Z = fftshift(fft2(input)) / sqrt(M * N);
    Y = abs(Z);
    phase = exp(1i * angle(Z));
    g = 2 * real(ifft2(ifftshift(phase .* (Y - target))) * sqrt(M * N));
end