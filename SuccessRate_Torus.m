% Run 1000 trials and compute success rate
% Success: final error < 1e-10
clear; close all; clc;

%% Configuration
numTrials = 100;
errorThreshold = 1e-10;
M = 11;
N = 11;

[X,Y] = meshgrid(linspace(-1,1,N), linspace(-1,1,M));
vx = [-0.43, -0.05, -0.10,  0.13,  0.05,  0.43,  0.18,  0.23, -0.22, -0.15];
vy = [ 0.00,  0.40,  0.05,  0.13, -0.38,  0.00,  0.43, -0.08, -0.3, -0.43];
in = inpolygon(X, Y, vx, vy);
aperture = double(in);

stepsize = 0.1;
maxIter  = 10000;
beta     = 0.0;

% Target
diffraction  = fftshift(fft2(aperture) / sqrt(M * N));
diffPattern  = abs(diffraction);
error_fn = @(input) norm(abs(fftshift(fft2(input) / sqrt(M * N))) - diffPattern, 'fro');

% Generate bases once
[vp, vc, S] = TorusDecompBasis(M, N);
numRotComp = length(S);

%% Run Trials
successCount = 0;
allErrors = zeros(numTrials, 1);

% Stage 1: Initialize
U = ifft2(ifftshift(diffPattern)) * sqrt(M * N);
U = real(U);
if sum(U, 'all') < 0
    U = -U;
end

% Stage 2: Torus Decomposition
U = reshape(U, [M * N, 1]);
proj = @(x, v) (v' * x) / (v' * v) * v;

x0 = 0;
for k = 1 : length(vp)
    x0 = x0 + proj(U, vp{k});
end

coeffRotX = zeros(numRotComp, 2);
for k = 1 : numRotComp
    coeffRotX(k, 1) = U' * vc{k, 1} / norm(vc{k, 1}) ^ 2;
    coeffRotX(k, 2) = U' * vc{k, 2} / norm(vc{k, 2}) ^ 2;
end

torus = @(t) torusRep(t, x0, coeffRotX, vc);

tic;
% Stage 3: Optimize
for trial = 1:numTrials
    fprintf('Trial %d/%d\n', trial, numTrials);
    
    tInitial = 2 * pi * rand(1, numRotComp);
    U = reshape(torus(tInitial), [M, N]);
    t = tInitial;
    grad = 0;
    
    for i = 1:maxIter
        U = reshape(U, [M * N, 1]);
        
        SexpStx = zeros(M * N, numRotComp);
        for k = 1 : numRotComp
            c  = cos(t(k));
            s  = sin(t(k));
            a0 =  coeffRotX(k, 2);
            b0 = -coeffRotX(k, 1);
            a =  a0 * c + b0 * s;
            b = -a0 * s + b0 * c;
            SexpStx(:,k) = a * vc{k, 1} + b * vc{k, 2};
        end
        
        % g    = - (U < 0)' * SexpStx;                             % positivity prior
        g    = 2 * (2 * U .^ 3 - 3 * U .^ 2 + U)' * SexpStx;    % 0,1-Prior
        grad = beta * grad + (1 - beta) * g;
        t    = t - stepsize * grad;
        
        U = torus(t);
        U = reshape(U, [M, N]);
    end
    
    % Stage 5: Finalization
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

%% Functions
function x = torusRep(t, x0, coeffRotX, vc)
    x = x0;
    for k = 1 : length(t)
        c  = cos(t(k));
        s  = sin(t(k));
        coeffk1 =  coeffRotX(k, 1) * c + coeffRotX(k, 2) * s;
        coeffk2 = -coeffRotX(k, 1) * s + coeffRotX(k, 2) * c;
        x = x + coeffk1 * vc{k, 1} + coeffk2 * vc{k, 2};
    end
end
