% TorusGD
% Run 1000 trials and compute noisy retrieval rate
% Success: final binary pixel error <= 1% up to trivial associates
clear; close all; clc;
rng(0);
algoName = "Torus-parameter GD on TorusDecompBasis (optimize angles t, noisy)";

%% Configuration
M = 101;
N = 101;

[X,Y] = meshgrid(linspace(-1,1,N), linspace(-1,1,M));

% Double slit
% aperture = (abs(abs(X) - 0.1) < 0.05) .* (abs(Y) < 0.2);

% Polygon
vx = [-0.30, -0.05, -0.10,  0.13,  0.05,  0.30,  0.18,  0.23, -0.22, -0.15];
vy = [ 0.00,  0.30,  0.05,  0.13, -0.21,  0.00,  0.35, -0.08, -0.3, -0.25];
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

% BadApple
% load('imgBadApple101.mat');
% aperture = zeros(M, N);
% aperture(33 : 66, 33 : 66) = imgBadApple;
% type = 'BadApple';

% Blobs
% r = 0.08;
% cx = [-0.3, -0.25, -0.3, -0.35, 0.25, 0.3, 0.25, -0.05, 0, 0.05, 0];           
% cy = [0.25, 0.3, 0.2, 0.25, 0.1, 0.15, 0.05,-0.25, -0.3, -0.25, -0.2];
% D2min = min( (X - reshape(cx,1,1,[])).^2 + (Y - reshape(cy,1,1,[])).^2 , [], 3);
% aperture = double(D2min <= r^2);

% Target
diffraction_clean  = fftshift(fft2(aperture) / sqrt(M * N));
diffPattern_clean  = abs(diffraction_clean);
cleanNorm = norm(diffPattern_clean, 'fro');
noiseLevel = 0.10;

% Trivial-associate evaluation
aperture_binary = double(aperture > 0.5);
fft_aperture = fft2(aperture_binary);
fft_aperture_ref = fft2(rot90(aperture_binary, 2));

% Functions
pixel_error_fn = @(input) best_assoc_pixel_error_fft(input, aperture_binary, fft_aperture, fft_aperture_ref, M, N);

% Generate bases once
[vp, vc] = TorusDecompBasis(M, N);
numRotComp = size(vc, 1);

%% Run Trials
stepsize = 0.05;
maxIter  = 10000;
beta     = 0.0;
numTrials = 1000;
pixelErrorThreshold = 0.001;
successCount = 0;

allPixelErrors = zeros(numTrials, 1);
allDataResiduals = zeros(numTrials, 1);
t0_all = 2 * pi * rand(numTrials, numRotComp);

tic;

for trial = 1:numTrials
    if mod(trial, 50) == 0
        fprintf('Trial %d/%d\n', trial, numTrials);
    end

    noise = randn(M, N);
    noise = noise / norm(noise, 'fro') * (noiseLevel * cleanNorm);
    diffPattern_obs = max(diffPattern_clean + noise, 0);

    data_residual_fn = @(input) ...
        norm(abs(fftshift(fft2(input) / sqrt(M * N))) - diffPattern_obs, 'fro') ...
        / norm(diffPattern_obs, 'fro');

    % Stage 1: Initialize from the noisy observation
    U = ifft2(ifftshift(diffPattern_obs)) * sqrt(M * N);
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

    % Stage 3: Optimize
    tInitial = t0_all(trial, :);
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

        g = 2 * (2 * U .^ 3 - 3 * U .^ 2 + U)' * SexpStx;
        % grad = beta * grad + (1 - beta) * g;
        t = t - stepsize * g;

        U = torus(t);
        U = reshape(U, [M, N]);

        % Early stopping
        if mod(i, 100) == 0 && i > 500
            U_temp = U .* (U > 0);
            U_temp = ones(size(U_temp)) .* (U_temp > 0.5);
            if pixel_error_fn(U_temp) <= pixelErrorThreshold
                break;
            end
        end
    end

    % Finalization
    U = U .* (U > 0);
    U = ones(size(U)) .* (U > 0.5);

    % Check success
    finalPixelError = pixel_error_fn(U);
    finalDataResidual = data_residual_fn(U);

    allPixelErrors(trial) = finalPixelError;
    allDataResiduals(trial) = finalDataResidual;

    if finalPixelError <= pixelErrorThreshold
        successCount = successCount + 1;
    end
end

totalTime = toc;

%% Results
fprintf('\n=== RESULTS ===\n');
fprintf('Algorithm: %s\n', algoName);
fprintf('Noise level: %.2f%%\n', noiseLevel * 100);
fprintf('Success: %d / %d\n', successCount, numTrials);
fprintf('Success rate: %.2f%%\n', (successCount/numTrials)*100);
fprintf('Mean pixel error: %.6f\n', mean(allPixelErrors));
fprintf('Std pixel error: %.6f\n', std(allPixelErrors));
fprintf('Min pixel error: %.6f\n', min(allPixelErrors));
fprintf('Max pixel error: %.6f\n', max(allPixelErrors));
fprintf('Mean data residual: %.6e\n', mean(allDataResiduals));
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

function err = best_assoc_pixel_error_fft(U, aperture, fft_ap, fft_ap_ref, M, N)
    U = double(U > 0.5);
    totalPixels = M * N;

    xcorr_normal = real(ifft2(fft2(U) .* conj(fft_ap)));
    xcorr_ref = real(ifft2(fft2(U) .* conj(fft_ap_ref)));

    sumU = sum(U(:));
    sumAp = sum(aperture(:));

    bestOverlap = max([max(xcorr_normal(:)), max(xcorr_ref(:))]);
    bestMismatch = sumU + sumAp - 2 * round(bestOverlap);
    bestMismatch = max(bestMismatch, 0);

    err = bestMismatch / totalPixels;
end