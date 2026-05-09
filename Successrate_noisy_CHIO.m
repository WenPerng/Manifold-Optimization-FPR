% CHIO with support
% Run 1000 trials and compute noisy retrieval rate
% Success: final binary pixel error <= 0.1% up to trivial associates
clear; close all; clc;
rng(0);
algoName = "CHIO (noisy)";

%% Configuration
M = 101;
N = 101;

[X,Y] = meshgrid(linspace(-1,1,N), linspace(-1,1,M));

% Double slit
% aperture = (abs(abs(X) - 0.1) < 0.05) .* (abs(Y) < 0.2);

% Polygon
% vx = [-0.30, -0.05, -0.10,  0.13,  0.05,  0.30,  0.18,  0.23, -0.22, -0.15];
% vy = [ 0.00,  0.30,  0.05,  0.13, -0.21,  0.00,  0.35, -0.08, -0.3, -0.25];
% in = inpolygon(X, Y, vx, vy);
% aperture = double(in);

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
% aperture(21 : 41, 21 : 41) = pad;

% 61x61 ICIP
% aperture = zeros(M, N);
% pad = [
%   1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 0 0;
%   1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 0;
%   0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 0 0 0 1 1 1;
%   0 0 0 0 1 1 0 0 0 0 0 1 1 1 0 0 0 0 1 1 1;
%   0 0 0 0 1 1 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0;
%   0 0 0 0 1 1 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0;
%   0 0 0 0 1 1 0 0 0 0 0 1 1 1 0 0 0 0 1 1 1;
%   0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 0 0 0 1 1 1;
%   1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 0;
%   1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1 0 0;
%   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
%   1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 0 0;
%   1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 0;
%   0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 1 1 1;
%   0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 1 1;
%   0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 1 1 1;
%   0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0;
%   0 0 0 0 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0;
%   0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0;
%   1 1 1 1 1 1 1 1 1 1 0 0 1 1 0 0 0 0 0 0 0;
%   1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 0 0 0 0 0 0;
% ];
% aperture(21 : 41, 21 : 41) = pad;

% BadApple
load('imgBadApple101.mat');
aperture = zeros(M, N);
aperture(33 : 66, 33 : 66) = imgBadApple;
type = 'BadApple';

% Blobs
% r = 0.08;
% cx = [-0.3, -0.25, -0.3, -0.35, 0.25, 0.3, 0.25, -0.05, 0, 0.05, 0];
% cy = [0.25, 0.3, 0.2, 0.25, 0.1, 0.15, 0.05, -0.25, -0.3, -0.25, -0.2];
% D2min = min((X - reshape(cx,1,1,[])).^2 + (Y - reshape(cy,1,1,[])).^2, [], 3);
% aperture = double(D2min <= r^2);

% Target
diffraction_clean = fftshift(fft2(aperture) / sqrt(M * N));
diffPattern_clean = abs(diffraction_clean);
cleanNorm = norm(diffPattern_clean, 'fro');
noiseLevel = 0.25;

% Trivial-associate evaluation
aperture_binary = double(aperture > 0.5);
fft_aperture = fft2(aperture_binary);
fft_aperture_ref = fft2(rot90(aperture_binary, 2));

% Functions
pixel_error_fn = @(input) best_assoc_pixel_error_fft( ...
    input, aperture_binary, fft_aperture, fft_aperture_ref, M, N);

%% Run Trials
maxIter = 10000;
numTrials = 1000;
pixelErrorThreshold = 0.001;
beta  = 0.7;
alpha = 0.4;
S = true(M, N);   % support
successCount = 0;

allPixelErrors = zeros(numTrials, 1);
allDataResiduals = zeros(numTrials, 1);
U0_all = rand(M, N, numTrials);

tic;

for trial = 1:numTrials
    if mod(trial, 50) == 0
        fprintf('Trial %d/%d\n', trial, numTrials);
    end

    % New noisy observation for each trial
    noise = randn(M, N);
    noise = noise / norm(noise, 'fro') * (noiseLevel * cleanNorm);
    diffPattern_obs = max(diffPattern_clean + noise, 0);

    PB = @(y) projb(y, diffPattern_obs);
    data_residual_fn = @(input) ...
        norm(abs(fftshift(fft2(input) / sqrt(M * N))) - diffPattern_obs, 'fro') ...
        / norm(diffPattern_obs, 'fro');

    % Stage 1: Random initialization
    U = U0_all(:, :, trial);

    % Stage 2: CHIO iteration
    for i = 1:maxIter
        U_old = U;

        % Fourier magnitude projection
        A = fftshift(fft2(U_old) / sqrt(M * N));
        A_proj = PB(A);
        U_prime = real(ifft2(ifftshift(A_proj) * sqrt(M * N)));

        % CHIO update
        U = U_old;

        idx1 = S & (U_prime >= alpha * U_old);
        U(idx1) = U_prime(idx1);

        idx2 = S & (U_prime >= 0) & (U_prime <= alpha * U_old);
        U(idx2) = U_old(idx2) - ((1 - alpha) / alpha) * U_prime(idx2);

        idx3 = ~(idx1 | idx2);
        U(idx3) = U_old(idx3) - beta * U_prime(idx3);

        U = min(max(U, 0), 1);

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
fprintf('Success rate: %.2f%%\n', (successCount / numTrials) * 100);
fprintf('Mean pixel error: %.6f\n', mean(allPixelErrors));
fprintf('Std pixel error: %.6f\n', std(allPixelErrors));
fprintf('Min pixel error: %.6f\n', min(allPixelErrors));
fprintf('Max pixel error: %.6f\n', max(allPixelErrors));
fprintf('Mean data residual: %.6e\n', mean(allDataResiduals));
fprintf('Time: %.2f seconds\n', totalTime);

%% Functions
function y = projb(y, diffPattern)
    y = diffPattern .* exp(1i * angle(y));
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