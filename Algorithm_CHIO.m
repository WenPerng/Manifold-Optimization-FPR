% Optimizing for M x N case
%--------------------------------------------------------------------------
% CHIO with support
clear; close all; clc;

%% Configurations
M = 271;
N = 271;

[X,Y] = meshgrid(linspace(-1,1,N), linspace(-1,1,M));

% Double slit
% aperture = (abs(abs(X) - 0.1) < 0.05) .* (abs(Y) < 0.2);

% Polygon
% vx = [-0.43, -0.05, -0.10,  0.13,  0.05,  0.43,  0.18,  0.23, -0.22, -0.15];
% vy = [ 0.00,  0.40,  0.05,  0.13, -0.38,  0.00,  0.43, -0.08, -0.3, -0.43];
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
% aperture(21 : 41, 21 : 41)=pad;

% Blobs
% r = 0.08;
% cx = [-0.3, -0.25, -0.3, -0.35, 0.25, 0.3, 0.25, -0.05, 0, 0.05, 0];           
% cy = [0.25, 0.3, 0.2, 0.25, 0.1, 0.15, 0.05,-0.25, -0.3, -0.25, -0.2];
% D2min = min( (X - reshape(cx,1,1,[])).^2 + (Y - reshape(cy,1,1,[])).^2 , [], 3);
% aperture = double(D2min <= r^2);

load('imgBadApple.mat');
aperture = zeros(M, N);
aperture(1:90,1:120) = imgBadApple;

maxIter = 10000;

%% Initialization =========================================================
% Target ------------------------------------------------------------------
diffraction  = fftshift(fft2(aperture) / sqrt(M * N));
diffPattern  = abs(diffraction);

% Analysis helper functions -----------------------------------------------
error   = @(input) norm(abs(fftshift(fft2(input) / sqrt(M * N))) - diffPattern, 'fro');
Err = zeros(1, maxIter);

PB = @(y) projb(y, diffPattern);

U = rand(M , N); % spatial domain
beta  = 0.7;
alpha = 0.4;

S = true(M, N);  % support

for i = 1:maxIter
    if mod(i, maxIter / 100) == 0
        disp(['Iteration: ', num2str(i), '/', num2str(maxIter)]);
    end

    U_old = U;

    % Fourier magnitude projection
    A      = fftshift(fft2(U_old) / sqrt(M*N));
    A_proj = PB(A);
    U_prime    = real(ifft2(ifftshift(A_proj) * sqrt(M*N)));

    % -------- CHIO update ----------
    U = U_old;

    idx1 = S & (U_prime >= alpha * U_old);
    U(idx1) = U_prime(idx1);

    idx2 = S & (U_prime >= 0) & (U_prime <= alpha * U_old);
    U(idx2) = U_old(idx2) - ((1 - alpha)/alpha) * U_prime(idx2);

    idx3 = ~(idx1 | idx2);
    U(idx3) = U_old(idx3) - beta * U_prime(idx3);
    % -------------------------------
    
    U = min(max(U,0),1);
    Err(i) = error(U);
end

%%
U = U .* (U > 0);
U = ones(size(U)) .* (U > 0.5);
Err = [Err, error(U)];

%%
% Plot convergence
subplot(2, 3, 1);
imshow(aperture, []);
colorbar;
title('Target Input Amplitude');
axis on;

diffraction_plot  = fftshift(fft2(aperture) / sqrt(M * N));
diffPattern_plot  = abs(diffraction_plot);

subplot(2, 3, 2);
imshow(diffPattern_plot, []);
colormap default;
colorbar;
title('Target Output Amplitude');
axis on;

subplot(2, 3, 3);
yyaxis left;
plot(0 : length(Err) - 1, Err, 'DisplayName', 'Err', 'LineWidth', 1.8);
grid on;
set(gca, 'YScale', 'log');
legend('Interpreter', 'latex');

subplot(2, 3, 4);
imshow(U, []);
colormap default;
colorbar;
title('Optimized Input Amplitude');
axis on;

a = fftshift(fft2(U)) / sqrt(M * N);

subplot(2, 3, 5);
imshow(abs(a), []);
colormap default;
colorbar;
title('Optimized Output Amplitude');
subtitle(['Total Error: ', num2str(Err(end))]);
axis on;

fontname('Times New Roman');
set(gcf, 'Color', [1, 1, 1]);

function y = projb(y, diffPattern)
    y = diffPattern .* exp(1i * angle(y));
end