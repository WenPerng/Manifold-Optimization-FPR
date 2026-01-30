% Dual Projected Gradient Descent
clear; close all; clc;

%% Configurations =========================================================
% Test Image
M = 21; % size of images
N = 21;

[X,Y] = meshgrid(linspace(-1,1,N), linspace(-1,1,M));

% aperture = rand(M, N);

% aperture = (abs(abs(X) - 0.1) < 0.05) .* (abs(Y) < 0.2);

% aperture = (abs(abs(X) - 0.05) < 0.01) .* (abs(abs(Y) - 0.05) < 0.01);

vx = [-0.43, -0.05, -0.10,  0.13,  0.05,  0.43,  0.18,  0.23, -0.22, -0.15];
vy = [ 0.00,  0.40,  0.05,  0.13, -0.38,  0.00,  0.43, -0.08, -0.3, -0.43];
in = inpolygon(X, Y, vx, vy);
aperture = double(in);

% aperture = ((abs(X) - 0.3) .^ 2 + Y .^ 2) <= 0.2 .^ 2;

% load('imgBadApple.mat');
% aperture = zeros(M, N);
% aperture(1:90,1:120) = imgBadApple;

% Gradient Descent Iteration ----------------------------------------------
stepsize = 0.05;
maxIter  = 2000;
beta     = 0.0;

%% Initializations ========================================================
% Target ------------------------------------------------------------------
diffraction  = fftshift(fft2(aperture) / sqrt(M * N));
diffPattern = abs(diffraction);

% Analysis helper functions -----------------------------------------------
error   = @(input) norm(abs(fftshift(fft2(input) / sqrt(M * N))) - diffPattern, 'fro');
sumNeg  = @(input) abs(sum(sum(input .* (input < 0))));
dist201 = @(input) norm(input .* (input - 1), 'fro');

Err      = [];
negVal   = [];
distTo01 = [];
gradSize = [];

% Gradient descent helper functions ---------------------------------------
g  = @(input) gradient(input);
PA = @(y) proja(y, M, N, diffPattern);

%% Gradient Descent
initialInput = PA(rand(M, N));

Err      = [Err,      error(initialInput)];
negVal   = [negVal,   sumNeg(initialInput)];
distTo01 = [distTo01, dist201(initialInput)];

U    = initialInput;
grad = 0;
tic;
for i = 1 : maxIter
    if mod(i, maxIter / 100) == 0
        disp(['Iteration: ', num2str(i), '/', num2str(maxIter)]);
    end
    
    grad = beta * grad + (1 - beta) * g(U);
    U = U - stepsize * grad;

    U = PA(U);
    
    Err      = [Err,      error(U)];
    negVal   = [negVal,   sumNeg(U)];
    distTo01 = [distTo01, dist201(U)];
    gradSize = [gradSize, norm(grad)];
end
toc;

%% Final Projection =======================================================
U = U .* (U > 0);
U = ones(size(U)) .* (U > 0.5);
Err = [Err, error(U)];

%% Plot Results ===========================================================
subplot(2, 3, 1);
imshow(aperture, []);
colorbar;
title('Target Signal');
axis on;

diffraction_plot  = fftshift(fft2(aperture) / sqrt(M * N));
diffPattern_plot  = abs(diffraction_plot);

subplot(2, 3, 2);
imshow(diffPattern_plot, []);
colormap default;
colorbar;
title('Target Fourier Amplitude');
axis on;

subplot(2, 3, 3);
yyaxis left;
plot(0 : length(distTo01) - 1, distTo01, 'DisplayName', 'Distance to 0,1', 'LineWidth', 1.8);
grid on;
set(gca, 'YScale', 'log');
yyaxis right;
plot(0 : length(negVal) - 1, negVal, 'DisplayName', 'Sum of Negative Components', 'LineWidth', 1.8);

set(gca, 'YScale', 'linear');
legend('Interpreter', 'latex');

subplot(2, 3, 4);
imshow(U, []);
colormap default;
colorbar;
title('Retrieved Signal');
axis on;

diffraction_out        = fftshift(fft2(U) / sqrt(M * N));
Diff_pattern_output    = abs(diffraction_out);

subplot(2, 3, 5);
imshow(Diff_pattern_output, []);
colormap default;
colorbar;
title('Optimized Fourier Amplitude');
subtitle(['Total Error: ', num2str(Err(end))]);
axis on;

subplot(2, 3, 6);
plot(1 : length(gradSize), gradSize, 'DisplayName', 'Gradient Size', 'LineWidth', 1.8);
grid on;
set(gca, 'YScale', 'log');
legend('Interpreter', 'latex');
axis on;

fontname('Times New Roman');
set(gcf, 'Color', [1, 1, 1]);
% set(gcf, 'Position', get(0, 'Screensize'));

%% Functions ==============================================================
function g = gradient(input)
    U = input;
    g = 2 * (2 * U .^ 3 - 3 * U .^ 2 + U);
end

function y = proja(y, M, N, diffPattern)
    x = fftshift(fft2(y) / sqrt(M * N));
    x = diffPattern .* exp(1i * angle(x));
    y = real(ifft2(ifftshift(x) * sqrt(M * N)));
end