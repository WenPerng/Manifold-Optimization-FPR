% Projected Gradient Descent
clear; close all; clc;

%% Initialization configurations
% Test Image
M = 21; % size of images
N = 21;

[X,Y] = meshgrid(linspace(-1,1,N), linspace(-1,1,M));

% aperture = (abs(abs(X) - 0.1) < 0.05) .* (abs(Y) < 0.2);

vx = [-0.43, -0.05, -0.10,  0.13,  0.05,  0.43,  0.18,  0.23, -0.22, -0.15];
vy = [ 0.00,  0.40,  0.05,  0.13, -0.38,  0.00,  0.43, -0.08, -0.3, -0.43];
in = inpolygon(X, Y, vx, vy);
aperture = double(in);

% aperture = (X .^ 2 + Y .^ 2) <= 0.2 .^ 2;

diffraction  = fftshift(fft2(aperture) / sqrt(M * N));
diffPattern = abs(diffraction);

% Initialization ----------------------------------------------------------
initialInput = rand(size(aperture));

% Gradient Descent --------------------------------------------------------
g = @(input) gradient(diffPattern, input);

% analysis
error = @(input) norm(abs(fftshift(fft2(input) / sqrt(M * N))) - diffPattern, 'fro');

%% Gradient Descent
stepsize = 0.05;
maxIter  = 10000;
beta     = 0.0;

Err    = zeros(maxIter + 1, 1);
Err(1) = error(initialInput);
grad   = zeros(size(initialInput));

U = initialInput;
tic;
for i = 1 : maxIter
    if mod(i, maxIter / 100) == 0
        disp(['Iteration: ', num2str(i), '/', num2str(maxIter)]);
    end

    grad = beta * grad + (1 - beta) * g(U);
    U = U - stepsize * grad;
    U = ones(size(U)) .* (U >= 1) + U .* (U < 1) .* (U > 0);    
    
    Err(i+1) = error(U);
end
toc;

%% Final Projection =======================================================
U = U .* (U > 0);
U = ones(size(U)) .* (U > 0.5);
% Err = [Err, error(U)];

%% Plot Results
figure;
subplot(2, 3, 1);
imshow(aperture, []);
colorbar;
title('Target Input Amplitude');
axis on;

diffraction  = fftshift(fft2(aperture) / sqrt(M * N));
diffPattern = abs(diffraction);
subplot(2, 3, 2);
imshow(diffPattern, []);
colormap default;
colorbar;
title('Target Output Amplitude');
axis on;

subplot(2, 3, 3);
plot(0 : maxIter, Err, 'DisplayName', ['$\mu = ', num2str(stepsize), '$'], 'LineWidth', 1.8);
title('Total Error', 'Interpreter', 'latex');
grid on;
legend('Interpreter', 'latex');
set(gca, 'YScale', 'log');

subplot(2, 3, 4);
imshow(U, []);
colormap default;
colorbar;
title('Optimized Input Amplitude');
axis on;

diffraction         = fftshift(fft2(U) / sqrt(M * N));
Diff_pattern_output = abs(diffraction);
subplot(2, 3, 5);
imshow(Diff_pattern_output, []);
colormap default;
colorbar;
title('Optimized Output Amplitude');
axis on;

subplot(2, 3, 6);
imshow(log(abs(diffPattern - Diff_pattern_output)), []);
colormap default;
colorbar;
title('Output Amplitude Difference (log)');
axis on;

fontname('Times New Roman');
set(gcf, 'Color', [1, 1, 1]);

%% Function
function g = gradient(target, input)
    M = size(target, 1);
    N = size(target, 2);
    Z = fftshift(fft2(input)) / sqrt(M *N);
    Y = abs(Z);
    
    phase = exp(1i * angle(Z));
    % delta = 1e-5;
    % eps   = 1e-8;
    % phase = (Z ./ (abs(Z) + eps)); % .^ double(Y < delta) + ...
            % (exp(1i * angle(Z))) .^ double(Y >= delta);
    
    g = 2 * real(ifft2(ifftshift(phase .* (Y - target))) * sqrt(M * N));
end
