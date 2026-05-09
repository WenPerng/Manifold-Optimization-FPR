% hessian_0131.m
clear; close all; clc;

% Load saved data
load('local_min_data.mat');
M = size(U, 1);
N = size(U, 2);

disp('=== Starting Full Hessian Eigenvalue Analysis ===');
disp(['Matrix size: ', num2str(M), ' x ', num2str(N), ' = ', num2str(M*N), ' variables']);
disp(['Hessian will be ', num2str(M*N), ' x ', num2str(M*N)]);
disp(' ');

% gradient
g  = @(input) gradient(input);

% torusdecomp U
[vp, vc] = TorusDecompBasis(M, N);
numRotComp = size(vc, 1);
U = reshape(U, [M * N, 1]);

coeffRotX = zeros(numRotComp, 2);
for k = 1 : numRotComp
    coeffRotX(k, 1) = U' * vc{k, 1} / norm(vc{k, 1}) ^ 2;
    coeffRotX(k, 2) = U' * vc{k, 2} / norm(vc{k, 2}) ^ 2;
end

% first derivative
SexpStx = zeros(M * N, numRotComp);
for k = 1 : numRotComp
    SexpStx(:, k) = -coeffRotX(k, 2) * vc{k, 1} + coeffRotX(k, 1) * vc{k, 2};
end

% Second derivative
SSexpStx = zeros(M * N, numRotComp);
for k = 1 : numRotComp
    SSexpStx(:, k) = -coeffRotX(k, 1) * vc{k, 1} - coeffRotX(k, 2) * vc{k, 2};
end

S = SexpStx;
D = SSexpStx;
grad = reshape(grad, [M * N, 1]);

% compute hessian
H = S' * diag( 12 * U .^ 2 - 12 * U + 2 ) * S + diag(D' * grad);
e = eig(H);
tol = 1e-10;

if all(e > tol)
    disp('global or local minimum');
else
    disp('saddle point');
end

% show image
subplot(1, 2, 1);
histogram(e, 1000);
grid on;

subplot(1, 2, 2);
U = reshape(U, [M , N]);
imshow(U, []); 
colormap default; 
colorbar; 
title('Optimized Input Amplitude'); 
axis on;

% function
function g = gradient(input)
    U = input;
    g = 2 * (2 * U .^ 3 - 3 * U .^ 2 + U);
end