function [vp, vc] = TorusDecompBasis(M, N)
% TorusDecompBasis
% Real-valued symmetry-adapted Fourier basis on Z_M x Z_N
%
% Outputs:
%   vp : cell (1x1 / 2x1 / 4x1), plane (self-conjugate) bases
%   vc : cell (dof x 2), {real part, imag part} of cylinder bases
%
% -------------------------------------------------------------

%% degrees of freedom & preallocation
numPlane = (1 + (mod(M,2)==0)) * (1 + (mod(N,2)==0));
dof      = (M*N - numPlane)/2;
vp       = cell(numPlane,1);
vc       = cell(dof,2);

%% 1D Fourier modes
w = exp(2*pi*1i/M);
z = exp(2*pi*1i/N);

u = cell(M,1);
v = cell(N,1);

for k = 0:M-1
    u{k+1} = w.^(k*(0:M-1));
end
for l = 0:N-1
    v{l+1} = z.^(l*(0:N-1));
end

%% decomposition
idx_p = 0;
idx_c = 0;

for k = 0:M-1
    for l = 0:N-1

        k2 = mod(-k, M);
        l2 = mod(-l, N);

        phi = u{k+1}.' * v{l+1};

        % -------- plane (self-conjugate) --------
        if k==k2 && l==l2
            idx_p = idx_p + 1;
            vp{idx_p} = reshape(real(phi), [M * N, 1]);

        % -------- cylinder (conjugate pairs) --------
        elseif (k < k2) || (k==k2 && l < l2)
            idx_c = idx_c + 1;

            vc{idx_c,1} = reshape(real(phi), [M * N, 1]);
            vc{idx_c,2} = reshape(imag(phi), [M * N, 1]);
        end
    end
end

%% safety check
assert(idx_c == dof, 'DOF mismatch: decomposition inconsistent');

end
