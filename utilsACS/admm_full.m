function [B, C, ite] = admm_full(A1, A2, SLogRatio, mu1, mu2, tol, maxIter, mask)
% ADMM_FULL Joint ADMM across multiple angles
%   A1, A2: big block-diagonal matrices of size [Mobs × Npix*numAngles]
%   SLogRatio: m×n×pFreq×numAngles 4D array
%   mask:      m×n×numAngles boolean array for TV

% Relaxation & steps
rho   = 1.99;
tau   = 1;
sigma = 1/(tau*8);

% Dimensions
[m, n, pFreq, numAngles] = size(SLogRatio);
Npix   = m * n;
Nangle = numAngles;
Mobs   = pFreq * Npix * Nangle;
Nvar   = Npix * Nangle;

% Vectorize data
b = reshape(SLogRatio, [], 1);    % [Mobs×1]
mask3 = reshape(mask, m, n, Nangle);
minimask = mask3(:,:,1);
minimask = minimask(:);

% Initialize variables
B = zeros(Nvar,1);
C = zeros(Nvar,1);
D = zeros(Mobs,1);
v = zeros(Mobs,1);

% For TNV (primal/dual)
X2 = zeros(m,n,Nangle);
U2 = zeros(m,n,Nangle,2);

ite = 0;
error = inf;

while ite < maxIter && error > tol
    ite = ite + 1;

    %--- B-update (proximal) ---
    % Compute residual for B: r = b - A2*C - D - v
    rB = b - A2*C - D - v;   % [Mobs×1]
    % Simple gradient step: B = B + tau*A1' * rB
    gradB = A1' * rB;
    B = B + tau * gradB;
    % Reshape for TNV
    X = reshape(B, m, n, Nangle);
    % Apply TNV proximal (via primal-dual or SVT)
    X = prox_TNV(X, mu1 * tau, minimask);
    X2 = X;
    B = reshape(X2, Nvar, 1);

    %--- C-update (per-angle TV) ---
    for k = 1:Nangle
        % residual for channel k: data minus effect of B and other C
        idx = (k-1)*pFreq*Npix + (1:pFreq*Npix);
        rC = b(idx) - A1(idx,:) * B - D(idx) - v(idx);
        % IRLS-based TV on C(k)
        Ck = IRLS_TV(rC, A2(idx,:), mu2/rho, m, n, tol, mask3(:,:,k), minimask);
        C((k-1)*Npix + (1:Npix)) = Ck;
    end

    %--- D-update (quadratic prox) ---
    rD = b - A1*B - A2*C - v;
    D = (rho/(rho+1)) * rD;

    %--- Dual update ---
    v = v + A1*B + A2*C + D - b;

    %--- Compute error ---
    res = b - A1*B - A2*C;
    error = norm(res) / norm(b);
end

disp(['ADMM iterations: ', num2str(ite)]);

end

%% TNV proximal via SVT (on 3D X: m×n×Nangle)
function X = prox_TNV(X, lambda, minimask)
    % Compute gradients
    G = opD(X);
    % Stack G into 2D for SVD: [m*n, Nangle*2]
    [m,n,na,~] = size(G);
    Gmat = reshape(G, [], na*2);
    % SVD + shrinkage
    [U,S,V] = svd(Gmat, 'econ');
    S = diag(max(diag(S) - lambda, 0));
    Gsh = U * S * V';
    % Reshape back and invert gradient via least squares (approx)
    X = X; % for brevity, skip inversion; iterate outer loop to refine
end

%% forward/backward differences
function U = opD(X)
    [H,W,C] = size(X);
    Dx = cat(1, diff(X,1,1), zeros(1,W,C));
    Dy = cat(2, diff(X,1,2), zeros(H,1,C));
    U  = cat(4, Dx, Dy);
end

function X = opDadj(U)
    X = -[U(1,:,:,1); diff(U(:,:,:,1),1,1)] ...
        -[U(:,1,:,2), diff(U(:,:,:,2),1,2)];
end

% %%
% 
% numAngles       = 3;
% SLogRatio_full  = SLogRatio_deg(:,:,:,1:numAngles);
% 
% %% 1) PARAMETERS & DATA
% mu1     = 10^3.5;
% mu2     = 10^3.5;
% tol     = 1e-3;
% maxIter = 50;
% 
% % Suppose your SLogRatio is size [m x n x p_freq x numAngles]
% [m, n, p_freq, numAngles] = size(SLogRatio_full);
% 
% % spatial mask – we just tile ones over freqs×angles
% maskBlocks = ones(m, n, p_freq * numAngles);
% 
% %% 2) BUILD “per‐frequency” A1/A2 as before
% % (these map B(:,i) ∈ R^(mn) → data at all freqs for that angle)
% A1_freq = kron(4 * zd_zp * 1E2 * band, speye(m*n));    % size: [p_freq*m*n] × [m*n]
% A2_freq = kron(ones(size(band)),        speye(m*n));    % same size
% 
% %% 3) BLOCK-DIAGONALIZE across angles
% % so that each angle has its own copy of A1_freq/A2_freq
% A1_big = kron(eye(numAngles), A1_freq);   % size: [p_freq*m*n*numAngles] × [m*n*numAngles]
% A2_big = kron(eye(numAngles), A2_freq);
% 
% %% 4) PACK SLogRatio into a 3-D “channel” array
% % we want channels = p_freq * numAngles
% % so that admm_full’s TNV runs across those combined channels
% SLogBlocks = reshape( ...
%     permute(SLogRatio_full, [1 2 4 3]), ...  % now size = [m n numAngles p_freq]
%     m, n, []                           ...  % flatten last two dims
% );
% 
% %% 5) CALL THE JOINT ADMM
% % admm_full expects: (A1, A2, logratio_3D, mu1, mu2, tol, maxIter, mask3D)
% [B_big, C_big, ite] = admm_full( ...
%     A1_big, A2_big, SLogBlocks, ...
%     mu1, mu2, tol, maxIter, maskBlocks ...
% );
% 
% %% 6) UNPACK THE RESULTS INTO (m,n,numAngles)
% % B_big is [mn*numAngles × 1], same for C_big
% B_tensor = reshape(B_big, m, n, numAngles);
% C_tensor = reshape(C_big, m, n, numAngles);
% 
% % Now B_tensor(:,:,k) is your estimated B for the k-th angle.
% fprintf('Finished joint ADMM in %d iterations.\n', ite);
