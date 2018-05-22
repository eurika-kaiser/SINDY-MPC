% from https://github.com/cwrowley/dmdtools/blob/master/python/scripts/total_dmd_example.py
clear all; clc; close all

n_rank = 2;  % True rank of the system
n = 250;     % Number of states
m = 1000;    % Number of snapshots
std = 5e+1;  % standard deviation of the noise


% True system
Alow = diag(exp([1i, 0.65i]));
data = zeros(n_rank, m+1);
data(:,1) = randn(n_rank,1) + 1i*randn(n_rank,1);

for ii=1:m
    data(:,ii+1) = Alow*data(:, ii);
end

Q = qr(randn(n, 2));
data = Q*data; % lift into higher-dimensional space
%data = [real(data), imag(data)];  % Split and stack real and image parts

% Add noise to the data
noisy_data = data + std*randn(size(data,1), size(data,2))+ 1i*std*randn(size(data,1), size(data,2));

% Create a new figure for output
fig = figure(1); hold on, box on
th = linspace(0, 2*pi, 101);
plot(cos(th), sin(th), '-', 'color','k')
plot(real(diag(Alow)), imag(diag(Alow)), 'ko')


% "Standard" DMD
% Note:  n_rank is doubled because we only deal with real numbers
X1 = data(:,1:end-1);
X2 = data(:,2:end);
r = n_rank;
[U, S, V] = svd(X1, 'econ');
U_r = U(:, 1:r); % truncate to rank-r
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);
Atilde = U_r' * X2 * V_r / S_r; % low-rank dynamics
[W_r, D] = eig(Atilde);
dmd_modes = X2 * V_r / S_r * W_r; % DMD modes

dmd_vals = diag(D); % discrete-time eigenvalues

% Plot the DMD eigenvalues
plot(real(dmd_vals), imag(dmd_vals), 'rv')

% "Total least squares" DMD a la M. Hemati
[~,~,V1] = svd([X1;X2]);
V1 = V1(:,1:r);
Xhat = X1*V1;
Yhat = X2*V1;
[U,S,V] = svd(Xhat,0);
Atilde = U'*Yhat*V*pinv(S);
[vtilde,evals] = eig(Atilde);
evals = diag(evals);
modes = U*vtilde;

plot(real(evals), imag(evals), 'bx')
