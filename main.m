%MAIN  Communication Systems - MATLAB homework 3: PLL with a limiting phase
%      detector, random processes through LTI filters, and ergodicity.

clear; close all; clc;
addpath(fileparts(mfilename('fullpath')));
rng(3);

%% Q1  Second-order PLL whose phase detector saturates at +-A
% e_d = sin(psi) clipped to [-A, A], loop filter Kp + Ki/s, VCO d(theta)/dt = K_v v.
fs = 20e3;  T = 1.0;  t = 0:1/fs:T-1/fs;
df = 15;                                         % input frequency step (Hz) at t = 0.1 s
theta_in = 2*pi*df*max(t - 0.1, 0);              % phase of the input relative to the free-running VCO
wn = 2*pi*8;  zeta = 0.707;                      % loop natural frequency and damping
Kp = 2*zeta*wn;  Ki = wn^2;                      % with the detector gain = 1 near lock
As = [0.2 0.35 0.5 0.9 1.0];
fig = figure('Visible', 'off', 'Position', [100 100 1100 420]);
lock = nan(size(As));
for k = 1:numel(As)
    A = As(k);
    th = 0;  integ = 0;  psi = zeros(size(t));
    for i = 1:numel(t)
        psi(i) = theta_in(i) - th;
        ed = min(max(sin(psi(i)), -A), A);       % limiting phase detector
        integ = integ + Ki*ed/fs;
        th = th + (Kp*ed + integ)/fs;            % VCO integrates its control signal
    end
    err = abs(mod(psi + pi, 2*pi) - pi);
    settled = find(err > 0.05, 1, 'last');
    if isempty(settled) || settled < numel(t), lock(k) = t(min(settled + 1, numel(t))) - 0.1; end
    if err(end) > 0.05, lock(k) = NaN; end
    subplot(1, 2, 1); plot(t, psi/(2*pi), 'LineWidth', 1.2); hold on;
end
subplot(1, 2, 1); grid on; xlabel('time (s)'); ylabel('phase error \psi / 2\pi (cycles)');
title(sprintf('Response to a %d Hz frequency step', df));
legend(arrayfun(@(A) sprintf('A = %.2f', A), As, 'UniformOutput', false), 'Location', 'northwest');
subplot(1, 2, 2); plot(As, lock, 'o-', 'LineWidth', 1.5); grid on;
xlabel('saturation level A'); ylabel('time to lock (s)'); title('Lock time grows as A shrinks');
save_figure(fig, 'q1_pll_lock');
fprintf('Q1  lock time vs. A: %s\n', strjoin(arrayfun(@(A, L) sprintf('A=%.2f: %.3f s', A, L), As, lock, 'UniformOutput', false), ', '));

%% Q2  Random process through h(t) = e^{-t} cos(28 pi t) u(t)
dt = 1e-3;  t2 = -20:dt:20-dt;  N = numel(t2);
f = (-N/2:N/2-1)/(N*dt);
H = 0.5*(1./(1 + 1j*2*pi*(f - 14)) + 1./(1 + 1j*2*pi*(f + 14)));   % exact transfer function
h = exp(-t2).*cos(28*pi*t2).*(t2 >= 0);
H_fft = fftshift(fft(ifftshift(h))) * dt;                            % numerical check
fprintf('Q2  max |H_fft - H_exact| = %.1e\n', max(abs(H_fft - H)));
inputs = {ones(size(f)), 'white, S_x = 1';  max(1 - abs(f/5), 0), 'colored, S_x = \Lambda(f/5)'};
fig = figure('Visible', 'off', 'Position', [100 100 1100 600]);
for k = 1:2
    Sx = inputs{k, 1};
    Sy = Sx .* abs(H).^2;
    Ry = real(fftshift(ifft(ifftshift(Sy)))) / dt;                   % R_y = IFT{S_y}
    if k == 1
        w0 = 28*pi;                                                      % exact h(t) * h(-t)
        Ry_exact = exp(-abs(t2))/4 .* (cos(w0*t2) + (cos(w0*t2) - w0*sin(w0*abs(t2)))/(1 + w0^2));
        fprintf('    white input: max |R_y - exact| = %.1e, R_y(0) = %.4f (= 1/4 + ...)\n', ...
                max(abs(Ry - Ry_exact)), Ry(N/2+1));
    end
    subplot(2, 2, 2*k-1); plot(f, Sy); xlim([-30 30]); grid on; xlabel('f (Hz)');
    title(['S_y(f) = S_x(f) |H(f)|^2, ' inputs{k, 2}]);
    subplot(2, 2, 2*k); plot(t2, Ry); xlim([-4 4]); grid on; xlabel('\tau (s)');
    title('R_y(\tau)');
end
save_figure(fig, 'q2_psd_autocorrelation');

%% Q3  X(t) = A cos(2 pi f0 t + theta): time averages vs. ensemble averages
A = 1;  f0 = 1;  t3 = linspace(0, 10, 10001);
for c = 1:2
    if c == 1, th = 2*pi*rand(20, 1);       range = '[0, 2\pi)';
    else,      th = -pi/4 + pi/2*rand(20, 1); range = '[-\pi/4, \pi/4]';
    end
    X = A*cos(2*pi*f0*t3 + th);
    time_mean = mean(mean(X, 2));  time_power = mean(mean(X.^2, 2));
    ens_mean_t = mean(X, 1);                                   % sample estimate across realisations
    if c == 1
        ens_mean_theory = zeros(size(t3));
    else
        ens_mean_theory = A*(2*sqrt(2)/pi)*cos(2*pi*f0*t3);    % E[cos(wt + theta)], theta ~ U[-pi/4, pi/4]
    end
    fprintf('Q3  theta ~ U%s: time averages <X> = %.3f, <X^2> = %.3f; ensemble mean %s\n', ...
            erase(range, '\'), time_mean, time_power, ...
            ternary(c == 1, '= 0 for all t -> ergodic in mean and power', ...
                    sprintf('= %.2f cos(2 pi t) (time-varying) -> not stationary, not ergodic', 2*sqrt(2)/pi)));
    fig = figure('Visible', 'off', 'Position', [100 100 1000 330]);
    subplot(1, 2, 1); plot(t3, X, 'Color', [0.2 0.4 0.8 0.3]); xlim([0 3]); grid on;
    title(['20 sample functions, \theta ~ U' range]); xlabel('t');
    subplot(1, 2, 2); plot(t3, ens_mean_t, t3, ens_mean_theory, 'k--', 'LineWidth', 1.2); xlim([0 3]); grid on;
    legend('ensemble average of 20 realisations', 'theory'); xlabel('t'); title('E[X(t)]');
    save_figure(fig, sprintf('q3_random_phase_%d', c));
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
