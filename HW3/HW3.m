clear all

% Define the time axis with increased resolution
t = -10:0.01:10;
% Define the frequency axis with increased resolution
f = -2:0.01:4;
% Generate the signal x
x = exp(1i*t.^2/10 - 1i*3*t).*((t>=-9)&(t<=1)) + exp(1i*t.^2/2 + 1i*6*t).*exp(-(t-4).^2/10);

tic                 % Timer Start
y = wdf(x, t, f);   % Call the wdf function to compute and plot WDF
toc                 % Timer End

function y = wdf(x, t, f)
    dt = mean(diff(t));
    df = mean(diff(f));
    
    % Calculate the number of samples in time and frequency
    N = floor(1 / (2 * dt * df));
    % Define indices for time and frequency
    n1 = floor(t(1) / dt);
    n2 = ceil(t(end) / dt);
    m1 = floor(f(1) / df);
    m2 = ceil(f(end) / df);
    m = mod((m1:m2) - 1, N) + 1;
    
    % Compute the dimensions of the WDF matrix
    Lt = ceil(t(end) / dt) - floor(t(1) / dt) + 1;
    Lf = ceil(f(end) / df) - floor(f(1) / df) + 1;
    y = zeros(Lf, Lt);
    
    % Loop over time indices
    for n = n1:n2
        % Calculate the lag parameter U
        U = max(0, min(n2 - n, n - n1));

        % Calculate the Q parameter
        Q = 2 * U;
        
        % Compute the auto-correlation function A
        A = x(1 - n1 + [n - U:n + U]) .* (x(1 - n1 + [n + U:-1:n - U])').';
        % Compute the FFT of A
        A1 = fft(A, N) * 2 * dt;
        % Calculate a1 parameter
        a1 = ceil(Q / N + 0.5) - 1;
        
        % Loop over a2 parameter
        for a2 = 2:a1
            % Perform inverse FFT and update A1
            A1 = ifft(A((a2 - 1) * N + 1:min(a2 * N, Q)), N) * sqrt(N / (2 * dt)) + A1;
        end
        % Compute the final WDF matrix entry
        y(:, n - n1 + 1) = (A1(m) .* exp(1i * 2 * pi / N * U * (m - 1))).';
    end
    
    % Plot the WDF
    image(t, f, abs(y) / max(max(abs(y))) * 400)
    colormap(gray(256))
    set(gca, 'Ydir', 'normal')
    set(gca, 'Fontsize', 12)
    xlabel('Time (Sec)', 'Fontsize', 12)
    ylabel('Frequency (Hz)', 'Fontsize', 12)
    title('Wigner Distribution Function', 'Fontsize', 12)
end