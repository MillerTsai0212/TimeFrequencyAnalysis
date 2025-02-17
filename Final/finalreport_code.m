clear all

%% Generate the original signal
object1_R = 1;
object2_R= 0.3;
x = linspace(0, 1, 1000);
y = object1_R * exp(-400*(x - 0.5).^2) + object2_R * exp(-200*(x - 0.2).^2);

%% Generate noise
snr_target = 10;
original_power = rms(y)^2;
noise_power = original_power / (10^(snr_target / 10));
noise = sqrt(noise_power) * randn(size(x));

%% Add noise to the original signal
y_noisy = y + noise;

%% Perform discrete wavelet transform
level = 5;
wavelet_type = 'sym5';
[coefficients, L] = dwtpa(y_noisy, level, wavelet_type);

%% Determine the threshold using SureHard thresholding
noise_std = std(noise); % Estimate of the standard deviation of the noise
n = numel(y);
threshold = noise_std * sqrt(2 * log(n));
coefficients_denoised = Softthresh(coefficients, threshold);

%% Perform inverse discrete wavelet transform to obtain the denoised signal
y_denoised = Reconstruct(coefficients_denoised, L, wavelet_type);

%% Plot the original signal, the noisy signal, and the denoised signal

% Plot the noisy signal
figure;
plot(x, y_noisy, 'b', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title(['Noisy signal SNR = ' num2str(snr(y, noise)) ' dB']);
grid on;
% Plot the original signal
figure;
plot(x, y, 'g', 'LineWidth', 2 );
xlabel('Time');
ylabel('Amplitude');
title('Original Signal');
grid on;


% Plot the denoised signal
figure;
plot(x, y_denoised, 'r--', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Amplitude');
title('Denoised Signal');
grid on;

% Compared the denoised signal andthe original signal
figure;
plot(x, y, 'g', 'LineWidth', 2, 'DisplayName', 'Original Signal');
hold on
plot(x, y_denoised, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Denoised Signal');
xlabel('Time');
ylabel('Amplitude');
legend;
grid on;

function [c, l] = dwtpa(x, level, LF_filter)

% Convert the filter name from string to character array
[Lo_f, Hi_f] = wfilters(LF_filter, 'd');

% Initialization
s = size(x);
x = x(:).'; % Convert the signal to a row vector
c = [];
l = zeros(1, level + 2, 'like', real(x([])));

l(end) = length(x);
for k = 1:level
    [x, d] = dwtc(x, Lo_f, Hi_f); % Decomposition
    c = [d, c]; % Store coefficients
    l(level + 2 - k) = length(d); % Store length
end

% Final low-frequency coefficients
c = [x, c];
l(1) = length(x);

if s(1) > 1
    c = c.'; 
    l = l';
end

end

function y = Softthresh(x, t)

    y = sign(x) .* max(abs(x) - t, 0);
    
end

function x = Reconstruct(c,l,varargin)

x = appcoef(c,l,varargin{:},0);

end

function [a, d] = dwtc(x, Lo_f, Hi_f)

% Check parameters for Extension and Shift.
DWT_Attribute = getappdata(0, 'DWT_Attribute');
Extension = DWT_Attribute.extMode; % Default: Extension.
Shift = DWT_Attribute.shift1D; % Default: Shift.

% Compute sizes and shape.
lf = length(Lo_f);
lx = length(x);

% Extend, Decompose & Extract coefficients.
first = 2 - Shift;
lenEXT = lf / 2;
last = 2 * ceil(lx / 2);

x_ex = wextend('1D', Extension, x, lenEXT);

% Compute LF coefficients.
z = wconv1(x_ex, Lo_f, 'valid');
a = z(first:2:last);

% Compute HF coefficients.
z = wconv1(x_ex, Hi_f, 'valid');
d = z(first:2:last);

end