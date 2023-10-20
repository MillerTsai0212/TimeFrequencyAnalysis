
clear all
dt = 0.05;
df = 0.05;
t1 = [0:dt:10-dt]; t2=[10:dt:20-dt]; t3=[20:dt:30];
t = [0:dt:30];
f = [-5:df:5];
x = [cos(2 * pi * t1), cos(6 * pi * t2), cos(4 * pi * t3)];
B = 1;
tic
y=recSTFT(x,t,f,B);
toc

function y = recSTFT(x, t, f, B)
%% Step 1 Calculate n0, m0, T, F, N, Q
C = 400;

% Calculate delta t and f
dt = t(2) - t(1);
df = f(2) - f(1);

% Calculate numbers of N points
N = round(1 / (dt * df));

% Calculate numbers of t points
nt = t ./ dt;
nt = round(nt);
nt0 = nt(1);
T = length(nt);

% Calculate numbers of f points
nf = f ./ df;
nf = round(nf); 
nf0 = nf(1);
F = length(nf);

% Calculate points in interval[-B , B]
Q = round(B/dt);

% Define parameter name
X = zeros(T,F);
t_value = zeros(1 , T);
f_value = zeros(1 , F);
x=[x,0];

%% Step 2 n = nt0
% time domain calculation
for n = 1 : T
%% Step 3 determine x1(q)
    % q tatol range
    q_range= zeros(1 , N);
    % q = 0 for 2Q < q <N
    q=[0 : 2 * Q];
    p = round(nt(n) - Q + q);
    p(p<1) = T + 1; 
    p(p>T) = T + 1;
    % Set x1 = 0 if x is out of range
    x1=[x(p), q_range]; 

%% Step 4 : X1(m) = FFT[x1(q)]
    X1=fft(x1 , N);

%% Step 5:Convert X1(m) into X(ndt, mdf)
% frequency domain calculation
    for m = 1 : F
        t_value(1, n) = nt(n) * dt;
        f_value(1, m) = nf(m) * df;
        X(n,m) = X1(1, mod(nf(m), N) + 1) * exp(1i * 2 * pi * (Q - nt(n)) * nf(m) / N) * dt;
    end

%% Step 6 : Set n = n + 1 and return to Step 3 until n = n0 + T - 1
    n = n + 1;
end
%% Step 7 : Draw plot
% Transpose X to let time on x-axis
y = transpose(X);
image(abs(y) / max(max(abs(y))) * C);
% Gray-level figure
colormap(gray(256));
% Let y-axis from small to large
set(gca, 'Ydir', 'normal');
% Set figure description
set(gca, 'Fontsize', 12)
xlabel('Time (Sec)','FontSize', 12)
ylabel('Frequency (Hz)','FontSize', 12)
title('recSTFT of x(t)', 'FontSize', 12)
end


