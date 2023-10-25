clear all

[a1, fs]=audioread('Chord.wav');

x=a1(:,1);

dtau=1/44100; 
dt=0.01; 
df=1; 
sgm=200;

tau = 0 : dtau : 1.6; 
t = 0: dt : max(tau); 
f = 20 : df : 1000;

tic
y=Gabor(x,tau,t,f,sgm);
toc

function y = Gabor(x, tau, t, f, sgm)
%% Step 1:Calculate n0, f0, tauo, T, F, Tau, N, Q
dtau = diff(tau(1:2));
dt = diff(t(1:2));
df = diff(f(1:2));

Tau = numel(tau);
F = numel(f);
T = numel(t);

N = fix(1 / (dtau * df));
Q = fix(1.9143 / (sqrt(sgm) * dtau));
S = fix(dt / dtau);

%% Step 2:n=c0
for time = 1 : T
    C = 400;
    %% Step3:Determine x1(q)
    
    q = (0 : 1 : N - 1)';
    q(q >= 2 * Q) = 0;
    
    x1 = zeros(N, 1);  

    for i = 1:N
        q_val = time * S - Q + q(i); % 算q值
        if q_val < 1 
            q_val = 1;  % 限制下界
        elseif q_val > Tau
            q_val = Tau; % 限制上界
        end
        x1(i) = x(q_val) * exp(-sgm * pi * ((Q - q(i)) * dtau) ^ 2);
    end

%% Step4:X1(m)=FFT(x1(q))
X1 = fft(x1, N);

%% Step5:ConvertX1(m) into X(ndt, mdf)
frequencies = 1:F;  % 頻率範圍
X(frequencies, time) = dtau * exp(j * 2 * pi * (Q - time * S) * frequencies / N) .* X1(frequencies)';

end
y=X;
image(t, f, abs(y) / max(max(abs(y))) * C);
colormap(gray(256));
set(gca, 'Ydir', 'normal');
end