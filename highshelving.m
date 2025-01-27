fs = 5000;           
fc = 100;            
G = 6;               
t = 0:1/fs:1;
[x,fs]=audioread("C:\Users\gaukh\OneDrive\Документы\MATLAB\edm.wav")
%x = sin(2*pi*100*t) + sin(2*pi*2000*t);
if size(x, 2) > 1
    x = mean(x, 2); 
end
t = (0:length(x)-1) / fs;



V0 = 10^(G / 20); 
K = tan(pi * fc / fs);

if G > 0
    % HF boost coefficients
    b0 = (V0 + sqrt(2*V0)*K + K^2) / (1 + sqrt(2)*K + K^2);
    b1 = (2 * (K^2 - V0)) / (1 + sqrt(2)*K + K^2);
    b2 = (V0 - sqrt(2*V0)*K + K^2) / (1 + sqrt(2)*K + K^2);
    a1 = (2 * (K^2 - 1)) / (1 + sqrt(2)*K + K^2);
    a2 = (1 - sqrt(2)*K + K^2) / (1 + sqrt(2)*K + K^2);
else
    % HF cut coefficients
    b0 = (1 + sqrt(2)*K + K^2) / (1 + sqrt(2)*K + V0*K^2);
    b1 = (2 * (K^2 - 1)) / (1 + sqrt(2)*K + V0*K^2);
    b2 = (1 - sqrt(2)*K + K^2) / (1 + sqrt(2)*K + V0*K^2);
    a1 = (2 * (K^2 - 1)) / (1 + sqrt(2)*K + V0*K^2);
    a2 = (1 - sqrt(2)*K + V0*K^2) / (1 + sqrt(2)*K + V0*K^2);
end

% Display coefficients
fprintf('b0: %.4f, b1: %.4f, b2: %.4f, a1: %.4f, a2: %.4f\n', b0, b1, b2, a1, a2);

x1 = 0; x2 = 0; 
y1 = 0; y2 = 0; 
y = zeros(size(x)); 

for n = 1:length(x)
    y(n) = b0 * x(n) + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
    x2 = x1;
    x1 = x(n);
    y2 = y1;
    y1 = y(n);
end

figure;
subplot(2, 1, 1);
plot(t, x);
title('Input Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(t, y);
title('Filtered Signal');
xlabel('Time (s)');
ylabel('Amplitude');

% Frequency
N = 2^nextpow2(length(x)); 
f = (0:N-1) * (fs / N); 
X = abs(fft(x, N));
Y = abs(fft(y, N));

figure;
subplot(2, 1, 1);
plot(f(1:N/2), X(1:N/2));
title('Input Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(2, 1, 2);
plot(f(1:N/2), Y(1:N/2));
title('Output Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

target_freqs = [100, 2000]; 
indices = round(target_freqs * N / fs); 

for i = 1:length(target_freqs)
    freq = target_freqs(i);
    idx = indices(i);
    boost = Y(idx) / X(idx); 
    boost_dB = 20 * log10(boost); 
    fprintf('Frequency: %d Hz, Boost: %.2f (linear), Gain: %.2f dB\n', freq, boost, boost_dB);
end
