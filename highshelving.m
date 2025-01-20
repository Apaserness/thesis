fs = 5000;           
fc = 100;            
Wc = 2 * fc / fs;    
G = -6;              

t = 0:1/fs:1; 
x = sin(2*pi*100*t) + sin(2*pi*2000*t); 
y = highshelving1(x, Wc, G);

function y = highshelving1(x, Wc, G)
    V = 10^(G/20); 
    H = V - 1;     
    c = 0.5 * ((tan(pi * Wc / 2) - V) / (tan(pi * Wc / 2) + V)); 
    fprintf('V: %.4f, H: %.4f, c: %.4f\n', V, H, c);
    xh = 0;
    y = zeros(size(x)); 
    
    for n = 1:length(x)
        xh_new = x(n) - c * xh; 
        ap_y = c * xh_new + xh; 
        xh = xh_new;            
        y(n) = V * x(n) + (1 - c) * ap_y; 
    end
end

N = length(x);
f = (0:N-1) * (fs / N); 
X = abs(fft(x));
Y = abs(fft(y));

figure;
subplot(2, 1, 1);
plot(f(1:floor(N/2)), X(1:floor(N/2)));
title('Input Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(2, 1, 2);
plot(f(1:floor(N/2)), Y(1:floor(N/2)));
title('Output Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Calculate boost
idx_100Hz = round(100 * N / fs);
idx_2000Hz = round(2000 * N / fs);
boost_100Hz = Y(idx_100Hz) / X(idx_100Hz);
boost_2000Hz = Y(idx_2000Hz) / X(idx_2000Hz);

fprintf('Boost at 100 Hz: %.2f\n', boost_100Hz);
fprintf('Boost at 2000 Hz: %.2f\n', boost_2000Hz);
