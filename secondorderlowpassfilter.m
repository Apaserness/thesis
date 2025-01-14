fs = 48000;          
fc = 1000;           
Q = 0.707 % butterworths quality parameter (the higher the Q then more resonance)

% Calculate Filter Coefficients (page 50)
K = tan(pi * fc / fs);       
b0 = K^2 * Q / (K^2 * Q + K + Q);
b1 = 2*K^2*Q / (K^2 * Q + K + Q);
b2 = K^2*Q / (K^2 * Q + K + Q); 
a1 = 2*Q * (K^2 - 1) / (K^2 * Q + K + Q);
a2 = (K^2 - K + Q) / (K^2 * Q + K + Q);
fprintf('Filter Coefficients:\nb0 = %.4f, b1 = %.4f, a1 = %.4f\n', b0, b1, a1);

filePath='C:\Users\gaukh\OneDrive\Документы\MATLAB\edm.wav';
[x,fs]=audioread("edm.wav");
%x=sin(2*pi)
fs
y = zeros(size(x));          % Pre-allocate filtered signal
xh = zeros(1);               % Initialize delayed sample for x
y = zeros(size(x)); % Pre-allocate filtered signal
for n = 3:length(x)
    % second-order difference equation
    y(n) = b0 * x(n) + b1 * x(n-1) + b2 * x(n-2) ...
           - a1 * y(n-1) - a2 * y(n-2);
end

audiowrite('second_order_lowpass_filtered_audio.wav', y, fs); 
%disp('"second_order_lowpass_filtered_audio.wav".');

% FFT
L = length(y);                
Y = fft(y);                   
P2 = abs(Y / L);              % Two-sided spectrum
P1 = P2(1:L/2+1);             % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);  % Adjust amplitude for single-sided
f = fs * (0:(L/2)) / L;       % Frequency vector

% Plot Original Signal and Frequency Response of Filtered Signal
time = (0:length(x)-1) / fs; % Time vector for plotting

figure;
% Plot original signal
subplot(2,1,1);
plot(x);
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

db=20*log10(P1);
db=db-max(db);
%f=(0:(length(db)-1))*(fs/length(db)); %so that the graph wouldnt go to -100

subplot(2,1,2);
plot(f, db); % Convert magnitude to dB
title('Frequency Response of Filtered Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
%xlim([0 100]);
%ylim([-60 0])
grid on;

%figure;
%freqz([b0,b1,b2],[a1,a2]);
%z=filter([b0,b1,b2],[a1,a2],x(:,1));
%sound(z, 48000)

figure;
lowpass(x,150,fs)
sound(x,48000)