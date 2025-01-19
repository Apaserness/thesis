fs = 5000;           % Sampling rate in Hz
fc = 100;            % Cutoff frequency in Hz
Wc = 2 * fc / fs;    % Normalized cutoff frequency
G = -6;                % Gain in dB for low-shelf boost

%filePath='C:\Users\gaukh\OneDrive\Документы\MATLAB\edm.wav';
%[x,fs]=audioread("edm.wav");

t = 0:1/fs:1; % 1 second duration
x = sin(2*pi*100*t) + sin(2*pi*2000*t);

%if size(x, 2) > 1
%    x = mean(x, 2); % Convert to mono by averaging channels
%end

y=highshelving1 (x,Wc,G);

function y = highshelving1 (x, Wc, G)
    V=10^(G/20); H=V - 1;
    if G>=0
        c=(tan(pi*Wc/2)-1) / (tan(pi*Wc/2)+1); %boost
    else
        c=(tan(pi*Wc/2)-V) / (tan(pi*Wc/2)+V); %cut
    end;
    xh=0;
    for n=1:length(x)
        xh_new=x(n)-c*xh;
        ap_y=c*xh_new+xh;
        xh=xh_new;
        y(n)=0.5*H*(x(n)-ap_y)+x(n); %change to plus for ls
    end;
end;

y = y / max(abs(y));

%sound(y, fs)
figure;
plot(y);

figure;

plot(x);
%freq response

% Frequency-domain analysis
N = length(x);                  % Length of the signal
f = (0:N-1) * (fs / N);         % Frequency axis (Hz)

% Compute FFT for input and output signals
X = abs(fft(x));                % Magnitude spectrum of input
Y = abs(fft(y));                % Magnitude spectrum of output

% Plot the spectra
figure;
subplot(2, 1, 1);
plot(f(1:N/2), X(1:N/2));       % Positive frequencies only
title('Input Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(2, 1, 2);
plot(f(1:N/2), Y(1:N/2));       % Positive frequencies only
title('Output Spectrum');
xlabel('Frequency (Hz)');
ylabel('Magnitude');


% Find indices corresponding to 100 Hz and 2000 Hz
idx_100Hz = round(100 * N / fs);
idx_2000Hz = round(2000 * N / fs);

% Compare magnitudes
boost_100Hz = Y(idx_100Hz) / X(idx_100Hz);
boost_2000Hz = Y(idx_2000Hz) / X(idx_2000Hz);

fprintf('Boost at 100 Hz: %.2f\n', boost_100Hz);
fprintf('Boost at 2000 Hz: %.2f\n', boost_2000Hz);
