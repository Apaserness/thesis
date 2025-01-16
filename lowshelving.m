fs = 48000;           % Sampling rate in Hz
f_c = 1000;            % Cutoff frequency in Hz
Wc = 2 * f_c / fs;    % Normalized cutoff frequency
G = 12;                % Gain in dB for low-shelf boost

filePath='C:\Users\gaukh\OneDrive\Документы\MATLAB\edm.wav';
[x,fs]=audioread("edm.wav");

if size(x, 2) > 1
    x = mean(x, 2); % Convert to mono by averaging channels
end

y=lowshelving1 (x,Wc,G);

function y = lowshelving1 (x, Wc, G)
    V0=10^(G/20); H0=V0 - 1;
    if G>=0
        c=(tan(pi*Wc/2)-1) / (tan(pi*Wc/2)+1); %boost
    else
        c=(tan(pi*Wc/2)-V0) / (tan(pi*Wc/2)+V0); %cut
    end;
    xh=0;
    for n=1:length(x)
        xh_new=x(n)-c*xh;
        ap_y=c*xh_new+xh;
        xh=xh_new;
        y(n)=0.5*H0*(x(n)+ap_y)+x(n); %change to minus for hs
    end;
end;
sound(y, fs)
figure;
plot(y);

figure;

plot(x);


