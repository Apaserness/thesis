[x, fs] = audioread('whitenoise.wav'); 
x = mean(x, 2); 

fc_low = 100;     % Low-shelf cutoff frequency
fc_high = 8000;   % High-shelf cutoff frequency
fc_peak = 1000;   % Peak filter center frequency

G_low = 10;        % Gain for low-shelf
G_high = -20;      % Gain for high-shelf
G_peak = 20;      % Gain for peak filter

Q = 2;            % Q factor for peak filter

y_low = lowshelving(x, fs, fc_low, G_low);
y_peak = peakfilter(y_low, fs, fc_peak, G_peak, Q);
y_eq = highshelving(y_peak, fs, fc_high, G_high);

sound(y_eq, fs);

function y = lowshelving(x, fs, fc, G)
    K = tan(pi * fc / fs);
    V0 = 10^(G/20);
    
    if G >= 0
        b0 = (1 + sqrt(2*V0)*K + V0*K^2) / (1 + sqrt(2)*K + K^2);
        b1 = (2 * (V0*K^2 - 1)) / (1 + sqrt(2)*K + K^2);
        b2 = (1 - sqrt(2*V0)*K + V0*K^2) / (1 + sqrt(2)*K + K^2);
        a1 = (2 * (K^2 - 1)) / (1 + sqrt(2)*K + K^2);
        a2 = (1 - sqrt(2)*K + K^2) / (1 + sqrt(2)*K + K^2);
    else
        b0 = (V0 * (1 + sqrt(2)*K + K^2)) / (V0 + sqrt(2*V0)*K + K^2);
        b1 = (2 * V0 * (K^2 - 1)) / (V0 + sqrt(2*V0)*K + K^2);
        b2 = (V0 * (1 - sqrt(2)*K + K^2)) / (V0 + sqrt(2*V0)*K + K^2);
        a1 = (2 * (K^2 - V0)) / (V0 + sqrt(2*V0)*K + K^2);
        a2 = (V0 - sqrt(2*V0)*K + K^2) / (V0 + sqrt(2*V0)*K + K^2);
    end

    N = length(x);
    y = zeros(N,1);
    
    for n = 3:N  
        y(n) = b0*x(n) + b1*x(n-1) + b2*x(n-2) - a1*y(n-1) - a2*y(n-2);
    end
end

function y = highshelving(x, fs, fc, G)
    K = tan(pi * fc / fs);
    V0 = 10^(G/20);

    if G >= 0
        b0 = (V0 + sqrt(2*V0)*K + K^2) / (1 + sqrt(2)*K + K^2);
        b1 = (2 * (K^2 - V0)) / (1 + sqrt(2)*K + K^2);
        b2 = (V0 - sqrt(2*V0)*K + K^2) / (1 + sqrt(2)*K + K^2);
        a1 = (2 * (K^2 - 1)) / (1 + sqrt(2)*K + K^2);
        a2 = (1 - sqrt(2)*K + K^2) / (1 + sqrt(2)*K + K^2);
    else
        b0 = (V0*(1 + sqrt(2)*K + K^2)) / (1 + sqrt(2*V0)*K + V0*K^2);
        b1 = (2 * V0 * (K^2 - 1)) / (1 + sqrt(2*V0)*K + V0*K^2);
        b2 = (V0*(1 - sqrt(2)*K + K^2)) / (1 + sqrt(2*V0)*K + V0*K^2);
        a1 = (2 * (V0 * K^2 - 1)) / (1 + sqrt(2*V0)*K + V0*K^2);
        a2 = (1 - sqrt(2*V0)*K + V0 * K^2) / (1 + sqrt(2*V0)*K + V0*K^2);
    end

    N = length(x);
    y = zeros(N,1);
    
    for n = 3:N  
        y(n) = b0*x(n) + b1*x(n-1) + b2*x(n-2) - a1*y(n-1) - a2*y(n-2);
    end
end

function y = peakfilter(x, fs, fc, G, Q)
    K = tan(pi * (fc / (fs/2)));
    V0 = 10^(G/20);

    if G >= 0
        b0 = (1 + (V0/Q)*K + K^2) / (1 + (1/Q)*K + K^2);
        b1 = (2*(K^2 - 1)) / (1 + (1/Q)*K + K^2);
        b2 = (1 - (V0/Q)*K + K^2) / (1 + (1/Q)*K + K^2);
        a1 = (2*(K^2 - 1)) / (1 + (1/Q)*K + K^2);
        a2 = (1 - (1/Q)*K + K^2) / (1 + (1/Q)*K + K^2);
    else
        b0 = (1 + (1/Q)*K + K^2) / (1 + (1/Q)*K + K^2);
        b1 = (2*(K^2 - 1)) / (1 + (1/Q)*K + K^2);
        b2 = (1 - (1/Q)*K + K^2) / (1 + (1/Q)*K + K^2);
        a1 = (2*(K^2 - 1)) / (1 + (1/(V0*Q))*K + K^2);
        a2 = (1 - (1/(V0*Q))*K + K^2) / (1 + (1/(V0*Q))*K + K^2);
    end

    N = length(x);
    y = zeros(N,1);
    
    for n = 3:N  
        y(n) = b0*x(n) + b1*x(n-1) + b2*x(n-2) - a1*y(n-1) - a2*y(n-2);
    end
end
freqz(y_eq)
