%function y = highshelving_second(x, fs, fc, G)
    [x, fs] = audioread('whitenoise.wav'); 
    
    %fs = 48000;
    fc = 2000; % Cutoff frequency
    G = 25;   % Gain in dB
    
    V0 = 10^(G/20);
    K = tan(pi * fc / fs);
    if G>=0
    
        b0 = (V0 + sqrt(2*V0)*K + K^2) / (1 + sqrt(2)*K + K^2);
        b1 = (2 * (K^2 - V0)) / (1 + sqrt(2)*K + K^2);
        b2 = (V0 - sqrt(2*V0)*K + K^2) / (1 + sqrt(2)*K + K^2);
    
        a0 = 1;
        a1 = (2 * (K^2 - 1)) / (1 + sqrt(2)*K + K^2);
        a2 = (1 - sqrt(2)*K + K^2) / (1 + sqrt(2)*K + K^2);
    else
        b0 = (V0*(1 + sqrt(2)*K + K^2)) / (1 + sqrt(2*V0)*K + V0*K^2);
        b1 = (2 * V0 * (K^2 - 1)) / (1 + sqrt(2*V0)*K + V0*K^2);
        b2 = (V0*(1 - sqrt(2)*K + K^2)) / (1 + sqrt(2*V0)*K + V0*K^2);
        a0 = 1;
        a1 = (2 * (V0 * K^2 - 1)) / (1 + sqrt(2*V0)*K + V0*K^2);
        a2 = (1 - sqrt(2*V0)*K + V0 * K^2) / (1 + sqrt(2*V0)*K + V0*K^2);
    end
    b = [b0, b1, b2];
    a = [a0, a1, a2];
    
    %freqz(b, a, 2048, fs);
    
    [H, w] = freqz(b, a, 1024, fs); 
    semilogx(w, 20*log10(abs(H)));


    z=fft(x);

    z = fftshift(z);
    N = length(x);
    f = [-N/2:N/2-1]/N;
    plot(f,abs(z))

    freqz(b,a)

        z=fft(x);

    z = fftshift(z);
    N = length(x);
    f = [-N/2:N/2-1]/N;
    plot(f,abs(z))
    freqz(b,a)
    
    N = length(x);
    y = zeros(N,1); 

    for n = 3:N
        y(n) = b0*x(n) + b1*x(n-1) + b2*x(n-2) - a1*y(n-1) - a2*y(n-2);
    end

    sound(y, fs);
    %N = length(x);
    %y = zeros(N,1);  % Ensure y is initialized properly

    % Apply filter loop
    %for n = 3:N  
    %    y(n) = b0*x(n) + b1*x(n-1) + b2*x(n-2) - a1*y(n-1) - a2*y(n-2);
    %end
   
