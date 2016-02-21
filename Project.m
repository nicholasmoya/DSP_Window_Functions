%% EE 480L Project

%{
Nicholas Moya
EE 480L
Prof. Blackstone
Nov 28, 2015

This code is the capstone project of the DSP lab. The main
purpose of this project is to survey the different types of window
functons used in filter design, namely

()Rectangle, ()Hanning, ()Hamming, ()Blackman

To do this, four window functions are written by me, named
RectangleWindow2, Hanning2, Hamming2, and Blackman2. Using these
functions, the code conducts this survey in four parts:

1. First, compare the different shapes of the windows
2. Second, compare the magnitudes of the user window functions
3. Third, compare the magnitudes of the Matlab window functions
4. Last, compare the windows as part of filter functions
   by generating a test signal, filtering it through each of the
   filters, and observing the output.
%}

%% FIR Filter

clear all
clc

% Linear Phase System
N = 55; % or N = 100
wl = 0.02*pi;
wu = 0.2*pi;
[h, nh] = FIRdesign(wl, wu, N);  % generate h[n]
figure(1)
subplot(2,1,1)
stem(nh, h)  % plot h[n]
title('h[n] for N = 55');

% Frequency response of h
[H, hf] = freqz(h);          % H is the frequency response of h
subplot(2,1,2)
stem(hf, H)
title('H[w] for N = 55');

% mag of H[w]
mag = 20*log(abs(H));       % mag is the magnitude of H
figure(2)
subplot(2,1,1)
semilogy(hf, mag);           % plots magnitude of H(w)
title('mag of h[n]');

% Phase of H[w]
[p, wp] = phasez(H);        % p is the phase of H
subplot(2,1,2);
plot(p, wp)                 % plots phase of H(w)
title('phase of h[n]');

%% Rectangle Window

clear all
clc

% Linear Phase System
N = 55; % or N = 100
wl = 0.02*pi;   % Lower frequency
wu = 0.2*pi;    % Upper frequency
[h, nh] = FIRdesign(wl, wu, N);  % generate h[n]

% Rectangular Window
[w, nw] = RectangularWindow2(N);
figure(3)
subplot(2,1,1)
stem(nw, w) % plot w[n]
title('w[n] for N = 55');

% Frequency response of w
[W, wf] = freqz(w);          % W(w) is the frequency response of h[n]
subplot(2,1,2)
stem(wf, W)
title('W[w] for N = 55');

% mag of W[w]
mag = 20*log(abs(W));       % mag is the magnitude of W
figure(4)
subplot(2,1,1)
semilogy(wf, mag);           % plots magnitude of W
title('Mag of W[w]');

% Phase of W[w]
[p, wp] = phasez(W);        % p is the phase of W
subplot(2,1,2);
plot(p, wp)                 % plots phase of W(w)
title('Phase of W[w]');

% Filter
filter = h.*w;  % filter = h[n]w[n]
figure(5)
subplot(2,1,1)
stem(nh, filter)    % plot filter
title('filter[n] = h[n]w[n] for N = 55');

% Frequency response of filter
[F, wfil] = freqz(filter);  % F is the frequency response of filter
subplot(2,1,2)
stem(wfil, F)
title('Filter[w] for N = 55');

% mag of W[w]
mag = 20*log(abs(F));       % mag is the magnitude of F
figure(6)
subplot(2,1,1)
semilogy(wfil, mag);           % plots magnitude of F
title('Mag of Filter[w]');

% Phase of W[w]
[p, wp] = phasez(F);        % p is the phase of F
subplot(2,1,2);
plot(p, wp)                 % plots phase of F
title('Phase of Filter[w]');

% signal
n = -49:50;                      % 100 samples
an = cos((0.1*pi)*n);
figure(7)
subplot(2,1,1)
stem(n, an)                     % plots a[n]
title('a[n] = cos[(0.1*pi)n]u[n]');
axis([-49 50 -1.5 1.5]);

% filtered signal
[ya, na] = convolve(filter, nh, an, n); % ya is the filtered signal
subplot(2,1,2)
stem(na, ya)    % plot ya
axis([-80 80 -1.5 1.5]);
title('ya[n]');

%% Hanning Window

clear all
clc

% Linear Phase System
N = 55; % or N = 100
wl = 0.02*pi;   % lower frequency
wu = 0.2*pi;    % upper frequency
[h, nh] = FIRdesign(wl, wu, N);  % generate h[n]

% Hanning Window
[w, nw] = Hanning2(N);
figure(8)
subplot(2,1,1)  % plot window
stem(nw, w) % plot w[n]
title('w[n] for N = 55');

% Frequency response of w
[W, wf] = freqz(w);          % W(w) is the frequency response of w[n]
subplot(2,1,2)
stem(wf, W)
title('W[w] for N = 55');

% mag of W[w]
mag = 20*log(abs(W));       % mag is the magnitude of W(w)
figure(9)
subplot(2,1,1)
semilogy(wf, mag);           % plots magnitude of W(w)
title('Mag of W[w]');

% Phase of W[w]
[p, wp] = phasez(W);        % p is the phase of W
subplot(2,1,2);
plot(p, wp)                 % plots phase of W(w)
title('Phase of W[w]');

% Filter
filter = h.*w;  % filter = h[n]w[n]
figure(10)
subplot(2,1,1)  % plot filter
stem(nh, filter)
title('filter[n] = h[n]w[n] for N = 55');

% Frequency response of filter
[F, wfil] = freqz(filter);  % F is the frequency response of filter
subplot(2,1,2)
stem(wfil, F)   % plot frequency response of filter
title('Filter[w] for N = 55');

% mag of W[w]
mag = 20*log(abs(F));       % mag is the magnitude of F
figure(11)
subplot(2,1,1)
semilogy(wfil, mag);           % plots magnitude of F(w)
title('Mag of Filter[w]');

% Phase of W[w]
[p, wp] = phasez(F);        % p is the phase of F
subplot(2,1,2);
plot(p, wp)                 % plots phase of F(w)
title('Phase of Filter[w]');

% signal
n = -49:50;                      % 100 samples
an = cos((0.1*pi)*n);   % test signal
figure(12)
subplot(2,1,1)
stem(n, an)                     % plots a[n]
title('a[n] = cos[(0.1*pi)n]u[n]');
axis([-49 50 -1.5 1.5]);

% filtered signal
[ya, na] = convolve(filter, nh, an, n); % ya is the filtered signal
subplot(2,1,2)
stem(na, ya)    % plot ya
axis([-80 80 -1.5 1.5]);
title('ya[n]');

%% Hamminging Window

clear all
clc

% Linear Phase System
N = 55; % or N = 100
wl = 0.02*pi;   % lower frequency
wu = 0.2*pi;    % upper frequency
[h, nh] = FIRdesign(wl, wu, N);  % generate h[n]

% Hamming Window
[w, nw] = Hamming2(N);
figure(13)
subplot(2,1,1)
stem(nw, w) % plot w[n]
title('w[n] for N = 55');

% Frequency response of w
[W, wf] = freqz(w);          % W(w) is the frequency response of w[n]
subplot(2,1,2)
stem(wf, W) % frequency response of w[n]
title('W[w] for N = 55');

% mag of W[w]
mag = 20*log(abs(W));       % mag is the magnitude of W(w)
figure(14)
subplot(2,1,1)
semilogy(wf, mag);           % plots magnitude of W(w)
title('Mag of W[w]');

% Phase of W[w]
[p, wp] = phasez(W);        % p is the phase of W
subplot(2,1,2);
plot(p, wp)                 % plots phase of W(w)
title('Phase of W[w]');

% Filter
filter = h.*w;  % filter = h[n]w[n]
figure(15)
subplot(2,1,1)
stem(nh, filter)
title('filter[n] = h[n]w[n] for N = 55');

% Frequency response of filter
[F, wfil] = freqz(filter);  % F is the frequency response of filter
subplot(2,1,2)
stem(wfil, F)
title('Filter[w] for N = 55');

% mag of W[w]
mag = 20*log(abs(F));       % mag is the magnitude of F
figure(16)
subplot(2,1,1)
semilogy(wfil, mag);           % plots magnitude of F(w)
title('Mag of Filter[w]');

% Phase of W[w]
[p, wp] = phasez(F);        % p is the phase of F
subplot(2,1,2);
plot(p, wp)                 % plots phase of F(w)
title('Phase of Filter[w]');

% signal
n = -49:50;                      % 100 samples
an = cos((0.1*pi)*n);   % test signal
figure(17)
subplot(2,1,1)
stem(n, an)                     % plots a[n]
title('a[n] = cos[(0.1*pi)n]u[n]');
axis([-49 50 -1.5 1.5]);

% filtered signal
[ya, na] = convolve(filter, nh, an, n); % ya is the filtered signal
subplot(2,1,2)
stem(na, ya)
axis([-80 80 -1.5 1.5]);
title('ya[n]');

%% Blackman Window

clear all
clc

% Linear Phase System
N = 55; % or N = 100
wl = 0.02*pi;   % lower frequency
wu = 0.2*pi;    % upper frequency
[h, nh] = FIRdesign(wl, wu, N);  % generate h[n]

% Blackman Window
[w, nw] = Blackman2(N);
figure(18)
subplot(2,1,1)
stem(nw, w) % plot w[n]
title('w[n] for N = 55');

% Frequency response of w
[W, wf] = freqz(w);          % W(w) is the frequency response of w[n]
subplot(2,1,2)
stem(wf, W)
title('W[w] for N = 55');

% mag of W[w]
mag = 20*log(abs(W));       % mag is the magnitude of W(w)
figure(19)
subplot(2,1,1)
semilogy(wf, mag);           % plots magnitude of W(w)
title('Mag of W[w]');

% Phase of W[w]
[p, wp] = phasez(W);        % p is the phase of W
subplot(2,1,2);
plot(p, wp)                 % plots phase of W(w)
title('Phase of W[w]');

% Filter
filter = h.*w;
figure(20)
subplot(2,1,1)
stem(nh, filter)
title('filter[n] = h[n]w[n] for N = 55');

% Frequency response of filter
[F, wfil] = freqz(filter);  % F is the frequency response of filter
subplot(2,1,2)
stem(wfil, F)
title('Filter[w] for N = 55');

% mag of W[w]
mag = 20*log(abs(F));       % mag is the magnitude of F
figure(21)
subplot(2,1,1)
semilogy(wfil, mag);           % plots magnitude of F(w)
title('Mag of Filter[w]');

% Phase of W[w]
[p, wp] = phasez(F);        % p is the phase of F
subplot(2,1,2);
plot(p, wp)                 % plots phase of F(w)
title('Phase of Filter[w]');

% signal
n = -49:50;                      % 100 samples
an = cos((0.1*pi)*n);   % test signal
figure(22)
subplot(2,1,1)
stem(n, an)                     % plots a[n]
title('a[n] = cos[(0.1*pi)n]u[n]');
axis([-49 50 -1.5 1.5]);

% filtered signal
[ya, na] = convolve(filter, nh, an, n); % ya is the filtered signal
subplot(2,1,2)
stem(na, ya)    % plot ya
axis([-80 80 -1.5 1.5]);
title('ya[n]');

%% Comparing Shape of Windows 

clear all
clc

N = 55; % or 100
n = -(N-1)/2:(N-1)/2;
w0 = ones(1,N); % Rectangle
w1 = 1/2+(1/2)*cos(2*pi*n/(N-1));   % Hanning
w2 = 0.54+0.46*cos(2*pi*n/(N-1));   % Hamming
w3 = 0.42+0.5*cos(2*pi*n/(N-1))+0.08*cos(4*pi*n/(N-1)); % Blackman

figure(23)
hold on
stem(n, w0)
stem(n, w1) % plot windows overlapping for comparison
stem(n, w2)
stem(n, w3)
hold off

title('Comparing Shape of Windows')
xlabel('Samples')
ylabel('Amplitude')
legend('w0 = Rectangle', 'w1 = Hanning', 'w2 = Hamming', 'w3 = Blackman')

%% Comparing Magnitude of W[w]

clear all
clc

N = 55; % or 100

% windows
[w0, nw0] = RectangularWindow2(N);
[w1, nw1] = Hanning2(N);
[w2, nw2] = Hamming2(N);
[w3, nw3] = Blackman2(N);

% frequency response of windows
[W0, wf0] = freqz(w0);
[W1, wf1] = freqz(w1);
[W2, wf2] = freqz(w2);
[W3, wf3] = freqz(w3);

% magnitude of frequency response of windows
mag0 = 20*log(abs(W0));
mag1 = 20*log(abs(W1));
mag2 = 20*log(abs(W2));
mag3 = 20*log(abs(W3));

% plot magnitude of windows overlapping for comparison
figure(24)
hold on
semilogy(wf0, mag0);
semilogy(wf1, mag1);
semilogy(wf2, mag2);
semilogy(wf3, mag3);
hold off

% legend
grid on
title('Comparing Magnitude of Windows')
xlabel('Samples')
ylabel('Magnitude (dB)')
legend('w0 = Rectangle', 'w1 = Hanning', 'w2 = Hamming', 'w3 = Blackman')

%% Matlab Rectangle Window Functions

clear all
clc

N = 55; % or 100

% windows
w0 = rectwin(N);
w1 = hann(N);
w2 = hamming(N);
w3 = blackman(N);

% frequency response of windows
[W0, wf0] = freqz(w0);
[W1, wf1] = freqz(w1);
[W2, wf2] = freqz(w2);
[W3, wf3] = freqz(w3);

% magnitude of frequency response of windows
mag0 = 20*log(abs(W0));
mag1 = 20*log(abs(W1));
mag2 = 20*log(abs(W2));
mag3 = 20*log(abs(W3));

% plot magnitude of windows overlapping for comparison
figure(25)
hold on
semilogy(wf0, mag0);
semilogy(wf1, mag1);
semilogy(wf2, mag2);
semilogy(wf3, mag3);
hold off

% legend
grid on
title('Comparing Magnitude of Windows')
xlabel('Samples')
ylabel('Magnitude (dB)')
legend('w0 = Rectangle', 'w1 = Hanning', 'w2 = Hamming', 'w3 = Blackman')


