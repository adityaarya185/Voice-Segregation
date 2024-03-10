%%
clc;close all;clear all;
% Load the audio file
audioFile = 'SAMPLE6.wav';
[y, Fs] = audioread(audioFile);
% Plot the waveform in the time domain
t = (0:length(y) - 1) / Fs;
figure;
subplot(2, 2, 1);
plot(t, y);
title('Audio Waveform');
xlabel('Time (s)');
ylabel('Amplitude');

% Objective 2: Effect of LPF filtering in Time Domain
% Design a low-pass filter using Butterworth
fc_lp = 3000; % Cutoff frequency for low-pass filter in Hz
order = 8; % Filter order
[b_lp, a_lp] = butter(order, fc_lp / (Fs / 2), 'low');

% Apply the low-pass filter to the audio using filter

y_lp = filter(b_lp, a_lp, y);

% Plot the filtered waveform in the time domain
subplot(2, 2, 2);
plot( t,y_lp);
title('Low-Pass Filtered Audio');
xlabel('Time (s)');
ylabel('Amplitude');

% Objective 3: Effect of HPF filtering in Time Domain
% Design a high-pass filter using Butterworth
fc_hp = 1000; % Cutoff frequency for high-pass filter in Hz
[b_hp, a_hp] = butter(order, fc_hp / (Fs / 2), 'high');

% Apply the high-pass filter to the audio using filter
y_hp = filter(b_hp, a_hp, y);

% Plot the filtered waveform in the time domain
subplot(2, 2, 3);
plot(t, y_hp);
title('High-Pass Filtered Audio');
xlabel('Time (s)');
ylabel('Amplitude');

% Objective 4: Segregate Voice and Un-voice parts based on frequency

% Design a low-pass filter for voice (below 6000 Hz)
fc_voice = 6000;  fc_unvoiced = 50; % Cutoff frequency for voice in Hz
f_pass=[fc_unvoiced,fc_voice];
[b_voice, a_voice] = butter(5, f_pass / (Fs / 2), 'bandpass');

% Design a high-pass filter for unvoiced (above 50 Hz)
fc_unvoiced = 50; % Cutoff frequency for unvoiced in Hz
[b_unvoiced, a_unvoiced] = butter(order, f_pass / (Fs / 2), 'stop');

% Apply the low-pass filter to the audio using filter to isolate voice
y_voice = filter(b_voice, a_voice, y);

% Apply the high-pass filter to the audio using filter to isolate unvoiced
y_unvoiced = filter(b_unvoiced, a_unvoiced, y);

% Plot the filtered voice and unvoiced signals in the time domain
subplot(2, 2, 4);
plot(t, y_voice);
title('Voiced and Unvoiced Segregation (Frequency-Based)');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Voiced', 'Unvoiced');


player1=audioplayer(y,Fs);
play(player1);
pause(length(y)/Fs);
player2=audioplayer(y_lp,Fs);
play(player2);    
pause(length(y_lp)/Fs);
player3=audioplayer(y_hp,Fs);
play(player3);
pause(length(y_hp)/Fs);
player4=audioplayer(y_voice,Fs);
play(player4);
pause(length(y_voice)/Fs);
player5=audioplayer(y_unvoiced,Fs);
play(player5);
%%
%Objective 1: Visual Analysis about the frequency Details
%Load the audio file
clc;close all;clear all;
[y, Fs] = audioread('SAMPLE6.wav');

% Calculate the time vector
t = (0:length(y) - 1) / Fs;

% Calculate the single-sided amplitude spectrum
N = length(y);
f = (0:(N-1)) * Fs / N;
Y = fft(y);
P = abs(Y);

% Plot the spectrum in the frequency domain
subplot(2, 2, 1);
plot(f, P);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Frequency Domain Spectrum');


% Objective 2: Effect of LPF filtering in Frequency Domain
% Load the audio file

% Design a low-pass filter in the frequency domain
cutoff_frequency = 6000;  % Adjust this cutoff frequency as needed

% Create the filter with the same size as Y
lpf = zeros(size(Y));
lpf(f<=cutoff_frequency,:) = 1;

% Apply the filter in the frequency domain
Y_lpf = Y .* lpf;
y_lpf=ifft(Y_lpf);
Y_lpf=abs(Y_lpf);
y_lpf=real(y_lpf);
% Compute the inverse FFT to get the filtered signal
sound(y_lpf,Fs);
pause(length(Y_lpf)/Fs);
% Take only the real part of the filtered signal (discard imaginary part)


% Plot the original and filtered signals

subplot(2, 2, 2);
plot(f,Y_lpf);
title('Original Signal');
xlabel('freq (s)');
ylabel('Amplitude');





% Objective 3: Effect of HPF filtering in Frequency Domain
% Load the audio file

% Design a high-pass filter in the frequency domain
cutoff_frequency = 1000;  % Adjust this cutoff frequency as needed

% Create the filter with the same size as Y
hpf = zeros(size(Y));
hpf(f >= cutoff_frequency,:) = 1;

% Apply the filter in the frequency domain
Y_hpf = Y .* hpf;

% Compute the inverse FFT to get the filtered signal
y_hpf = ifft(Y_hpf);

% Take only the real part of the filtered signal (discard imaginary part)
y_hpf = real(y_hpf);


Y_hpf=abs(Y_hpf);
subplot(2, 2, 3);
plot(f,Y_hpf);
title('High-Pass Filtered Signal (in Frequency Domain)');
xlabel('Freq (s)');
ylabel('Amplitude');
sound(y_hpf,Fs);







