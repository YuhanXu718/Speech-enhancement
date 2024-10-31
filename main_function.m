clear; clc; close all;

% load and process file: clean bubble stationary 
[clean, fs_clean] = audioread('clean_speech.wav');
[noise, fs_noise] = audioread('babble_noise.wav');
[noise_sta, fs_noise_sta] = audioread('stationary speech-shaped noise.wav');
[clean, noise, noise_sta, fs] = load_audio(fs_clean, fs_noise, fs_noise_sta,clean, noise, noise_sta);

% generate noisy siganl
desired_SNR = 10;
[noisy, noisy_sta] = generate_noisy(clean, noise, noise_sta, desired_SNR);

audiowrite('babble_clean.wav', noisy, fs);
audiowrite('stationary_clean.wav', noisy_sta, fs);

% compute stft of noisy signal
window_length = 0.025; % 25 ms 
window_overlap = 0.015; % 15 ms
nfft = 512; % 256 1024
nfft_sta = 512;
[S_noisy, F, T, window, noverlap] = win_stft(noisy, fs, window_length, window_overlap, nfft);
[S_noisy_sta, F_sta, T_sta, window_sta, noverlap_sta] = win_stft(noisy_sta, fs, window_length, window_overlap, nfft_sta);

%% Method 1 estimate nois PSD using the first 0.5 seconds of the pure noise segment
noise_duration = 0.5;
noise_est_m1 = m1_pure_noise(noisy, fs, noise_duration, window_length, window_overlap, nfft);
noise_est_m1_sta = m1_pure_noise(noisy_sta, fs, noise_duration, window_length, window_overlap, nfft_sta);
% apply wiener
alpha = 1; 
gain_floor = 0.1;
smoothing_length = 1;
[denoised_m1, H_smoothed_m1] = wiener(S_noisy, noisy, fs, window, noverlap, nfft, noise_est_m1, alpha, gain_floor, smoothing_length);
[denoised_m1_sta, H_smoothed_m1_sta] = wiener(S_noisy_sta, noisy_sta, fs, window_sta, noverlap_sta, nfft_sta, noise_est_m1_sta, alpha, gain_floor, smoothing_length);


%% Method 2 estimate nois PSD using VADs from noisy signal
energy_threshold = median(sum(abs(S_noisy).^2, 1)) * 0.5;  
% threshold, using energy for every fram to estimat and distinguish noise and speech (0.5 may change? see effect)
noise_est_m2 = m2_vads(S_noisy, energy_threshold);

energy_threshold_sta = median(sum(abs(S_noisy_sta).^2, 1)) * 0.5;
noise_est_m2_sta = m2_vads(S_noisy_sta, energy_threshold_sta);

% apply wiener
alpha = 1; % stationary set to 1 % nonstationary bigger
gain_floor = 0.1;
smoothing_length = 1; %  stationary 1 2 nonstation 3 5
[denoised_m2, H_smoothed_m2] = wiener(S_noisy, noisy, fs, window, noverlap, nfft, noise_est_m2, alpha, gain_floor, smoothing_length);
[denoised_m2_sta, H_smoothed_m2_sta] = wiener(S_noisy_sta, noisy_sta, fs, window_sta, noverlap_sta, nfft_sta, noise_est_m2_sta, alpha, gain_floor, smoothing_length);

%% Method 3 estimate nois PSD using welch
beta = 1.3; % 0.05 0.5 0.2
gain_floor = 0.1;
[denoised_m3, H_smoothed_m3, freq_welch, noisy_psd, noise_psd] = m3_welch_wiener(noisy, noise_duration, fs, window, noverlap, nfft, beta, gain_floor);
[denoised_m3_sta, H_smoothed_m3_sta, freq_welch_sta, noisy_psd_sta, noise_psd_sta] = m3_welch_wiener(noisy_sta, noise_duration, fs, window_sta, noverlap_sta, nfft_sta, beta, gain_floor); %, beta, gain_floor
figure;
plot(freq_welch, 10*log10(noisy_psd), 'Color', [0 0 0.5], 'LineWidth', 1.8); % Plot noisy signal PSD in dB
hold on;
plot(freq_welch, 10*log10(noise_psd), 'Color', [0.3 0.75 0.93], 'LineWidth', 1.5); % Plot noise-only PSD in dB
hold on;
plot(freq_welch, 10*log10(noisy_psd_sta), 'Color', [0.5 0 0], 'LineWidth', 1.8); % Plot noisy signal PSD in dB
hold on;
plot(freq_welch, 10*log10(noise_psd_sta), 'Color', [1 0.4 0.4], 'LineWidth', 1.5); % Plot noise-only PSD in dB
title('PSD Comparison: Noisy Signal vs Noise-Only');
xlabel('Frequency/Hz');
ylabel('Power/Frequency dB/Hz');
legend('Noisy Signal', 'Noise-Only', 'Noisy Signal sta', 'Noise-Only sta');
grid on;
%% Evaluation

% snr
noise_only_noisy = noisy - clean;
snr_before = 10 * log10(sum(clean.^2) / sum(noise_only_noisy.^2));

% m1
noise_only_m1 = denoised_m1 - clean;
snr_after_m1 = 10 * log10(sum(clean.^2) / sum(noise_only_m1.^2));

noise_only_m1_sta = denoised_m1_sta - clean;
snr_after_m1_sta = 10 * log10(sum(clean.^2) / sum(noise_only_m1_sta.^2));

% m2
noise_only_m2 = denoised_m2 - clean;
snr_after_m2 = 10 * log10(sum(clean.^2) / sum(noise_only_m2.^2));

noise_only_m2_sta = denoised_m2_sta - clean;
snr_after_m2_sta = 10 * log10(sum(clean.^2) / sum(noise_only_m2_sta.^2));

% m3
noise_only_m3 = denoised_m3 - clean;
snr_after_m3 = 10 * log10(sum(clean.^2) / sum(noise_only_m3.^2));

noise_only_m3_sta = denoised_m3_sta - clean;
snr_after_m3_sta = 10 * log10(sum(clean.^2) / sum(noise_only_m3_sta.^2));

fprintf('Before enhancement SNR: %.2f dB\n', snr_before);
fprintf('Babble noisy signal after Method1 enhancement SNR: %.2f dB\n', snr_after_m1);
fprintf('Babble noisy signal after Method2 enhancement SNR: %.2f dB\n', snr_after_m2);
fprintf('Babble noisy signal after Method3 enhancement SNR: %.2f dB\n', snr_after_m3);

fprintf('Stationary noisy signal after Method1 enhancement SNR: %.2f dB\n', snr_after_m1_sta);
fprintf('Stationary noisy signal after Method2 enhancement SNR: %.2f dB\n', snr_after_m2_sta);
fprintf('Stationary noisy signal after Method3 enhancement SNR: %.2f dB\n', snr_after_m3_sta);

% SNR improvement
snr_improvement_m1 = snr_after_m1 - snr_before;
snr_improvement_m2 = snr_after_m2 - snr_before;
snr_improvement_m3 = snr_after_m3 - snr_before;

snr_improvement_m1_sta = snr_after_m1_sta - snr_before;
snr_improvement_m2_sta = snr_after_m2_sta - snr_before;
snr_improvement_m3_sta = snr_after_m3_sta - snr_before;

fprintf('SNR improvement of Babble noisy signal via Method 1: %.2f dB\n', snr_improvement_m1);
fprintf('SNR improvement of Babble noisy signal via Method 2: %.2f dB\n', snr_improvement_m2);
fprintf('SNR improvement of Babble noisy signal via Method 3: %.2f dB\n', snr_improvement_m3);
fprintf('SNR improvement of Stationary noisy signal via Method 1: %.2f dB\n', snr_improvement_m1_sta);
fprintf('SNR improvement of Stationary noisy signal via Method 2: %.2f dB\n', snr_improvement_m2_sta);
fprintf('SNR improvement of Stationary noisy signal via Method 3: %.2f dB\n', snr_improvement_m3_sta);


% Audio file
denoised_m1 = denoised_m1 / max(abs(denoised_m1));
denoised_m2 = denoised_m2 / max(abs(denoised_m2));
denoised_m3 = denoised_m3 / max(abs(denoised_m3));
denoised_m1_sta = denoised_m1_sta / max(abs(denoised_m1_sta));
denoised_m2_sta = denoised_m2_sta / max(abs(denoised_m2_sta));
denoised_m3_sta = denoised_m3_sta / max(abs(denoised_m3_sta));
audiowrite('denoised_method1_babble.wav', denoised_m1, fs);
audiowrite('denoised_method2_babble.wav', denoised_m2, fs);
audiowrite('denoised_method3_babble.wav', denoised_m3, fs);
audiowrite('denoised_method1_stationary.wav', denoised_m1_sta, fs);
audiowrite('denoised_method2_stationary.wav', denoised_m2_sta, fs);
audiowrite('denoised_method3_stationary.wav', denoised_m3_sta, fs);

addpath PESQ
% PESQ
pesq_before = pesq('clean_speech.wav','babble_clean.wav');
pesq_after_m1 = pesq('clean_speech.wav','denoised_method1_babble.wav');
pesq_after_m2 = pesq('clean_speech.wav','denoised_method2_babble.wav');
pesq_after_m3 = pesq('clean_speech.wav','denoised_method3_babble.wav');
disp(['PESQ score before enhancement Bubble:', num2str(pesq_before)]);
disp(['PESQ score after Method 1 Babble:', num2str(pesq_after_m1)]);
disp(['PESQ score after Method 2 Babble:', num2str(pesq_after_m2)]);
disp(['PESQ score after Method 3 Babble:', num2str(pesq_after_m3)]);


pesq_before_sta = pesq('clean_speech.wav','stationary_clean.wav');
pesq_after_sta_m1 = pesq('clean_speech.wav','denoised_method1_stationary.wav');
pesq_after_sta_m2 = pesq('clean_speech.wav','denoised_method2_stationary.wav');
pesq_after_sta_m3 = pesq('clean_speech.wav','denoised_method3_stationary.wav');
disp(['PESQ score before enhancement Stationary:', num2str(pesq_before_sta)]);
disp(['PESQ score after Method 1 Stationary:', num2str(pesq_after_sta_m1)]);
disp(['PESQ score after Method 2 Stationary:', num2str(pesq_after_sta_m2)]);
disp(['PESQ score after Method 3 Stationary:', num2str(pesq_after_sta_m3)]);


% stoi
stoi_before = stoi(clean,noisy,fs);
stoi_after_m1 = stoi(clean,denoised_m1,fs);
stoi_after_m2 = stoi(clean,denoised_m2,fs);
stoi_after_m3 = stoi(clean,denoised_m3,fs);
disp(['STOI score before enhancement Babble:', num2str(stoi_before)]);
disp(['STOI score after Method 1 Babble:', num2str(stoi_after_m1)]);
disp(['STOI score after Method 2 Babble:', num2str(stoi_after_m2)]);
disp(['STOI score after Method 3 Babble:', num2str(stoi_after_m3)]);

stoi_before_sta = stoi(clean,noisy_sta,fs);
stoi_after_sta_m1 = stoi(clean,denoised_m1_sta,fs);
stoi_after_sta_m2 = stoi(clean,denoised_m2_sta,fs);
stoi_after_sta_m3 = stoi(clean,denoised_m3_sta,fs);
disp(['STOI score before enhancement Stationary:', num2str(stoi_before_sta)]);
disp(['STOI score after Method 1 Stationary:', num2str(stoi_after_sta_m1)]);
disp(['STOI score after Method 2 Stationary:', num2str(stoi_after_sta_m2)]);
disp(['STOI score after Method 3 Stationary:', num2str(stoi_after_sta_m3)]);

% Visualization

t = (0:length(noisy) - 1) / fs;


figure;
subplot(4,2,1);
plot(t, clean);
title('Clean signal');
xlabel('time/s');
ylabel('Amplitude');
grid on;
subplot(4,2,2);
spectrogram(clean,window_length * fs,window_overlap * fs,nfft,fs,'yaxis');
subplot(4,2,3);
plot(t, denoised_m1);
title('Method 1 enhanced bubble noisy');
xlabel('time/s');
ylabel('Amplitude');
grid on;
subplot(4,2,4);
spectrogram(denoised_m1,window_length * fs,window_overlap * fs,nfft,fs,'yaxis');

subplot(4,2,5);
plot(t, denoised_m2);
title('Method 2 enhanced bubble noisy');
xlabel('time/s');
ylabel('Amplitude');
grid on;
subplot(4,2,6);
spectrogram(denoised_m2,window_length * fs,window_overlap * fs,nfft,fs,'yaxis');

subplot(4,2,7);
plot(t, denoised_m3);
title('Method 3 enhanced bubble noisy');
xlabel('time/s');
ylabel('Amplitude');
grid on;
subplot(4,2,8);
spectrogram(denoised_m3,window_length * fs,window_overlap * fs,nfft,fs,'yaxis');


% stationary
figure;
subplot(4,2,1);
plot(t, clean);
title('Clean signal');
xlabel('time/s');
ylabel('Amplitude');
grid on;
subplot(4,2,2);
spectrogram(clean,window_length * fs,window_overlap * fs,nfft,fs,'yaxis');
subplot(4,2,3);
plot(t, denoised_m1_sta);
title('Method 1 enhanced stationary noisy');
xlabel('time/s');
ylabel('Amplitude');
grid on;
subplot(4,2,4);
spectrogram(denoised_m1_sta,window_length * fs,window_overlap * fs,nfft,fs,'yaxis');
subplot(4,2,5);
plot(t, denoised_m2_sta);
title('Method 2 enhanced stationary noisy');
xlabel('time/s');
ylabel('Amplitude');
grid on;
subplot(4,2,6);
spectrogram(denoised_m2_sta,window_length * fs,window_overlap * fs,nfft,fs,'yaxis');
subplot(4,2,7);
plot(t, denoised_m3_sta);
title('Method 3 enhanced stationary noisy');
xlabel('time/s');
ylabel('Amplitude');
grid on;
subplot(4,2,8);
spectrogram(denoised_m3_sta,window_length * fs,window_overlap * fs,nfft,fs,'yaxis');
% Waveform comparasion
figure;
plot(t,denoised_m1 , 'k', 'LineWidth', 1.5);
hold on;
plot(t, denoised_m1_sta, 'r');
plot(t, clean, 'y');
title('Comparasion between bubble and stationary denoised signal Method 1');
xlabel('time/s');
ylabel('Amplitude');
legend('bubble', 'stationary', 'clean');
grid on;

figure;
plot(t,denoised_m2 , 'k', 'LineWidth', 1.5);
hold on;
plot(t, denoised_m2_sta, 'r');
hold on;
plot(t, clean, 'y');
title('Comparasion between bubble and stationary denoised signal Method 2');
xlabel('time/s');
ylabel('Amplitude');
legend('bubble', 'stationary', 'clean');
grid on;

figure;
plot(t, denoised_m3, 'k', 'LineWidth', 1.5);
hold on;
plot(t, denoised_m3_sta, 'r');
plot(t, clean, 'y');
title('Comparasion between bubble and stationary denoised signal Method 3');
xlabel('time/s');
ylabel('Amplitude');
legend('bubble', 'stationary', 'clean');
grid on;



% frame_idx = 50; 
% figure;
% plot(F, 20*log10(H_smoothed_m1(:,frame_idx) + 1e-10), 'y');
% hold on;
% % plot(F, 20*log10(H_smoothed_m2(:,frame_idx) + 1e-10), 'r');
% % hold on;
% % plot(F, 20*log10(H_smoothed_m3(:,frame_idx) + 1e-10), 'g');
% % title(['Wiener Filter Gain Comparison (Frame ', num2str(frame_idx), ')']);
% % xlabel('Frequency (Hz)');
% % ylabel('Gain (dB)');
% % legend('Method 1', 'Method 2', 'Method 3');
% % grid on;


