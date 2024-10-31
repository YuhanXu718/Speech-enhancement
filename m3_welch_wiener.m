function [denoised, H_matrix, freq_welch, noisy_psd, noise_psd] = m3_welch_wiener(noisy, noise_duration, fs, window, noverlap, nfft, beta, gain_floor)
    % noisy signal psd
    [noisy_psd, freq_welch] = pwelch(noisy, window, noverlap, nfft, fs);
    
    % noise psd
    noise_samples = min(round(noise_duration * fs), length(noisy));
    noise_only = noisy(1:noise_samples);
    [noise_psd, ~] = pwelch(noise_only, window, noverlap, nfft, fs);
%     figure;
    % plot(freq_welch, 10*log10(noisy_psd), 'b');
    % hold on;
    % plot(freq_welch, 10*log10(noise_psd), 'r');
    % title('psd');
    % xlabel('f');
    % grid on;
    
    % frequency resolution
    delta_f = fs / nfft;
    
    % Power estimate in frequency component
    noisy_est = noisy_psd * delta_f;
    noise_est = noise_psd * delta_f;
    
    % smooth noise power spectrum
    noise_est_smoothed = movmean(noise_est, 5);% 5
    
    % estimate signal power spectrum
    signal_est = max(noisy_est - beta * noise_est_smoothed, 1e-10);
    
    % claculate gain of wiener
    H_wiener = signal_est ./ (signal_est + noise_est_smoothed + eps);
    
    H_wiener = max(min(H_wiener, 1), gain_floor); % limit the gain between [gain_floor, 1]
    
    % calculate the stft of noisy signal
    [S_noisy, F, T] = stft(noisy, fs, 'Window', window, 'OverlapLength', noverlap, 'FFTLength', nfft);

    
    % interpolate wiener gain to match stft frequency
    H_interp = interp1(abs(freq_welch), H_wiener, abs(F), 'pchip', 'extrap');
%     disp(freq_welch(1:5)); 
%     disp(F(1:5));
%     disp(size(H_interp)); 
%     disp(size(S_noisy));

    % apply gain to all signal
    
    H_matrix = repmat(H_interp, 1, size(S_noisy, 2)); %H_interp

     % apply wiener
    S_denoised = H_matrix .* S_noisy;

     % istft
    denoised = istft(S_denoised, fs, 'Window', window, 'OverlapLength', noverlap, 'FFTLength', nfft);
    
    if length(denoised) < length(noisy)
        denoised = [denoised; zeros(length(noisy) - length(denoised), 1)];
    elseif length(denoised) > length(noisy)
        denoised = denoised(1:length(noisy));
        
    end
     denoised = real(denoised);
%Plot the noisy signal and noise-only PSD
% figure;
% plot(freq_welch, 10*log10(noisy_psd), 'b', 'LineWidth', 1.5); % Plot noisy signal PSD in dB
% hold on;
% plot(freq_welch, 10*log10(noise_psd), 'r', 'LineWidth', 1.5); % Plot noise-only PSD in dB
% title('PSD Comparison: Noisy Signal vs Noise-Only');
% xlabel('Frequency/Hz');
% ylabel('Power/Frequency dB/Hz');
% legend('Noisy Signal', 'Noise-Only');
% grid on;
% % Plot the Wiener filter gain curve (H_wiener)
% figure;
% plot(freq_welch, 20*log10(H_wiener), 'g', 'LineWidth', 1.5); % Plot Wiener gain in dB
% title('Wiener Filter Gain Curve');
% xlabel('Frequency (Hz)');
% ylabel('Gain (dB)');
% legend('Wiener Gain');

end


