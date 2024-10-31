function noise_est = m1_pure_noise(noisy, fs, noise_duration, window_length, window_overlap, nfft)
    % pure noise segment
    noise_samples = min(round(noise_duration * fs), length(noisy));
    noise_only = noisy(1:noise_samples);
    
    % spectrum in different segment
    [S_noise, ~, ~,~, ~] = win_stft(noise_only, fs, window_length, window_overlap, nfft);
   
    % noise PSD
    noise_est = mean(abs(S_noise).^2, 2);
end
