function [denoised, H_smoothed] = wiener(S_noisy, noisy, fs, window, noverlap, nfft, noise_est, alpha, gain_floor, smoothing_length)
    % noise estimation matrix
    noise_est_matrix = repmat(noise_est, 1, size(S_noisy, 2));
    
    % estimate signal power spectrum 
    signal_est = max(abs(S_noisy).^2 - alpha * noise_est_matrix, 1e-10);
    
    % wiener gain
    H = signal_est ./ (signal_est + noise_est_matrix + eps);
    
    % limit gain
    H = max(min(H, 1), 0);
    
    % smoothing moving average
    H_smoothed = movmean(H, smoothing_length, 2);
    H_smoothed = max(min(H_smoothed, 1), 0);
    
    % limit gain via gain floor
    H_smoothed = max(H_smoothed, gain_floor);
    
    % apply wiener
    S_denoised = H_smoothed .* S_noisy;
    
    % istft
    denoised = istft(S_denoised, fs, 'Window', window, 'OverlapLength', noverlap, 'FFTLength', nfft);
    
    % adjust length of denoised signal to original length
    if length(denoised) < length(noisy)
        denoised = [denoised; zeros(length(noisy) - length(denoised), 1)];
    elseif length(denoised) > length(noisy)
        denoised = denoised(1:length(noisy));
    end
    denoised = real(denoised);
end
