function [noisy, noisy_sta] = generate_noisy(clean, noise, noise_sta, desired_SNR)
    % calculate power of signal and noise
    clean_power = sum(clean.^2) / length(clean);
    noise_power = sum(noise.^2) / length(noise);
    noise_power_sta = sum(noise_sta.^2) / length(noise_sta);
    
    % scaling factor: to meet the required snr
    scaling_factor = sqrt(clean_power / (10^(desired_SNR/10) * noise_power));
    scaling_factor_sta = sqrt(clean_power / (10^(desired_SNR/10) * noise_power_sta));
    
    % generate noisy signal 
    noisy = clean + noise * scaling_factor;
    noisy_sta = clean + noise_sta * scaling_factor_sta;
end
