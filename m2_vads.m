function noise_est = m2_vads(S, energy_threshold)
  
    % select noise frame  
    frame_energy = sum(abs(S).^2, 1);
    noise_frames = frame_energy < energy_threshold;
    
    % check if write
    if sum(noise_frames) < 1
        error('No frames');
    end
    
    % Extract all frames marked as noise
    S_noise = S(:, noise_frames);
    
    % noise PSD
    noise_est = mean(abs(S_noise).^2, 2);
end


