function [clean, noise, noise_sta, fs] = load_audio(fs_clean, fs_noise, fs_noise_sta,clean, noise, noise_sta)
    
    % check sampling rate
    if fs_clean ~= fs_noise || fs_clean ~= fs_noise_sta
        error('sampling rate is different');
    end
    fs = fs_clean;
    
    % check the length and set to the same
    if length(noise) < length(clean)
        num_repeats = ceil(length(clean) / length(noise)); % #repeat
        noise = repmat(noise, num_repeats, 1);    % repeat
        noise = noise(1:length(clean));               % cut to the same length
    else
        noise = noise(1:length(clean));
    end

    if length(noise_sta) < length(clean)
        num_repeatss = ceil(length(clean) / length(noise_sta)); % #repeat, up 
        noise_sta = repmat(noise_sta, num_repeatss, 1);    % repeat
        noise_sta = noise_sta(1:length(clean));               % cut to the same length
    else
        noise_sta = noise_sta(1:length(clean));
    end

end

