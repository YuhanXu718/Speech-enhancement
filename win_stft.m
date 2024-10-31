function [S, F, T, window, noverlap] = win_stft(signal, fs, window_length, window_overlap, nfft)
    window = hamming(round(window_length * fs), 'periodic');
    noverlap = round(window_overlap * fs);
    [S, F, T] = stft(signal, fs, 'Window', window, 'OverlapLength', noverlap, 'FFTLength', nfft);
end
