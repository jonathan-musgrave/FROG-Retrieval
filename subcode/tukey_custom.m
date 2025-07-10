function w = tukey_custom(N, alpha, n, m)
    % N: Total length of the output window
    % alpha: Shape parameter of the Tukey window (0 -> rectangular, 1 -> Hann)
    % n: Index where the window is centered
    % m: The number of zeros on each side of the window (total taper = 2*m)

    % Error checking
    if m < 0 || 2*m >= N
        error('Invalid value for m. Ensure 2*m < N.');
    end

    % Calculate the width of the non-zero window (N - 2*m)
    width = N - 2*m;

    % Create a Tukey window of width (N - 2*m)
    tukey_win = tukeywin(width, alpha);

    % Create the full window with N length and zero padding
    w = zeros(1, N);

    % Calculate the start and end indices for the non-zero portion
    start_idx = max(1, n - floor(width/2));
    end_idx = min(N, start_idx + width - 1);

    % Insert the Tukey window into the full window
    w(start_idx:end_idx) = tukey_win;
end
