function filtered_signal = edgeAdjustedAverageFilter2(input_signal, window_size)
    % Check if window size is even
    if mod(window_size, 2) == 0
        error('Window size must be odd');
    end
    
    % Initialize filtered signal
    filtered_signal = zeros(size(input_signal));
    len = length(input_signal);
    half_window_size = (window_size - 1) / 2;

    % Compute the filtered signal
    for i = 1:len
        if i <= half_window_size
            % Near the start edge
            
            filtered_signal(i) = mean(input_signal(1:(2*i-1)));
        elseif i > len - half_window_size
            % Near the end edge
            current_window_size = len-i;
            filtered_signal(i) = mean(input_signal((len-2*current_window_size):len));
        else
            % General case
            filtered_signal(i) = mean(input_signal(i - half_window_size:i + half_window_size));
        end
    end
end
