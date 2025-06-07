% Define two sequences
seq1 = 'AGCTTAGCTA';
seq2 = 'CGTTAGCTAG';

% Define sliding window size
windowSize = 3;

% Initialize the dot plot matrix
len1 = length(seq1);
len2 = length(seq2);
dotPlotMatrix = zeros(len1, len2);

% Sliding window comparison
for i = 1:len1 - windowSize + 1
    for j = 1:len2 - windowSize + 1
        % Extract window segments
        windowSeq1 = seq1(i:i + windowSize - 1);
        windowSeq2 = seq2(j:j + windowSize - 1);
        
        % Check if windows match
        if strcmp(windowSeq1, windowSeq2)
            dotPlotMatrix(i, j) = 1;
        end
    end
end

% Plot the dot plot
figure;
spy(dotPlotMatrix, 'k'); % 'k' for black dots
xlabel('Sequence 2');
ylabel('Sequence 1');
title(['Dot Plot with Sliding Window (Window Size: ', num2str(windowSize), ')']);
