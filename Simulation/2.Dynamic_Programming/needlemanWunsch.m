function [alignedSeq1, alignedSeq2, scoreMatrix, tracebackMatrix] = needlemanWunsch(seq1, seq2, match, mismatch, gap)
    % Needleman-Wunsch Algorithm for Sequence Alignment
    % Inputs:
    %   seq1: First sequence (string or character array)
    %   seq2: Second sequence (string or character array)
    %   match: Score for a match
    %   mismatch: Penalty for a mismatch
    %   gap: Penalty for a gap
    % Outputs:
    %   alignedSeq1: Aligned sequence 1
    %   alignedSeq2: Aligned sequence 2
    %   scoreMatrix: Scoring matrix
    %   tracebackMatrix: Traceback matrix
    % Convert sequences to character arrays if they are strings
    seq1 = char(seq1);
    seq2 = char(seq2);
    % Initialize score and traceback matrices
    len1 = length(seq1);
    len2 = length(seq2);
    scoreMatrix = zeros(len1 + 1, len2 + 1);
    tracebackMatrix = zeros(len1 + 1, len2 + 1); % 1: diagonal, 2: up, 3: left
    % Fill the first row and column with gap penalties
    for i = 1:len1 + 1
        scoreMatrix(i, 1) = (i - 1) * gap;
        tracebackMatrix(i, 1) = 2; % Up
    end
    for j = 1:len2 + 1
        scoreMatrix(1, j) = (j - 1) * gap;
        tracebackMatrix(1, j) = 3; % Left
    end
    tracebackMatrix(1, 1) = 0; % No direction at the origin
    % Fill the rest of the score matrix
    for i = 2:len1 + 1
        for j = 2:len2 + 1
            matchScore = scoreMatrix(i - 1, j - 1) + (seq1(i - 1) == seq2(j - 1)) * match + ...
                         (seq1(i - 1) ~= seq2(j - 1)) * mismatch;
            gapSeq1 = scoreMatrix(i - 1, j) + gap;
            gapSeq2 = scoreMatrix(i, j - 1) + gap;
            % Choose the maximum score and determine the traceback direction
            [scoreMatrix(i, j), tracebackMatrix(i, j)] = max([matchScore, gapSeq1, gapSeq2]);
        end
    end
    % Traceback to find the optimal alignment
    alignedSeq1 = '';
    alignedSeq2 = '';
    i = len1 + 1;
    j = len2 + 1;
    while i > 1 || j > 1
        if tracebackMatrix(i, j) == 1 % Diagonal
            alignedSeq1 = [seq1(i - 1), alignedSeq1];
            alignedSeq2 = [seq2(j - 1), alignedSeq2];
            i = i - 1;
            j = j - 1;
        elseif tracebackMatrix(i, j) == 2 % Up
            alignedSeq1 = [seq1(i - 1), alignedSeq1];
            alignedSeq2 = ['-', alignedSeq2];
            i = i - 1;
        elseif tracebackMatrix(i, j) == 3 % Left
            alignedSeq1 = ['-', alignedSeq1];
            alignedSeq2 = [seq2(j - 1), alignedSeq2];
            j = j - 1;
        end
    end
    % Display results
    disp('Score Matrix:');
    disp(scoreMatrix);
    disp('Aligned Sequences:');
    disp(alignedSeq1);
    disp(alignedSeq2);
    disp('traceback Matrix:');
    disp(tracebackMatrix);
end