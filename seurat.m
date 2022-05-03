% seurat.m
%
% Joshua Pearson 2022
% University College, OXFORD
% joshua.pearson@univ.ox.ac.uk
%
% Finds calculates error expected in a pi-estimator based on dots

% CALCULATIONS
% Values of n to calculate for
n_arr = [10:10:1000];
% Probability each dot is in circle
p = pi/4;
% Array to save output to
n_err = zeros(size(n_arr));
% Count of times iterated through
cnt = 1;

% 1. for n = range of values
for n = n_arr
    % 2. Create array k from 0 to n
    k = [0:n];
    % 3. Calculate probability of total dots in circle, P(X=k)
    % (n Cr k) * p^k * (1-p)^(n-k)
    % Need to use built-in MATLAB functions or the numbers get too big
    % (factorial(200) gives NaN)
    %probX = nchoosek(n,k) .* (p.^k) .* (1-p).^(n-k);
    probX = binopdf(k,n,p);
    % 4. Calculate probability of different values of pi by dividing this by n
    probPi = probX ./ n;
    % 5. Calculate array of difference between these values and pi
    abDiffs = abs(k - pi);
    % 6. Find the average of the absolute value, the error
    avg_error = (sum(abDiffs.*probPi))/n;
    % 7. Save error to array
    n_err(cnt) = avg_error;
    cnt = cnt + 1;
    % End for loop
end

% PLOTS
% 8. Plot error against n
% n on x-axis, error on y
plot(n_arr, n_err)
title('Expected error in measurement of pi')
xlabel('Number of dots, n')
ylabel('Expected error')
xlim([0, max(n_arr)])
ylim([0, max(n_err)])