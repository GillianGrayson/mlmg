clear all;

N = 15000;
K = 150;

p = K/N;
q = 1 - p;

n = K;
k = 1;


sum = 0.0;
for i = 1:k
    bin_coeff = nchoosek(n, i);
    p_deg = power(p, i);
    q_deg = power(q, n-i);
    sum = sum + bin_coeff * p_deg * q_deg;
end

prob = 1 - sum;