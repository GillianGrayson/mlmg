import scipy.stats as ss

n = 150
p = 150.0 / 15000.0
max_bets = 1

hh = ss.binom(n, p)

total_p = 0
for k in range(1, max_bets + 1):  # DO NOT FORGET THAT THE LAST INDEX IS NOT USED
    tmp = hh.pmf(k)
    total_p += hh.pmf(k)

print(1-total_p)