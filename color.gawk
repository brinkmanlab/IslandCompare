# For an infinite discreet series n we want to map to a range [0, 1] such that the values are as distant (distinct)
# as possible without recalculating the predecessors.
# This can be done in O(1) time via recursive subdivision of the range.
#   n  = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17 ...
# d(n) = 0, 1, 0.5, 0.25, 0.75, 0.125, 0.375 ...
#      = 0/1, 1/1, 1/2, 1/4, 3/4, 1/8, 3/8, 5/8, 7/8, 1/16, 3/16, ...
# You can see the numerator and denominator are increasing with a logaritmic pattern.
# Because we are recursively subdividing by 2 it is ceil(log_2(n)).
# d(n) = g(n)/2^f(n)
# g(0) = 0
# f(0) = 0
#
# f(n) = ceil(log_2(n))
#      = 0, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5
#2^f(n)= 1, 1, 2, 4, 4, 8, 8, 8, 8, 16, 16, 16, 16, 16, 16, 16, 16, 32
# 
# The numerator is a little more complicated as it is the integer remainder of log_2(n) or n - 2^floor(log_2(n))
# n - 2^(floor(log_2(n))) = 0, 0, 0, 1, 0, 1, 0, 1, 2, 0, 1, 2, 3, ...
# This series is (x-1)/2 away from the series we want so we adjust by *2+1
# g(n) = (n - 2^(floor(log_2(n)))) * 2 + 1
#      = 0, 1, 1, 1, 3, 1, 3, 5, 7, 1, 3, 5, 7, 9, 11, 13, 15, 1, ...
#
# It is worth noting that floor(log_2()) can be optimised using bitwise operations and intermediate values can be
# exploited to avoid the 2^(log_2()) sillyness.

function generate_distinct_color(n){
  n -= 1; # Clusters are 1-indexed, shift to 0-index
  if (n > 0) {
    #calculate log2 n TODO implement using integer log2 as it is more efficient
    log2n=log(n)/log(2);
    floor=int(log2n);
    ceil=(log2n>floor?floor+1:floor);
    h=((n - 2^floor) * 2 + 1) / (2^ceil);
    h *= 0.9; # Shift color range to exclude circular h values
  } else {
    h=0;
  }
  s = 1;
  v = 1;
  # The following maps into a hsv color gradient
  # Taken from https://gist.github.com/mjackson/5311256
  i = int(h * 6);
  f = h * 6 - i;
  p = v * (1 - s);
  q = v * (1 - f * s);
  t = v * (1 - (1 - f) * s);

  switch (i % 6) {
    case 0: r = v; g = t; b = p; break;
    case 1: r = q; g = v; b = p; break;
    case 2: r = p; g = v; b = t; break;
    case 3: r = p; g = q; b = v; break;
    case 4: r = t; g = p; b = v; break;
    case 5: r = v; g = p; b = q; break;
  }

  return sprintf("#%02X%02X%02X", r * 255, g * 255, b * 255)
}
