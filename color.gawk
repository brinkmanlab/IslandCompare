function generate_distinct_color(n){
  n -= 1;
  #calculate log2 n
  log2n=log(n)/log(2)
  floor=int(log2n)
  ceil=(log2n>floor?floor+1:floor)
  h=((n - 2^ceil) * 2 + 1) / ceil
  h *= 0.9; # Shift color range to exclude circular h values
  s = 1;
  v = 1;
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