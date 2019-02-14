function generate_distinct_color(n){
  return colors[n]
#  n -= 1;
#  h = #TODO generate distinct values
#  h *= 0.9; # Shift color range to exclude circular h values
#  s = 1;
#  v = 1;
#  # Taken from https://gist.github.com/mjackson/5311256
#  i = int(h * 6);
#  f = h * 6 - i;
#  p = v * (1 - s);
#  q = v * (1 - f * s);
#  t = v * (1 - (1 - f) * s);
#
#  switch (i % 6) {
#    case 0: r = v; g = t; b = p; break;
#    case 1: r = q; g = v; b = p; break;
#    case 2: r = p; g = v; b = t; break;
#    case 3: r = p; g = q; b = v; break;
#    case 4: r = t; g = p; b = v; break;
#    case 5: r = v; g = p; b = q; break;
#  }
#
#  return sprintf("#%02X%02X%02X", r * 255, g * 255, b * 255)
}
BEGIN { 
    split("#e6194b #3cb44b #ffe119 #4363d8 #f58231 #911eb4 #46f0f0 #f032e6 #bcf60c #fabebe #008080 #e6beff #9a6324 #fffac8 #800000 #aaffc3 #808000 #ffd8b1 #000075 #808080 #ffffff #000000", colors, " ");
    print "##gff-version 3";
}
(NF > 1) { for (i=1; i<=NF; i++) print gensub(/([^:]+):([0-9]+)-([0-9]+)/, "\\1\tMCL\tgenomic_island\t\\2\t\\3\t.\t.\t.\tcluster=", 1, $i) NR ";color=" generate_distinct_color(NR); }
