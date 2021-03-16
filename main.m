function defineW(a, b, n)
  w = linspace(0,0,n);
  for i = 0:n/2
    w(i) = ((b-a)/(2*n))*i;
    w(n-1-i) = w(i);
  endfor
endfunction

