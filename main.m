function defineW(a, b, n)
  w = linspace(0,0,n);
  for i = 0:n/2
    w(i) = ((b-a)/(2*n))*i;
    w(n-1-i) = w(i);
  endfor
endfunction

function defineT(a, b, n, w)
  t = linspace(0,0,n);
  for i = 0:n/2
    t(i) = a + i*w(i)/2;
    t(n-1-i) = (a + b) - t(i);
  endfor
endfunction


