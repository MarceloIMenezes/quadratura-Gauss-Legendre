function w = defineW(a, b, n)
  w = linspace(0,0,n);
  for i = 0:n/2
    w(i) = ((b-a)/(2*n))*i;
    w(n-1-i) = w(i);
  endfor
endfunction

function t = defineT(a, b, n, w)
  t = linspace(0,0,n);
  for i = 0:n/2
    t(i) = a + i*w(i)/2;
    t(n-1-i) = (a + b) - t(i);
  endfor
  if (rem(n,2) != 0)
    t(n/2) = (a+b)/2;
  endif
endfunction

function soma = calculaSoma(n, w, t)
  soma = linspace(0, 0, 2*n);
  for j = 0:2*n
    for i = 0:n
      soma(j) += w(i)*(t(0)**j);
    endfor
  endfor
endfunction

function g = calculaIntegralX(a, b, m)
  g = linspace(0, 0, 2*n);
  for j = 0:2*n
    f = @(x) x**j;
    g(j) = calculaIntegralTrapezio(a, b, m, f);
  endfor
endfunction

function I = calculaIntegralTrapezio(a, b, m, f)
  h = (b-a)/m;
  x = a;
  I = 0;
  for i = 0:(m+1)
    x = a + i*h;
    if (i == 0 || i == m)
      c = 1;
    else
      c = 2;
    endif
    I += c * f(x); 
  endfor
  I = (h/2)*I;
endfunction
  


