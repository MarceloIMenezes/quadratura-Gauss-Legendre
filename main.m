function [w] = defineW(a, b, n)
  w = linspace(0,0,n);
  for i = 1:n/2
    w(i) = ((b-a)/(2*n))*i;
    w(n+1-i) = w(i);
  endfor
endfunction

function [t] = defineT(a, b, w)
  n = lenght(w);
  t = linspace(0,0,n);
  for i = 1:n/2
    t(i) = a + i*w(i)/2;
    t(n+1-i) = (a + b) - t(i);
  endfor
  if (rem(n,2) != 0)
    t(n/2) = (a+b)/2;
  endif
endfunction

function [soma] = calculaSoma(w, t)
  n = lenght(w);
  soma = linspace(0, 0, 2*n);
  for j = 1:2*n
    for i = 1:n
      soma(j) += w(i)*(t(0)**j-1);
    endfor
  endfor
endfunction

function [g] = calculaIntegral(a, b, m, f)
  g = linspace(0, 0, 2*n);
  for j = 1:2*n
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

function [f] = calculaF(w, t, g)
  n = lenght(w);
  f = linspace(0, 0, 2*n);
  [soma] = calculaSoma(w, t);
  for j = 1:2*n
    f(j) = soma(i) - g(i);
  endfor
endfunction

function [J] = montaJacobi(g, w, t, f, tol)
  n = lenght(w);
  J = zeros(n, 2*n);
  for l = 1:n
    for c = 1:n
      w(c) += tol;
      waux = linspace(0, 0, 2*n);
      [waux] = calculaF(w, t, g);
      J(l,c) = (waux(l) - f(l))/tol;
      w(c) -= tol;
      t(c) += tol;
      taux = linspace(0, 0, 2*n);
      [taux] = calculaF(w, t, g);
      J(l, c + n) = (taux(l) - f(l))/tol;
      t(c) -= tol;
    endfor
  endfor
endfunction

function [w, t] = metodoNewton(F, a, b, m, tol)
    w = defineW(a, b, m);
    t = defineT(a, b, w);
   [g] = calculaIntegral(a, b, m, f);
   e = power(10, -8);
   while (true)
     [f] = calculaF(w, t, g);
     [J] = montaJacobi(g, w, t, f, e);
     s = J \ (-f);
     wprox = linspace(0, 0, n);
     tprox = linspace(0, 0, n);
     for i = 1:2*n
       if (i <= n)
        wprox(i) = s(i) + w(i);
        w(i) = wprox(i);
       else
        tprox(i) = s(i) + t(i);
        t(i) = tprox(i);
       endif
     endfor
     if (norm(s,inf) <= tol)
       break;
     endif
   endwhile
   h = (b-a)/m;
   save dados.txt F a b m tol w t e h
endfunction
  
  
