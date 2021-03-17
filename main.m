1
format long;

  [w, t] = metodoNewton(-1, 1, 7, 10e-8);
  #[g] = calculaIntegral(-1, 1, 7, 1000);
  #f = calculaF(w, t, g);
  #soma = calculaSoma(w, t);
  #[J] = montaJacobi(g, w, t, f, e);
  #s = J\(-f);
  F = @(x) (exp(-x + 1));
  r = quadraturaGaussiana(w, t, F);
  save respostas.txt w t F r;

function [w] = defineW(a, b, n)
  w = double(linspace(0,0,n));
  for i = 1:ceil(n/2)
    w(i) = ((b-a)/(2*n))*i;
    w(n+1-i) = w(i);
  endfor
endfunction

function [t] = defineT(a, b, w)
  n = length(w);
  t = double(linspace(0,0,n));
  for i = 1:ceil(n/2)
    t(i) = a + i*w(i)/2;
    t(n+1-i) = (a + b) - t(i);
  endfor
  if (rem(n,2) != 0)
    t(ceil(n/2)) = (a+b)/2;
  endif
endfunction

function [soma] = calculaSoma(w, t)
  n = length(w);
  soma = double(linspace(0, 0, 2*n));
  for j = 1:2*n
    for i = 1:n
      soma(j) += (w(i)).*((t(i)).^(j-1));
    endfor
  endfor
endfunction

function [g] = calculaIntegral(a, b, n, m)
  g = double(linspace(0, 0, 2*n));
  for j = 1:2*n
    y = @(x) x^(j-1);
    g(j) = calculaIntegralTrapezio(a, b, m, y);
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
  n = length(w);
  f = double(zeros(2*n, 1));
  [soma] = calculaSoma(w, t);
  for j = 1:2*n
    f(j) = soma(j) - g(j);
  endfor
endfunction

function [J] = montaJacobi(g, w, t, f, tol)
  n = length(w);
  J = double(zeros(2*n, 2*n));
  waux = double(linspace(0,0,n));
  taux = double(linspace(0,0,n));
  for l = 1:2*n
    for c = 1:(n)
      waux = w;
      waux(c) += tol;
      [fwaux] = calculaF(waux, t, g);
      J(l,c) = (fwaux(l) - f(l))/tol;
      taux = t;
      taux(c) += tol;
      [ftaux] = calculaF(w, taux, g);
      J(l, c + n) = (ftaux(l) - f(l))/tol;
    endfor
  endfor
endfunction

function [w, t] = metodoNewton(a, b, n, tol)
    w = defineW(a, b, n);
    t = defineT(a, b, w);
   [g] = calculaIntegral(a, b, n, 1);
   e = 10e-8;
   s = linspace(0, 0, 2*n);
   while (true)
     [f] = calculaF(w, t, g);
     [J] = montaJacobi(g, w, t, f, e);
     s = J\(-f);
     wprox = linspace(0, 0, n);
     tprox = linspace(0, 0, n);
     for i = 1:2*n
       if i <= n
          wprox(i) = s(i) + w(i);
          w(i) = wprox(i);
        else
          tprox(i-n) = s(i) + t(i-n);
          t(i-n) = tprox(i-n);
        endif
     endfor
     if (norm(s,inf) <= tol)
       break;
     endif
   endwhile
   h = (b-a)/1000;
   save dados.txt g a b n tol e h s;
endfunction

function [r] = quadraturaGaussiana(w, t, F)
    n = length(w);
    r = 0;
    for i = 1:n,
       r += w(i) * F(t(i));
    endfor
endfunction
