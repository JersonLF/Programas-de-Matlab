function [m] = Matriz_mij( pmEin )

x1 = pmEin(1);     x2 = pmEin(3);    x3 = pmEin(5);   x4 = pmEin(7);
y1 = pmEin(2);     y2 = pmEin(4);    y3 = pmEin(6);   y4 = pmEin(8);
gam= pmEin(9);     g  = pmEin(10);


Ae = (x2-x1)*(y4-y1);    % Area del elemento finito i 
Ve = Ae;                 % Volumen del elemento finito i
we = gam*Ve;             % Peso del elemento finito i
me = we/g;               % Masa del elemento finito i 
mn = me/4.*ones(1,8);    % Masa tributaria por nodo
m = diag(mn);

end
 



