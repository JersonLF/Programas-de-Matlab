function [k] = Matriz_kij( pkEin )

x1 = pkEin(1);     x2 = pkEin(3);    x3 = pkEin(5);   x4 = pkEin(7);
y1 = pkEin(2);     y2 = pkEin(4);    y3 = pkEin(6);   y4 = pkEin(8);
E  = pkEin(9);     nu = pkEin(10);

w1 = 1;
w2 = 1;

Xi  = 1/sqrt(3).*[-1 1]'; 
Eta = 1/sqrt(3).*[-1 1]'; 

ip = length(Xi);     % Número de puntos de integración

D = E/( (1+nu)*(1-2*nu) ) .* [1-nu nu 0; nu 1-nu 0 ;0 0 (1-2*nu)/2];          % Modelo constitutivo de deformaciones planas. 


k = zeros(8, 8);
for i = 1:ip  
    for j = 1:ip 
        
         if i == 2; 
         Xi = 1/sqrt(3).*[1,-1]';
         end

    J11 = 1/4.*( -(1-Eta(i,1))*x1 + (1-Eta(i,1))*x2 + (1+Eta(i,1))*x3 - (1+Eta(i,1))*x4);
    J12 = 1/4.*( -(1-Eta(i,1))*y1 + (1-Eta(i,1))*y2 + (1+Eta(i,1))*y3 - (1+Eta(i,1))*y4);
    J21 = 1/4.*( -(1-Xi(j,1))*x1 - (1+Xi(j,1))*x2 + (1+Xi(j,1))*x3 + (1-Xi(j,1))*x4);
    J22 = 1/4.*( -(1-Xi(j,1))*y1 - (1+Xi(j,1))*y2 + (1+Xi(j,1))*y3 + (1-Xi(j,1))*y4);
    
    J = [J11 J12; J21 J22];
    det(J);
    A = 1/det(J) .* [ J22 -J12 0 0 ...
                    ; 0 0 -J21 J11 ...
                    ; -J21 J11 J22 -J12 ];
        
        
    G = 1/4.*[ -1+Eta(i,1) 0 1-Eta(i,1) 0 1+Eta(i,1) 0 -1-Eta(i,1) 0 ...
             ; -1+Xi(j,1) 0 -1-Xi(j,1) 0 1+Xi(j,1) 0 1-Xi(j,1) 0 ...
             ; 0 -1+Eta(i,1) 0 1-Eta(i,1) 0 1+Eta(i,1) 0 -1-Eta(i,1) ...
             ; 0 -1+Xi(j,1) 0 -1-Xi(j,1) 0 1+Xi(j,1) 0 1-Xi(j,1)];
             
    B = A*G;
        
    kij = B.' * D * B * det(J);
    k = k + kij;
  
    end

end

end




