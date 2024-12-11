function Plot(Cnx, Cnz, Ccz, nnzb, nex, nen, nnx, nnzm, ntex, nezm, ncn, nn, en, pnb, pnc, dx, dz, nc, X, Z, cID, nID, eID)

Cnz = abs(Cnz);
Ccz = abs(Ccz);
dz  = abs(dz);
Z   = abs(Z); 


% #1
figure;                            % Grafica la malla generada, las condiciones de borde y los ID de los nodos

title('Modelo del depósito de suelo discretizado y capas de suelo')
ylabel('Profundidad [ft]'); xlabel('Ancho [ft]');
set(gca,'FontSize',13); 
set(gcf,'color','w');

for i = 1:3
   hold on; grid on; grid minor

   cx = Cnx( pnb(1,i):pnb(1,i+1) );         % El cx obtiene las coordenadas en x de todos los nodos del bloque i
   
   xv = ones(nnzb(1,i),nex(1,i)+1).*cx;     % Las siguientes 4 líneas son para graficar la malla. 
   yv = - Cnz(nen(1,i)+1:end); 
   [r, c] = size(xv);
   xh = xv(1,:);
   yh = ones(1, c).*yv;
   
   plot( xv, yv, '-o', 'Color','b', 'LineWidth',0.5, 'MarkerSize', min(0.4*max(dx), 4) );              % Verticales
   for k = 1:r
     plot( xh, yh(k, :), '-o', 'Color','b', 'LineWidth',0.5, 'MarkerSize', min(0.4*max(dx), 4));      % horizontales
   end
   plot( Cnx, -ones(1,nnx).*Cnz(nnzm,1), '^', 'Color', 'g', 'LineWidth',2, 'MarkerSize', min(0.6*max(dx), 7));        % Apoyos tipo pin 0.7*mean(dx)
   if i == 1
     plot( ones(1,r-1).*Cnx(1,i), yv(1:end-1), '-o', 'Color', 'g', 'LineWidth',2, 'MarkerSize', min(0.6*max(dx), 7));      % Apoyos tipo rollers  %0.7*mean(dx)
   end
   if i == 3
     plot( ones(1,r-1).*Cnx(1,end), yv(1:end-1), '-o', 'Color', 'g', 'LineWidth',2, 'MarkerSize', min(0.6*max(dx), 7) );    % Apoyos tipo rollers 0.7*mean(dx) 
   end
   xlim([-5 sum(X)+5]); ylim([-max(Z)-5 5]); axis equal
   end
%return
  col = ['b','g','r','c','m','y','k','g'];
  col = repmat(col,1,round(nc/8)+1);
   
 hold on

for j = 1:3
  for i = 1:nc
      
  cx = Cnx(pnb(1,j) : pnb(1,j+1))  .* nn( pnc(i,1), pnb(1,j) : pnb(1,j+1));
  cz = [Ccz(i,1), Ccz(i+1,1)];
     if any(cx(2:end)==0) == 1
         cx = 0.*cx;
         cz = 0.*[Ccz(i,1), Ccz(i+1,1)];
     end
     
     if ncn(1,j) >= i
     cx = 0.*cx;
     cz = 0.*[Ccz(i,1), Ccz(i+1,1)];
     end
     
  xBox = [cx(1,1) cx(1, end) cx(1, end) cx(1,1)];
  yBox = -[cz(1, 2) cz(1, 2) cz(1, 1) cz(1, 1)];

  hold on
  fill(xBox, yBox, col(i), 'FaceAlpha', 0.25)
  
  if j == 1
     txt = ['Capa ',  num2str(cID(i,1))];         % Texto para agregar el número de capa 
     text( sum(X)+1, - Ccz(i+1, 1) + 0.5*( Ccz(i+1, 1) - Ccz(i, 1) ), txt, 'FontSize', 8, 'fontweight', 'bold')   % Muestra el número de capa
  end
  
  end
  
end



% #2
figure;                            % Grafica el ID de los nodos

title('Modelo del depósito de suelo discretizado con número de nodos')
ylabel('Profundidad [ft]'); xlabel('Ancho [ft]');
set(gca,'FontSize',13);
set(gcf,'color','w');

for i = 1:3
   hold on; grid on; grid minor

   cx = Cnx( pnb(1,i):pnb(1,i+1) );         % El cx obtiene las coordenadas en x de todos los nodos del bloque i
   
   xv = ones(nnzb(1,i),nex(1,i)+1).*cx;    % Las siguientes 4 líneas son para graficar la malla. 
   yv = - Cnz(nen(1,i)+1:end); 
   [r, c] = size(xv);
   xh = xv(1,:);
   yh = ones(1, c).*yv;
   
   plot( xv, yv, '-o', 'Color','b', 'LineWidth',0.5, 'MarkerSize',min(0.4*max(dx), 4));              % Verticales
   for k = 1:r
     plot( xh, yh(k, :), '-o', 'Color','b', 'LineWidth',0.5, 'MarkerSize',min(0.4*max(dx), 4));      % horizontales
   end
   plot( Cnx, -ones(1,nnx).*Cnz(nnzm,1), '^', 'Color', 'g', 'LineWidth',2, 'MarkerSize',min(0.6*max(dx), 7));        % Apoyos tipo pin
   if i == 1
     plot( ones(1,r-1).*Cnx(1,i), yv(1:end-1), '-o', 'Color', 'g', 'LineWidth',2, 'MarkerSize',min(0.6*max(dx), 7));      % Apoyos tipo rollers
   end
   if i == 3
     plot( ones(1,r-1).*Cnx(1,end), yv(1:end-1), '-o', 'Color', 'g', 'LineWidth',2, 'MarkerSize',min(0.6*max(dx), 7));    % Apoyos tipo rollers
   end
   xlim([-5 sum(X)+5]); ylim([-max(Z)-5 5]); axis equal
   end
%return
  col = ['b','g','r','c','m','y','k','g'];
  col = repmat(col,1,round(nc/8)+1);
   
 hold on

for j = 1:3
  for i = 1:nc
  
  nnb = nn( pnc(i,1), pnb(1,j) : pnb(1,j+1));      % Obtiene el vector de nodos nulos por bloque para cada capa    
  cx = Cnx(pnb(1,j) : pnb(1,j+1))  .* nnb ;        
  cz = [Ccz(i,1), Ccz(i+1,1)];
    
     if any(cx(2:end)==0) == 1
         cx = 0.*cx;
         cz = 0.*[Ccz(i,1), Ccz(i+1,1)];
     end
     
     if ncn(1,j) >= i
     cx = 0.*cx;
     cz = 0.*[Ccz(i,1), Ccz(i+1,1)];
     end
     
  xBox = [cx(1,1) cx(1, end) cx(1, end) cx(1,1)];
  yBox = -[cz(1, 2) cz(1, 2) cz(1, 1) cz(1, 1)];

  hold on
  fill(xBox, yBox, col(i), 'FaceAlpha', 0.25)
  
  end
 
  end

hold on
for i = 1:nnx
    for j = 1:nnzm
lbl = nID(j,i).*nn(j,i);    % obtiene el ID de los nodos por bloque  
      if lbl ~= 0
       txt = num2str(lbl);         % Texto para agregar el número de nodo
       text( Cnx(1,i)+0.05*min(dx), -Cnz(j,1)-0.1*min(dz), txt, 'FontSize', 3, 'fontweight', 'bold')   % Muestra el número de nodo 2*mean(dx)
      end
    end
end
  


% #3
figure;                            % Grafica la malla generada, las condiciones de borde y los ID de los nodos

title('Modelo del depósito de suelo discretizado con números de elementos')
ylabel('Profundidad [ft]'); xlabel('Ancho [ft]');
set(gca,'FontSize',13);
set(gcf,'color','w');

for i = 1:3
   hold on; grid on; grid minor

   cx = Cnx( pnb(1,i):pnb(1,i+1) );         % El cx obtiene las coordenadas en x de todos los nodos del bloque i
   
   xv = ones(nnzb(1,i),nex(1,i)+1).*cx;    % Las siguientes 4 líneas son para graficar la malla. 
   yv = - Cnz(nen(1,i)+1:end); 
   [r, c] = size(xv);
   xh = xv(1,:);
   yh = ones(1, c).*yv;
   
   plot( xv, yv, '-o', 'Color','b', 'LineWidth',0.5, 'MarkerSize', 0.4*mean(dx) );              % Verticales
   for k = 1:r
     plot( xh, yh(k, :), '-o', 'Color','b', 'LineWidth',0.5, 'MarkerSize', 0.4*mean(dx) );      % horizontales
   end
   plot( Cnx, -ones(1,nnx).*Cnz(nnzm,1), '^', 'Color', 'g', 'LineWidth',2, 'MarkerSize', 0.6*mean(dx) );        % Apoyos tipo pin
   if i == 1
     plot( ones(1,r-1).*Cnx(1,i), yv(1:end-1), '-o', 'Color', 'g', 'LineWidth',2, 'MarkerSize', 0.6*mean(dx) );      % Apoyos tipo rollers
   end
   if i == 3
     plot( ones(1,r-1).*Cnx(1,end), yv(1:end-1), '-o', 'Color', 'g', 'LineWidth',2, 'MarkerSize', 0.6*mean(dx));    % Apoyos tipo rollers
   end
   xlim([-5 sum(X)+5]); ylim([-max(Z)-5 5]); axis equal
   end
%return
  col = ['b','g','r','c','m','y','k','g'];
  col = repmat(col,1,round(nc/8)+1);
   
 hold on

for j = 1:3
  for i = 1:nc
  
  nnb = nn( pnc(i,1), pnb(1,j) : pnb(1,j+1));      % Obtiene el vector de nodos nulos por bloque para cada capa    
  cx = Cnx(pnb(1,j) : pnb(1,j+1))  .* nnb ;        
  cz = [Ccz(i,1), Ccz(i+1,1)];
    
     if any(cx(2:end)==0) == 1
         cx = 0.*cx;
         cz = 0.*[Ccz(i,1), Ccz(i+1,1)];
     end
     
     if ncn(1,j) >= i
     cx = 0.*cx;
     cz = 0.*[Ccz(i,1), Ccz(i+1,1)];
     end
     
  xBox = [cx(1,1) cx(1, end) cx(1, end) cx(1,1)];
  yBox = -[cz(1, 2) cz(1, 2) cz(1, 1) cz(1, 1)];

  hold on
  fill(xBox, yBox, col(i), 'FaceAlpha', 0.25)
  
  end
 
  end

hold on
for i = 1:ntex
    for j = 1:nezm
lbl = eID(j,i).*en(j,i);    % obtiene el ID de los nodos por bloque  
      if lbl ~= 0
       txt = num2str(lbl);         % Texto para agregar el número de capa 
       text( Cnx(1,i) + 0.4*( Cnx(1,i+1) - Cnx(1,i) ), -Cnz(j,1) + 0.5*(Cnz(j,1)-Cnz(j+1,1)), txt, 'FontSize', 3, 'fontweight', 'bold')   % Muestra el número de capa 2*mean(dx)
      end
    end
end

end