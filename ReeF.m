%-------------------------------- REDES.m -----------------------------%
% ------------------- Versión inicial: 14-Julio-2022 ---------------------%
% Programa para calcular la respuesta sísmica de un depósito de suelo en 2D
% formado por "n" capas. El programa utiliza elementos finitos tipo Q4 y el
% modelo de deformaciones planas. El sismo de entrada puede ser aplicado en
% la base o en la superficie.
% --------------- Última actualización: 10-Diciembre-2024 -----------------%
% Nota: Todos los datos ingresados deben tener unidades consistentes

clc; clearvars; close all; clear all
tic

%---------------------------- Datos de entrada ---------------------------% 
%-------------------------- Definición del modelo ------------------------%


nc  = 2;                    % Número de capas de suelo.
hc  = [30 30];              % Espesor de cada capa empezando desde la superficie.
zi  = [0.05 0.05];          % Razones de amortiguamiento de las capas de suelo.
nec = [15 15];              % Define el número de EF en cada capa.
nex = [10 10 10];           % Define el número de EF a lo largo de cada bloque.
c   = [600 1200];           % Velocidad de ondas de corte por capa.
gam = [0.10 0.11];          % Pesos unitarios de los suelos.
nu  = [0.30 0.30];          % Razones de Poisson de los suelos.
ncb = [2 2 2];              % Número de capas por bloque.
X   = [20 20 20];           % Ancho de los bloques (la suma es igual al ancho del depósito).
g   = 32.2;                 % Aceleración de la gravedad.       
sv = 'NO';                  % SI o NO: Define si se guardan o no las Funciones de Transferencia en la superficie.
sr = 'NO';                  % SI o NO: Define si se guardan los registros en la superficie.
dr = [0 20];                % Define las distancias (desde el extremo derecho del bloque 1) a las que se desea calcular la respuesta.
op = 'C';                   % Define la operación a realizar 'C': Convolución, 'D': Deconvolución.        
PGAs = 0.0;                 % Máxima aceleracion para escalar el sismo de entrada, 0 si no se desea escalar



%--------------------------- Datos de entrada ----------------------------%
% ------- Definición de cargas sismícas y vector de frecuencias y --------%
% --------------- direcciones para leer y guardar archivos ---------------%

nom  = 'Eco2Hx';	        % Nombre del archivo con el acelerograma
dt   = 0.005;		        % Intervalo de tiempo de muestreo: seg
omm = 300;                  % Frecuencia máxima para graficar

f_in  = 'C:\Users\Jerson LF\OneDrive - University of Puerto Rico\Desktop\UPRM\Tesis\Registros sismicos y espectros\Registros Eco-Electrica\Registros originales';    % Carpeta con registro de aceleraciones de entrada
f_out = 'C:\Users\Jerson LF\OneDrive - University of Puerto Rico\Desktop\UPRM\Tesis\Registros sismicos y espectros\Registros Eco-Electrica\Registros originales';                                                    % Carpeta para guardar los archivos de salida
addpath(f_out); addpath(f_in);


%------------ Genera vectores con las propiedades de cada capa -----------%
 
Gci = c.^2 .* gam./g;         % Módulo de corte de las capas de suelo 
Eci = Gci .* 2.*(1+nu);       % Módulo de elasticidad de las capas.
E = (1+2*1i.*zi).*Eci;        % Cálcula los módulos de elasticidad complejos


%-------- Verifica que no haya ninguna inconsistencia en el modelo -------%

if 6*nc ~= length(hc)+length(nec)+length(E)+length(gam)+length(nu)+length(zi)
    disp('Hay una inconsistencia entre el número de capas definidas y las propiedades de las mismas')
    return
end

if any(ncb == 0) == 1 || max(ncb)>nc || max(ncb)<nc
    disp('Hay una inconsistencia en la definición del número de capas por bloque')
    return
end

%---------- Genera vectores con las propiedades de los elementos ---------%

ncn = nc - ncb;                          % Calcula el número de capas nulas (a borrar) en cada bloque

for i = 1:3                          
     nez(i) = sum(nec(ncn(i)+1: end));   % Número de elementos por bloque horizontal
     Z(i)   = sum(hc(ncn(i)+1: end));    % Altura de los bloques horizontales
end

neo = sum(nec)*sum(nex);     % Calcula el número total de elementos del modelo (incluyendo nulos)
nen = sum(nec)-nez;          % Calcula el número de elementos nulos por bloque horizontal
dxb = X./nex;                % Dimensión horizontal de los elementos de cada bloque: ft
dec = hc./nec;               % Dimensión vertical de los elementos por capa: ft
nnx = sum(nex)+1;            % Número de nodos totales en la dirección x
nnzb = nez+1;                % Número de nodos por bloque en la dirección z


%-- Guarda las propiedades de los EF del modelo en la dirección vertical--%

dz  = [];       % Altura del elemento i
Ei  = [];       % Módulo de elasticidad del elemento i
nui = [];       % Razón de Poisson del elemento i
gmi = [];       % Peso volumetrico del elemento i

for i = 1:nc
   dv = dec(i).*ones(1,nec(1,i));   
   dz = horzcat(dz,dv);       
   
   Ec = E(i)*ones(1,nec(1,i));
   Ei = horzcat(Ei,Ec);       
   
   nuc = nu(i)*ones(1,nec(1,i));
   nui = horzcat(nui,nuc);   
   
   gmc = gam(i)*ones(1,nec(1,i));
   gmi = horzcat(gmi,gmc);
end

%------------------------ Convierte filas a columnas ---------------------%

dz = dz';       % Vector con dimensiones de los elementos en la dirección vertical
Ei = Ei';       % Vector con módulos de elasticidad de los elementos en la dirección vertical
nui = nui';     % Vector con razones de Poisson de los elementos en la dirección vertical
gmi = gmi';     % Vector con densidades de los elementos en la dirección vertical


%------------------ Guarda el ancho de los EF por bloque  ----------------%

dxb1 = dxb(1,1) .* ones([1, nex(1,1) ]);   % Vector con ancho de los elementos del bloque 1
dxb2 = dxb(1,2) .* ones([1, nex(1,2) ]);   % Vector con ancho de los elementos del bloque 2
dxb3 = dxb(1,3) .* ones([1, nex(1,3) ]);   % Vector con ancho de los elementos del bloque 3
dx = [dxb1 dxb2 dxb3];                     % Vector con ancho de los elementos del modelo completo

%------------- Calcula y guarda la coordenada de los elementos -----------%

Cnx = 0;                                % Coordenada incial en x
for j = 2:nnx
    Cnx = [ Cnx, sum(dx(1, 1:j-1))];    % Guarda las coordenadas de todos los EF (nodos) en la dirección horizontal
end

Cnz = 0;                                % Coordenada incial en z
for j = 2:max(nnzb)
    Cnz = [ Cnz; sum(dz(1:j-1, 1))];    % Guarda las coordenadas de todos los EF (nodos) en la dirección vertical
end

Ccz = 0;                                % Coordenada incial en z
for j = 2:nc+1
    Ccz = [ Ccz; sum(hc(1, 1:j-1))];    % Guarda las coordenadas de las capas en la dirección vertical
end

nnzm = max(nnzb);                       % Número de nodos máximo en la dirección z 
nezm = max(nez);                        % Número de elementos máximo en la dirección z 
ntex = sum(nex);                        % Número de elementos en la dirección x 

nIDt = reshape( 1:nnzm*nnx, nnzm, nnx );      % Genera una matriz con el orden de los nodos
eIDt = reshape( 1:nezm*ntex, nezm, ntex );    % Genera una matriz con el orden de los EF
cID = reshape( 1:nc, nc, 1 );                 % Genera una matriz con el orden de las capas


%-------- Loop para guardar propiedades generales de los elementos -------%

nnpet = zeros(ntex * nezm, 4);        % Matriz para guardar los 4 nodos de cada elemento, ordenados en filas  
GdLet = zeros(ntex * nezm, 8);        % Matriz para guardar los 8 GdL por elemento, ordenados en filas 
pkEi   =  zeros(ntex * nezm, 10);     % Matriz para guardar las propiedades de rigidez del elemento i
pkEi   =  zeros(ntex * nezm, 10);     % Matriz para guardar las propiedades de masa del elemento i
fi = 1;                               % Para rastrear el índice de la fila en GdLe
nz = 0;

% En este loop se definen las propiedades de todos los elementos finitos (conectividades, materiales y geometría).

for j = 1:ntex  % j = 1:nex
     for i = 1:nezm   % for i = 1:nez
        
       k = (j - 1) * nezm + i;
       pkEi(k, :) = [Cnx(j), -Cnz(i+1), Cnx(j+1), -Cnz(i+1), Cnx(j+1), -Cnz(i), Cnx(j), -Cnz(i), Ei(i), nui(i)];     % propiedad del elemento k
       pmEi(k, :) = [Cnx(j), -Cnz(i+1), Cnx(j+1), -Cnz(i+1), Cnx(j+1), -Cnz(i), Cnx(j), -Cnz(i), gmi(i), g];         % propiedad del elemento k
       Elem = nezm*(j-1)+i;                % Sigue la secuencia de los elementos
       [row,col] = find(eIDt==Elem);       % Obtiene la fila y la columna del elemento Elem en la matriz eID

       nodos = [ nnzm*(col-1)+row     nnzm*(col-1)+row+nnzm...
                nnzm*(col-1)+row+1+nnzm     nnzm*(col-1)+row+1 ];

        nodos = [ nnzm*(col-1)+row+1        nnzm*(col-1)+row+1+nnzm...
                  nnzm*(col-1)+row+nnzm     nnzm*(col-1)+row];
            
       % Obtiene los grados de libertad (x,y) globales desde el nodo 1 al 4 del elemento i 
       GdLj  = [ nodos(1,1)*2-1 nodos(1,1)*2        nodos(1,2)*2-1  nodos(1,2)*2 ...
                 nodos(1,3)*2-1 nodos(1,3)*2        nodos(1,4)*2-1  nodos(1,4)*2 ];  
       
       nnpet(fi, :) = nodos;      % Guarda los nodos de cada elemento
       GdLet(fi, :) = GdLj;       % Guarda los GdL de cada elemento

     % Incrementar el índice de resultados
        fi = fi + 1;

     end
end


%------------------- Identifica nodos y elementos nulos ------------------%

[en, nn] = Elementos_nulos(nez, nex, nen);  % Obtiene matrices con ceros y unos dependiendo si el elemento o nodo ha sido eliminado 
  
ven =  reshape(en, [], 1);                  % Convierte la matriz de nodos nulos a un vector 
vnn =  reshape(nn, [], 1);                  % Convierte la matriz de nodos nulos a un vector 
nte = sum(ven);                             % Número total de elementos (excluyendo los nulos)
ntn = sum(vnn);                             % Calcula el número total de nodos del modelo
enn = ven(ven ~= 0);                        % Vector con elementos no-nulos (para re-enumerar)
nnn = vnn(vnn ~= 0);                        % Vector con nodos no-nulos (para re-enumerar)
nGdL   = 2*ntn;                             % Número de grados de libertad.


%----------- Elimina propiedades el elementos nulos y reenumera ----------%

pkEi = pkEi.*ven;          % Elimina las propiedades de los elementos nulos.
fz = any(pkEi , 2);        % Elimina las filas que únicamente tienen ceros.
pkEin  = pkEi (fz , :);    % Guarda las propiedades de los elementos existentes.

pmEi = pmEi.*ven;          % Elimina las propiedades de los elementos nulos.
fz = any(pmEi , 2);        % Elimina las filas que únicamente tienen ceros.
pmEin  = pmEi (fz , :);    % Guarda las propiedades de los elementos existentes.


%-------------------------- Reenumerando los nodos -----------------------%

nIDt = nIDt.*nn; 
u1 = sort(unique(nIDt(nIDt ~= 0)));
nID =  0.*(nIDt);

for i = 1:length(u1)
    nID(nIDt == u1(i)) = i;
end

fz = any(nID, 2); 
nID = nID(fz , :);     % Redefine el número de los elementos


%------------------------ Reenumerando los elementos ---------------------%

eIDt = eIDt.*en; 
u1 = sort(unique(eIDt(eIDt ~= 0)));
eID =  0.*(eIDt);

for i = 1:length(u1)
    eID(eIDt == u1(i)) = i;
end

fz = any(eID, 2); 
eID = eID(fz , :);      % Redefine el número de los elementos
eIDb = eID(end,:);      % Obtiene los números de los elementos en la base



%----------- Reenumerando los números de nodos de los elementos ----------%

nnpet = nnpet.*ven;
u1 = sort(unique(nnpet(nnpet ~= 0)));
nnpe =  0.*(nnpet);

for i = 1:length(u1)
    nnpe(nnpet == u1(i)) = i;
end

fz = any(nnpe, 2); 
nnpe = nnpe(fz , :);     % Redefine el número de nodos por elemento



%------------------- Reenumerando los GdL de los elementos ---------------%

GdLet = GdLet.*ven;
u1 = sort(unique(GdLet(GdLet ~= 0)));
GdLe = 0.*GdLet;

for i = 1:length(u1)
    GdLe(GdLet == u1(i)) = i;
end

fz = any(GdLe, 2);
GdLe = GdLe(fz , :);            % Redefine el número de GdL por elemento
nGdL   = 2*ntn;                 % Recalcula el número de grados de libertad
  
Kg = zeros(nGdL, nGdL);         % Matriz de ceros para ensamblar la matriz de rigidez global    
Mg = zeros(nGdL, nGdL);         % Matriz de ceros para ensamblar la matriz de masa global


%----------- Calcula y ensambla las matrices de masa y rigidez -----------%

for j = 1:nte  % j = 1:nex 
   if ~ismember(j, eIDb)                          % Este if es para no ensamblar los elementos de la base
     [k] = Matriz_kij( pkEin(j,:) );              % Calcula la matriz de rigidez del elemento j
     [m] = Matriz_mij( pmEin(j,:) );              % Calcula la matriz de rigidez del elemento j
     GdLj = GdLe(j, :);                   
   
      for f = 1:8                                 % Filas de la matriz de rigidez del elemento i   
        fi = GdLj(1,f);                           % Fila i de la matriz de rigidez i
          for q = 1:8                             % Columnas de la matriz de rigidez del elemento i 
            ci = GdLj(1,q);                       % Columna i de la matriz de rigidez i
            Kg(fi,ci) = Kg(fi,ci) + k(f,q);       % Ensambla la matriz de rigidez K
            Mg(fi,ci) = Mg(fi,ci) + m(f,q);       % Ensambla la matriz de masa M
          end                         
      end
   end
end     
   
disp(['Matrices de rigidez y de masa calculadas: ', num2str(nte) ]); disp(' ');


cni = nonzeros( nID(1:end-1,1) );           % Columna con nodos del extremo izquierdo del modelo
cnd = nonzeros( nID(1:end-1,end) );         % Columna con nodos del extremo derecho del modelo
cnb = nonzeros( nID(end,:) )';              % Fila con los nodos de la base 
cnv = nonzeros( vertcat(cni, cnd) );        % Columna que contiene los nodos de los extremos del modelo



%--------------------- Obtiene los G.d.L restringidos --------------------%

GdLry = vertcat(cnv*2, 2*cnb');             % Grados de libertad restringidos verticalmente 
GdLry = unique(GdLry);                      % Elimina los repetidos GdL repetidos           
GdLrx = 2*cnb' - 1;                         % Grados de libertad restringidos horizontalmente
GdLrt = sort( vertcat(GdLry, GdLrx) );      % Grados de libertad restringidos (totales)



%----------- Aplica condiciones de borde a las matrices M y K ------------%


Kg(:, GdLrt) = [];                     % Elimina los GdL restringidos
Kg(GdLrt,:) = [];                      % Elimina los GdL restringidos

Mg(:, GdLrt) = [];                     % Elimina los GdL restringidos
Mg(GdLrt,:) = [];                      % Elimina los GdL restringidos


GdLn = GdLe;                           % Esto es para reenumerar los GdL

for i = 1:length(GdLrt)
    GdLn(GdLn == GdLrt(i)) = 0;
end


%------------------- Reenumerando los GdL de los elementos ---------------%

u1 = sort(unique(GdLn(GdLn ~= 0)));
GdLe = 0.*GdLn;

for i = 1:length(u1)
    GdLe(GdLn == u1(i)) = i;
end

clear GdLn 
 


GdLe_b = GdLe(eIDb, :);    % Obtiene los G.d.L de los elementos de la base
pkEib  = pkEin(eIDb, :);   % Obtiene las propiedades de los elementos de la base
pmEib  = pmEin(eIDb, :);   % Obtiene las propiedades de los elementos de la base


%-------- Calcula y ensambla las propiedades de masa y rigidez de --------%
%------------------------ los elementos de la base -----------------------%   

for j = 1:ntex 
   [k] = Matriz_kij( pkEib(j,:) );       % Calcula la matriz de rigidez del elemento j
   [m] = Matriz_mij( pmEib(j,:) );       % Calcula la matriz de rigidez del elemento j
   GdLj = GdLe_b(j, :);                   
   
     for f = 1:8                                 % Filas de la matriz de rigidez del elemento i   
        fi = GdLj(1,f);                          % Fila i de la matriz de rigidez i
          if fi ~= 0
          
            for q = 1:8                          % Columnas de la matriz de rigidez del elemento i 
              ci = GdLj(1,q);                    % Columna i de la matriz de rigidez i
               
              if ci ~=0
                 Kg(fi,ci) = Kg(fi,ci) + k(f,q);       % Ensambla la matriz de rigidez K
                 Mg(fi,ci) = Mg(fi,ci) + m(f,q);       % Ensambla la matriz de masa M
              end
            
            end
          end
     end
end   

[nGdL, nGdL] = size(Kg); 

for i = 1:size(eID, 2)             % Obtiene el número de los elementos 
    col = eID(:, i);               % en la superficie
    idx = find(col ~= 0, 1);  
    if ~isempty(idx)
        esup(i) = col(idx);  
    end
end

GdLx = vertcat( GdLe(:, 1), GdLe(:, 3), GdLe(:, 5), GdLe(:, 7) );     % Obtiene los GdL horizontales
GdLx = sort( nonzeros( unique (GdLx) ) );
GdLy = vertcat( GdLe(:, 2), GdLe(:, 4), GdLe(:, 6), GdLe(:, 8) );     % Obtiene los GdL verticales
GdLy = sort( nonzeros( unique (GdLy) ) );
GdLe_s = GdLe(esup, :);                                               % Obtiene los GdL de los elementos en la superficie 
GdLx_s = vertcat( GdLe_s(:, 5) , GdLe_s(:, 7) );
GdLx_s = unique(GdLx_s)';

%---------------------------- Grafica el modelo --------------------------%

for i = 1:3
  seb(i) = sum(nex(1:i))+1;         % Calcula la suma de los elementos por bloque
end
seb = horzcat(1, seb);


if nc == 1                         
    i = 1;
    sec = sum(nec(1:i))+1;
else
  for i = 1:nc-1
    sec(i) = sum(nec(1:i))+1;       % Calcula el elemento en la parte superior de cada capa 
  end
end

sec = vertcat(1, sec');

%----- Grafica el perfil del depósito de suelos y la malla generada ------%

Plot(Cnx, Cnz, Ccz, nnzb, nex, nen, nnx, nnzm, ntex, nezm, ncn, nn, en, seb, sec, dx, dz, nc, X, Z, cID, nID, eID)

Kg = conj(Kg);
Mgs = sparse(Mg);
Kgs = sparse(Kg);
r = ones(nGdL,1);        % Genera el vector de coeficientes de influencia
r(GdLy) = 0;             % Elimina los GdL restringidos
r = sparse(r);

%---------------- Genera un vector de frecuencias discretas --------------%

[om, dw, xgi, omc, ncut1] = F_frec_dis(nom, dt, omm, PGAs, g);
lw = length(om);

%----------------- Calcula las funciones de transferencia ----------------%

disp(' ');
disp('Calculando las Funciones de Transferencia'); disp('  ');


B =  - Mgs * r;
for i = 1:lw;
     if om(i) <= omc 
        A = Kgs - ( om(i).^2 .* Mgs );
        H(:,i) = -om(i)^2 .* (A \ B) + 1;
     
        if i == round(0.2*length( om(1:ncut1) ) )
         disp('20% calculadas'); disp('  ')
        end
             
        if i == round(0.4*length( om(1:ncut1) ) )
         disp('40% calculadas'); disp('     ')
        end  

        if i == round(0.6*length( om(1:ncut1) ) )
         disp('60% calculadas'); disp('     ')
        end  
     
        if i == round(0.8*length( om(1:ncut1) ) )
         disp('80% calculadas'); disp('     ')
        end 

     else
      H(:,i) = 0;  
     end
     
end

eIDb1 = eID(:, 1:nex(1,1) ) ;      % Obtiene todos los elementos del bloque 1 
eIDb1 = eIDb1(any(eIDb1, 2), :);   % Elimina filas de ceros si las hay
eIDb1 = eIDb1(1,:);                % Obtiene los elementos superficiales del bloque 1 

GdL_pl = zeros(nex(1,1), 2);

for i = 1:nex(1,1)                 % Este loop obtiene los GdL superficiales del bloque 1
  GdL_pl(i, :) = [ GdLe( eIDb1(1, i), 5 ), GdLe( eIDb1(1, i), 7 ) ];      % Los G.d.L 5 y 7 son horizontales
end 


GdL_pl   = unique(GdL_pl);          % Elimina los GdL repetidos; estos son los GdL para graficar
dri = abs( X(1,1) - dr);
[~, pos] = ismember(dri, Cnx);      % Obtiene los grados de libertad en el intervalo especificado en ig
pos      = unique(pos);             % Si los puntos especificados en dr no coinciden con los nodos, dará un vector de ceros

if sum(pos) == 0;
    disp('Verificar que los puntos en que se está calculado la respuesta coincida con los nodos de la superficie')
end

GdL_pl   = GdL_pl(pos(1:end));     % Obtiene los G.d.L para graficar las F.T.
out_pos = Cnx(pos);
nr = length(GdL_pl);               % Obtiene el número de grados de libertad 


for j = 1 : nr
        dis(j) = X(1,1) - out_pos(j);  % Obtiene las distancias de los puntos obtenidos medidas desde la base del tanque
end

switch op
    case 'C'
         H = H;

    case 'D'
        oj = find(om == omc);
        H = [ 1./H(:, 1:oj), H(:, oj+1:end) ];

    otherwise
        disp('Redefinir apropiadamente la variable "op" ')
end 


for j = 1:nr
    figure;            % Grafica el ID de los nodos
    hold on; grid on; grid minor
    Hj = H(GdL_pl(j), :);
    ab = abs(Hj);
    plot(om', ab, '-o', 'Color', 'b', 'LineWidth', 1.5, 'MarkerSize', 0.5); hold on
    title(['Función de Transferencia a: ', num2str(dis(j)), 'ft de la base'])
    ylabel('H(\omega)'); xlabel('Frecuencia (rad/seg)');
    xlim([0 omm]);

    % Ajustar la letra y visibilidad de las líneas de la cuadrícula
    set(gca,'FontSize',20); set(gcf,'color','w');  
    ax = gca;                       % Obtener el objeto de los ejes actuales
    ax.GridLineStyle = '-';         % Establecer el estilo de línea de la cuadrícula principal (puede ser '-', '--', etc.)
    ax.MinorGridLineStyle = ':';    % Establecer el estilo de línea para la cuadrícula menor (puede ser ':', '--', etc.)
    ax.GridAlpha = 1;               % Cuadrícula principal completamente opaca
    ax.MinorGridAlpha = 1;          % Cuadrícula menor con algo de transparencia


    [xg, xgo, N] = F_convolucion(nom, dt, Hj, omm, dis(j), PGAs, g, j);
    Espectros(xg, xgo, dt, N, nom, dis(j)); % calcula espectros de respuesta en superficie y en la base

    dis_val = dis(j);  % Valor correcto de dis


    switch sr
        case 'SI'
            % Esto es para generar el nombre de los archivos de salida.
            k = 'Registro de acc obtenido a ';
            k2 = ' ft de la base';
            Txt = strcat(k, num2str(dis_val), k2, '.txt');
            name = char(Txt);
            % Elimina espacios adicionales
            name = strtrim(name);
            
            out = xgo';
            QQ = fullfile(f_out, name);  % Utiliza fullfile para crear la ruta
            
            % Verificar si la carpeta de destino existe
            if ~exist(f_out, 'dir')
                warning('La carpeta de destino no existe: %s', f_out);
            else
                try
                    % Abre el archivo para escribir el comentario y los datos
                    fileID = fopen(QQ, 'w');
                    if fileID == -1
                        error('No se pudo abrir el archivo: %s', QQ);
                    end
                    % Escribe el comentario al inicio del archivo
                    fprintf(fileID, '%% dt = %.3f sec\n', dt);
                    % Escribe los datos en formato ASCII
                    fclose(fileID);
                    save(QQ, 'out', '-ascii', '-append');
                catch ME
                    warning('No se pudo guardar el archivo: %s\nError: %s', QQ, ME.message);
                end
            end

        case 'NO'
            if j == nr
                disp('No se guardó el registro de aceleraciones obtenido');
            end
        otherwise
            if j == nr
                disp('No se guardó el registro de aceleraciones obtenido');
            end
     end


   switch sv
     case 'SI'
        
        % Esto es para generar el nombre de los archivos de salida.
        k    = {'F de T. a '};
        k2    = {' ft de la base'};
        Txt = strcat(k, num2str(dis_val), k2, '.txt');
        name = char(Txt);

        % Elimina espacios adicionales
        name = strtrim(name);
 
        out = [ om', ab' ];
        QQ = fullfile(f_out, name);  % Utiliza fullfile para crear la ruta
        
        % Verificar si la carpeta de destino existe
        if ~exist(f_out, 'dir')
                warning('La carpeta de destino no existe: %s', f_out);
            else
                try
                    % Abre el archivo para escribir el comentario y los datos
                    fileID = fopen(QQ, 'w');
                    if fileID == -1
                        error('No se pudo abrir el archivo: %s', QQ);
                    end
                    % Escribe el comentario al inicio del archivo
                    fprintf(fileID, '%% om (rad/seg), |H(om)|\n');
                    % Escribe los datos en formato ASCII
                    fclose(fileID);
                    save(QQ, 'out', '-ascii', '-append');
                catch ME
                    warning('No se pudo guardar el archivo: %s\nError: %s', QQ, ME.message);
                end
            end
       
       case 'NO'
         if j == nr
           disp('No se guardaron  las funciones de transferencia')
         end
       otherwise
         if j == nr
           disp('No se guardaron  las funciones de transferencia')
         end
  end    

end


disp('  '); disp('********************** Resultados del programa *********************** '); disp('  ')
disp([' Dimensiones del modelo: L = ', num2str(sum(X)), ', H = ', num2str(sum(hc))  ]); disp(' ')
disp([' Número total de capas: ', num2str(nc) ]); disp(' ')
disp([' Espesor de cada capa: ', num2str(hc) ]); disp(' ')
disp([' Pesos unitarios de los suelos: ', num2str(gam) ]); disp(' ')
disp([' Velocidad de ondas de corte por capa: ', num2str(c) ]); disp(' ')
disp([' Módulo de corte por capa: ', num2str(Gci) ]); disp(' ')
disp([' Módulo de elasticidad por capa: ', num2str(Eci) ]); disp(' ')
disp([' Razón de amortiguamiento por capa: ', num2str(zi) ]); disp(' ')
disp([' Razones de Poisson por capa: ', num2str(nu) ]); disp(' ')
disp([' Puntos en que se calculó la respuesta medidos desde el extremo izquierdo: ', num2str(out_pos) ]); disp(' ')

toc
