function [w, dw, xgi, omc, ncut1] = F_frec_dis(nom, dt, omm, PGAs, g)


npt  = 1;                                     % porciento de total de ptos p/graficar en tiempo
npw  = 1;                                     % porciento de total de ptos p/graficar en frec.
wcut = omm;                                   % frecuencia de corte para filtrar el registro y las F. de transferencia


% ----------------- Lectura y gráfico del registro sísmico --------------------

terr = load ([nom,'.at2']);                 % lee archivo con el registro sísmico
[nr,nc]  = size(terr); 	                    % nro. de filas y columnas del registro
nt       = nr*nc;						    % número original de puntos del registro
xgi(1:nt) = terr';                          % copia el acelerograma en un vector

amx = max(abs(xgi));                        % PGA del registro original
if PGAs ~= 0                                % averigua si hay que escalar el registro
        xgi = PGAs/amx*xgi; 	            % escala el registro original (adimensional)
    else
        PGAs = amx;                         % define el PGA si no se escala
end
%xg  = [xg*g, zeros(1,100)]; 	            % agrega ceros al final del registro
xgi  =  xgi*g;                              % no agrega ceros al final del registro
N   = length(xgi);	      		            % nuevo número de ptos. del registro

if N ~= 2*fix(N/2)                          % verifica que N sea un número par
   N     = N + 1;                           % agrega un punto si N es impar
   xgi(N) = 0;                              % agrega un cero al registro
end
tf  = (N-1) * dt;	     			        % duración final del registro: seg
t   = 0 : dt : tf;				            % vector con los tiempos de muestreo



% -------------- Transformada de Fourier discreta del registro --------------------

T   = N * dt;                               % constante de tiempo: seg 
wny = pi/dt;				                % frecuencia de Nyquist: rad/s
dw  = 2*pi / T;                             % intervalo de frecuencias: rad/s
w   = 0: dw: N/2*dw;                        % vector con frecs. discretas
Xg  = fft(xgi);                             % Transf. de Fourier del registro


ncut1 = round(wcut/dw)+5;
%ncut2 = round((2*wny-wcut)/dw);
ncut2 = round(wny/dw);
Xg(ncut1:ncut2) = 0;
omc = w(ncut1);


Np  = round(npw*N/2);                        % número de frecuencias para graficar 
xgi = ifft(Xg);                              % aceleración en la superficie filtrada
xgi = real(xgi/g);                           % recupera la parte real y divide por g
xgi = detrend(xgi,2);                        % corrección de línea de base de acelerog. superf.

amx = max(abs(xgi));                         % PGA del registro en superficie filtrado

 end
