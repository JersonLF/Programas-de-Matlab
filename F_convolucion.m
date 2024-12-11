function [xgi, xgo, N] = F_convolucion(nom, dt, H, omm, nsupl, PGAs, g, j)

[rH, cH] = size(H);
npt  = 1;                                     % porciento de total de ptos p/graficar en tiempo
npw  = 1;                                     % porciento de total de ptos p/graficar en frec.
wcut = omm;                                   % frecuencia de corte para filtrar el registro


% ----------------- Lectura y gr�fico del registro s�smico --------------------

terr = load ([nom,'.at2']);                 % lee archivo con el registro s�smico
[nr,nc]  = size(terr); 	                    % nro. de filas y columnas del registro
nt       = nr*nc;						    % n�mero original de puntos del registro
xgi(1:nt) = terr';                          % copia el acelerograma en un vector

amx = max(abs(xgi));                        % PGA del registro original
if PGAs ~= 0                                % averigua si hay que escalar el registro
        xgi = PGAs/amx*xgi; 	            % escala el registro original (adimensional)
    else
        PGAs = amx;                         % define el PGA si no se escala
end
%xg  = [xg*g, zeros(1,100)]; 	            % agrega ceros al final del registro
xgi  =  xgi*g;                              % no agrega ceros al final del registro
N   = length(xgi);	      		            % nuevo n�mero de ptos. del registro

if N ~= 2*fix(N/2)                          % verifica que N sea un n�mero par
   N     = N + 1;                           % agrega un punto si N es impar
   xgi(N) = 0;                              % agrega un cero al registro
end
tf  = (N-1) * dt;	     			        % duraci�n final del registro: seg
t   = 0 : dt : tf;				            % vector con los tiempos de muestreo

if j == 1
np = round(npt*nt);                         % nr. de ptos. p/graficar acelerograma
figure; plot( t(1:np),xgi(1:np)/g,'k', t(1:np),zeros(1,np) ); axis tight; 
xlabel('Tiempo (seg)'); ylabel('Aceleraci�n (g)'); 
title(['Acelerograma de entrada en el modelo: ',nom]);
legend(['PGA = ',num2str(amx,3),'g'])
grid on; grid minor;

% Ajustar la letra y visibilidad de las l�neas de la cuadr�cula
set(gca,'FontSize',20); set(gcf,'color','w');  
ax = gca;                       % Obtener el objeto de los ejes actuales
ax.GridLineStyle = '-';         % Establecer el estilo de l�nea de la cuadr�cula principal (puede ser '-', '--', etc.)
ax.MinorGridLineStyle = ':';    % Establecer el estilo de l�nea para la cuadr�cula menor (puede ser ':', '--', etc.)
ax.GridAlpha = 1;               % Cuadr�cula principal completamente opaca
ax.MinorGridAlpha = 1;          % Cuadr�cula menor con algo de transparencia

end



% -------------- Transformada de Fourier discreta del registro --------------------

T   = N * dt;                               % constante de tiempo: seg 
wny = pi/dt;				                % frecuencia de Nyquist: rad/s
dw  = 2*pi / T;                             % intervalo de frecuencias: rad/s
w   = 0: dw: N/2*dw;                        % vector con frecs. discretas
Xg  = fft(xgi);                             % Transf. de Fourier del registro


% ------------- C�lculo de las frecuencias y periodos naturales --------------

Fd1 = diff(abs(H),1) / dw;                        % 1ra derivada de |H|
Fd2 = diff(abs(H),2) / dw^2;                      % 2da derivada de |H|
nw = length(w);                                   % nro.de frecuencias discretas
Z   = Fd1(1:nw-2).*Fd1(2:nw-1);
I   = find(Z < 0 & Fd2 < 0);                      % �ndices de las frecs. nat.
wn  = w(I);                                       % frecuencias naturales: rad/s
fn  = wn/(2*pi);                                  % frecuencias naturales: ciclos/s
Tn = 2*pi./wn;                                    % periodos naturales: s

disp(' ')
disp(['Per�odos naturales a ', num2str(nsupl), ' ft de la base (seg): ', num2str(Tn(1:4)) ]);
disp(' ')

ncut1 = round(wcut/dw)+5;
%ncut2 = round((2*wny-wcut)/dw);
ncut2 = round(wny/dw);
Xg(ncut1:ncut2) = 0;

Np  = round(npw*N/2);                       % n�mero de frecuencias para graficar 

if j == 1
figure; plot( w(1:Np),abs(Xg(1:Np)) ); grid on; axis tight
title(['Espectro de Fourier del registro de entrada: ',nom])
xlabel('Frecuencia \omega (rad/seg)'); ylabel('Amplitud')
xlim([0, omm]);
grid on; grid minor;

% Ajustar la letra y visibilidad de las l�neas de la cuadr�cula
set(gca,'FontSize',20); set(gcf,'color','w');  
ax = gca;                       % Obtener el objeto de los ejes actuales
ax.GridLineStyle = '-';         % Establecer el estilo de l�nea de la cuadr�cula principal (puede ser '-', '--', etc.)
ax.MinorGridLineStyle = ':';    % Establecer el estilo de l�nea para la cuadr�cula menor (puede ser ':', '--', etc.)
ax.GridAlpha = 1;               % Cuadr�cula principal completamente opaca
ax.MinorGridAlpha = 1;          % Cuadr�cula menor con algo de transparencia
end

xgi = ifft(Xg);                              % aceleraci�n en la superficie filtrada
xgi = real(xgi)/g;                           % recupera la parte real y divide por g
xgi = detrend(xgi,2);                        % correcci�n de l�nea de base de acelerog. superf.
amx = max(abs(xgi));                         % PGA del registro en superficie filtrado



%- ------------ T. de Fourier de la aceleraci�n absoluta del output -------------

N1  = N/2 + 1;
for i = 1:rH
    H(i, N1+1:N)  = conj( H(i, N/2:-1:2) ); % Genera los valores para n > N/2
end

 
Xs =  H .* Xg;


xsi = ifft( Xs );                       % aceleraci�n en la roca en funci�n del tiempo
xsi = real(xsi)/g;                      % recupera la parte real y divide por g
xgo = detrend(xsi,'linear');            % correcci�n de l�nea de base de acelerograma en superficie
PGAs= max(abs( xsi ));                  % m�xima aceleraci�n en la superficie en g
 

figure; plot( t,xgo,'r'); grid on; axis tight
title(['Aceleraciones obtenidas a: ', num2str(nsupl), ' ft de la base' ]);
xlabel('Tiempo (seg)'); ylabel('Aceleraci�n (g)');
legend(['PGA = ',num2str(PGAs,3),' g'])
grid on; grid minor;

% Ajustar la letra y visibilidad de las l�neas de la cuadr�cula
set(gca,'FontSize',20); set(gcf,'color','w');  
ax = gca;                       % Obtener el objeto de los ejes actuales
ax.GridLineStyle = '-';         % Establecer el estilo de l�nea de la cuadr�cula principal (puede ser '-', '--', etc.)
ax.MinorGridLineStyle = ':';    % Establecer el estilo de l�nea para la cuadr�cula menor (puede ser ':', '--', etc.)
ax.GridAlpha = 1;               % Cuadr�cula principal completamente opaca
ax.MinorGridAlpha = 1;          % Cuadr�cula menor con algo de transparencia



% Grafica las Funciones de Transferencia    
% figure;   
% plot( w(1:Np), imag(H(1:Np)), 'Color', 'r', 'LineWidth',1.5, 'MarkerSize', 0.5); hold on
% plot( w(1:Np), real(H(1:Np)), 'Color', 'b', 'LineWidth',1.5, 'MarkerSize', 0.5); hold off
% axis tight; grid on; grid minor;
% title('Funci�n de Transferencia F_r_s(\omega) entre roca y superficie')
% xlabel('Frecuencia \omega [rad/seg]'); ylabel('Amplitud')
% legend(': Parte imaginaria', ': Parte real')
% set(gca,'FontSize',13);
% set(gcf,'color','w');
% xlim([ 0 omm]);
end   


