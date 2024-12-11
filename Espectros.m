function [PSAs] = Espectros(xg,as,dt,nt,nom, pos)
% Programa "function" para calcular y graficar los espectros de respuesta de % 
% seudoaceleraci�n en la roca y en la superficie de un dep�sito de suelo.    %
% Se usa la soluci�n recursiva de la integral de Duhamel. Los dos registros  %
% de aceleraciones (xg, as) se pasan como argumentos al programa Espectros.m %
%                                                                            %
% ------------- Revisado en: 8 - enero - 2018- Luis E. Su�rez --------------%

zi   = 0.05;                               % raz�n de amortiguamiento
Tmax = 3.0;                                % periodo m�ximo del espectro
dT   = 0.0005;                               % incremento en los periodos
Ti   = 0.0005;                              % periodo inicial del espectro
T    = Ti: dT: Tmax;					   % vector con periodos
nper = length(T);						   % nro. de periodos del espectro

PSAr = zeros(1,nper);				       % Espectro de seudoacel. en la base
PSAs = zeros(1,nper);				       % Espectro de seudoacel. en superf.
for j = 1 : nper
   om     = 2*pi/T(j);
   ur     = duhamel(om,zi,1,dt,nt,0,0,-xg);
   us     = duhamel(om,zi,1,dt,nt,0,0,-as);
   PSAr(j) = om^2 * max( abs(ur) );
   PSAs(j) = om^2 * max( abs(us) );
end

Sm = 1.05*max( [PSAr,PSAs] );

figure; plot(T,PSAr,'-', T,PSAs, 'LineWidth', 2); axis ([Ti Tmax 0 Sm]); 
xlabel('Per�odo (seg)'); ylabel('Sa (g)');
title(['Espectros de aceleraciones del terremoto de ',nom, ', a ', num2str(pos) ' ft de la base'])
legend(': entrada',': salida'); grid on; grid minor;

% Ajustar la letra y visibilidad de las l�neas de la cuadr�cula
set(gca,'FontSize',20); set(gcf,'color','w');  
ax = gca;                       % Obtener el objeto de los ejes actuales
ax.GridLineStyle = '-';         % Establecer el estilo de l�nea de la cuadr�cula principal (puede ser '-', '--', etc.)
ax.MinorGridLineStyle = ':';    % Establecer el estilo de l�nea para la cuadr�cula menor (puede ser ':', '--', etc.)
ax.GridAlpha = 1;               % Cuadr�cula principal completamente opaca
ax.MinorGridAlpha = 1;          % Cuadr�cula menor con algo de transparencia




