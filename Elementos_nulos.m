function [en, nn] = Elementos_nulos(nez, nex, nen)


b = cell(1, 3);  % Esta celda es para guardar resultados.

for i = 1:3
    if nez(1, i) == 0
        b{i} = zeros(nez(1, i) + 1, nex(1, i) + 1);
    else
        b{i} = ones(nez(1, i) + 1, nex(1, i) + 1);
    end
end

x = cell(1, 3);  % Esta celda es para guardar resultados.

for i = 1:3
    if nen(1, i) == 0
        x{i} = [];
    else
        x{i} = zeros(nen(1, i), nex(1, i) + 1);
    end
    x{i} = [x{i}; b{i}];
end

for i = 1:2
    if sum(x{i}(:, end)) >= sum(x{i+1}(:, 1))
        x{i+1} = x{i+1}(:, 2:end);
    else
        x{i} = x{i}(:, 1:end-1);
    end
end

nn = cat(2, x{:}); 

enb = zeros(max(nez), 3);
for i = 1:3
    enb(:, i) = [zeros(1, nen(1, i)) ones(1, nez(1, i))]';  % Define los elementos nulos de cada bloque
end

c = cell(1, 3); 

for i = 1:3
    c{i} = repmat(enb(:, i), 1, nex(1, i));
end

en = cat(2, c{:});  


end





