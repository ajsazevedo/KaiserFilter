%   PROJETO 3: FILTRO FIR PELO METODO DA JANELA

% Requisitos do projeto:
%Frequencia de amostragem
fsamp = 2*pi;

% atenuacao minima >= 50db
mags = [1 0];
devs = [0.1 0.00315];

% regiao de transicao delta < 0.1pi
% frequencia de corte = pi/2 (a -6db)
fcutsa = (pi/2)-(pi*0.04);
fcutsb = (pi/2)+(pi*0.04);
fcuts = [fcutsa fcutsb];

% Usando a janela de Kaiser

[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);

% Forma direta
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'scale');

% Forma cascata (seções de segunda ordem)
[z,p,k] = tf2zpk(hh,1);
sos = zp2sos(z,p,k);

% filtro na forma direta
fvtool(hh);
% filtro em cascata em seções de segunda ordem
fvtool(sos);

% arredondando para 5 casas decimais
hh5 = round(hh,5,'significant');
fvtool(hh5);
sos5 = round(sos,5,'significant');
fvtool(sos5);

% arredondando para 4 casas decimais
hh4 = round(hh,4,'significant');
fvtool(hh4);
sos4 = round(sos,4,'significant');
fvtool(sos4);

% arredondando para 3 casas decimais
hh3 = round(hh,3,'significant');
fvtool(hh3);
sos3 = round(sos,3,'significant');
fvtool(sos3);

% arredondando para 2 casas decimais
hh2 = round(hh,2,'significant');
fvtool(hh2);
sos2 = round(sos,2,'significant');
fvtool(sos2);

% arredondando para 1 casa decimal
hh1 = round(hh,1,'significant');
fvtool(hh1);
sos1 = round(sos,1,'significant');
fvtool(sos1);

% filtro com z^-1 => -z^-1; passa alta
sys1 = zeros(1);
for i = 1:length(hh)
    sys1(i) = hh(i)*(-1)^(i-1);
end
fvtool(sys1);

% filtro com z^-1 => z^2; rejeita faixa
sys2 = zeros(1);
for i = 1:length(hh)
    sys2((i*2)) = hh(i);
end
fvtool(sys2);

% filtro com z^-1 => -z^-2; passa faixa
sys3 = zeros(1);
for i = 1:length(hh)
    sys3((i*2)-1) = hh(i)*(-1)^(i-1);
end
fvtool(sys3);