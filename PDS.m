%   PROJETO 3: FILTRO FIR PELO METODO DA JANELA

% Requisitos do projeto:
%Frequencia de amostragem
% W = wT = 2pif/fs = 2pi/N
% Para evitar aliasing a FA deve ser maior que a frequencia maxima do
% sinal, dessa forma, a maior FA eh 2pi.
fsamp = 2*pi;

% regiao de transicao Deltaw < 0.1pi; usamos Deltaw = 0.08
% frequencia de corte: wc = pi/2 (a -6db, padrao)
% aqui observamos o inicio e o final da faixa de transicao
% fcutsa: wc - ws
% fcutsb: wc + wp
% wp - ws = Deltaw
fcutsa = (pi/2) - (pi*0.04);
fcutsb = (pi/2) + (pi*0.04);
fcuts = [fcutsa fcutsb];

% atenuacao minima >= 50db
% A = -20 log dta => Escolhemos A = 50db
% Oscilacao das bandas => Pela formula acima, inferimos que dta = 0.0031622776602 
% M = A-8 / 2,285*(Deltaw) => M = 73.134656563015362169508950346232?
% Como A = 50, B = 0.5842(A - 21)^0.4 + 0.07886(A - 21) => B = 4.5335141209812482327179764338239?
% define os valores dos ripples
devs = [0.1 0.0031622776602];
% define o filtro como passa baixa
mags = [1 0];

% Usando a janela de Kaiser
% Nessa funcao obtemos os parametros para o Filtro de Kaiser:
% M = 74
% Beta = 4.5335
% Wn = 0.5 - Frequencia de corte (0.5*pi)
% ftype = low (Low Pass Filter)
[M, Wn, Beta, ftype] = kaiserord(fcuts, mags, devs, fsamp);

% Janela de Kaiser na forma direta
% Usando os parametros calculados a mao nao muda muita coisa nao
hh = fir1(74, Wn, ftype, kaiser(74+1, 4.5335141209812482327179764338239), 'scale');

% Forma cascata (seções de segunda ordem)
% Discrete-time transfer function to zero-pole conversion. => Converte a
% funcao filtro para uma matriz de zeros/polos/ganhos
[z, p, k] = tf2zpk(hh, 1);
% Zero-pole-gain to second-order sections model conversion
% Converte a matriz de zeros/polos/ganhos para um modelo de secoes de
% segunda ordem
sos = zp2sos(z, p, k);

% Filtro ideal
% na forma direta
fvtool(hh);
% em cascata em seções de segunda ordem
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

% filtro com z^-1 => -z^-2; passa faixa
sys2 = zeros(1);
for i = 1:length(hh)
    sys2((i*2)-1) = hh(i)*(-1)^(i-1);
end
fvtool(sys2);