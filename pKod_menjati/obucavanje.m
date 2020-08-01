% Projekat za Petnicu - Detekcija Photoshop-a na slikama
% Obucavanje klasifikatora

% Ova skripta se koristi za obucavanje klasifikatora. Podrazumeva se
% sledece:
% - Slike se nalaze u folderu baza
% - Originalne slike su oznacene kao oi.jpg gde je i broj slike
% - Modifikovane slike su oznacene kao mi.jpg gde je i broj slike
% - Slika mi.jpg je dobijena koriscenjem liqufy ili blur alata u
% photoshop-u od slike oi.jpg
% - Svaka slika je jedna tekstura (za pocetak jedna, pa ako algoritam radi
% moze da se poboljsa dodavanjem segmentacije slike po teksturama

%

close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N - broj slika u bazi
% No - broj slika koje cine obucavajuci skup

N = 70; % U konacnom projektu naravno vise slika
No = 60; % No je oko 0.8N

folder = 'baza4\';

% Za pocetak, citamo svaku sliku, i za celu sliku koja predstavlja teksturu
% odredimo vrednosti obelezja. 
TextureParameters = zeros(12, N);
L1 = 10; L2 = 25;
for i = 1:N
    % Ucitaj modifikovanu sliku 
    name = ['s' num2str(i,'%d')]; 
    img = rgb2gray(imread(strcat(folder,name), 'jpg'));
    
    % Racunamo parametre teksture
    TextureParameters(:, i) = parametriTeksture(img, L1, L2); 
end

% Sve slike su 300x500x3 dimenzija, prema tome ima 300x500x3xN piksela

% Na koliko blokova se deli slika L = 10, znaci 10x10 blokova

% Svaka slika ima L*L blokova po svakoj dimenziji
X = zeros(12, N * L1 * L2); % za svaku sliku po 100 kolona
Y = zeros(1, N * L1 * L2); % 0 ili 1

% Ucitaj izlaz za svaki blok
for i =1:N
    %Ucitaj binarizovanu razliku sliku
    name_bw = ['os' num2str(i,'%d')]; 
    bw = (imread(strcat(folder,name_bw), 'jpg'));
    
    razlika = zaokruziRazliku(bw, L1,L2, 0.95);
    
    Y(:, (i-1)*L1*L2+1 : i*L1*L2) = razlika(:);
end

for i=1:N
    % Ucitaj modifikovanu sliku 
    name = ['s' num2str(i,'%d')]; 
    img = rgb2gray(imread(strcat(folder,name), 'jpg'));
    
    % Nadji obelezja i stavi u X.
    o = obelezja(img, TextureParameters(:, i), L1, L2);
    X(:, (i-1)*L1*L2+1 : i*L1*L2) = o;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prvih No slika su obucavajuce, znaci L*L*No vrednosti vektora X su
% obucavajuci skup, ostatak sluzi za testiranje

X_O = X(:, 1:No*L1*L2);
X_T = X(:, No*L1*L2+1 : end);

Y_O = Y(:, 1:No*L1*L2);
Y_T = Y(:, No*L1*L2+1 : end);

Klasa1 = (X(:, Y == 1))';
Klasa2 = (X(:, Y == 0))';

% Srednja vrednost svake klase, kovarijaciona matrica i verovatnoca
% pojavljivanja klase.
% 1. Klasa - modifikovani blokovi
M1 = mean(Klasa1)'; 
S1 = cov(Klasa1)';
P1 = sum(Y)/(N*L1*L2);

% 2. Klasa - nemodifikovani blokovi
M2 = mean(Klasa2)'; 
S2 = cov(Klasa2)';
P2 = 1- P1;
% Pripremi matrice potrebne za redukciju dimenzija.


Sw = P1*S1 + P2*S2; % Within class scatter matrix.
M0 = P1*M1 + P2*M2; % Combined expected value vector.
Sb = P1*(M1-M0)*(M1-M0)'+P2*(M2-M0)*(M2-M0)'; % Between class scatter matrix.
Sm = Sw + Sb; % Mixed matrix.

% Take eigenvalues. 
S1 = Sb; S2 = Sw;
S = S2^(-1)*S1;
[F Lambda] = eig(S); %vraca dijagonalnu matricu D na kojoj se nalaze sop.v
                     %Eigenvalues and eigenvectors.

% Take 2 eigenvectors for which lambda is maximum.
[l, ind] = sort(diag(Lambda), 'descend');
A = [F(:, ind(1)), F(:, ind(2))];

% Plot the information perserved by taking first n dimensions.
l_sum = sum(l);
I = zeros(size(l));
m = zeros(size(l));
for i=1:length(I)
    I(i) = sum(l(1:i))/l_sum;
    m(i) = i;
end

figure(3),
plot(m, I*100), title('I(m) %'),
xticks(m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Klasa1 = Klasa1';
Klasa2 = Klasa2';

Z1 = zeros(2, size(Klasa1, 2));
Z2 = zeros(2, size(Klasa2, 2));

%Transform input samples
for i = 1: size(Klasa1,2)
    Z1(:, i) = A' * Klasa1(:, i);
end
%Transform input samples
for i = 1: size(Klasa2,2)
    Z2(:, i) = A' * Klasa2(:, i);
end

figure(1)
plot(Z1(1, :), Z1(2, :), 'sb', Z2(1, :), Z2(2, :), 'hg');
