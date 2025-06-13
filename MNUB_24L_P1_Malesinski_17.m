clear
close all
clc

%% Zadanie 1

% ładowanie pliku z danymi
load('MNUB_24L_P1_dane17.mat');

% inicjowanie pustych macierzy
A = zeros(3, 3);
b = zeros(3, 165);

% obliczanie współczynników macierzy A
for r = 1:3
    A(r,1) = 2*(-R(r,1)+R(4,1));
    A(r,2) = 2*(-R(r,2)+R(4,2));
    A(r,3) = 2*(-R(r,3)+R(4,3));
end

% obliczanie współczynników macierzy b
for n = 1:165
    for r = 1:3
        b(r,n) = D(n,r)^2-D(n,4)^2-R(r,1)^2+R(4,1)^2-R(r,2)^2+R(4,2)^2-R(r,3)^2+R(4,3)^2;
    end
end

%% Zadanie 2

% Rozkład LU macierzy A
[L, U] = lu(A);

% inicjowanie pustej macierzy
x = zeros(3, 165);

% obliczanie rozwiązań równania metodą LU
for n = 1:165
    x(:,n) = U\(L\b(:,n));
end

% Odtwarzanie trajektorii ruchu osoby w 3D
figure;
plot3(x(1,:), x(2,:), x(3,:), 'k', 'LineWidth', 3); % Trajektoria ruchu osoby
hold on;
% Zaznaczanie położeń R1, R2, R3 i R4 jako trójkąty
scatter3(R([1, 3, 2, 4],1), R([1, 3, 2, 4],2), R([1, 3, 2, 4],3), 150, 'b', 'filled', 'Marker', 'v');
text(R(:,1), R(:,2), R(:,3), "R" + num2str((1:length(R))'), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right'); % Numeracja punktów
xlabel('x[m]');
ylabel('y[m]');
zlabel('z[m]');
title('Rozmieszczenie czujników radarowych i trówymiarowa trajektoria ruchu osoby');
grid on;
hold off;

% Odtwarzanie trajektorii ruchu osoby na płaszczyźnie 2D
figure;
plot(t, x(1,:));% Trajektoria ruchu osoby na płaszczyźnie x
hold on
plot(t, x(2,:));% Trajektoria ruchu osoby na płaszczyźnie y
hold on
plot(t, x(3,:));% Trajektoria ruchu osoby na płaszczyźnie z
h_legend = legend('x', 'y', 'z');
set(h_legend, 'FontSize', 10);
set(h_legend, 'Position', [0.7,0.5,0.2,0.2]);
title('Wartości współrzędnych środka ciężkości osoby w trakcie ruchu');
xlabel('Czas [s]');
ylabel('x, y, z [m]');
grid on;
hold off;

%% Zadanie 3

% Zdefiniowanie odchylenia standardowego
sigma_d = sort([10^(-3) 2*10^(-3) 3*10^(-3) 4*10^(-3) 5*10^(-3) 6*10^(-3) 7*10^(-3) 8*10^(-3) 9*10^(-3) 10^(-2) 2*10^(-2) 3*10^(-2) 4*10^(-2) 5*10^(-2) 6*10^(-2) 7*10^(-2) 8*10^(-2) 9*10^(-2) 10^(-1)]);

% Zaburzenie multiplikatywne
D_zaburzone = ones(165,4,19);
for i = 1:19
    D_zaburzone(:,:,i)=(1+randn(size(D))*sigma_d(i)).*D;
end

% Wyznaczenie maksymalnego zagagregowanego błędu
x_zaburzone = zeros(3, 165, 19);
b_zaburzone = zeros(3, 165, 19);

for i = 1:19
    for n = 1:165
        for r = 1:3
            b_zaburzone(r,n,i) = D_zaburzone(n,r,i)^2-D_zaburzone(n,4,i)^2-R(r,1)^2+R(4,1)^2-R(r,2)^2+R(4,2)^2-R(r,3)^2+R(4,3)^2;
        end
    end
   for n = 1:165
            x_zaburzone(:,n, i) = U\(L\b_zaburzone(:,n,i));
   end    
end

% Obliczanie zagregowanego błędu
zagregowany_blad=zeros(1, 19);

for i = 1:19
        zagregowany_blad(1,i)= norm(x_zaburzone(:,:,i)-x, 'fro')/norm(x, 'fro');
end

% Obliczanie zagregowanego błędu dla wektora b
zagregowany_blad_b=zeros(1, 19);

for i = 1:19
        zagregowany_blad_b(1,i)= norm(b_zaburzone(:,:,i)-b, 'fro')/norm(b, 'fro');
end

oszacowany_zagregowany_blad=zeros(1,19);

% Szacowanie zagrogowanego błędu na podstawie liczby kondycyjnej A i
% zagregowanego błędu wektora b
for i = 1:19
    oszacowany_zagregowany_blad(1,i)=cond(A)*zagregowany_blad_b(1,i);
end

%Rysowanie wykresu
figure;
loglog(sigma_d, zagregowany_blad);
hold on;
loglog(sigma_d, oszacowany_zagregowany_blad);
hold off;
h_legend = legend('wyznaczona wartość błędu', 'oszacowana wartość błędu', 'Location', 'best');
set(h_legend, 'FontSize', 10);
xlabel('Sigma');
ylabel('Maksymalny zagregowany błąd względny');
title('Zależność wyznaczonego oraz oszacowanego maksymalnego błędu względnego w zależności od wartości sigmy');
grid on;

%% Zadanie 4
% Zdefiniowanie odchylenia standardowego
sigma_z = linspace(-0.1, 0.1, 150);

% Zaburzenie wartości R
R_zaburzone_1 = zeros(4,3,150);
R_zaburzone_2 = zeros(4,3,150);
R_zaburzone_3 = zeros(4,3,150);
R_zaburzone_4 = zeros(4,3,150);

for i=1:150
    R_zaburzone_1(:,:,i)=R;
    R_zaburzone_2(:,:,i)=R;
    R_zaburzone_3(:,:,i)=R;
    R_zaburzone_4(:,:,i)=R;
    R_zaburzone_1(1,3,i)=R(1,3)*(1+sigma_z(i));
    R_zaburzone_2(2,3,i)=R(2,3)*(1+sigma_z(i));
    R_zaburzone_3(3,3,i)=R(3,3)*(1+sigma_z(i));
    R_zaburzone_4(4,3,i)=R(4,3)*(1+sigma_z(i));
end

b_R_1=zeros(3, 165, 150);
b_R_2=zeros(3, 165, 150);
b_R_3=zeros(3, 165, 150);
b_R_4=zeros(3, 165, 150);
A_R_1=zeros(3,3,150);
A_R_2=zeros(3,3,150);
A_R_3=zeros(3,3,150);
A_R_4=zeros(3,3,150);
x_R_1=zeros(3,165,150);
x_R_2=zeros(3,165,150);
x_R_3=zeros(3,165,150);
x_R_4=zeros(3,165,150);

% Obliczanie współczynników macierzy A_R
for i = 1:150
    for r = 1:3
        A_R_1(r,1, i) = 2*(-R_zaburzone_1(r,1,i)+R_zaburzone_1(4,1,i));
        A_R_1(r,2, i) = 2*(-R_zaburzone_1(r,2,i)+R_zaburzone_1(4,2,i));
        A_R_1(r,3, i) = 2*(-R_zaburzone_1(r,3,i)+R_zaburzone_1(4,3,i));

        A_R_2(r,1, i) = 2*(-R_zaburzone_2(r,1,i)+R_zaburzone_2(4,1,i));
        A_R_2(r,2, i) = 2*(-R_zaburzone_2(r,2,i)+R_zaburzone_2(4,2,i));
        A_R_2(r,3, i) = 2*(-R_zaburzone_2(r,3,i)+R_zaburzone_2(4,3,i));

        A_R_3(r,1, i) = 2*(-R_zaburzone_3(r,1,i)+R_zaburzone_3(4,1,i));
        A_R_3(r,2, i) = 2*(-R_zaburzone_3(r,2,i)+R_zaburzone_3(4,2,i));
        A_R_3(r,3, i) = 2*(-R_zaburzone_3(r,3,i)+R_zaburzone_3(4,3,i));

        A_R_4(r,1, i) = 2*(-R_zaburzone_4(r,1,i)+R_zaburzone_4(4,1,i));
        A_R_4(r,2, i) = 2*(-R_zaburzone_4(r,2,i)+R_zaburzone_4(4,2,i));
        A_R_4(r,3, i) = 2*(-R_zaburzone_4(r,3,i)+R_zaburzone_4(4,3,i));
    end
end

% Obliczanie współczynników macierzy b_R
for i = 1:150
    for n = 1:165
        for r = 1:3
            b_R_1(r,n,i) = D(n,r)^2-D(n,4)^2-R_zaburzone_1(r,1,i)^2+R_zaburzone_1(4,1,i)^2-R_zaburzone_1(r,2,i)^2+R_zaburzone_1(4,2,i)^2-R_zaburzone_1(r,3,i)^2+R_zaburzone_1(4,3,i)^2;
            b_R_2(r,n,i) = D(n,r)^2-D(n,4)^2-R_zaburzone_2(r,1,i)^2+R_zaburzone_2(4,1,i)^2-R_zaburzone_2(r,2,i)^2+R_zaburzone_2(4,2,i)^2-R_zaburzone_2(r,3,i)^2+R_zaburzone_2(4,3,i)^2;
            b_R_3(r,n,i) = D(n,r)^2-D(n,4)^2-R_zaburzone_3(r,1,i)^2+R_zaburzone_3(4,1,i)^2-R_zaburzone_3(r,2,i)^2+R_zaburzone_3(4,2,i)^2-R_zaburzone_3(r,3,i)^2+R_zaburzone_3(4,3,i)^2;
            b_R_4(r,n,i) = D(n,r)^2-D(n,4)^2-R_zaburzone_4(r,1,i)^2+R_zaburzone_4(4,1,i)^2-R_zaburzone_4(r,2,i)^2+R_zaburzone_4(4,2,i)^2-R_zaburzone_4(r,3,i)^2+R_zaburzone_4(4,3,i)^2;
        end
    end
end

% Obliczanie rozwiązań równania metodą LU
for i=1:150
    [L, U] = lu(A_R_1(:,:,i));
    for n = 1:165
        x_R_1(:,n,i) = U\(L\b_R_1(:,n,i));
    end

    [L, U] = lu(A_R_2(:,:,i));
    for n = 1:165
        x_R_2(:,n,i) = U\(L\b_R_2(:,n,i));
    end

    [L, U] = lu(A_R_3(:,:,i));
    for n = 1:165
        x_R_3(:,n,i) = U\(L\b_R_3(:,n,i));
    end

    [L, U] = lu(A_R_4(:,:,i));
    for n = 1:165
        x_R_4(:,n,i) = U\(L\b_R_4(:,n,i));
    end
end

% Obliczanie zagregowanego błędu
zagregowany_blad_R=zeros(1, 150, 4);

for i = 1:150
    zagregowany_blad_R(1,i,1)= norm(x_R_1(:,:,i)-x, 'fro')/norm(x, 'fro');
    zagregowany_blad_R(1,i,2)= norm(x_R_2(:,:,i)-x, 'fro')/norm(x, 'fro');
    zagregowany_blad_R(1,i,3)= norm(x_R_3(:,:,i)-x, 'fro')/norm(x, 'fro');
    zagregowany_blad_R(1,i,4)= norm(x_R_4(:,:,i)-x, 'fro')/norm(x, 'fro');
end

% Wyznaczanie wartości wyznacznika macierzy A oraz jej liczby kondycyjnej
wyznacznik_macierzy_A = zeros(1,150,4);
liczba_kondycyjna_A = zeros(1,150,4);

for i=1:150
    wyznacznik_macierzy_A(1,i,1)=det(A_R_1(:,:,i));
    liczba_kondycyjna_A(1,i,1)=cond(A_R_1(:,:,i));

    wyznacznik_macierzy_A(1,i,2)=det(A_R_2(:,:,i));
    liczba_kondycyjna_A(1,i,2)=cond(A_R_2(:,:,i));

    wyznacznik_macierzy_A(1,i,3)=det(A_R_3(:,:,i));
    liczba_kondycyjna_A(1,i,3)=cond(A_R_3(:,:,i));

    wyznacznik_macierzy_A(1,i,4)=det(A_R_4(:,:,i));
    liczba_kondycyjna_A(1,i,4)=cond(A_R_4(:,:,i));
end

%Rysowanie wykresów
figure;
grid on;
plot(sigma_z, zagregowany_blad_R(:,:,1));
hold on;
plot(sigma_z, zagregowany_blad_R(:,:,2));
plot(sigma_z, zagregowany_blad_R(:,:,3));
plot(sigma_z, zagregowany_blad_R(:,:,4));
hold off;
xlabel('Sigma');
ylabel('Maksymalny zagregowany błąd względny');
h_legend = legend('zaburzony z1', 'zaburzony z2',  'zaburzony z3', 'zaburzony z4','Location', 'best');
set(h_legend, 'FontSize', 10);
title('Zależność sigmy od zagregowanego błędu')

figure;
grid on;
plot(sigma_z, wyznacznik_macierzy_A(:,:,1));
hold on;
plot(sigma_z, wyznacznik_macierzy_A(:,:,2));
plot(sigma_z, wyznacznik_macierzy_A(:,:,3));
plot(sigma_z, wyznacznik_macierzy_A(:,:,4));
hold off;
xlabel('Sigma');
ylabel('Wyznacznik macierzy A');
h_legend = legend('zaburzony z1', 'zaburzony z2',  'zaburzony z3', 'zaburzony z4','Location', 'best');
set(h_legend, 'FontSize', 10);
title('Zależność sigmy od wyznacznika macierzy')

figure;
grid on;
plot(sigma_z, liczba_kondycyjna_A(:,:,1));
hold on;
plot(sigma_z, liczba_kondycyjna_A(:,:,2));
plot(sigma_z, liczba_kondycyjna_A(:,:,3));
plot(sigma_z, liczba_kondycyjna_A(:,:,4));
hold off;
xlabel('Sigma');
ylabel('Liczba kondycyjna A');
h_legend = legend('zaburzony z1', 'zaburzony z2',  'zaburzony z3', 'zaburzony z4','Location', 'best');
set(h_legend, 'FontSize', 10);
title('Zależność sigmy od liczby kondycyjnej macierzy')