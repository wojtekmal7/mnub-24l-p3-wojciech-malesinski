clear
close all
clc

%% Zad1

% Parametry
p1 = 13; % Współczynnik wzrostu dla populacji x
p2 = 0.13; % Współczynnik spadku populacji x spowodowany interakcją z y
p3 = 0.05; % Współczynnik wzrostu populacji y spowodowany interakcją z x
p4 = 11; % Współczynnik spadku dla populacji y

% Warunki początkowe
x0 = 470; % Początkowa wartość populacji x
y0 = 20; % Początkowa wartość populacji y

% Układ równań Lotki-Volterry zdefiniowany jako funkcja anonimowa
lotka_volterra = @(t, z) [p1 * z(1) - p2 * z(1) * z(2); % Równanie dla x
                          p3 * z(1) * z(2) - p4 * z(2)]; % Równanie dla y

% Zakres czasu
tspan = [0 1]; % Od 0 do 1

% Rozwiązanie za pomocą ode45
[t, z] = ode45(lotka_volterra, tspan, [x0 y0]);

% Wykres rozwiązania
figure;
plot(t, z(:, 1), '-r', t, z(:, 2), '-b');
xlabel('Czas t');
ylabel('Populacje');
legend('x(t)', 'y(t)');
title('Rozwiązanie układu równań Lotki-Volterry za pomocą ode45');
grid on;

%% Zad2 (a) Otwarta metoda Eulera

h = 0.005; % Krok całkowania
N = 1/h; % Liczba kroków
t = 0:h:1; % Wektor czasu
x = zeros(1, N+1); % Inicjalizacja wektora populacji x
y = zeros(1, N+1); % Inicjalizacja wektora populacji y
x(1) = x0; % Początkowa wartość x
y(1) = y0; % Początkowa wartość y

% Pętla dla otwartej metody Eulera
for n = 1:N
    x(n+1) = x(n) + h * (p1 * x(n) - p2 * x(n) * y(n)); % Aktualizacja x
    y(n+1) = y(n) + h * (p3 * x(n) * y(n) - p4 * y(n)); % Aktualizacja y
end

% Zapis wyników dla otwartej metody Eulera
y_a = y;
t_a = t;

% Rysowanie wykresu dla otwartej metody Eulera
figure;
plot(t, x, '-r', t, y, '-b');
xlabel('Czas t');
ylabel('Populacje');
legend('x(t)', 'y(t)');
title('Rozwiązanie układu równań Lotki-Volterry metodą Eulera otwartego');
grid on;

%% Zad2 (b) Zamknięta metoda Eulera

% Przestrzeń czasowa
N = length(t);

% Inicjalizacja wektorów
x = zeros(N, 1);
y = zeros(N, 1);
x(1) = x0;
y(1) = y0;

% Rozwiązanie za pomocą zamkniętej metody Eulera
for n = 1:N-1
    % Aktualne wartości
    x_current = x(n);
    y_current = y(n);

    % Funkcja dla fsolve
    fun = @(z) [
        z(1) - x_current - h * (p1 * z(1) - p2 * z(1) * z(2));
        z(2) - y_current - h * (p3 * z(1) * z(2) - p4 * z(2))
    ];

    % Początkowe przybliżenie
    z0 = [x_current; y_current];

    % Rozwiązanie układu równań
    z_next = fsolve(fun, z0, optimset('Display', 'off'));

    % Zapisanie wyników
    x(n + 1) = z_next(1);
    y(n + 1) = z_next(2);
end

y_b = y;
t_b = t;

% Rysowanie wyników
figure;
plot(t, x, 'r', t, y, 'b');
legend('x(t)', 'y(t)');
xlabel('Czas t');
ylabel('Populacja');
title('Zamknięta metoda Eulera dla układu Lotki-Volterry');
grid on;

%% Zad2 (c) Metoda Adamsa-Bashfortha rzędu 2

N = 1/h; % Liczba kroków
x = zeros(1, N+1); % Inicjalizacja wektora populacji x
y = zeros(1, N+1); % Inicjalizacja wektora populacji y
x(1) = x0; % Początkowa wartość x
y(1) = y0; % Początkowa wartość y

% Metoda Eulera na pierwszym kroku
x(2) = x(1) + h * (p1 * x(1) - p2 * x(1) * y(1));
y(2) = y(1) + h * (p3 * x(1) * y(1) - p4 * y(1));

% Pętla dla metody Adamsa-Bashfortha rzędu 2
for n = 2:N
    fx1 = p1 * x(n) - p2 * x(n) * y(n); % Obliczanie fx1
    fy1 = p3 * x(n) * y(n) - p4 * y(n); % Obliczanie fy1
    fx2 = p1 * x(n-1) - p2 * x(n-1) * y(n-1); % Obliczanie fx2
    fy2 = p3 * x(n-1) * y(n-1) - p4 * y(n-1); % Obliczanie fy2
    x(n+1) = x(n) + h/2 * (3 * fx1 - fx2); % Aktualizacja x
    y(n+1) = y(n) + h/2 * (3 * fy1 - fy2); % Aktualizacja y
end

% Zapis wyników dla metody Adamsa-Bashfortha rzędu 2
y_c = y;
t_c = t;

% Rysowanie wykresu dla metody Adamsa-Bashfortha rzędu 2
figure;
plot(t, x, '-r', t, y, '-b');
xlabel('Czas t');
ylabel('Populacje');
legend('x(t)', 'y(t)');
title('Rozwiązanie układu równań Lotki-Volterry metodą Adamsa-Bashfortha rzędu 2');
grid on;

%% Zad2 (d) Dwukrokowa metoda trapezów

x = zeros(1, N+1); % Inicjalizacja wektora populacji x
y = zeros(1, N+1); % Inicjalizacja wektora populacji y
x(1) = x0; % Początkowa wartość x
y(1) = y0; % Początkowa wartość y

% Metoda Eulera na pierwszym kroku
fx1 = p1 * x(1) - p2 * x(1) * y(1);
fy1 = p3 * x(1) * y(1) - p4 * y(1);
x(2) = x(1) + h * fx1;
y(2) = y(1) + h * fy1;

for n = 2:N
    % Predykcja
    fx1 = p1 * x(n) - p2 * x(n) * y(n); % Obliczanie fx1
    fy1 = p3 * x(n) * y(n) - p4 * y(n); % Obliczanie fy1
    x_predict = x(n) + h * fx1; % Predykcja x
    y_predict = y(n) + h * fy1; % Predykcja y
    
    % Korekcja
    fx2 = p1 * x_predict - p2 * x_predict * y_predict; % Obliczanie fx2
    fy2 = p3 * x_predict * y_predict - p4 * y_predict; % Obliczanie fy2
    x(n+1) = x(n) + h/2 * (fx1 + fx2); % Aktualizacja x
    y(n+1) = y(n) + h/2 * (fy1 + fy2); % Aktualizacja y
end

% Zapis wyników dla dwukrokowej metody trapezów
y_d = y;
t_d = t;

% Rysowanie wykresu dla dwukrokowej metody trapezów
figure;
plot(t, x, '-r', t, y, '-b');
xlabel('Czas t');
ylabel('Populacje');
legend('x(t)', 'y(t)');
% title('Rozwiązanie układu równań Lotki-Volterry dwukrokową metodą trapezów');
grid on;

%% Zad3

% Parametry tolerancji dla funkcji ode45
RelTol = 1e-8; % Tolerancja względna
AbsTol = 1e-12; % Tolerancja bezwzględna

% Inicjalizacja macierzy błędów dla każdej metody
errors = zeros(4, 1); % Macierz przechowująca błędy dla każdej z 4 metod

% Dla każdej metody oblicz estymaty y_i
for method = 1:4
    switch method
        case 1 % Otwarta metoda Eulera
            y = y_a; % Populacje y dla metody otwartej Eulera
            t = t_a; % Czas t dla metody otwartej Eulera
        case 2 % Zamknięta metoda Eulera
            y = y_b; % Populacje y dla metody zamkniętej Eulera
            t = t_b; % Czas t dla metody zamkniętej Eulera
        case 3 % Metoda Adamsa-Bashfortha rzędu 2
            y = y_c; % Populacje y dla metody Adamsa-Bashfortha
            t = t_c; % Czas t dla metody Adamsa-Bashfortha
        case 4 % Dwukrotna metoda trapezów
            y = y_d; % Populacje y dla dwukrotnej metody trapezów
            t = t_d; % Czas t dla dwukrotnej metody trapezów
    end
    
    % Rozwiązanie referencyjne za pomocą funkcji ode45 z dokładnymi parametrami
    [t_ref, z_ref] = ode45(lotka_volterra, tspan, [x0 y0], odeset('RelTol', RelTol, 'AbsTol', AbsTol));
    
    y_ref = z_ref(:,2);

    % Interpolacja estymowanych wartości na siatkę czasową t_ref
    y_interp = interp1(t, y, t_ref, 'linear');

    % Obliczanie zagregowanego błędu względnego
    error = sqrt(sum((y_interp - y_ref).^2)) / sqrt(sum(y_ref.^2));
    
    % Zapisanie błędu dla danej metody
    errors(method) = error;
end

% Wyświetlenie błędów dla każdej z metod
disp('Zagregowany błąd względny dla każdej metody od 2a do 2d:');
disp(errors);

%% Zad4

% Wygenerowanie wektora h
k = 10000:-100:100;
h_values = 1 ./ k;

% Sprawdzanie, czy każde 1/h jest liczbą całkowitą
valid_indices = mod(1 ./ h_values, 1) == 0;

% Wybieranie tylko tych wartości h, dla których 1/h jest liczbą całkowitą
h_values = h_values(valid_indices);

y_a = cell(1, length(h_values));
y_b = cell(1, length(h_values));
y_c = cell(1, length(h_values));
y_d = cell(1, length(h_values));
t_a = cell(1, length(h_values));
t_b = cell(1, length(h_values));
t_c = cell(1, length(h_values));
t_d = cell(1, length(h_values));

for i = 1:length(h_values)

    h = h_values(i);

    % (a) Otwarta metoda Eulera

    N = 1/h;
    t = 0:h:1;
    x = zeros(1, N+1);
    y = zeros(1, N+1);
    x(1) = x0;
    y(1) = y0;

    for n = 1:N
        x(n+1) = x(n) + h * (p1 * x(n) - p2 * x(n) * y(n));
        y(n+1) = y(n) + h * (p3 * x(n) * y(n) - p4 * y(n));
    end

    y_a{i} = y;
    t_a{i} = t;

    % (b) Zamknięta metoda Eulera

    % Przestrzeń czasowa
    N = length(t);
    
    % Inicjalizacja wektorów
    x = zeros(N, 1);
    y = zeros(N, 1);
    x(1) = x0;
    y(1) = y0;

    % Rozwiązanie za pomocą zamkniętej metody Eulera
    for n = 1:N-1
        % Aktualne wartości
        x_current = x(n);
        y_current = y(n);
    
        % Funkcja dla fsolve
        fun = @(z) [
            z(1) - x_current - h * (p1 * z(1) - p2 * z(1) * z(2));
            z(2) - y_current - h * (p3 * z(1) * z(2) - p4 * z(2))
        ];
    
        % Początkowe przybliżenie
        z0 = [x_current; y_current];
    
        % Rozwiązanie układu równań
        z_next = fsolve(fun, z0, optimset('Display', 'off'));
    
        % Zapisanie wyników
        x(n + 1) = z_next(1);
        y(n + 1) = z_next(2);
    end

    y_b{i} = y;
    t_b{i} = t;

    % (c) Metodą Adamsa-Bashfortha rzędu 2

    N = 1/h;
    x = zeros(1, N+1);
    y = zeros(1, N+1);
    x(1) = x0;
    y(1) = y0;

    x(2) = x(1) + h * (p1 * x(1) - p2 * x(1) * y(1));
    y(2) = y(1) + h * (p3 * x(1) * y(1) - p4 * y(1));

    for n = 2:N
        fx1 = p1 * x(n) - p2 * x(n) * y(n);
        fy1 = p3 * x(n) * y(n) - p4 * y(n);
        fx2 = p1 * x(n-1) - p2 * x(n-1) * y(n-1);
        fy2 = p3 * x(n-1) * y(n-1) - p4 * y(n-1);
        x(n+1) = x(n) + h/2 * (3 * fx1 - fx2);
        y(n+1) = y(n) + h/2 * (3 * fy1 - fy2);
    end

    y_c{i} = y;
    t_c{i} = t;

    % (d) Dwukrokową metodą trapezów

    x = zeros(1, N+1);
    y = zeros(1, N+1);
    x(1) = x0;
    y(1) = y0;

    % Metoda Eulera na pierwszym kroku
    fx1 = p1 * x(1) - p2 * x(1) * y(1);
    fy1 = p3 * x(1) * y(1) - p4 * y(1);
    x(2) = x(1) + h * fx1;
    y(2) = y(1) + h * fy1;
    
    for n = 2:N
        % Predykcja
        fx1 = p1 * x(n) - p2 * x(n) * y(n); % Obliczanie fx1
        fy1 = p3 * x(n) * y(n) - p4 * y(n); % Obliczanie fy1
        x_predict = x(n) + h * fx1; % Predykcja x
        y_predict = y(n) + h * fy1; % Predykcja y
        
        % Korekcja
        fx2 = p1 * x_predict - p2 * x_predict * y_predict; % Obliczanie fx2
        fy2 = p3 * x_predict * y_predict - p4 * y_predict; % Obliczanie fy2
        x(n+1) = x(n) + h/2 * (fx1 + fx2); % Aktualizacja x
        y(n+1) = y(n) + h/2 * (fy1 + fy2); % Aktualizacja y
    end

    y_d{i} = y;
    t_d{i} = t;

end

% Inicjalizacja błędów dla każdej metody
errors_a = zeros(1, length(h_values));
errors_b = zeros(1, length(h_values));
errors_c = zeros(1, length(h_values));
errors_d = zeros(1, length(h_values));

for i = 1:length(h_values)
    % Dla każdej metody oblicz estymaty y_i
    for method = 1:4
        switch method
            case 1 % Otwarta metoda Eulera
                y=y_a{i};
                t=t_a{i};
            case 2 % Zamknięta metoda Eulera
                y=y_b{i};
                t=t_b{i};
            case 3 % Metoda Adamsa-Bashfortha rzędu 2
                y=y_c{i};
                t=t_c{i};
            case 4 % Dwukrotna metoda trapezów
                y=y_d{i};
                t=t_d{i};
        end
    
        % Rozwiązanie referencyjne za pomocą funkcji ode45
        [t_ref, z_ref] = ode45(@(t, y) [p1*y(1) - p2*y(1)*y(2); p3*y(1)*y(2) - p4*y(2)], [0 1], [x0; y0], odeset('RelTol', RelTol, 'AbsTol', AbsTol));
        
        y_ref = z_ref(:,2);

        % Interpolacja estymowanych wartości na t_ref
        y_interp = interp1(t, y, t_ref, 'linear');
    
        % Obliczanie błędu względnego
        error = sqrt(sum((y_ref - y_interp).^2)) / sqrt(sum(y_ref.^2));

        switch method
            case 1 % Otwarta metoda Eulera
               errors_a(i) = error;
            case 2 % Zamknięta metoda Eulera
               errors_b(i) = error;
            case 3 % Metoda Adamsa-Bashfortha rzędu 2
               errors_c(i) = error;
            case 4 % Dwukrotna metoda trapezów
               errors_d(i) = error;
        end
    end
end

figure;
loglog(h_values, errors_a, 'Color', [0, 0.4470, 0.7410]);
hold on;
loglog(h_values, errors_b, 'Color', [0.8500, 0.3250, 0.0980]);
loglog(h_values, errors_c, 'Color', [0.4660, 0.6740, 0.1880]);
loglog(h_values, errors_d, 'Color', [0.4940, 0.1840, 0.5560]);
hold off;

legend({'Otwarta metoda Eulera', 'Zamknięta metoda Eulera', 'Metoda Adamsa-Bashfortha rzędu 2', 'Dwukrokowa metoda trapezów'}, 'Location', 'best');
xlabel('Długość kroku całkowania h');
ylabel('Zagregowany błąd względny y');
% title('Zależność długości kroku całkowania od zagregowanego błędu względnego y');
grid on;
%% Zad5
% Wczytanie danych
data = readmatrix('MNUB_24L_P3_dane17.csv');
t_data = data(:, 1);
x_data = data(:, 2);
y_data = data(:, 3);

% Początkowe wartości parametrów
p0 = [p1, p2, p3, p4];

% Początkowe warunki
x0 = x_data(1);
y0 = y_data(1);

% Definiowanie funkcji kosztu
cost_function = @(p) calculate_cost(p, t_data, x_data, y_data, [x0, y0]);

% Opcje dla fminsearch (opcjonalnie)
options = optimset('Display', 'iter', 'PlotFcns', @optimplotfval);

% Minimalizacja kosztu
p_opt = fminsearch(cost_function, p0, options);

% Wyświetlenie optymalnych wartości parametrów
fprintf('Optymalne wartości parametrów: p1 = %f, p2 = %f, p3 = %f, p4 = %f\n', p_opt);

function cost = calculate_cost(p, t_data, x_data, y_data, initial_conditions)
    % Rozwiązanie układu równań różniczkowych
    [~, result] = ode45(@(t, z) lotka_volterra5(t, z, p), t_data, initial_conditions);
    
    % Interpolacja wyników
    x_interp = interp1(t_data, result(:, 1), t_data);
    y_interp = interp1(t_data, result(:, 2), t_data);
    
    % Obliczanie wartości funkcji kosztu
    cost = sum((x_data - x_interp).^2) + sum((y_data - y_interp).^2);
end

function dzdt = lotka_volterra5(~, z, p)
    dzdt = [p(1) * z(1) - p(2) * z(1) * z(2);
            p(3) * z(1) * z(2) - p(4) * z(2)];
end
