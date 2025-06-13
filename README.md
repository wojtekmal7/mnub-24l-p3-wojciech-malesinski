
# Analiza układu równań Lotki-Volterry

Ten projekt zawiera implementację różnych metod numerycznych do rozwiązania układu równań różniczkowych opisujących populacje drapieżnika i ofiary (model Lotki-Volterry). 

## Zawartość

- **Zadanie 1**: Rozwiązanie za pomocą `ode45`.
- **Zadanie 2**: Rozwiązania metodami:
  - Euler otwarty (jawny)
  - Euler zamknięty (niejawny)
  - Adams-Bashforth rzędu 2
  - Dwukrokowa metoda trapezów
- **Zadanie 3**: Porównanie błędów względnych każdej metody względem rozwiązania referencyjnego.
- **Zadanie 4**: Analiza wpływu kroku całkowania na błąd względny metod z Zadania 2.
- **Zadanie 5**: Dopasowanie parametrów modelu do danych pomiarowych z pliku `MNUB_24L_P3_dane17.csv` za pomocą optymalizacji (minimalizacja funkcji kosztu).

## Uruchamianie

Otwórz plik `.m` w środowisku MATLAB i uruchom całość. Dane są wczytywane z pliku `MNUB_24L_P3_dane17.csv` — upewnij się, że plik znajduje się w tym samym katalogu co skrypt.

## Dane wejściowe

- `MNUB_24L_P3_dane17.csv`: zawiera dane empiryczne do dopasowania modelu w postaci trzech kolumn:
  1. Czas
  2. Populacja x (ofiary)
  3. Populacja y (drapieżnika)

## Wymagane pakiety

Skrypt korzysta z:
- `ode45` (solver ODE)
- `fsolve` (rozwiązywanie układów nieliniowych)
- `fminsearch` (optymalizacja nieliniowa)

Upewnij się, że posiadasz pakiet Optimization Toolbox.

## Autorzy

- Wojciech Malesiński

