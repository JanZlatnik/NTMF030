# Zápočtový problém z NTMF030 Kvantová teorie rozptylu

### Stručný úvod
  - samotný výpočetní program je napsán v Rustu
  - pro vytvoření grafů je využíván Python, konkrétně knihovna `matplotlib`
  - Mathematica notebook [`Well.nb`](Well.nb) obsahuje velmi jednoduchý program pro orientační analytické ověření výsledků pro potenciálovou jámu
  - soubor [`QST I.pdf`](QST%20I.pdf) obsahuje řešení zápočtového problému s komentářem k jednotlivým bodům

### Spuštění programu
  - složka [`target/release`](target/release/) obsahuje aplikaci [`Potential_scattering.exe`](target/release/Potential_scattering.exe), která provádí výpočet
  - ke svému spuštění vyžaduje aplikace přítomnost [`settings.toml`](settings.toml) ve stejném adresáři, tento soubor obsahuje parametry výpočtu (hmotnost, počítané parciální vlny, síly potenciálové jámy,...) a výpočetní nastavení (délka gridu, spacing gridu,...)
  - výsledky výpočtů jsou automaticky ukládány do složek odpovídajících jednotlivým úlohám ([`Task 1`](Task%201/),[`Task 2`](Task%202/) a [`Test`](Test/) )

### Zdrojový kód
  - zdrojové kódy pro aplikaci [`Potential_scattering.exe`](target/release/Potential_scattering.exe) se nachází ve složce [`src`](src/)
  - modul [`core.rs`](src/core.rs) -> obsahuje základní fyzikální konstanty a inicializace nastavení
  - modul [`math.rs`](src/math.rs) -> zde jsou naprogramovány základní metody jako tvorba gridu, integrace, Riccati-Besselovy funkce
  - modul [`scattering.rs`](src/scattering.rs) -> zde je implementována metoda pro výpočet fázových posunutí a jejich následné převedení na účinné průřezy
  - modul [`bound_states.rs`](src/bound_states.rs) -> implementuje metodu pro výpočet vázaných stavů pomocí nelezení nul matching-funkce bisekcí
  - modul [`main.rs`](src/main.rs) -> spouští jednotlivé metody a ukládá data do souborů
  - soubor [`Cargo.toml`](Cargo.toml) -> obsahuje informace o použitých knihovnách, využity jsou knihovny pro práci se soubory, časem, numerické vlastnosti a knihovna pro numerickou algebru

### Tvorba grafů
  - skripty pro tvorbu grafů jsou napsány v Pythonu a spolu s grafy se nachází ve složce [`scripts`](scripts/)
  - jedná se převážně o tu stejnou šablonu pouze několikrát zkopírovanou pro snadnější úpravu vzniklých grafů
    
