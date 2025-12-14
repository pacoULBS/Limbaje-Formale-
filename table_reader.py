# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# lr1_table_pretty_print.py

# Citește un result.csv (formatul generat de lr1_table_generator_csv.py) și afișează
# în consolă un tabel "human-readable" frumos al tabelelor ACTION / GOTO,
# cu aliniere și wrapping/indentare pentru celule lungi.

# Usage:
#   python3 lr1_table_pretty_print.py --file result.csv
#   python3 lr1_table_pretty_print.py --file result.csv --wrap 40   # limitează lățimea coloanelor la 40 caractere
#   python3 lr1_table_pretty_print.py --file result.csv --no-wrap  # fără wrapping, colanele vor fi largi după nevoie

# Observații:
#  - CSV-ul trebuie să aibă prima coloană "state", urmată de terminale (inclusiv $)
#    și apoi neterminale (GOTO).
#  - Scriptul tratează celulele goale corect și afișează r[A->rhs], sN, acc sau numere GOTO.
# """
# from __future__ import annotations
# import csv
# import argparse
# import shutil
# import textwrap
# import sys
# from typing import List

# def read_csv_rows(path: str) -> (List[str], List[List[str]]):
#     with open(path, newline='', encoding='utf-8') as f:
#         reader = csv.reader(f)
#         rows = [r for r in reader]
#     if not rows:
#         raise RuntimeError("Fișier CSV gol")
#     header = rows[0]
#     data = rows[1:]
#     return header, data

# def compute_column_widths(header: List[str], data: List[List[str]], wrap: bool, wrap_limit: int|None, pad: int=1):
#     cols = len(header)
#     # lungimi initiale bazate pe conținut
#     maxlens = [len(h) for h in header]
#     for r in data:
#         # dacă r are mai puține coloane, completează cu stringuri goale
#         if len(r) < cols:
#             r += [''] * (cols - len(r))
#         for i, cell in enumerate(r[:cols]):
#             l = len(cell or '')
#             if l > maxlens[i]:
#                 maxlens[i] = l
#     # aplicare wrap limit dacă e cerut (nu extinde coloana peste wrap_limit)
#     col_widths = list(maxlens)
#     if wrap and wrap_limit and wrap_limit > 0:
#         for i, w in enumerate(col_widths):
#             if w > wrap_limit:
#                 col_widths[i] = wrap_limit
#     # adăugăm padding stânga/dreapta
#     col_widths = [w + 2*pad for w in col_widths]
#     return col_widths

# def wrap_cell_lines(cell: str, width: int) -> List[str]:
#     # width include padding: scădem 2 pentru padding implicit (1 stânga, 1 dreapta)
#     inner_width = max(1, width - 2)
#     if cell is None:
#         cell = ''
#     wrapped = textwrap.wrap(cell, inner_width) or ['']
#     # adaugăm spații laterale (1 spatiu stg/dreapta)
#     return [f" {line.ljust(inner_width)} " for line in wrapped]

# def print_table(header: List[str], data: List[List[str]], col_widths: List[int], pad: int=1):
#     cols = len(header)
#     # pregătim linii împachetate pentru header (fără wrap, header rămâne pe o singură linie)
#     header_cells = []
#     for i, h in enumerate(header):
#         w = col_widths[i]
#         inner_w = w - 2*pad
#         cell = h.center(inner_w)
#         header_cells.append(f"{' '*pad}{cell}{' '*pad}")
#     sep_col = " | "
#     header_line = sep_col.join(header_cells)
#     total_width = sum(col_widths) + len(sep_col) * (cols - 1)
#     print(header_line)
#     print("-" * total_width)

#     # pentru fiecare rând, construim lista de linii per celulă (după wrap) și apoi le afișăm rând cu rând
#     for r in data:
#         if len(r) < cols:
#             r += [''] * (cols - len(r))
#         cell_lines_per_col = []
#         for i in range(cols):
#             w = col_widths[i]
#             cell = r[i] if i < len(r) else ''
#             # aplicăm wrap folosind inner width (width minus paddings)
#             lines = wrap_cell_lines(cell, w)
#             cell_lines_per_col.append(lines)
#         # înălțimea rândului este maxul liniilor din coloane
#         row_height = max(len(lines) for lines in cell_lines_per_col)
#         for line_idx in range(row_height):
#             out_cells = []
#             for i in range(cols):
#                 lines = cell_lines_per_col[i]
#                 w = col_widths[i]
#                 # dacă nu există linie la indexul respectiv -> spațiu gol padding
#                 if line_idx < len(lines):
#                     out = lines[line_idx]
#                 else:
#                     inner_w = w - 2*pad
#                     out = ' ' * pad + ' ' * inner_w + ' ' * pad
#                 out_cells.append(out)
#             print(sep_col.join(out_cells))
#         # separator sub fiecare stare (opțional: liniă fină)
#         print("-" * total_width)

# def terminal_width() -> int:
#     try:
#         return shutil.get_terminal_size().columns
#     except Exception:
#         return 80

# def main():
#     parser = argparse.ArgumentParser(description="Afișează frumos result.csv (ACTION/GOTO).")
#     parser.add_argument("--file", "-f", help="Fișier CSV (implicit: result.csv)", default="result.csv")
#     parser.add_argument("--wrap", type=int, help="Lățimea maximă (în caractere) pentru coloane înainte de wrapping. 0 = fără wrap", default=30)
#     parser.add_argument("--no-wrap", dest='no_wrap', help="Dezactivează wrapping (echivalent cu --wrap 0)", action='store_true')
#     parser.add_argument("--pad", type=int, help="Număr spații padding stânga/dreapta în celule", default=1)
#     args = parser.parse_args()

#     csv_file = args.file
#     try:
#         header, data = read_csv_rows(csv_file)
#     except FileNotFoundError:
#         print(f"Fișierul '{csv_file}' nu a fost găsit.", file=sys.stderr)
#         sys.exit(2)
#     except Exception as e:
#         print(f"Eroare la citirea CSV: {e}", file=sys.stderr)
#         sys.exit(2)

#     # determinăm dacă facem wrap
#     do_wrap = not args.no_wrap and (args.wrap and args.wrap > 0)
#     wrap_limit = args.wrap if do_wrap else None

#     col_widths = compute_column_widths(header, data, wrap=do_wrap, wrap_limit=wrap_limit, pad=args.pad)

#     # Dacă total width e mai lat decât terminalul și nu avem wrap activ, afișăm un mic avert:
#     total_width = sum(col_widths) + (3 * (len(header) - 1))  # estimare sep " | "
#     termw = terminal_width()
#     if not do_wrap and total_width > termw:
#         print(f"Atentie: tabelul are {total_width} cols iar terminalul {termw}. Poți rula cu --wrap N pentru a limita lățimea coloanelor.\n")

#     print_table(header, data, col_widths, pad=args.pad)

# if __name__ == "__main__":
#     main()








#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
table_reader.py (updated)

Acest script poate citi atât:
 - formatul vechi (result.csv generat ca tabel state x symbols: prima col "state", următoarele coloane terminale + nonterminale și fiecare rând = stare)
 - noul format action_table.csv (fiecare rând = simbol, coloane = stare 0..N-1)

Dacă primește action_table.csv (simbol-per-rând), transformă în forma state-per-rând
și afișează tabelul frumos (human-readable), exact ca înainte.

Usage:
  python3 table_reader.py --file action_table.csv
  python3 table_reader.py --file result.csv   # dacă result.csv este în format state-per-row
  python3 table_reader.py --wrap 40
"""
from __future__ import annotations
import csv
import argparse
import shutil
import textwrap
import sys
from typing import List

def read_all_rows(path: str) -> List[List[str]]:
    with open(path, newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        rows = [r for r in reader]
    return rows

def detect_and_transform(rows: List[List[str]]) -> (List[str], List[List[str]]):
    """
    Detect format and, if needed, transform symbol-per-row -> state-per-row.

    Returns (header, data) where header is list of column names and data is list of rows (each row list).
    The returned format matches what the pretty-printer expects: header[0] == 'state', data rows start with state index.
    """
    if not rows:
        raise RuntimeError("Fișier CSV gol")

    # Heuristic: if first row first cell is 'state' (case-insensitive), assume old format (state-per-row, header present)
    first = rows[0]
    if first and first[0].strip().lower() == 'state' and len(rows) > 1:
        header = first
        data = rows[1:]
        return header, data

    # Another heuristic: if every row's length equals and first column of each row is a symbol (non-numeric),
    # and there are multiple rows, treat as symbol-per-row (action_table)
    # symbol-per-row layout: each row = [symbol, val_state0, val_state1, ...]
    # We'll transform to state-per-row: header = ['state', sym1, sym2, ...], rows: [0, cell(sym1,state0), cell(sym2,state0), ...], ...
    row0_len = len(rows[0])
    # Validate consistent column count
    consistent = all(len(r) == row0_len for r in rows)
    if consistent and row0_len >= 2:
        symbols = [r[0] for r in rows]
        num_states = row0_len - 1
        header = ['state'] + symbols
        data = []
        for s in range(num_states):
            row = [str(s)]
            for r in rows:
                # r is a symbol row: r[1 + s] is the cell for state s
                val = ""
                if 1 + s < len(r):
                    val = r[1 + s]
                row.append(val)
            data.append(row)
        return header, data

    # Fallback: maybe the CSV is already state-per-row but without header.
    # We try to detect: if first column values are numeric (state indices) across many rows, assume state-per-row without header
    all_first_numeric = all((len(r) > 0 and r[0].strip().isdigit()) for r in rows)
    if all_first_numeric:
        # create header as generic: 'state', col1..colN-1
        cols = len(rows[0])
        header = ['state'] + [f"col{i}" for i in range(1, cols)]
        data = rows
        return header, data

    # If we cannot detect, raise error
    raise RuntimeError("Format CSV necunoscut. Aștept fie 'state' header, fie action_table (simbol-per-rând), fie state-per-rând fără header.")

def compute_column_widths(header: List[str], data: List[List[str]], wrap: bool, wrap_limit: int|None, pad: int=1):
    cols = len(header)
    maxlens = [len(h) for h in header]
    for r in data:
        if len(r) < cols:
            r += [''] * (cols - len(r))
        for i, cell in enumerate(r[:cols]):
            l = len(cell or '')
            if l > maxlens[i]:
                maxlens[i] = l
    col_widths = list(maxlens)
    if wrap and wrap_limit and wrap_limit > 0:
        for i, w in enumerate(col_widths):
            if w > wrap_limit:
                col_widths[i] = wrap_limit
    col_widths = [w + 2*pad for w in col_widths]
    return col_widths

def wrap_cell_lines(cell: str, width: int, pad: int):
    inner_width = max(1, width - 2*pad)
    if cell is None:
        cell = ''
    wrapped = textwrap.wrap(cell, inner_width) or ['']
    return [f"{' '*pad}{line.ljust(inner_width)}{' '*pad}" for line in wrapped]

def print_table(header: List[str], data: List[List[str]], col_widths: List[int], pad: int=1):
    cols = len(header)
    header_cells = []
    for i, h in enumerate(header):
        w = col_widths[i]
        inner_w = w - 2*pad
        cell = h.center(inner_w)
        header_cells.append(f"{' '*pad}{cell}{' '*pad}")
    sep_col = " | "
    header_line = sep_col.join(header_cells)
    total_width = sum(col_widths) + len(sep_col) * (cols - 1)
    print(header_line)
    print("-" * total_width)

    for r in data:
        if len(r) < cols:
            r += [''] * (cols - len(r))
        cell_lines_per_col = []
        for i in range(cols):
            w = col_widths[i]
            cell = r[i] if i < len(r) else ''
            lines = wrap_cell_lines(cell, w, pad)
            cell_lines_per_col.append(lines)
        row_height = max(len(lines) for lines in cell_lines_per_col)
        for line_idx in range(row_height):
            out_cells = []
            for i in range(cols):
                lines = cell_lines_per_col[i]
                w = col_widths[i]
                if line_idx < len(lines):
                    out = lines[line_idx]
                else:
                    inner_w = w - 2*pad
                    out = ' ' * pad + ' ' * inner_w + ' ' * pad
                out_cells.append(out)
            print(sep_col.join(out_cells))
        print("-" * total_width)

def terminal_width() -> int:
    try:
        return shutil.get_terminal_size().columns
    except Exception:
        return 80

def main():
    parser = argparse.ArgumentParser(description="Afișează frumos action_table.csv sau result.csv (state x symbols).")
    parser.add_argument("--file", "-f", help="Fișier CSV (implicit: action_table.csv)", default="action_table.csv")
    parser.add_argument("--wrap", type=int, help="Lățimea maximă (în caractere) pentru coloane înainte de wrapping. 0 = fără wrap", default=30)
    parser.add_argument("--no-wrap", dest='no_wrap', help="Dezactivează wrapping (echivalent cu --wrap 0)", action='store_true')
    parser.add_argument("--pad", type=int, help="Număr spații padding stânga/dreapta în celule", default=1)
    args = parser.parse_args()

    csv_file = args.file
    try:
        rows = read_all_rows(csv_file)
    except FileNotFoundError:
        print(f"Fișierul '{csv_file}' nu a fost găsit.", file=sys.stderr)
        sys.exit(2)
    except Exception as e:
        print(f"Eroare la citirea CSV: {e}", file=sys.stderr)
        sys.exit(2)

    try:
        header, data = detect_and_transform(rows)
    except RuntimeError as e:
        print(f"Eroare: {e}", file=sys.stderr)
        sys.exit(2)

    do_wrap = not args.no_wrap and (args.wrap and args.wrap > 0)
    wrap_limit = args.wrap if do_wrap else None

    col_widths = compute_column_widths(header, data, wrap=do_wrap, wrap_limit=wrap_limit, pad=args.pad)

    total_width = sum(col_widths) + (3 * (len(header) - 1))
    termw = terminal_width()
    if not do_wrap and total_width > termw:
        print(f"Atentie: tabelul are {total_width} cols iar terminalul {termw}. Poți rula cu --wrap N pentru a limita lățimea coloanelor.\n")

    print_table(header, data, col_widths, pad=args.pad)

if __name__ == "__main__":
    main()