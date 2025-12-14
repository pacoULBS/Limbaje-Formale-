from collections import defaultdict, deque
import csv
import argparse
import pprint
import sys

ENDMARK = '$'
EPS = 'ε'  # simbol folosit intern pentru epsilon


# -------------------------
# Parsare gramatică din text
# -------------------------
def parse_grammar(text):
    """
    Parsează o gramatică exprimată liniar (ex: "E -> E + T | T").

    Returnează un dictionar: Nonterminal -> listă de producții (fiecare producție este o listă de simboluri).
    - Productiile vide se reprezintă ca listă vidă.
    """
    G = defaultdict(list)
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '->' not in line:
            raise ValueError(f"Lipsea '->' pe linia: {line}")
        left, right = line.split('->', 1)
        A = left.strip()
        alts = [alt.strip() for alt in right.split('|')]
        for alt in alts:
            if alt == '' or alt in ('ε', 'eps', 'empty'):
                rhs = []
            else:
                rhs = alt.split()
            G[A].append(rhs)
    return dict(G)


# ---------------------------------------
# Determinare terminale și non-terminale
# ---------------------------------------
def compute_terminals_and_nonterminals(G):
    """
    Returnează două mulțimi: (terminale, nonterminale)
    """
    nonterms = set(G.keys())
    terms = set()
    for A, prods in G.items():
        for rhs in prods:
            for sym in rhs:
                if sym == EPS:
                    continue
                if sym not in nonterms:
                    terms.add(sym)
    return terms, nonterms


# -----------------------------
# Calcul FIRST pentru simboluri
# -----------------------------
def compute_first_sets(G, terminals, nonterms):
    """
    Calculează FIRST pentru toate simbolurile (folosit la construcția itemilor LR(1)).
    RETURN: dictionar simbol -> set(terminale sau EPS)
    """
    FIRST = {t: {t} for t in terminals}
    for N in nonterms:
        FIRST[N] = set()
    FIRST[EPS] = {EPS}

    changed = True
    while changed:
        changed = False
        for A, prods in G.items():
            for rhs in prods:
                if not rhs:
                    if EPS not in FIRST[A]:
                        FIRST[A].add(EPS)
                        changed = True
                    continue
                add_eps = True
                for sym in rhs:
                    before = len(FIRST[A])
                    if sym not in FIRST:
                        FIRST[sym] = {sym}
                    FIRST[A].update(x for x in FIRST[sym] if x != EPS)
                    after = len(FIRST[A])
                    if after > before:
                        changed = True
                    if EPS in FIRST[sym]:
                        add_eps = True
                        continue
                    else:
                        add_eps = False
                        break
                if add_eps:
                    if EPS not in FIRST[A]:
                        FIRST[A].add(EPS)
                        changed = True
    return FIRST


def first_of_sequence(seq, FIRST):
    """
    FIRST pentru o secvență de simboluri (folosit în closure).
    Returnează un set de terminale (sau EPS dacă secvența deriva epsilon).
    """
    if not seq:
        return {EPS}
    res = set()
    for sym in seq:
        if sym not in FIRST:
            res.add(sym)
            return res
        res.update(x for x in FIRST[sym] if x != EPS)
        if EPS in FIRST[sym]:
            continue
        else:
            break
    else:
        res.add(EPS)
    return res


# -----------------------
# Reprezentare item LR(1)
# -----------------------
def make_item(A, rhs, dot, la):
    # Un item e un tuplu: (A, rhs_tuple, dot_pos, lookahead)
    return (A, tuple(rhs), dot, la)


def item_to_str(item):
    # Reprezentare textuală pentru afișare
    A, rhs, dot, la = item
    rhs_list = list(rhs)
    rhs_with_dot = rhs_list[:]
    rhs_with_dot.insert(dot, '•')
    rhs_s = ' '.join(rhs_with_dot).replace(' •', '•')
    return f"[{A} -> {rhs_s}, {la}]"


# -------------------------
# closure și goto (LR(1))
# -------------------------
def closure(items, G, FIRST):
    """
    Construiește închiderea LR(1) pentru un set de itemi.
    - items: set de itemi LR(1)
    - G: gramatică
    - FIRST: seturi FIRST calculate
    """
    I = set(items)
    added = True
    while added:
        added = False
        new_items = set()
        for item in I:
            A, rhs, dot, la = item
            if dot < len(rhs):
                B = rhs[dot]
                if B in G:
                    beta = list(rhs[dot+1:])
                    for prod in G[B]:
                        seq = beta + [la]
                        first_seq = first_of_sequence(seq, FIRST)
                        for b in first_seq:
                            if b == EPS:
                                continue
                            new_item = make_item(B, prod, 0, b)
                            if new_item not in I:
                                new_items.add(new_item)
        if new_items:
            I |= new_items
            added = True
    return frozenset(I)


def goto(I, X, G, FIRST):
    """
    GOTO pentru colecția de itemi: mută dot-ul peste simbolul X și apoi face closure.
    """
    moved = set()
    for item in I:
        A, rhs, dot, la = item
        if dot < len(rhs) and rhs[dot] == X:
            moved.add(make_item(A, list(rhs), dot+1, la))
    if not moved:
        return frozenset()
    return closure(moved, G, FIRST)


# ----------------------------------------
# Construire colecție canonică LR(1)
# ----------------------------------------
def canonical_LR1_collection(G):
    """
    Construiește stările (seturi de itemi) și tranzițiile pentru automatul LR(1).
    Returnează: states (list of frozenset(items)), transitions dict((state_id, symbol) -> state_id),
               G_aug (gramatică augmentată), S_prime (start augmentat), FIRST
    """
    starts = list(G.keys())
    if not starts:
        raise RuntimeError("Gramatică goală")
    start = starts[0]
    # Augmentăm gramatica: S' -> S
    G_aug = {k: [list(p) for p in v] for k, v in G.items()}
    S_prime = start + "'"
    i = 1
    while S_prime in G_aug:
        S_prime = f"{start}'{i}"
        i += 1
    G_aug[S_prime] = [[start]]

    terminals, nonterms = compute_terminals_and_nonterminals(G_aug)
    FIRST = compute_first_sets(G_aug, terminals, set(G_aug.keys()))

    init_item = make_item(S_prime, [start], 0, ENDMARK)
    I0 = closure({init_item}, G_aug, FIRST)
    states = [I0]
    state_ids = {I0: 0}
    transitions = dict()

    q = deque([I0])
    while q:
        I = q.popleft()
        sid = state_ids[I]
        symbols = set()
        for item in I:
            A, rhs, dot, la = item
            if dot < len(rhs):
                symbols.add(rhs[dot])
        for X in symbols:
            J = goto(I, X, G_aug, FIRST)
            if not J:
                continue
            if J not in state_ids:
                state_ids[J] = len(states)
                states.append(J)
                q.append(J)
            transitions[(sid, X)] = state_ids[J]
    return states, transitions, G_aug, S_prime, FIRST


# --------------------------------------
# Construire tabele ACTION și GOTO
# --------------------------------------
def build_parsing_table(states, transitions, G_aug, S_prime, FIRST):
    """
    Construiește structurile ACTION și GOTO:
    - ACTION: dict (state, terminal) -> ("shift", s) / ("reduce", (A, rhs)) / ("accept",)
    - GOTO: dict (state, nonterminal) -> state
    """
    ACTION = dict()
    GOTO = dict()
    conflicts = []

    prods = []
    for A, plist in G_aug.items():
        for rhs in plist:
            prods.append((A, tuple(rhs)))

    terminals, nonterms = compute_terminals_and_nonterminals(G_aug)
    terminals = set(terminals)
    nonterms = set(G_aug.keys())

    for sid, I in enumerate(states):
        for item in I:
            A, rhs, dot, la = item
            if dot < len(rhs):
                a = rhs[dot]
                # dacă simbolul după dot este terminal => shift
                if a not in nonterms:
                    if (sid, a) in transitions:
                        t = transitions[(sid, a)]
                        key = (sid, a)
                        entry = ("shift", t)
                        if key in ACTION and ACTION[key] != entry:
                            conflicts.append(("conflict", sid, a, ACTION[key], entry))
                        else:
                            ACTION[key] = entry
            else:
                # dot la dreapta => reducere (sau accept dacă e S' -> S . , $)
                if A == S_prime:
                    if la == ENDMARK:
                        ACTION[(sid, ENDMARK)] = ("accept",)
                else:
                    key = (sid, la)
                    entry = ("reduce", (A, tuple(rhs)))
                    if key in ACTION and ACTION[key] != entry:
                        conflicts.append(("conflict", sid, la, ACTION[key], entry))
                    else:
                        ACTION[key] = entry
        # GOTO pentru neterminale
        for B in nonterms:
            if (sid, B) in transitions:
                GOTO[(sid, B)] = transitions[(sid, B)]
    return ACTION, GOTO, conflicts


# -------------------------
# Formatare celulă ACTION
# (folosit pentru printare)
# -------------------------
def format_action_cell_for_print(v):
    if not v:
        return ""
    if v[0] == "shift":
        return f"s{v[1]}"
    elif v[0] == "reduce":
        A, rhs = v[1]
        rhs_s = " ".join(rhs) if rhs else EPS
        return f"r[{A}->{rhs_s}]"
    elif v[0] == "accept":
        return "acc"
    else:
        return str(v)


# ---------------------------------------------------
# Export în formatele cerute de LRParser2 și table_reader
# ---------------------------------------------------
def export_action_and_prod_tables(action, goto, terminals, nonterms, num_states, prods_list, S_prime, filename_action="action_table.csv", filename_prod="result.csv"):
    """
    Scrie două fișiere:
      - action_table.csv : fiecare rând = simbol, coloane = stare 0..N-1
          - Terminal rows: folosesc prefixul 'd' pentru shift (ex: d3), 'r<nr>' pentru reduceri (ex: r2), 'acc' pentru accept.
          - Nonterminal rows: conțin numere GOTO (starea) sau gol.
      - result.csv : lista producțiilor numerotate (prod_num, LHS, RHS)
          - Notă: numerotarea produselor exclude producția augmentată S' -> S.
    Returnează dicționarul prod_index (mapping (A, rhs) -> prod_num ca string).
    """
    # Construim numerotarea producțiilor (începând de la 1) - LRParser2 se așteaptă la chei ca string
    prod_index = {}
    prod_rows = []
    idx = 1
    for A, rhs in prods_list:
        # sărim producțiile augmentate (S' -> S)
        if A == S_prime:
            continue
        key = (A, tuple(rhs))
        if key not in prod_index:
            prod_index[key] = str(idx)
            rhs_s = " ".join(rhs) if rhs else ""
            # reținem rândul pentru result.csv: prod_num, LHS, RHS
            prod_rows.append((str(idx), A, rhs_s))
            idx += 1

    # Simboluri pentru rânduri în action_table: terminale + $ și apoi nonterminale (GOTO)
    ter_list = sorted(list(terminals) + [ENDMARK])
    non_list = sorted(list(nonterms))

    # Scriem action_table.csv (fără antet, LRParser2 se așteaptă la rânduri de tip simbol,...)
    with open(filename_action, "w", newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        # Terminale (inclusiv '$') - fiecare rând începe cu simbolul
        for t in ter_list:
            row = [t]
            for s in range(num_states):
                cell = ""
                entry = action.get((s, t), "")
                if entry:
                    if entry[0] == "shift":
                        # prefix 'd' pentru compatibilitate cu LRParser2 (care caută 'd...' la shift)
                        cell = "d" + str(entry[1])
                    elif entry[0] == "reduce":
                        A, rhs = entry[1]
                        prod_num = prod_index.get((A, tuple(rhs)))
                        if prod_num:
                            cell = "r" + prod_num
                        else:
                            cell = "r?"  # fallback, arată că există o reducere dar fără număr
                    elif entry[0] == "accept":
                        cell = "acc"
                    else:
                        cell = str(entry)
                row.append(cell)
            writer.writerow(row)
        # Neterminale -> rânduri GOTO (valoarea celulei e numărul stării sau gol)
        for N in non_list:
            row = [N]
            for s in range(num_states):
                dst = goto.get((s, N), "")
                row.append(str(dst) if dst != "" else "")
            writer.writerow(row)

    # Scriem result.csv (producții), rând per producție: prod_num, LHS, RHS
    with open(filename_prod, "w", newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        for pr in prod_rows:
            writer.writerow(pr)

    print(f"Am creat fișierele: {filename_action} și {filename_prod}")
    return prod_index


# Exemplu de gramatică folosit când nu este furnizat fișier
def example_grammar_text():
    return """
S -> E
E -> E + T | T
T -> T * F | F
F -> ( E ) | id
""".strip()


# -------------------------
# Entry point (main)
# -------------------------
def main():
    parser = argparse.ArgumentParser(description="Generator tabele LR(1) și export action_table.csv + result.csv.")
    parser.add_argument("--file", "-f", help="Fișier text cu gramatica", default=None)
    parser.add_argument("--action", help="Nume fișier action table output", default="action_table.csv")
    parser.add_argument("--prod", help="Nume fișier productions output", default="result.csv")
    args = parser.parse_args()

    if args.file:
        try:
            with open(args.file, encoding='utf-8') as fh:
                text = fh.read()
        except FileNotFoundError:
            print(f"Eroare: fișierul '{args.file}' nu a fost găsit.", file=sys.stderr)
            sys.exit(2)
    else:
        # Dacă nu s-a dat fișier, folosim o gramatică exemplu pentru demonstrație.
        print("Nu s-a dat fișier; folosesc exemplu clasic (E -> E + T ...).")
        text = example_grammar_text()

    # Parsăm gramatică și afișăm sumar (utile pentru debugging)
    G = parse_grammar(text)
    print("Gramatică parsată (sumar):")
    pprint.pprint(G)
    print()

    # Construim colecția LR(1)
    states, transitions, G_aug, S_prime, FIRST = canonical_LR1_collection(G)
    print(f"Gramatică augmentată start: {S_prime}")
    print(f"Număr stări: {len(states)}")

    # Construim tabelele ACTION și GOTO
    ACTION, GOTO, conflicts = build_parsing_table(states, transitions, G_aug, S_prime, FIRST)

    # Construim lista de producții (folosită pentru numerotare și pentru result.csv)
    prods_list = []
    for A, plist in G_aug.items():
        for rhs in plist:
            prods_list.append((A, list(rhs)))

    terminals, nonterms = compute_terminals_and_nonterminals(G_aug)

    # Exportăm action_table.csv și result.csv (compatibil cu LRParser2)
    prod_index = export_action_and_prod_tables(ACTION, GOTO, terminals, nonterms, len(states), prods_list, S_prime, filename_action=args.action, filename_prod=args.prod)

    # Afișăm conflictele (dacă există) pentru debugging
    if conflicts:
        print("Conflicte detectate:")
        for c in conflicts:
            print(" ", c)
    else:
        print("Niciun conflict detectat.")


if __name__ == "__main__":
    main()