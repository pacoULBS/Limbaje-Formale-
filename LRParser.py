import csv


def read_action_table():
    action_table = {}
    with open('tema 2\\action_table.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            symbol = row[0]  # simbol (id, +, etc.)
            values = row[1:]  # restul valorilor din rand
            action_table[symbol] = values
    return action_table

def read_prod():
    prod = {}
    with open('tema 2\\productions.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            symbol = row[0]
            values = row[1:]
            prod[symbol] = values
    return prod

action_table = read_action_table()
print("Action Table:")
for key, values in action_table.items():
    print(f"{key}: {','.join(values)}")

prod = read_prod()
print("\nProductions:")
for key, values in prod.items():
    print(f"{key}: {','.join(values)}")

def parse_input(input_string):
    state_stack = [0]  # initializam stiva de stari cu 0
    token_stack = ['$'] #initializam stiva de tokeni cu simbolul de start
    input_tokens = input_string.split() + ['$']  # adaugam simbolul de sfarsit
    pointer = 0  # pointer pentru input_tokens

    while True:
        current_state = state_stack[-1] # starea curenta
        current_token = input_tokens[pointer] # tokenul curent

        try:
            action = action_table.get(current_token, [])[current_state] # actiunea din tabel
        except IndexError: # daca starea nu exista in tabel
            print("Input rejected: Invalid state or token.") # invalidam inputul
            return False

        if action.startswith('d'):  # actiune de deplasare
            next_state = int(action[1:]) # extragem starea urmatoare
            state_stack.append(next_state) # adaugam starea urmatoare in stiva
            token_stack.append(current_token) #adaugam simbolul curent
            pointer += 1 # consumam tokenul curent
            print(f"Shift: Move to state {next_state}, consume '{current_token}'") 
        elif action.startswith('r'):  # actiune de reductie
            prod_number = action[1:] # extragem numarul productiei
            entry = prod.get(prod_number)   # obtinem productia corespunzatoare
            if not entry: # daca productia nu exista
                print(f"Input rejected: Unknown production {prod_number}.") #invalidam inputul
                return False

            lhs = entry[0] # partea stanga a productiei
            rhs = entry[1].strip() if len(entry) > 1 else '' # partea dreapta a productiei
            rhs_tokens = rhs.split() if rhs else [] # tokenii din partea dreapta
            rhs_length = len(rhs_tokens) # lungimea partii drepte

            for _ in range(rhs_length): # eliminam elementele din stiva conform lungimii partii drepte
                state_stack.pop()
                token_stack.pop()
            
            token_stack.append(lhs) #punem pe stiva simbolul la care s-a ajuns (din reductie)

            try: 
                goto_cell = action_table.get(lhs, [])[state_stack[-1]] # obtinem starea din tabela de salt
                goto_state = int(goto_cell) # convertim la int
            except (IndexError, ValueError): # daca starea nu exista sau nu e valida
                print("Input rejected: Invalid goto for reduction.") #invalidam inputul
                return False

            state_stack.append(goto_state) # adaugam starea de salt in stiva
            rhs_display = rhs if rhs else 'ε' # afisam ε daca partea dreapta e vida
            print(f"Reduce: Using production {prod_number}: {lhs} -> {rhs_display}, goto state {goto_state}") 
        elif action in ('acc', 'accept'):  # acceptam inputul
            print("Input accepted.")
            return True
        else:
            print("Input rejected.") # invalidam inputul
            return False


def testGood():
    if (
        parse_input("id * id + id")==True
     ):
        return True
    else:
        return False

def testBad():
    if(
        parse_input("id * * + id")==False
    ):
        return True
    else:
        return False

def runTests():
    if(testGood()):
        if(testBad()):
            return print("Test passed for string: id * id + id \n Tests passed for bad string: id * * + id")
        else: 
            return print("Test passed for string: id * id + id \n Test failed for bad string: id * * + id")
    else:
        return print("Test failed for string: id * id + id \n Test failed for bad string: id * * + id")
    

#run tests:
# runTests()


# fara input de la tastatura:
# parse_input("id + id * id")

#input de la tastatura:
parse_input(input("Introduceti un sir de terminale: "))