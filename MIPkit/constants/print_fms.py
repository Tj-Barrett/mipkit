from MIPkit.constants.constants import encode_fms, fm2smiles

def print_smiles():
    fm2sm = fm2smiles()

    for fm in fm2sm:
        fm = fm
        sm = fm2sm[fm]

        print(fm.ljust(7)+" : "+sm.rjust(30))

def print_fms():
    enfm = encode_fms()

    for fm in enfm:
        fm = fm
        en = enfm[fm]

        print(fm.ljust(7)+" : "+en.rjust(4))
