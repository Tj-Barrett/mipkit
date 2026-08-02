from MIPkit.constants.constants import res_list, atom2PT, res_atoms_pdb, res_atoms_gmx
from MIPkit.constants.constants import encode_fms, decode_fms, gmx_backbone
from datetime import datetime
from MIPkit.utils.read_pdb import read_pdb
import os

def print_backbone(mols, filename, conect=False, encode=True, amino='gmx', header=True):
    fpdb = open(f"{filename}", "w+")
    if header:
        fpdb.write("HEADER :  BMBT MIPkit \n")
        fpdb.write("AUTHOR :  Tj Barrett (tbar@tf.uni-kiel.de, barrett.t@northeastern.edu) \n")
        fpdb.write("TITLE  :  Molecule print from BMBT MIPkit\n")

    if amino == 'gmx':
        resatoms = res_atoms_gmx()
    elif amino == 'pdb':
        resatoms = res_atoms_pdb()
    residues = set(res_list())
    encodefm = encode_fms()
    decodefm = set(decode_fms())
    backbone = gmx_backbone()

    symbol = atom2PT()

    atmnum = 0

    for _, mol in enumerate(mols):
        n = 0
        offset = atmnum

        for _i, atm in enumerate(mol["atoms"]):
            # print(f'{mol["name"]} : {_i} : {atm}')

            # ATOM or HETATM label
            if mol["name"] in residues:
                a = "ATOM".ljust(6)  # atom#6s
            else:
                a = "HETATM".ljust(6)  # atom#6s

            # Atom Number
            b = str(atmnum).rjust(5)  # aomnum#5d

            # Atom Type
            if mol["name"] in residues:
                if mol["name"] in resatoms and atm != 'H' and n < len(resatoms[mol["name"]]):
                    new = resatoms[mol["name"]][n]
                    n+= 1
                elif mol["name"] in resatoms and atm != 'H' and n == len(resatoms[mol["name"]]):
                    # COOH Terminal Unit, rename Oxygen
                    new = 'OH'
                elif atm in symbol:
                    new = symbol[atm]
                else:
                    new = atm
                    
                c = new.center(4)  # atomname$#4s
            elif atm in symbol:
                new = symbol[atm]
                c = new.center(4)  # atomname$#4s
            else:
                new = atm
                c = atm.center(4)  # atomname$#4s

            # Residue Name
            if mol["name"] in residues:
                d = mol["name"].ljust(3)  # resname#1s
            elif mol["name"] in decodefm:
                d = mol["name"].ljust(3)  # resname#1s
            else:
                if encode:
                    efm = encodefm.get(mol["name"], mol["name"])
                else:
                    efm = mol["name"]
                d = efm.ljust(3)  # resname#1s

            # Chain
            try:
                e = mol['chain'].rjust(1)  # Astring
            except (KeyError, AttributeError) as ex:
                e = "A".rjust(1)  # Astring

            # Residue Number
            f = str(mol["number"]).rjust(4)  # resnum

            # Coordinates
            g = str("%8.3f" % (float(mol["coords"][_i][0]))).rjust(8)  # x
            h = str("%8.3f" % (float(mol["coords"][_i][1]))).rjust(8)  # y
            i = str("%8.3f" % (float(mol["coords"][_i][2]))).rjust(8)  # z\

            # Extra
            j = "".rjust(6)  # str('%6.2f'%(float(j[9]))).rjust(6)#occ
            k = "".rjust(6)  # str('%6.2f'%(float(j[10]))).ljust(6)#temp

            # Element Name
            if atm in symbol:
                l = symbol[atm].rjust(12)  # elname
            else:
                l = atm.rjust(12)  # elname

            if new in ['N', 'CA', 'C', 'O']:
                fpdb.write(
                    "%s%s %s %s %s%s    %s%s%s%s%s%s\n"
                    % (a, b, c, d, e, f, g, h, i, j, k, l)
                )
            # fpdb.write("%s%s %s %s %s%s    %s%s%s\n" % (a,b,c,d,e,f,g,h,i))
            atmnum += 1

    fpdb.write("ENDMDL \n")
    # fpdb.write("TER \n")
    fpdb.close()


def print_heavies(mols, filename, conect=False, encode=True, amino='gmx', header=True):
    fpdb = open(f"{filename}", "w+")
    if header:
        fpdb.write("HEADER :  BMBT MIPkit \n")
        fpdb.write("AUTHOR :  Tj Barrett (tbar@tf.uni-kiel.de, barrett.t@northeastern.edu) \n")
        fpdb.write("TITLE  :  Molecule print from BMBT MIPkit\n")

    if amino == 'gmx':
        resatoms = res_atoms_gmx()
    elif amino == 'pdb':
        resatoms = res_atoms_pdb()
    residues = set(res_list())
    encodefm = encode_fms()
    decodefm = set(decode_fms())

    symbol = atom2PT()

    atmnum = 0

    for _, mol in enumerate(mols):
        n = 0
        offset = atmnum

        for _i, atm in enumerate(mol["atoms"]):
            # print(f'{mol["name"]} : {_i} : {atm}')

            # ATOM or HETATM label
            if mol["name"] in residues:
                a = "ATOM".ljust(6)  # atom#6s
            else:
                a = "HETATM".ljust(6)  # atom#6s

            # Atom Number
            b = str(atmnum).rjust(5)  # aomnum#5d

            # Atom Type
            if mol["name"] in residues:
                if mol["name"] in resatoms and atm != 'H' and n < len(resatoms[mol["name"]]):
                    new = resatoms[mol["name"]][n]
                    n+= 1
                elif mol["name"] in resatoms and atm != 'H' and n == len(resatoms[mol["name"]]):
                    # COOH Terminal Unit, rename Oxygen
                    new = 'OH'
                elif atm in symbol:
                    new = symbol[atm]
                else:
                    new = atm
                    
                c = new.center(4)  # atomname$#4s
            elif atm in symbol:
                new = symbol[atm]
                c = new.center(4)  # atomname$#4s
            else:
                c = atm.center(4)  # atomname$#4s

            # Residue Name
            if mol["name"] in residues:
                d = mol["name"].ljust(3)  # resname#1s
            elif mol["name"] in decodefm:
                d = mol["name"].ljust(3)  # resname#1s
            else:
                if encode:
                    efm = encodefm.get(mol["name"], mol["name"])
                else:
                    efm = mol["name"]
                d = efm.ljust(3)  # resname#1s

            # Chain
            try:
                e = mol['chain'].rjust(1)  # Astring
            except (KeyError, AttributeError) as ex:
                e = "A".rjust(1)  # Astring

            # Residue Number
            f = str(mol["number"]).rjust(4)  # resnum

            # Coordinates
            g = str("%8.3f" % (float(mol["coords"][_i][0]))).rjust(8)  # x
            h = str("%8.3f" % (float(mol["coords"][_i][1]))).rjust(8)  # y
            i = str("%8.3f" % (float(mol["coords"][_i][2]))).rjust(8)  # z\

            # Extra
            j = "".rjust(6)  # str('%6.2f'%(float(j[9]))).rjust(6)#occ
            k = "".rjust(6)  # str('%6.2f'%(float(j[10]))).ljust(6)#temp

            # Element Name
            if atm in symbol:
                l = symbol[atm].rjust(12)  # elname
            else:
                l = atm.rjust(12)  # elname
            
            if atm != 'H':
                fpdb.write(
                    "%s%s %s %s %s%s    %s%s%s%s%s%s\n"
                    % (a, b, c, d, e, f, g, h, i, j, k, l)
                )
            # fpdb.write("%s%s %s %s %s%s    %s%s%s\n" % (a,b,c,d,e,f,g,h,i))
            atmnum += 1

    fpdb.write("ENDMDL \n")
    # fpdb.write("TER \n")
    fpdb.close()