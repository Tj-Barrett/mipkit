import os, re
from MIPkit.constants.constants import decode_fms, exclusions
from MIPkit.initiators.initiators_yaml import init_attackers

def isfloat(n):
    try:
        float(n)
        return True
    except ValueError:
        return False


def read_pdb(filename, conect=False, decode=False, hydrogens=True, initiators=False):
    molecules = []

    with open(filename, "r") as f:
        lines = f.readlines()

    decodefm = decode_fms()

    molecule = {"name": "", "number": "", "atoms": [], "coords": [], "bonds": {}}

    new = True
    curres = ""
    curnum = ""
    atmnum = 0
    molnum = 0

    step_pts = []

    #PDB files have standard layouts
    #  6s5d 4s 3s 1s4d    8.3f8.3f8.3f6s6s12s
    # "%s%s %s %s %s%s    %s%s%s%s%s%s
    # PDB_ints = [(0, 6), (6, 11),
    #             (12, 16), (17, 20),
    #             (21, 22), (22, 26),
    #             (30, 38), (38, 46),
    #             (46, 54), (54, 60),
    #             (60, 66), (66, 78)]

    for line_idx, line in enumerate(lines):
        _sl = line.split()

        HETATM = line[0:6].strip()
        if HETATM == "ATOM" or HETATM == "HETATM":

            # _ml = []
            # for (a,b) in PDB_ints:
            #     _ml.append(line[a:b])

            res = line[17:20].strip()         #_ml[3].strip()
            atm = line[12:16].strip()         #_ml[2].strip()
            anum = line[6:11].strip()         #_ml[1].strip()
            molnum = line[22:26].strip()      #_ml[5].strip()
            x = float(line[30:38])            #float(_ml[6])
            y = float(line[38:46])            #float(_ml[7])
            z = float(line[46:54])            #float(_ml[8])
            forces = False

            # first instance
            if curres != res and new:
                curres = res
                curnum = molnum
                new = False

                # part of the same residue
            if curres == res and curnum == molnum:
                if molecule["name"] == "":
                    if res in decodefm and decode:
                        molecule["name"] = decodefm[res]
                    else:
                        molecule["name"] = res
                        molecule["number"] = molnum
                # skip if it is a hydrogen and strip any numbers if hetatm
                if HETATM == "HETATM":
                    atm = re.sub(r'\d+', '', atm)
                if atm == "H" and hydrogens == False:
                    continue
                molecule["atoms"].append(atm)
                step_pts.append((x, y, z))

            elif curres == res and curnum != molnum:
                molecule["coords"] = step_pts

                molecules.append(molecule)
                molecule = {
                    "name": "",
                    "number": "",
                    "atoms": [],
                    "coords": [],
                    "bonds": {},
                }

                curres = res
                curnum = molnum
                if molecule["name"] == "":
                    if res in decodefm and decode:
                        molecule["name"] = decodefm[res]
                    else:
                        molecule["name"] = res
                        molecule["number"] = molnum

                # skip if it is a hydrogen and strip any numbers if hetatm
                if HETATM == "HETATM":
                    atm = re.sub(r'\d+', '', atm)
                if atm == "H" and hydrogens == False:
                    continue
                molecule["atoms"].append(atm)

                step_pts = []
                step_pts.append((x, y, z))

            elif curres != res:
                molecule["coords"] = step_pts

                molecules.append(molecule)
                molecule = {
                    "name": "",
                    "number": "",
                    "atoms": [],
                    "coords": [],
                    "bonds": {},
                }

                curres = res
                curnum = molnum
                if molecule["name"] == "":
                    if res in decodefm and decode:
                        molecule["name"] = decodefm[res]
                    else:
                        molecule["name"] = res
                        molecule["number"] = molnum
                # skip if it is a hydrogen and strip any numbers if hetatm
                if HETATM == "HETATM":
                    atm = re.sub(r'\d+', '', atm)
                if atm == "H" and hydrogens == False:
                    continue
                molecule["atoms"].append(atm)

                step_pts = []
                step_pts.append((x, y, z))

        elif _sl and _sl[0] == "CONECT" and conect:
            molecule['bonds'][int(_sl[1])] = []
            bonded_atms = _sl[2:]
            for ba in bonded_atms:
                molecule['bonds'][int(_sl[1])].append(int(ba))

        elif len(_sl) == 3 or (_sl and _sl[0] in ['TER', 'END', 'ENDMDL']) or line_idx == len(lines) - 1 :
            if len(step_pts) > 0:
                molecule["coords"] = step_pts

                molecules.append(molecule)
            
            molecule = {
                    "name": "",
                    "number": "",
                    "atoms": [],
                    "coords": [],
                    "bonds": {},
                }

            step_pts = []

            continue
    
    _molecules = molecules
    molecules = []
    inits = exclusions()
    for i, mol in enumerate(_molecules):
        if mol["name"] not in inits or initiators:
            molecule = {
                "name": mol["name"],
                "number": i + 1,
                "atoms": mol["atoms"],
                "coords": mol["coords"],
                "bonds": mol["bonds"],
            }
            molecules.append(molecule)

    return molecules


def clean_pdb(filename):
    exclude = ["MASTER", "ENDMDL", "MODEL", "TITLE", "COMPND", "AUTHOR", "HEADER", "REMARK", "CRYST1"]

    with open(filename, "r") as f:
        lines = f.readlines()

    with open(filename, "w") as f:
        for line_idx, line in enumerate(lines):
            _sl = line.split()
            if _sl and _sl[0] in exclude:
                continue
            else:
                f.write(line)
