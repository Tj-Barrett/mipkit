
################################################################
#
# Acpype
#
################################################################

def read_itp(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    itp = {
        "atoms": [],
        "pairs": [],
        "bonds": [],
        "angles": [],
        "propers": [],
        "impropers": [],
    }

    atomtypes = []

    btypes = False
    batoms = False
    bpairs = False
    bbonds = False
    bangle = False
    bprope = False
    bimpro = False

    for i, _l in enumerate(lines):
        _sl = _l.split()

        if len(_sl) > 0:
            if _sl[0] == ";":
                continue

            if _sl[0] == "[":
                if _sl[1] == "atomtypes":
                    btypes = True
                elif _sl[1] == "moleculetype":
                    btypes = False
                elif _sl[1] == "atoms":
                    batoms = True
                elif _sl[1] == "bonds":
                    batoms = False
                    bbonds = True
                elif _sl[1] == "pairs":
                    bbonds = False
                    bpairs = True
                elif _sl[1] == "angles":
                    bpairs = False
                    bangle = True
                elif _sl[1] == "dihedrals" and len(_sl) > 4 and _sl[4] == "propers":
                    bangle = False
                    bprope = True
                elif _sl[1] == "dihedrals" and len(_sl) > 4 and _sl[4] == "impropers":
                    bprope = False
                    bimpro = True
                elif _sl[1] == "dihedrals" and len(_sl) <= 4:
                    # no inline "propers"/"impropers" comment; assume propers,
                    # matching a bare "[ dihedrals ]" header
                    bangle = False
                    bprope = True
                continue

            if btypes:
                atomtypes.append(_l)

            if batoms:
                atype = _sl[1]
                atom = _sl[4]
                cgnr = _sl[5]
                charge = _sl[6]
                mass = _sl[7]

                itp["atoms"].append([atype, atom, cgnr, charge, mass])

            if bpairs:
                ai = int(_sl[0])
                aj = int(_sl[1])
                funct = _sl[2]

                itp["pairs"].append([ai, aj, funct])

            if bbonds:
                ai = int(_sl[0])
                aj = int(_sl[1])
                funct = _sl[2]
                r = _sl[3]
                k = _sl[4]

                itp["bonds"].append([ai, aj, funct, r, k])

            if bangle:
                ai = int(_sl[0])
                aj = int(_sl[1])
                ak = int(_sl[2])
                funct = _sl[3]
                theta = _sl[4]
                cth = _sl[5]

                itp["angles"].append([ai, aj, ak, funct, theta, cth])

            if bprope:
                i = int(_sl[0])
                j = int(_sl[1])
                k = int(_sl[2])
                l = int(_sl[3])
                funct = _sl[4]
                phase = _sl[5]
                kd = _sl[6]
                pn = _sl[7]

                itp["propers"].append([i, j, k, l, funct, phase, kd, pn])

            if bimpro:
                i = int(_sl[0])
                j = int(_sl[1])
                k = int(_sl[2])
                l = int(_sl[3])
                funct = _sl[4]
                phase = _sl[5]
                kd = _sl[6]
                pn = _sl[7]

                itp["impropers"].append([i, j, k, l, funct, phase, kd, pn])

    return itp, atomtypes


def read_atomtypes(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    itp = {
        "atoms": [],
        "pairs": [],
        "bonds": [],
        "angles": [],
        "propers": [],
        "impropers": [],
    }

    atomtypes = []

    btypes = False
    batoms = False
    bpairs = False
    bbonds = False
    bangle = False
    bprope = False
    bimpro = False

    for i, _l in enumerate(lines):
        _sl = _l.split()

        if len(_sl) > 0:
            if _sl[0] == ";":
                continue

            if _sl[0] == "[":
                if _sl[1] == "atomtypes":
                    btypes = True
                elif _sl[1] == "moleculetype":
                    btypes = False
                continue

            if btypes:
                atomtypes.append(_l)

    return atomtypes


def read_atoms(filename):
    molecules = []
    # print(filename)
    with open(filename, "r") as f:
        lines = f.readlines()

    atoms = []

    btypes = False
    batoms = False
    bpairs = False
    bbonds = False
    bangle = False
    bprope = False
    bimpro = False

    for i, _l in enumerate(lines):
        _sl = _l.split()

        if len(_sl) > 0:
            if _sl[0] == ";":
                continue

            if _sl[1] == "atoms":
                batoms = True
                continue
            elif _sl[1] == "bonds":
                batoms = False
            elif _sl[1] == "pairs":
                batoms = False
            if batoms:
                atype = _sl[1]
                atom = _sl[4]
                cgnr = _sl[5]
                charge = _sl[6]
                mass = _sl[7]

                atoms.append([atype, atom, cgnr, charge, mass])

    return atoms


def merge_itps(filenames, outname, residues, molname="CLX"):
    if len(residues) < len(filenames):
        raise ValueError(
            f"merge_itps: residues list ({len(residues)} entries) is shorter "
            f"than the number of itp files being merged ({len(filenames)})."
        )
    merged = {
        "atoms": [],
        "pairs": [],
        "bonds": [],
        "angles": [],
        "propers": [],
        "impropers": [],
    }

    atomtypes = []

    buffer = 0
    resname = []
    resnumb = []
    rn = 0
    for _i, filename in enumerate(filenames):
        rn += 1
        itp, atmtyp = read_itp(filename)

        atomtypes.extend(atmtyp)

        merged["atoms"].extend(itp["atoms"])

        for _j in range(len(itp["atoms"])):
            resname.append(residues[_i])
            resnumb.append(rn)

        for pair in itp["pairs"]:
            ai = pair[0] + buffer
            aj = pair[1] + buffer
            funct = pair[2]

            merged["pairs"].append([ai, aj, funct])

        for bond in itp["bonds"]:
            ai = bond[0] + buffer
            aj = bond[1] + buffer
            funct = bond[2]
            r = bond[3]
            k = bond[4]

            merged["bonds"].append([ai, aj, funct, r, k])

        for angle in itp["angles"]:
            ai = angle[0] + buffer
            aj = angle[1] + buffer
            ak = angle[2] + buffer
            funct = angle[3]
            theta = angle[4]
            cth = angle[5]

            merged["angles"].append([ai, aj, ak, funct, theta, cth])

        for proper in itp["propers"]:
            i = proper[0] + buffer
            j = proper[1] + buffer
            k = proper[2] + buffer
            l = proper[3] + buffer
            funct = proper[4]
            phase = proper[5]
            kd = proper[6]
            pn = proper[7]

            merged["propers"].append([i, j, k, l, funct, phase, kd, pn])

        for proper in itp["impropers"]:
            i = proper[0] + buffer
            j = proper[1] + buffer
            k = proper[2] + buffer
            l = proper[3] + buffer
            funct = proper[4]
            phase = proper[5]
            kd = proper[6]
            pn = proper[7]

            merged["impropers"].append([i, j, k, l, funct, phase, kd, pn])

        buffer = len(merged["atoms"])

    file = open(outname, "w+")

    atmset = list(set(atomtypes))
    file.write(
        "[ atomtypes ]\n;name   bond_type     mass     charge   ptype   sigma         epsilon       Amb\n"
    )
    for atm in atmset:
        if atm[0] != ";":
            file.write(atm)

    file.write("\n[ moleculetype ]\n; name  nrexcl\n")
    # resset = list(set(residues))
    resset = [molname]
    for res in resset:
        a = res.rjust(5)
        b = "3".rjust(8)
        file.write("%s%s\n" % (a, b))

    file.write(
        "\n[ atoms ]\n;   nr  type  resi  res  atom  cgnr     charge      mass       \n"
    )
    for i, atom in enumerate(merged["atoms"]):
        a = str(i + 1).rjust(6)
        b = atom[0].rjust(5)
        c = str(resnumb[i]).rjust(6)
        d = resname[i].rjust(6)
        e = atom[1][0].rjust(6)
        f = str(atom[2]).rjust(5)
        g = str(atom[3]).rjust(13)
        h = str(atom[4]).rjust(13)
        file.write("%s%s%s%s%s%s%s%s\n" % (a, b, c, d, e, f, g, h))

    file.write("\n[ pairs ]\n;   ai     aj    funct\n")
    for i, pair in enumerate(merged["pairs"]):
        a = str(pair[0]).rjust(6)
        b = str(pair[1]).rjust(7)
        c = pair[2].rjust(7)
        file.write("%s%s%s\n" % (a, b, c))

    file.write("\n[ bonds ]\n;   ai     aj funct   r             k\n")
    for i, bond in enumerate(merged["bonds"]):
        a = str(bond[0]).rjust(6)
        b = str(bond[1]).rjust(7)
        c = bond[2].rjust(4)
        d = bond[3].rjust(14)
        e = bond[4].rjust(14)
        file.write("%s%s%s%s%s\n" % (a, b, c, d, e))

    file.write("\n[ angles ]\n;   ai     aj     ak    funct   theta         cth\n")
    for i, angle in enumerate(merged["angles"]):
        a = str(angle[0]).rjust(6)
        b = str(angle[1]).rjust(6)
        c = str(angle[2]).rjust(7)
        d = angle[3].rjust(7)
        e = angle[4].rjust(14)
        f = angle[5].rjust(14)
        file.write("%s%s%s%s%s%s\n" % (a, b, c, d, e, f))

    file.write(
        "\n[ dihedrals ] ; propers\n;    i      j      k      l   func   phase     kd      pn\n"
    )
    for i, dihedral in enumerate(merged["propers"]):
        a = str(dihedral[0]).rjust(6)
        b = str(dihedral[1]).rjust(7)
        c = str(dihedral[2]).rjust(7)
        d = str(dihedral[3]).rjust(7)
        e = dihedral[4].rjust(7)
        f = dihedral[5].rjust(9)
        g = dihedral[6].rjust(10)
        h = dihedral[7].rjust(4)
        file.write("%s%s%s%s%s%s%s%s\n" % (a, b, c, d, e, f, g, h))

    file.write(
        "\n[ dihedrals ] ; impropers\n;    i      j      k      l   func   phase     kd      pn\n"
    )
    for i, dihedral in enumerate(merged["impropers"]):
        a = str(dihedral[0]).rjust(6)
        b = str(dihedral[1]).rjust(7)
        c = str(dihedral[2]).rjust(7)
        d = str(dihedral[3]).rjust(7)
        e = dihedral[4].rjust(7)
        f = dihedral[5].rjust(9)
        g = dihedral[6].rjust(10)
        h = dihedral[7].rjust(4)
        file.write("%s%s%s%s%s%s%s%s\n" % (a, b, c, d, e, f, g, h))

    file.close()


def merge_itps_no_acpype(filenames, outname, residues):
    merged = {
        "atoms": [],
        "pairs": [],
        "bonds": [],
        "angles": [],
        "propers": [],
        "impropers": [],
    }

    atomtypes = []

    buffer = 0
    resname = []
    resnumb = []
    rn = 0
    for _i, filename in enumerate(filenames):
        rn += 1
        itp, atmtyp = read_itp(filename)

        atomtypes.extend(atmtyp)

        merged["atoms"].extend(itp["atoms"])

        for _j in range(len(itp["atoms"])):
            resname.append(residues[_i])
            resnumb.append(rn)

        for pair in itp["pairs"]:
            ai = pair[0] + buffer
            aj = pair[1] + buffer
            funct = pair[2]

            merged["pairs"].append([ai, aj, funct])

        for bond in itp["bonds"]:
            ai = bond[0] + buffer
            aj = bond[1] + buffer
            funct = bond[2]

            merged["bonds"].append([ai, aj, funct])

        for angle in itp["angles"]:
            ai = angle[0] + buffer
            aj = angle[1] + buffer
            ak = angle[2] + buffer
            funct = angle[3]

            merged["angles"].append([ai, aj, ak, funct])

        for proper in itp["propers"]:
            i = proper[0] + buffer
            j = proper[1] + buffer
            k = proper[2] + buffer
            l = proper[3] + buffer
            funct = proper[4]

            merged["propers"].append([i, j, k, l, funct])

        for proper in itp["impropers"]:
            i = proper[0] + buffer
            j = proper[1] + buffer
            k = proper[2] + buffer
            l = proper[3] + buffer
            funct = proper[4]

            merged["impropers"].append([i, j, k, l, funct])

        buffer = len(merged["atoms"])

    file = open(outname, "w+")

    file.write("\n[ moleculetype ]\n; name  nrexcl\n")
    resset = ["CLX"]
    for res in resset:
        a = res.rjust(5)
        b = "3".rjust(8)
        file.write("%s%s\n" % (a, b))

    file.write(
        "\n[ atoms ]\n;   nr  type  resi  res  atom  cgnr     charge      mass       \n"
    )
    for i, atom in enumerate(merged["atoms"]):
        a = str(i + 1).rjust(6)
        b = atom[0].rjust(5)
        c = str(resnumb[i]).rjust(6)
        d = resname[i].rjust(6)
        e = atom[1][0].rjust(6)
        f = str(atom[2]).rjust(5)
        g = str(atom[3]).rjust(13)
        h = str(atom[4]).rjust(13)
        file.write("%s%s%s%s%s%s%s%s\n" % (a, b, c, d, e, f, g, h))

    file.write("\n[ pairs ]\n;   ai     aj    funct\n")
    for i, pair in enumerate(merged["pairs"]):
        a = str(pair[0]).rjust(6)
        b = str(pair[1]).rjust(7)
        c = pair[2].rjust(7)
        file.write("%s%s%s\n" % (a, b, c))

    file.write("\n[ bonds ]\n;   ai     aj funct   r             k\n")
    for i, bond in enumerate(merged["bonds"]):
        a = str(bond[0]).rjust(6)
        b = str(bond[1]).rjust(7)
        c = bond[2].rjust(4)
        file.write("%s%s%s\n" % (a, b, c))

    file.write("\n[ angles ]\n;   ai     aj     ak    funct   theta         cth\n")
    for i, angle in enumerate(merged["angles"]):
        a = str(angle[0]).rjust(6)
        b = str(angle[1]).rjust(6)
        c = str(angle[2]).rjust(7)
        d = angle[3].rjust(7)
        file.write("%s%s%s%s\n" % (a, b, c, d))

    file.write(
        "\n[ dihedrals ] ; propers\n;    i      j      k      l   func   phase     kd      pn\n"
    )
    for i, dihedral in enumerate(merged["propers"]):
        a = str(dihedral[0]).rjust(6)
        b = str(dihedral[1]).rjust(7)
        c = str(dihedral[2]).rjust(7)
        d = str(dihedral[3]).rjust(7)
        e = dihedral[4].rjust(7)
        file.write("%s%s%s%s%s\n" % (a, b, c, d, e))

    file.write(
        "\n[ dihedrals ] ; impropers\n;    i      j      k      l   func   phase     kd      pn\n"
    )
    for i, dihedral in enumerate(merged["impropers"]):
        a = str(dihedral[0]).rjust(6)
        b = str(dihedral[1]).rjust(7)
        c = str(dihedral[2]).rjust(7)
        d = str(dihedral[3]).rjust(7)
        e = dihedral[4].rjust(7)
        file.write("%s%s%s%s%s\n" % (a, b, c, d, e))

    file.close()


def merge_atomtypes(itpfrom, itpinto, cleaned):
    atomtypes = []

    for _i, filename in enumerate([itpfrom, itpinto]):
        atmtyp = read_atomtypes(filename)
        atomtypes.extend(atmtyp)

    # ----------------------------------------------------
    # Merge Files
    # ----------------------------------------------------
    with open(itpinto, "r") as fin:
        itplines = fin.readlines()

    file = open(itpinto, "w+")

    atmset = list(set(atomtypes))
    file.write(
        "[ atomtypes ]\n;name   bond_type     mass     charge   ptype   sigma         epsilon       Amb\n"
    )
    for atm in atmset:
        if atm[0] != ";":
            file.write(atm)

    file.write("\n")

    btype = True
    for line in itplines:
        sl = line.split()

        if len(sl) > 2 and sl[0] == "[" and sl[1] == "moleculetype":
            btype = False

        if btype:
            continue
        else:
            file.write(line)

    file.close()

    del itplines
    # ----------------------------------------------------
    # Merge Files
    # ----------------------------------------------------
    with open(itpfrom, "r") as fin:
        itplines = fin.readlines()

    file = open(cleaned, "w+")
    btype = True
    for line in itplines:
        sl = line.split()

        if len(sl) > 2 and sl[0] == "[" and sl[1] == "moleculetype":
            btype = False

        if btype:
            continue
        else:
            file.write(line)

    file.close()


def replace_residues(itpfrom, itpinto, residues, molname):
    atoms = read_atoms(itpfrom)

    if len(residues) < len(atoms):
        raise ValueError(
            f"replace_residues: residues list ({len(residues)} entries) is shorter "
            f"than the atom list parsed from {itpfrom} ({len(atoms)} atoms)."
        )

    # ----------------------------------------------------
    # Merge Files
    # ----------------------------------------------------
    with open(itpfrom, "r") as fin:
        itplines = fin.readlines()

    file = open(itpinto, "w+")

    btype = False
    moltype = False
    atom_write = True
    mol_write = True
    for line in itplines:
        sl = line.split()

        if len(sl) > 2 and sl[0] == "[" and sl[1] == "moleculetype":
            moltype = True

        if moltype and mol_write:
            file.write("\n[ moleculetype ]\n; name  nrexcl\n")
            # resset = list(set(residues))
            resset = [molname]
            for res in resset:
                a = res.rjust(5)
                b = "3".rjust(8)
                file.write("%s%s\n" % (a, b))

            mol_write = False

        if len(sl) > 2 and sl[0] == "[" and sl[1] == "atoms":
            btype = True
            moltype = False

        if len(sl) > 2 and sl[0] == "[" and sl[1] == "bonds":
            btype = False

        if btype and atom_write:
            resnumb = 1
            file.write(
                "\n[ atoms ]\n;   nr  type  resi  res  atom  cgnr     charge      mass       \n"
            )
            for i, atom in enumerate(atoms):
                a = str(i + 1).rjust(6)
                b = atom[0].rjust(5)

                if i > 0 and residues[i] != residues[i - 1]:
                    resnumb += 1
                c = str(resnumb).rjust(6)
                d = residues[i].rjust(6)
                e = atom[1][0].rjust(6)
                f = str(atom[2]).rjust(5)
                g = str(atom[3]).rjust(13)
                h = str(atom[4]).rjust(13)
                file.write("%s%s%s%s%s%s%s%s\n" % (a, b, c, d, e, f, g, h))

            atom_write = False
            file.write("\n")

        elif not btype and not moltype:
            file.write(line)

    file.close()


################################################################
#
# Gromacs and or acpype
#
################################################################
def read_top(filename):
    molecules = []
    # print(filename)
    with open(filename, "r") as f:
        lines = f.readlines()

    itp = {
        "atoms": [],
        "pairs": [],
        "bonds": [],
        "angles": [],
        "propers": [],
        "impropers": [],
    }

    atomtypes = []

    btypes = False
    batoms = False
    bpairs = False
    bbonds = False
    bangle = False
    bprope = False
    bimpro = False

    for i, _l in enumerate(lines):
        _sl = _l.split()

        if len(_sl) > 0:
            if _sl[0] == ";":
                continue

            if _sl[0] == "#ifdef":
                bprope = False
                continue

            if _sl[0] == "[":
                if _sl[1] == "atomtypes":
                    btypes = True
                elif _sl[1] == "moleculetype":
                    btypes = False
                elif _sl[1] == "atoms":
                    batoms = True
                elif _sl[1] == "bonds":
                    batoms = False
                    bbonds = True
                elif _sl[1] == "pairs":
                    bbonds = False
                    bpairs = True
                elif _sl[1] == "angles":
                    bpairs = False
                    bangle = True
                elif _sl[1] == "dihedrals":
                    bangle = False
                    bprope = True
                continue

            if btypes:
                atomtypes.append(_l)

            if batoms:
                atype = _sl[1]
                atom = _sl[4]
                cgnr = _sl[5]
                charge = _sl[6]
                mass = _sl[7]

                itp["atoms"].append([atype, atom, cgnr, charge, mass])

            if bpairs:
                ai = int(_sl[0])
                aj = int(_sl[1])
                funct = _sl[2]

                itp["pairs"].append([ai, aj, funct])

            if bbonds:
                ai = int(_sl[0])
                aj = int(_sl[1])
                funct = _sl[2]

                itp["bonds"].append([ai, aj, funct])

            if bangle:
                ai = int(_sl[0])
                aj = int(_sl[1])
                ak = int(_sl[2])
                funct = _sl[3]

                itp["angles"].append([ai, aj, ak, funct])

            if bprope:
                i = int(_sl[0])
                j = int(_sl[1])
                k = int(_sl[2])
                l = int(_sl[3])
                funct = _sl[4]

                itp["propers"].append([i, j, k, l, funct])

    return itp, atomtypes


def convert_itp(filenames, outname, residues):
    merged = {
        "atoms": [],
        "pairs": [],
        "bonds": [],
        "angles": [],
        "propers": [],
        "impropers": [],
    }

    atomtypes = []

    buffer = 0
    resname = []
    resnumb = []
    rn = 0
    for _i, filename in enumerate(filenames):
        rn += 1
        itp, atmtyp = read_top(filename)

        atomtypes.extend(atmtyp)

        merged["atoms"].extend(itp["atoms"])

        for _j in range(len(itp["atoms"])):
            resname.append(residues[_i])
            resnumb.append(rn)

        for pair in itp["pairs"]:
            ai = pair[0] + buffer
            aj = pair[1] + buffer
            funct = pair[2]

            merged["pairs"].append([ai, aj, funct])

        for bond in itp["bonds"]:
            ai = bond[0] + buffer
            aj = bond[1] + buffer
            funct = bond[2]

            merged["bonds"].append([ai, aj, funct])

        for angle in itp["angles"]:
            ai = angle[0] + buffer
            aj = angle[1] + buffer
            ak = angle[2] + buffer
            funct = angle[3]

            merged["angles"].append([ai, aj, ak, funct])

        for proper in itp["propers"]:
            i = proper[0] + buffer
            j = proper[1] + buffer
            k = proper[2] + buffer
            l = proper[3] + buffer
            funct = proper[4]

            merged["propers"].append([i, j, k, l, funct])

        for proper in itp["impropers"]:
            i = proper[0] + buffer
            j = proper[1] + buffer
            k = proper[2] + buffer
            l = proper[3] + buffer
            funct = proper[4]

            merged["impropers"].append([i, j, k, l, funct])

        buffer = len(merged["atoms"])

    file = open(outname, "w+")

    atmset = list(set(atomtypes))
    file.write(
        "[ atomtypes ]\n;name   bond_type     mass     charge   ptype   sigma         epsilon       Amb\n"
    )
    for atm in atmset:
        if atm[0] != ";":
            file.write(atm)

    file.write("\n[ moleculetype ]\n; name  nrexcl\n")
    # resset = list(set(residues))
    resset = ["Protein"]
    for res in resset:
        a = res.rjust(5)
        b = "3".rjust(8)
        file.write("%s%s\n" % (a, b))

    file.write(
        "\n[ atoms ]\n;   nr  type  resi  res  atom  cgnr     charge      mass       \n"
    )
    for i, atom in enumerate(merged["atoms"]):
        a = str(i + 1).rjust(6)
        b = atom[0].rjust(5)
        c = str(resnumb[i]).rjust(6)
        d = resname[i].rjust(6)
        e = atom[1][0].rjust(6)
        f = str(atom[2]).rjust(5)
        g = str(atom[3]).rjust(13)
        h = str(atom[4]).rjust(13)
        file.write("%s%s%s%s%s%s%s%s\n" % (a, b, c, d, e, f, g, h))

    file.write("\n[ pairs ]\n;   ai     aj    funct\n")
    for i, pair in enumerate(merged["pairs"]):
        a = str(pair[0]).rjust(6)
        b = str(pair[1]).rjust(7)
        c = pair[2].rjust(7)
        file.write("%s%s%s\n" % (a, b, c))

    file.write("\n[ bonds ]\n;   ai     aj funct   r             k\n")
    for i, bond in enumerate(merged["bonds"]):
        a = str(bond[0]).rjust(6)
        b = str(bond[1]).rjust(7)
        c = bond[2].rjust(4)
        file.write("%s%s%s\n" % (a, b, c))

    file.write("\n[ angles ]\n;   ai     aj     ak    funct   theta         cth\n")
    for i, angle in enumerate(merged["angles"]):
        a = str(angle[0]).rjust(6)
        b = str(angle[1]).rjust(6)
        c = str(angle[2]).rjust(7)
        d = angle[3].rjust(7)
        file.write("%s%s%s%s\n" % (a, b, c, d))

    file.write(
        "\n[ dihedrals ] ; propers\n;    i      j      k      l   func   phase     kd      pn\n"
    )
    for i, dihedral in enumerate(merged["propers"]):
        a = str(dihedral[0]).rjust(6)
        b = str(dihedral[1]).rjust(7)
        c = str(dihedral[2]).rjust(7)
        d = str(dihedral[3]).rjust(7)
        e = dihedral[4].rjust(7)
        file.write("%s%s%s%s%s\n" % (a, b, c, d, e))

    file.close()
