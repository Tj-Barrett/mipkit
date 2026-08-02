

def read_ndx(filename):
    molecules = []
    # print(filename)
    with open(filename, "r") as f:
        lines = f.readlines()

    groups = []
    names = []

    for i, _l in enumerate(lines):
        _sl = _l.split()
        if len(_sl) > 0 and _sl[0] == "[":
            names.append(_sl[1])
            continue
        elif len(_sl) > 0:
            if _sl[0] == ";":
                continue
            indx = []
            for _s in _sl:
                try:
                    indx.append(int(_s))
                except ValueError:
                    # a comment or other non-index token on this line
                    break

            groups.append(indx)

    return groups, names

def clean_ndx(filename, outname, include=False, exclude=False):
    with open(filename, "r") as f:
        lines = f.readlines()

    new_include = []
    if include:
        for i in include:
            new_include.append(i.upper())
            new_include.append(i)

    new_exclude = []
    if exclude:
        for i in exclude:
            new_exclude.append(i.upper())
            new_exclude.append(i)

    clean = False
    with open(outname, "w+") as f:
        for i, _l in enumerate(lines):
            _sl = _l.split()
            if len(_sl) > 1:
                if include:
                    if _sl[1] in new_include:
                        clean = True
                    elif _sl[0] == "[":
                        clean = False
                if exclude:
                    if _sl[1] in new_exclude:
                        clean = False
                    elif _sl[0] == "[":
                        clean = True

                if clean == True:
                    f.write(_l)

def merge_ndx_groups(filename, outname, mergename="FMs"):
    with open(filename, "r") as f:
        lines = f.readlines()

    big_group = []
    for i, _l in enumerate(lines):
        _sl = _l.split()
        if len(_sl) > 0 and _sl[0] == "[":
            continue
        elif len(_sl) > 0:
            if _sl[0] == ";":
                continue
            for _s in _sl:
                try:
                    big_group.append(int(_s))
                except ValueError:
                    # a comment or other non-index token on this line
                    break

    with open(outname, "w+") as f:
        f.write(f'[ {mergename} ]\n')
        n = 0
        for idx in big_group:
            if n < 15:
                f.write(f'{idx} ')
                n += 1
            else:
                n = 0
                f.write(f'{idx} \n')
        f.write('\n')