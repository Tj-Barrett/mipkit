import os, glob

def strip_string(filename, pre, post):
    f1 = filename.replace(pre,"")
    f2 = f1.replace(post,"")
    return int(f2)

def wildcard(filename):
    matches = []

    if "*" not in filename:
        print('No wildcard ("*") found in filename. Use quotes around the string and include a "*".')
        exit()

    for infile in glob.glob(filename):
        matches.append(infile)

    if len(matches) == 0:
        print('No wildcard files found. Check wildcard name and use quotes for the string.')
        exit()
    else:
        basenames = []
        filepath = ""
        for fname in matches:
            basenames.append(os.path.basename(fname))
            if filepath == "":
                filepath = os.path.split(fname)[0]

        basefilename = os.path.basename(filename)
        wildcard_loc = basefilename.find("*")

        wpost = basefilename[wildcard_loc+1:]
        wpre = basefilename[:wildcard_loc]

        basenames.sort(key=lambda x: strip_string(x,wpre,wpost))

        ordered_files = []
        for base in basenames:
            ordered_files.append(filepath+"/"+base)

    return ordered_files