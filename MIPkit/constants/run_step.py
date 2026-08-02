from math import floor

def alphanumeric():
    return '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'

def run_step(i):
    a = alphanumeric()
    b = a
    if i > len(a):
        j = floor(i/len(a))
        k = i - j*len(a)
    else:
        j = 0
        k = i

    return 'X'+a[k]+b[j]