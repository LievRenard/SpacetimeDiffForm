import numpy as np

def delta(i,j):
    if i == j: return 1
    else: return 0

def E4(i, j, k, l):
    return np.linalg.det(np.array([
        [delta(0,i), delta(0,j), delta(0,k), delta(0,l)],
        [delta(1,i), delta(1,j), delta(1,k), delta(1,l)],
        [delta(2,i), delta(2,j), delta(2,k), delta(2,l)],
        [delta(3,i), delta(3,j), delta(3,k), delta(3,l)]
    ]))

def how_permut(v: list) -> int:
    if len(v) == 4:
        return int(E4(v[0], v[1], v[2], v[3]))
    if len(v) == 3:
        vsor = np.sort(v)
        return int(np.linalg.det(np.array([
            [delta(vsor[0],v[0]), delta(vsor[0],v[1]), delta(vsor[0],v[2])],
            [delta(vsor[1],v[0]), delta(vsor[1],v[1]), delta(vsor[1],v[2])],
            [delta(vsor[2],v[0]), delta(vsor[2],v[1]), delta(vsor[2],v[2])]
        ])))
    if len(v) == 2:
        if v[0] < v[1]: return 1
        else: return -1
    if len(v) <= 1:
        return 1

'''def Hodge(d):
    signcheck = np.array(d).tolist()
    res = []; eta = 1
    for i in range(4):
        if not i in d: signcheck.append(i), res.append(i)
        elif i != 0: eta *= -1
    return (eta*E4(signcheck[0], signcheck[1], signcheck[2], signcheck[3]), res)
'''

class KForm:
    def __init__(self, sign: int, val: str, diff: list):
        self.sign = np.sign(sign)
        self.val = val
        self.diff = diff

    def sort(self) -> None:
        self.sign = self.sign * how_permut(self.diff)
        self.val = self.val
        self.diff.sort

    def get_sorted(self):
        return KForm(self.sign * how_permut(self.diff), self.val, np.sort(self.diff).tolist())

    def __repr__(self) -> str:
        text = ""
        diff_sym = ["dt", "dx", "dy", "dz"]
        if self.sign == -1: text = text + "-"
        text = text + self.val + " "
        for i in self.diff:
            text = text + diff_sym[i] + "∧"
        return text[:-1]

    def Hodge(self):
        signcheck = np.array(self.diff).tolist()
        res = []; eta = 1
        for i in range(4):
            if not i in self.diff : signcheck.append(i), res.append(i)
            elif i != 0: eta *= -1
        return KForm(self.sign * int(eta*E4(signcheck[0], signcheck[1], signcheck[2], signcheck[3])), self.val, res).get_sorted()

    def __invert__(self):
        return self.Hodge()

    def wedge(self, other):
        for i in self.diff:
            if i in other.diff: return KForm(1, "0", [])
        return KForm(self.sign*other.sign, self.val+"*"+other.val, self.diff+other.diff).get_sorted()
    
    def __xor__(self, other):
        return self.wedge(other)
    
    def __rxor__(self, other):
        return self.wedge(other)
        
Ex = KForm(1, "Ex", [0,1])
epx = KForm(1, "epx", [0,1])
print(~epx^Ex)