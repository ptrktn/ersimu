#!/usr/bin/env python3

"""
ERSimu: A software for exploring chemical kinetics and dynamics in elementary 
reaction systems.
"""

__version__ = "0.0.2"

import sys
import string
import argparse

# FIXME SymPy will be needed for deriving Jacobian (implementation pending)
# LaTeX output requires sympy
# For some reason (?) 'sympy' needs to be imported before 're'

config = {}
config["has_sympy"] = True

try:
    from sympy import *
except ImportError:
    config["has_sympy"] = False

import re
import os
import getopt

config["verbose"] = False
config["octave"] = False
config["latex"] = False
config["scipy"] = False
config["run"] = False
config["name"] = "simulation"


def dbg(msg):
    if config["verbose"]:
        sys.stderr.write(f"DEBUG: {msg}\n")


def err(msg, estatus=1):
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(estatus)


def get_arg_parser():
    parser = argparse.ArgumentParser(
        prog=os.path.basename(__file__), description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--latex", action="store_true", dest="latex", default=False, help="generate LaTeX output")
    parser.add_argument("--lsodac", action="store_true", dest="lsodac", default=True, help="generate LSODA C output (default)")
    parser.add_argument("--name", type=str, dest="name", default="simulation", help="output name")
    parser.add_argument("--octave", action="store_true", dest="octave", default=False, help="generate GNU Octave output")
    parser.add_argument("--run", action="store_true", dest="run", default=False, help="run the simulation")
    parser.add_argument("--scipy", action="store_true", dest="scipy", default=False, help="generate SciPy output")
    parser.add_argument("--verbose", action="store_true", dest="verbose", default=False, help="be verbose")
    parser.add_argument("inputfile", metavar="FILE", type=str, help="input file name")

    return parser


class R:
    def __init__(self):
        self.text = None
        self.i = 0
        self.reactants = None
        self.products = None
        self.frate = None
        self.rrate = None
        # FIXME variable rates need to be renamed
        self.rates = [0, 0]
        self.k = ["kf", "kr"]
        self.s = {}
        self.kf = 1
        self.kr = 1

    def reaction(self, r, n, ers):
        self.i = n

        if " <==> " in r:
            self.rates = [1, 1]
            cs = " <==> "
        elif " ==> " in r:
            self.rates = [1, 0]
            cs = " ==> "
        else:
            raise Exception("Type not found: \" <==> \" or \" ==> \" (%d)" % n)

        self.text = " ".join(" ".join(r.split()).split(" "))

        if "@" in r:
            [rct, rts] = map(str.strip, r.split("@"))
            if 1 == len(rts.split()) and " ==> " == cs:
                self.frate = rts.split()[0]
            elif 2 == len(rts.split()) and " <==> " == cs:
                [self.frate, self.rrate] = rts.split()
            else:
                raise Exception("Malformed rate input: \"%s\" (%d)" % (r, n))
        elif "|" in r:
            [rct, rts] = map(str.strip, r.split("|"))
            if 1 == len(rts.split()) and " ==> " == cs:
                self.kf = rts.split()[0]
            elif 2 == len(rts.split()) and " <==> " == cs:
                [self.kf, self.kr] = rts.split()
            else:
                raise Exception("Malformed constant rate input: \"%s\" (%d)"
                                % (r, n))
        else:
            rct = r.strip()

        [self.reactants, self.products] = rct.split(cs)

        for c in rct.replace(cs, "+").replace("*", "+").split("+"):
            c = c.strip()
            for d in c.split():
                if not(d in ers.constants):
                    self.s[d] = 1 + self.s.get(d, 0)

    def kinet(self, d=True):
        if d:
            a = self.reactants
            b = self.frate
            j = 0
        else:
            a = self.products
            b = self.rrate
            j = 1
            
        v = False
        
        if self.rates[j] > 0:
            if b:
                v = " * ".join(b.split("*"))
            else:
                v = self.k[j] + str(self.i) + " * " + " * ".join(a.split("+"))

        return v
    
    def species(self):
        return self.s

    def smc(self, d, y):
        c = 0
        
        if d:
            r = self.reactants
        else:
            r = self.products
        
        for i in map(str.strip, r.split("+")):
            if y == i:
                c += 1
                
        return c


class ERSimu:
    def __init__(self):
        self.reactions = []
        self.x = []
        self.excess = []
        self.xdot = []
        self.xdot_raw = []
        self.system_species = [] # FIXME dict?
        self.initial = {}
        self.constants = {}
        self.simulation = {}
        self.simulation["T_POINTS"] = 10
        self.latex = {}
        self.bibitem = []
        self.title = "Simulation"
        self.author = None
        self.keywords = []
        self.kf = {}
        self.kr = {}
        self.kinet = {}
        self.kinet_keys = []
        self.lsode_atol = "1E-6"
        self.lsode_rtol = "1E-6"
        self.n_reactions = None
        self.n_equations = None


    def load(self, fname):
        fp = open(fname)
        n = 0
        nr = 1
        r = R()
        for line in fp:
            n += 1
            line = line.strip()

            if config["verbose"]:
                dbg(f"LINE {n:05d} '{line}'")

            if len(line) > 0 and '#' != line[0]:
                if re.search(r'^EXCESS\s', line):
                    if len(line.split()) > 1:
                        for a in line.split()[1:]:
                            if not a in self.excess:
                                self.excess.append(a)
                    continue
                elif re.search(r'^INITIAL\s', line):
                    if len(line.split()) > 2:
                        a, val = line.split()[1:3]
                        self.initial[a] = val
                    continue
                elif re.search(r'^CONSTANT\s', line):
                    if len(line.split()) > 2:
                        a, val = line.split()[1:3]
                        self.constants[a] = val
                    continue
                elif re.search(r'^SIMULATION\s', line):
                    vals = line.split()
                    if len(vals) > 2:
                        if "T_END" == vals[1]:
                            self.simulation["T_END"] = vals[2]
                        elif "T_POINTS" == vals[1]:
                            self.simulation["T_POINTS"] = vals[2]
                        elif "MAXIMUM_STEP_SIZE" == vals[1]:
                            self.simulation["MAXIMUM_STEP_SIZE"] = vals[2]
                        elif "INITIAL_STEP_SIZE" == vals[1]:
                            self.simulation["INITIAL_STEP_SIZE"] = vals[2]
                        else:
                            raise Exception("Invalid SIMULATION argument: %s"
                                            % vals[1])
                        continue
                elif re.search(r'^LATEX\s', line):
                    vals = line.split()
                    if len(vals) > 2:
                        self.latex[vals[1]] = vals[2]
                    continue
                elif re.search(r'^AUTHOR\s', line):
                    vals = line.split()
                    if len(vals) > 2:
                        self.author = " ".join(line.split()[1:])
                    continue
                elif re.search(r'^TITLE\s', line):
                    vals = line.split()
                    if len(vals) > 2:
                        self.title = " ".join(line.split()[1:])
                    continue
                elif re.search(r'^KEYWORD\s', line):
                    vals = line.split()
                    if len(vals) > 2:
                        self.keywords.append(" ".join(line.split()[1:]))
                    continue
                elif re.search(r'^BIBITEM\s', line):
                    vals = line.split()
                    if len(vals) > 2:
                        self.bibitem.append(" ".join(line.split()[1:]))
                    continue

                r.reaction(line, nr, self)
                self.reactions.append(r)
                r = R()
                nr += 1

        fp.close()

        for a in self.excess:
            if a not in self.initial:
                raise Exception("Missing initial value for: %s" % a)

        for p in ["T_END", "T_POINTS"]:
            if p not in self.simulation:
                raise Exception("Missing SIMULATIONS parameter: %s" % p)

        for r in self.reactions:
            for s in r.species():
                i = s.strip()
                if not(i in self.system_species):
                    self.system_species.append(i)
    
        self.n_reactions = len(self.reactions)
        self.n_equations = len(self.xdot_raw)


def xsmc(r, x):
    c = 0
    for i in map(str.strip, r.split("+")):
        if x == i:
            c += 1
    return c


def proc_rspcs():
    for r in rctns:
        for s in r.species():
            i = s.strip()
            if not(i in rspcs):
                rspcs.append(i)


def symbols2(x):
    r = []
    a = "+-*"

    for i in " ".join(x.split()).split(" "):
        i = i.strip()
        if i is None:
            pass
        elif i.isdigit() or i in a:
            r.append(i)
        else:
            r.append("__%s__" % i)

    return " ".join(r)

                
def symbols(x):
    r = []

    for i in " ".join(x.split()).split(" "):
        i = i.strip()
        if i is None:
            pass
        elif i.isdigit() or "+" == i or "-" == i or "*" == i:
            r.append(i)
        else:
            r.append("Symbol('__%s__')" % i)

    return " ".join(r)


def pretty_kinet(s, y, k):
    if 1 == y:
        return "%s %s" % (s, k)
    return "%s %d * %s" % (s, y, k)


def proc_r(ers):

    for s in ers.system_species:
        if s in ers.excess or s in ers.constants:
            continue
        ers.x.append(s)

    for r in ers.reactions:
        k = r.kinet(True)
        if k:
            v = "fkinet%d" % r.i
            ers.kinet[v] = subst_x(ers, str(symbols2(k)))
            ers.kinet_keys.append(v)
        k = r.kinet(False)
        if k:
            v = "rkinet%d" % r.i
            ers.kinet[v] = subst_x(ers, str(symbols2(k)))
            ers.kinet_keys.append(v)

    for s in ers.system_species:
        
        if s in ers.excess or s in ers.constants:
            continue
        
        a = []
        b = []
        for r in ers.reactions:
            # A + B -> C + D
            k = r.kinet(True)
            z = "fkinet%d" % r.i
            y = r.smc(True, s)
            if y > 0 and k:
                a.append("- %d * %s" % (y, k))
                # b.append("- %d * %s" % (y, z))
                b.append(pretty_kinet("-", y, z))
            y = r.smc(False, s)
            if y > 0 and k:
                a.append("+ %d * %s" % (y, k))
                # b.append("+ %d * %s" % (y, z))
                b.append(pretty_kinet("+", y, z))
                
            # A + B <- C + D
            k = r.kinet(False)
            z = "rkinet%d" % r.i
            y = r.smc(False, s)
            if y > 0 and k:
                a.append("- %d * %s" % (y, k))
                # b.append("- %d * %s" % (y, z))
                b.append(pretty_kinet("-", y, z))
            y = r.smc(True, s)
            if y > 0 and k:
                a.append("+ %d * %s" % (y, k))
                # b.append("+ %d * %s" % (y, z))
                b.append(pretty_kinet("+", y, z))

        ers.xdot_raw.append(" ".join(b))

        if config["has_sympy"]:
            global df
            cmd = "df = %s" % symbols(" ".join(a))
            exec(cmd, globals())
            dbg(df)
            if config["octave"]:
                ers.xdot.append(subst_x(ers, str(df)))
            if config["latex"]:
                ers.xdot.append(str(df))


def subst_x(ers, e):
    c = list(ers.constants.keys())
    
    for s in ers.system_species + c:
        s1 = "__%s__" % s
        
        if s in ers.excess or s in c:
            s2 = " %s " % s
        else:
            i = 1 + ers.x.index(s)
            s2 = " x(%i) " % i
            
        e = e.replace(s1, s2)
        
    e = re.sub(r'__(k[f,r])(\d+)__', r' \1(\2) ', e)

    return " ".join(e.split())


def latex_sub(ers, s, eq=True):
    for j in ers.latex:
        if s == j:
            if eq:
                return ers.latex[j]
            else:
                return "$" + ers.latex[j] + "$"

    return s


def latex_sub2(ers, s):
    # FIXME constants, excess
    x = s
    x = re.sub(r'\*\*(\d{1,})', r'^{\1}', s)
    x = re.sub(r'(__kf)(\d{1,})(__)', r'k_{\2}', x)
    x = re.sub(r'(__kr)(\d{1,})(__)', r'k_{-\2}', x)
    x = re.sub(r'(__)([a-zA-Z0-9]{1,})(__)', r'[\\mathrm{\2}]', x)
    x = re.sub(r'(\d{1,})(\*)', r'\1', x)
    x = re.sub(r'\*', r'', x)

    if len(ers.constants):
        for a in ers.constants:
            if a in ers.latex:
                x = x.replace("[\\mathrm{" + a + "}]", latex_sub(ers, a))
            else:
                x = x.replace("[\\mathrm{" + a + "}]", "\\mathrm{%s}" % a)

    for a in ers.latex:
        l = "{[\\mathrm{%s}]}" % ers.latex[a]
        x = x.replace("[\\mathrm{" + a + "}]", l)

    return x


def latex_reaction_fmt(ers, s):
    res = []
    seen = {}
    
    for i in s.split():
        x = i

        if seen.get(x):
            continue

        if "+" == x:
            continue

        seen[x] = True
        n = 0

        for j in s.split():
            if x == j:
                n += 1

        if x not in ers.latex:
            if n > 1:
                res.append(str(n) + " " + x)
            else:
                res.append(x)                    
            continue
                
        for j in ers.latex:
            if i == j:
                x = "$\\mathrm{%s}$" % ers.latex[j]
                if n > 1:
                    res.append(str(n) + " " + x)
                else:
                    res.append(x)                    
                break

    return " + ".join(res)


def chembi(i):
    return "$\\mathop{\\kern0pt \\rightleftharpoons}\\limits^{k_{%d}}_{k_{-%d}}$" % (i, i)


def chemuni(i):
    return "$\\mathop{\\kern0pt \\rightarrow}\\limits^{k_{%d}}$" % i


def eqnarray_rhs(s, label):
    m = 3
    n = 0
    r = []
    a = s.split()
    c = 0
    
    for i in a:
        c += 1
        if i != "+" and i != "-":
            n += 1
        r.append(i)
        if n == m and c < len(a):
            r.append("\\nonumber \\\\\n")
            r.append(" & & ")
            n = 0

    r.append("\\label{rate:%s}" % label)

    return " ".join(r)

# FIXME unit to separate function
# Reaction rate coefficient formatter
def latex_rc(s, a):
    notcounted = ["H2O"]
    # FIXME add constants to notcounted
    n = 0

    for i in a.split():
        if "+" != i and i not in notcounted:
            n += 1

    x = re.sub(r'E([\+\-.]*\d{1,})', r' $\\times 10^{\1}$ ', str(s))

    if 0 == n:
        # FIXME Better croak here?
        print("WARNING: Suspicious values latex_rc(%s, %s)" % (s, a))
        return x
    elif 1 == n:
        return "%s s$^{-1}$" % x
    
    return "%s mol$^{-%d}$dm$^{%d}$s$^{-1}$" % (x, (n - 1), (3 * (n - 1)))


def cmp_rate(a, b):
    isa = bool(re.match(r'__k\d{1,}__', a.lower()))
    isb = bool(re.match(r'__k\d{1,}__', b.lower()))

    if isa and isb:
        r = cmp(a, b)
    elif isa:
        r = -1
    else:
        r = 1

    return r

def cmp_dx(a, b):
    # Numeric values
    isa = bool(re.match(r'\d{1,}', a.lower()))
    isb = bool(re.match(r'\d{1,}', b.lower()))

    if isa:
        return -1
    elif isb:
        return 1

    # kf/kr
    isa = bool(re.match(r'__k(f|r)*\d{1,}__', a.lower()))
    isb = bool(re.match(r'__k(f|r)*\d{1,}__', b.lower()))

    if isa and isb:
        r = cmp(a, b)
    elif isa:
        r = -1
    else:
        r = 1

    return r


def rate_sorted(v):
    r = v.split("*")
    r.sort()
    return "*".join(r)


def dx_sorted(v):
    r = []

    # Defensive hacks to ensure expected formatting
    v = re.sub(r'^-', r'- ', v)
    v = re.sub(r'^\+', r'+ ', v)

    for i in v.split():
        if "+" == i or "-" == i:
            r.append(i)
        else:
            # Hack to preserve '**2'
            i = i.replace("**2", "^2")
            t = i.split("*")
            t.sort()
            r.append("*".join(t).replace("^2", "**2"))

    return " ".join(r)


def latex_vrate(ers, v):
    global tmp_vrate
    exec("tmp_vrate = str(%s)" % symbols(" * ".join(v.split("*"))), globals())
    dbg(f"VRATE = {tmp_vrate}")
    r = rate_sorted(tmp_vrate)
    r = latex_sub2(ers, r)
    return r


def latex_output(ers, fbase, src):
    fname = "%s.tex" % fbase
    n = len(ers.reactions)
    neq = len(ers.xdot_raw)
    
    try:
        fp = open(fname, "w")
    except:
        raise

    fp.write("%%%% START SRC %s\n" % os.path.basename(src))

    with open(src) as f:
        for l in f:
            fp.write("%%%% %s" % l)

    fp.write("%%%% END SRC %s\n" % os.path.basename(src))

    fp.write("\n")
    fp.write("\\documentclass[12pt]{article}\n"
             "\\usepackage{tabularx}\n"
             "\\usepackage{adjustbox}\n"
             "\\usepackage{graphicx}\n"
             "\\usepackage[margin=0.75in]{geometry}\n"
             "\\begin{document}\n"
             "\\title{%s}\n"
             "\\date{\\today}\n"
             "\\author{%s}\n"
             "\\maketitle\n\n" % (ers.title, ers.author))
             
    fp.write("\\section{Abstract}\n"
             "\\IfFileExists{abstract.tex}{\n\\input{abstract}\n}{\n"
             "Results from simulations of a %d--reaction, %d--species\n"
             "chemical reaction system model are reported.\n}\n\n"
             % (n, neq))

    fp.write("\\section{Introduction}\n\n"
             "\\IfFileExists{introduction.tex}{\n\\input{introduction}"
             "\n}{\n%%%% Tabula rasa\n}\n\n")

    fp.write("\\section{Results}\n\n"
             "\\IfFileExists{results.tex}{\n\\input{results}"
             "\n}{\n"
             "The %d reaction system model in Table~\\ref{table:reactions}\n"
             "was converted to ordinary differential equations programmatically.\n"
             "Then the resulting %d ordinary differention equations were\n"
             "integrated by LSODA, which can handle both non--stiff\n"
             "and stiff initial value problems~\\cite{lsoda}.\n"
             "Numerical parameters were: relative tolerance %s,\n"
             "absolute tolerance %s and maximum allowed time step %s.\n"
             "\n}\n\n" % (n, neq, ers.lsode_rtol, ers.lsode_atol,
                          ers.simulation.get("MAXIMUM_STEP_SIZE", "unlimited")))
    
    sx = 0
    for r in ers.reactions:
        if len(r.text) > sx:
            sx = len(r.text)

    # Elementary reaction table
    fp.write("\\begin{table}\n\\label{table:reactions}\n"
             "\\begin{adjustbox}{width=\\textwidth}\n"
             "\\begin{tabular}{lrclll}\n")

    fp.write(" & \multicolumn{3}{c}{Reactions} & $k_{forward}$ & $k_{reverse}$\\\\\n\\hline\n")
    sx += 3
    fmt = "  %%-%ds (R%%d)\n" % sx
    nr = 0
    ne = 0
    nekinet = []

    for r in ers.reactions:
        nr += 1
        fp.write("%%")
        fp.write(fmt % (r.text, r.i))
        fp.write("(R%d) & " % r.i)
        fp.write("%s & " % latex_reaction_fmt(ers, r.reactants))

        if [1, 0] == r.rates:
            fp.write("%s & " % chemuni(r.i))
        else:
            fp.write("%s & " % chembi(r.i))
        
        fp.write("%s " % latex_reaction_fmt(ers, r.products))

        # Forward reaction rates (or exceptions)
        fp.write(" & ")
        if r.rates[0] and not(r.frate):
            fp.write("%s" % latex_rc(r.kf, r.reactants))
        elif r.rates[0] and r.frate:
            ne += 1
            nekinet.append("$v_{%d} = %s$" % (r.i, latex_vrate(ers, r.frate)))
            fp.write(" ($v_{%d}$) " % r.i)

        # Reverse reactions rates (or exceptions)
        fp.write(" & ")
        if r.rates[1] and not(r.rrate):
            fp.write(" %s " % latex_rc(r.kr, r.products))
        elif r.rates[1] and r.rrate:
            ne += 1
            nekinet.append("$v_{-%d}$ = %s" % (r.i, r.rrate))
            fp.write(" ($v_{-%d}$) " % r.i)

        fp.write("\\\\\n")
    
    fp.write("\n")
        
    fp.write("\\end{tabular}\n"
             "\\end{adjustbox}\n"
             "\\caption{Elementary reactions and rate constants of the model system.\n")

    if len(nekinet):
        fp.write("Exceptions in mass action kinetic model: %s.\n"
                 % ", ".join(nekinet))

    if len(ers.constants):
        c = []
        for a in ers.constants:
            c.append("%s = %s" % (latex_sub(ers, a, False), ers.constants[a]))
        fp.write("Constants: %s." % ", ".join(c))
        
    fp.write("}\n"
             "\\label{tab1e:reactions}\n"
             "\\end{table}\n\n")

    fp.write("\\begin{eqnarray}\n")

    i = 0
    for dx in ers.xdot:
        dx = dx_sorted(dx)
        dx = latex_sub2(ers, dx)
        dx = eqnarray_rhs(dx, ers.x[i])
        fp.write("%% /* %s */\n" % ers.x[i])
        fp.write("\\frac{d [\\mathrm{%s}]}{dt} & = & %s" %
                 (latex_sub(ers, ers.x[i]), dx))
        if i != (neq - 1):
            fp.write("\\\\")
        fp.write("\n")
        i += 1
    
    fp.write("\\end{eqnarray}\n\n")
    
    fp.write("\\section{Conclusions}\n\n"
             "\\IfFileExists{conclusions.tex}{\n\\input{conclusions}"
             "\n}{\n%%%% Tabula rasa\n}\n\n")

    fp.write("\\section{Acknowledgments}\n\n"
             "\\IfFileExists{acknowledgments.tex}{\n\\input{acknowledgments}"
             "\n}{\n%%%% Tabula rasa\n}\n\n")

    fp.write("\\section{Keywords}\n\n")
    if len(ers.keywords):
        fp.write(", ".join(ers.keywords))
    fp.write("\n\n")

    fp.write("\\begin{thebibliography}{99}\n"
             "\\bibitem{lsoda} A.C. Hindmarsh, {\\em ODEPACK, A Systematized Collection of ODE Solvers}, in {\\em Scientific Computing}, R.S. Stepleman et al. (Eds.), North--Holland, Amsterdam, {\\bf 1983}, pp. 55-64.\n")
    for i in ers.bibitem:
        fp.write("\\bibitem{%s} %s\n"
                 % (i.split()[0], " ".join(i.split()[1:])))
    fp.write("\\end{thebibliography}\n\n")

    fp.write("\\end{document}\n")
    
    i = 0
    fp.write("\n")
    for v in ers.x:
        fp.write("%%")
        fp.write("  x_%d = %s\n" % (1 + i, v))
        i += 1

    fp.write("\n")

    if len(ers.excess):
        fp.write("%%%% Species in excess: %s\n" % " ".join(ers.excess))
        for a in ers.excess:
            if ers.initial.get(a):
                fp.write("%%%% #define %s (%s)\n" % (a, ers.initial[a]))
        fp.write("\n")

    fp.close()

    try:
        os.chmod(fname, 0o644)
    except:
        pass


def lsoda_c_output(ers, fbase):
    fname = "%s.c" % fbase
    n = len(ers.reactions)
    neq = len(ers.xdot_raw)

    with open("lsoda.c") as fp:
        lsodac = fp.read()
    
    try:
        fp = open(fname, "w")
        fp.write(lsodac)
        fp.write("#ifndef __ERHELPER_H__\n#define __ERHELPER_H__\n\n/*\n")
        
        sx = 0
        for r in ers.reactions:
            if len(r.text) > sx:
                sx = len(r.text)

        sx += 3
        fmt = "  %%-%ds (R%%d)\n" % sx
        for r in ers.reactions:
            fp.write(fmt % (r.text, r.i))

        i = 0
        fp.write("\n")
        for v in ers.x:
            fp.write("  x(%d) %s\n" % (1 + i, v))
            i += 1

        fp.write("\n  Plot command:\n\n xplot.sh -t '%s' FILE"
                 % " ".join(ers.x))
            
        fp.write("\n*/\n")

        defs = ("\n#define x(i) (x[i-1])\n"
                "#define x0(i) (x[i])\n"
                "#define xdot(i) (xdot[i-1])\n"
                "#define kf(i) kf[i]\n"
                "#define kr(i) kr[i]\n")
        
        fp.write("%s" % defs)

        if len(ers.excess):
            fp.write("\n/* Species in excess %s */\n" % " ".join(ers.excess))
            for a in ers.excess:
                if ers.initial.get(a):
                    fp.write("#define %s (%s)\n" % (a, ers.initial[a]))

        if len(ers.constants):
            fp.write("\n/* Constants %s */\n" % " ".join(ers.constants))
            for a in ers.constants:
                fp.write("#define %s (%s)\n" % (a, ers.constants[a]))

        fp.write("\n#define NEQ %d\n" % neq)
        fp.write("#define N_REACTIONS %d\n" % n)
        fp.write("#define T_END %s\n" % ers.simulation["T_END"])
        fp.write("#define T_DELTA (1/ (double) %s)\n" % ers.simulation["T_POINTS"])
        fp.write("#define LSODE_ATOL %s\n" % ers.lsode_atol)
        fp.write("#define LSODE_RTOL %s\n" % ers.lsode_rtol)
        fp.write("#define LSODE_HMAX %s\n" %
                 ers.simulation.get("MAXIMUM_STEP_SIZE", "0.0"))
        fp.write("#define LSODE_H0 %s\n" %
                 ers.simulation.get("INITIAL_STEP_SIZE", "0.0"))

        fp.write("\nstatic double kf[N_REACTIONS+1], kr[N_REACTIONS+1];\n")
        fp.write("int get_neq(void) { return NEQ; }\n")
        fp.write("\nvoid ersimu_init(double *x, double *abstol, double *reltol, double *t_end, double *dt, double *h0, double *hmax)\n{\n")

        i = 0
        while i <= neq:
            fp.write("\tabstol[%d] = LSODE_ATOL;\n" % i)
            fp.write("\treltol[%d] = LSODE_RTOL;\n" % i)
            fp.write("\tx0(%d) = 0;\n" % i)
            i += 1
        
        fp.write("\n\t/* initial conditions */\n")
        for a in ers.x:
            if a in ers.initial:
                i = 1 + ers.x.index(a)
                fp.write("\t/* %s */\n\tx0(%d) = %s ;\n" % (a, i, ers.initial[a]))

        fp.write("\n\t/* forward reaction rates */\n")
        for r in ers.reactions:
            if r.rates[0]:
                fp.write("\tkf(%d) = %s ;\n" % (r.i, r.kf))
        fp.write("\n\t/* reverse reaction rates */\n")
        for r in ers.reactions:
            if r.rates[1]:
                fp.write("\tkr(%d) = %s ;\n" % (r.i, r.kr))

        fp.write("\n\t*t_end = T_END;\n\t*dt = T_DELTA;\n"
                 "\t*h0 = LSODE_H0;\n\t*hmax = LSODE_HMAX;\n"
                 "\n}\n")
                
        fp.write("\nvoid fex(double t, double *x, double *xdot, void *data)\n{\n")

        for i in ers.kinet_keys:
            fp.write("\tdouble %s = %s ;\n" % (i, ers.kinet[i]))
        fp.write("\n")

        # FIXME zeros
        i = 0
        for dx in ers.xdot_raw:
            fp.write("\t/* %s */\n" % ers.x[i])
            fp.write("\txdot(%d) = %s ; \n" % (i + 1, ers.xdot_raw[i]))
            i += 1

        fp.write("\n}\n\n#endif\n")

    except:
        raise
    finally:
        fp.close()

    try:
        os.chmod(fname, 0o644)
    except:
        pass


def octave_output(ers, fbase):
    fname = "%s.m" % fbase
    mname = "%s.dat" % fbase
    n = len(ers.reactions)

    try:
        fp = open(fname, "w")

        fp.write("#!/usr/bin/octave -qf\n")
        fp.write("# -*-octave-*-\n")

        dump_model_in_comments(ers, fp, "#")

        fp.write("more off\n")
        fp.write("global kf kr ;\n")

        if len(ers.excess):
            fp.write("\n# Species in excess\nglobal %s ;\n" % " ".join(ers.excess))
            for a in ers.excess:
                if ers.initial.get(a):
                    fp.write("%s = %s ;\n" % (a, ers.initial[a]))

        if len(ers.constants):
            fp.write("\n# Constants\nglobal %s ;\n" % " ".join(ers.constants))
            for a in ers.constants:
                fp.write("%s = %s ;\n" % (a, ers.constants[a]))

        fp.write("\n# forward reaction rates\nkf = zeros(%d, 1) ;\n" % n)
        for r in ers.reactions:
            if r.rates[0]:
                fp.write("kf(%d) = %s ;\n" % (r.i, r.kf))

        fp.write("\n# reverse reaction rates\nkr = zeros(%d, 1) ;\n" % n)
        for r in ers.reactions:
            if r.rates[1]:
                fp.write("kr(%d) = %s ;\n" % (r.i, r.kr))

        fp.write("\n# initial conditions\nx0 = zeros(%d, 1) ;\n" % len(ers.xdot))

        # FIXME Jacobian

        for a in ers.x:
            if a in ers.initial:
                i = 1 + ers.x.index(a)
                fp.write("# %s\nx0(%d) = %s ;\n" % (a, i, ers.initial[a]))
        fp.write("\nfunction xdot = f (x, t)\n")
        fp.write("    global kf kr ;\n")

        if len(ers.excess):
            fp.write("    global %s ;\n" % " ".join(ers.excess))
        if len(ers.constants):
            fp.write("    global %s ;\n" % " ".join(ers.constants))

        fp.write("    xdot = zeros(%d, 1) ;\n" % len(ers.xdot))
        i = 0
        for dx in ers.xdot:
            fp.write("    xdot(%d) = %s ; \n" % (i + 1, ers.xdot[i]))
            i += 1
        
        fp.write("endfunction\n")

        fp.write("\nlsode_options(\"integration method\", \"stiff\") ;\n")

        if "MAXIMUM_STEP_SIZE" in ers.simulation:
            fp.write("lsode_options(\"maximum step size\", %s) ;\n" %
                     ers.simulation["MAXIMUM_STEP_SIZE"])

        fp.write("lsode_options(\"absolute tolerance\", %s) ;\n" % ers.lsode_atol)
        fp.write("lsode_options(\"relative tolerance\", %s) ;\n" % ers.lsode_rtol)
        fp.write("\nt_end = %s ;\n" % ers.simulation["T_END"])
        fp.write("t_points = %s ;\n" % ers.simulation["T_POINTS"])
        fp.write("t = linspace(0, t_end, t_end * t_points) ;\n")
        fp.write("tstart = cputime ;\n")
        fp.write("  [y, istate, msg] = lsode (\"f\", x0, t) ;\n")
        fp.write("lsode_options\n")
        fp.write("istate\n")
        fp.write("msg\n")
        fp.write("printf('Total CPU time: %f seconds\\n', cputime - tstart) ;\n")
        fp.write("if istate != 2\n    exit(1) ;\nendif\n")
        fp.write("mat = [rot90(t, -1), y] ;\n")
        fp.write("save %s mat ;\n" % mname)
  
    except:
        raise
    finally:
        fp.close()

    try:
        os.chmod(fname, 0o755)
    except:
        pass


def dump_model_in_comments(ers, fp, comment_string="#"):
    sx = 0
    for r in ers.reactions:
        if len(r.text) > sx:
            sx = len(r.text)

    sx += 3
    fmt = "%s %%-%ds (R%%d)\n" % (comment_string, sx)
    for r in ers.reactions:
        fp.write(fmt % (r.text, r.i))


def scipy_output(ers, fbase):
    fname = f"{fbase}.py"
    n = len(ers.reactions)

    try:
        fp = open(fname, "w")

        fp.write("#!/usr/bin/env python3\n")
        fp.write("# -*-python-*-\n")
        fp.write("""
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
import math
\n""")

        dump_model_in_comments(ers, fp, "#")

        if len(ers.excess):
            fp.write("\n# Species in excess\n#global %s ;\n" % " ".join(ers.excess))
            for a in ers.excess:
                if ers.initial.get(a):
                    fp.write("%s = %s ;\n" % (a, ers.initial[a]))

        if len(ers.constants):
            fp.write("\n# Constants\n#global %s ;\n" % " ".join(ers.constants))
            for a in ers.constants:
                fp.write("%s = %s ;\n" % (a, ers.constants[a]))

        fp.write("\n# forward reaction rates\nkf = np.zeros(%d) ;\n" % n)
        for r in ers.reactions:
            if r.rates[0]:
                fp.write("kf[%d] = %s ;\n" % (r.i - 1, r.kf))

        fp.write("\n# reverse reaction rates\nkr = np.zeros(%d) ;\n" % n)
        for r in ers.reactions:
            if r.rates[1]:
                fp.write("kr[%d] = %s ;\n" % (r.i - 1, r.kr))

        fp.write("\n# initial conditions\nx0 = np.zeros(%d) ;\n" % len(ers.xdot_raw))
        for a in ers.x:
            if a in ers.initial:
                i = ers.x.index(a)
                fp.write("# %s\nx0[%d] = %s ;\n" % (a, i, ers.initial[a]))

        fp.write("\ndef f(t, x):\n")
        fp.write("    xdot = np.zeros((%d, 1))\n" % len(ers.xdot_raw))

        for i in ers.kinet_keys:
            e = re.sub(r'([x,kf,kr])\((\d+)\)', lambda z: z.groups()[0] + "[" + str(int(z.groups()[1]) - 1) + "]", ers.kinet[i])
            fp.write("    %s = %s\n" % (i, e))
        fp.write("\n")

        i = 0
        for dx in ers.xdot_raw:
            fp.write("    xdot[%d] = %s\n" % (i, ers.xdot_raw[i]))
            i += 1

        fp.write("\n    return xdot\n\n\n")

        fp.write("max_step_ =  %s\n" % ers.simulation.get("MAXIMUM_STEP_SIZE", "np.inf"))
        fp.write("atol_ = %s\n" % ers.lsode_atol)
        fp.write("rtol_ = %s\n" % ers.lsode_rtol)
        fp.write("t_end = %s\n" % ers.simulation["T_END"])
        fp.write("t_points = %s\n" % ers.simulation["T_POINTS"])
        fp.write("t = np.linspace(0, t_end, t_end * t_points)\n")
        fp.write("sol = solve_ivp(f, (0, t_end), x0, t_eval=t, "
                 "method='LSODA', dense_output=True, vectorized=True, "
                 "max_step=max_step_, atol=atol_, rtol=rtol_)\n")

        n = len(ers.x)
        fp.write(f"fig, ax = plt.subplots({n}, 1, figsize=(15, 8))\n")
        for i in range(n):
            subplot = f"[{i}]" if n > 1 else ""
            fp.write(f"ax{subplot}.plot(sol.t.T, sol.y[{i}].T)\n"
                     f"ax{subplot}.set_ylabel('{ers.x[i]}')\n")
        fp.write(f"plt.show()\nplt.savefig('{fbase}.pdf', bbox_inches='tight')\n")

    except:
        raise
    finally:
        fp.close()

    try:
        os.chmod(fname, 0o755)
    except:
        pass


def validate_input():
    for a in excess:
        if a not in initial:
            raise Exception("Missing initial value for: %s" % a)

    for p in ["T_END", "T_POINTS"]:
        if p not in simulation:
            raise Exception("Missing SIMULATIONS parameter: %s" % p)


def main(argv):
    opts = get_arg_parser().parse_args()

    if opts.latex or opts.octave:
        if not(config["has_sympy"]):
            err("sympy is not available")

    if "ersimu" == opts.name or "lsoda" == opts.name:
        err("bad name")
    elif "-" in opts.name:
        err("name can not contain dashes")

    config["verbose"] = opts.verbose
    config["octave"] = opts.octave
    config["latex"] = opts.latex
    config["name"] = opts.name
    fname = opts.inputfile

    ers = ERSimu()
    ers.load(fname)
    proc_r(ers) # FIXME move to ERSimu
    
    for r in ers.reactions:
        dbg(f"reaction = {r.i}")
        dbg(f"reactants = {r.reactants}")
        dbg(f"products = {r.products}")
        dbg(f"rates = {r.rates}")
        dbg(f"fkinet = {r.kinet()}")
        dbg(f"rkinet = {r.kinet(False)}")
        dbg(f"species = {r.species()}")
    dbg("SYSTEM SPECIES %s" % ers.system_species)
    dbg("X %s" % ers.x)
    dbg("EXCESS %s" % ers.excess)
    dbg("INITIAL %s" % ers.initial)
    dbg("CONSTANTs %s" % ers.constants)

    for k in ers.kinet:
        dbg("KINET %s = %s" % (k, ers.kinet[k]))

    i = 0
    for dx in ers.xdot_raw:
        if config["latex"] or config["octave"]:
            dbg("xdot(%d) = %s ; " % (i + 1, ers.xdot[i]))
        dbg("xdot_raw(%d) = %s ; " % (i + 1, ers.xdot_raw[i]))
        i += 1

    for ltx in ers.latex:
        dbg("LATEX %s = %s" % (ltx, ers.latex[ltx]))

    if opts.scipy:
        scipy_output(ers, opts.name)
        if opts.run:
            exec(f"import {opts.name}")
    elif opts.octave:
        octave_output(ers, opts.name)
    elif opts.latex:
        latex_output(ers, opts.name, fname)
    else:
        lsoda_c_output(ers, opts.name)


if "__main__" == __name__:
    main(sys.argv[1:])
