#!/usr/bin/env python

import sys
import string

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

config["verbose"] = 0
config["opt_octave"] = False
config["opt_latex"] = False
config["title"] = None
config["author"] = None
config["bibitem"] = []
config["keywords"] = []

rctns = list(())
x = list(())
excess = list(())
xdot = list(())
xdot_raw = list(())
rspcs = []
initial = {}
constant = {}
simulation = {}
latex = {}
kf = {}
kr = {}
kinet = {}
kinet_keys = []
lsode_atol = "1E-6"
lsode_rtol = "1E-6"

def dbg(msg):
    if config["verbose"] > 0:
        print("DEBUG: %s" % str(msg))


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

    def reaction(self, r, n):
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
            [rct, rts] = map(string.strip, r.split("@"))
            if 1 == len(rts.split()) and " ==> " == cs:
                self.frate = rts.split()[0]
            elif 2 == len(rts.split()) and " <==> " == cs:
                [self.frate, self.rrate] = rts.split()
            else:
                raise Exception("Malformed rate input: \"%s\" (%d)" % (r, n))
        elif "|" in r:
            [rct, rts] = map(string.strip, r.split("|"))
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
                if not(d in constant):
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
        
        for i in map(string.strip, r.split("+")):
            if y == i:
                c += 1
                
        return c

    
def xsmc(r, x):
    c = 0
    for i in map(string.strip, r.split("+")):
        if x == i:
            c += 1
    return c


def read_r(fname):
    with open(fname) as f:
        n = 0
        nr = 1
        r = R()
        for l in f:
            n +=1
            line = l.strip()

            if config["verbose"] > 1:
                dbg("LINE %05d: %s" % (n, line))

            if len(line) > 0 and '#' != line[0]:
                if re.search(r'^EXCESS\s', line):
                    if len(line.split()) > 1:
                        for a in line.split()[1:]:
                            if not a in excess:
                                excess.append(a)
                    continue
                elif re.search(r'^INITIAL\s', line):
                    if len(line.split()) > 2:
                        a, val = line.split()[1:3]
                        initial[a] = val
                    continue
                elif re.search(r'^CONSTANT\s', line):
                    if len(line.split()) > 2:
                        a, val = line.split()[1:3]
                        constant[a] = val
                    continue
                elif re.search(r'^SIMULATION\s', line):
                    vals = line.split()
                    if len(vals) > 2:
                        if "T_END" == vals[1]:
                            simulation["T_END"] = vals[2]
                        elif "T_POINTS" == vals[1]:
                            simulation["T_POINTS"] = vals[2]
                        elif "MAXIMUM_STEP_SIZE" == vals[1]:
                            simulation["MAXIMUM_STEP_SIZE"] = vals[2]
                        elif "INITIAL_STEP_SIZE" == vals[1]:
                            simulation["INITIAL_STEP_SIZE"] = vals[2]
                        else:
                            raise Exception("Invalid SIMULATION argument: %s"
                                            % vals[1])
                        continue
                elif re.search(r'^LATEX\s', line):
                    vals = line.split()
                    if len(vals) > 2:
                        latex[vals[1]] = vals[2]
                    continue
                elif re.search(r'^AUTHOR\s', line):
                    vals = line.split()
                    if len(vals) > 2:
                        config["author"] = " ".join(line.split()[1:])
                    continue
                elif re.search(r'^TITLE\s', line):
                    vals = line.split()
                    if len(vals) > 2:
                        config["title"] = " ".join(line.split()[1:])
                    continue
                elif re.search(r'^KEYWORD\s', line):
                    vals = line.split()
                    if len(vals) > 2:
                        config["keywords"].append(" ".join(line.split()[1:]))
                    continue
                elif re.search(r'^BIBITEM\s', line):
                    vals = line.split()
                    if len(vals) > 2:
                        config["bibitem"].append(" ".join(line.split()[1:]))
                    continue

                r.reaction(line, nr)
                rctns.append(r)
                r = R()
                nr += 1


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


def proc_r():

    for s in rspcs:
        if s in excess or s in constant:
            continue
        x.append(s)

    for r in rctns:
        k = r.kinet(True)
        if k:
            v = "fkinet%d" % r.i
            kinet[v] = subst_x(str(symbols2(k)))
            kinet_keys.append(v)
        k = r.kinet(False)
        if k:
            v = "rkinet%d" % r.i
            kinet[v] = subst_x(str(symbols2(k)))
            kinet_keys.append(v)

    for s in rspcs:
        
        if s in excess or s in constant:
            continue
        
        a = []
        b = []
        for r in rctns:
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

        xdot_raw.append(" ".join(b))

        if config["has_sympy"]:
            exec("df = %s" % symbols(" ".join(a)))
            dbg(df)
            if config["opt_octave"]:
                xdot.append(subst_x(str(df)))
            if config["opt_latex"]:
                xdot.append(str(df))



def subst_x(e):
    c = constant.keys()
    
    for s in rspcs + c:
        s1 = "__%s__" % s
        
        if s in excess or s in c:
            s2 = " %s " % s
        else:
            i = 1 + x.index(s)
            s2 = " x(%i) " % i
            
        e = e.replace(s1, s2)
        
    e = re.sub(r'__(k[f,r])(\d+)__', r' \1(\2) ', e)

    return " ".join(e.split())


def latex_sub(s, eq=True):
    for j in latex:
        if s == j:
            if eq:
                return latex[j]
            else:
                return "$" + latex[j] + "$"

    return s


def latex_sub2(s):
    # FIXME constants, excess
    x = s
    x = re.sub(r'\*\*(\d{1,})', r'^{\1}', s)
    x = re.sub(r'(__kf)(\d{1,})(__)', r'k_{\2}', x)
    x = re.sub(r'(__kr)(\d{1,})(__)', r'k_{-\2}', x)
    x = re.sub(r'(__)([a-zA-Z0-9]{1,})(__)', r'[\\mathrm{\2}]', x)
    x = re.sub(r'(\d{1,})(\*)', r'\1', x)
    x = re.sub(r'\*', r'', x)

    if len(constant):
        for a in constant:
            if a in latex:
                x = x.replace("[\\mathrm{" + a + "}]", latex_sub(a))
            else:
                x = x.replace("[\\mathrm{" + a + "}]", "\\mathrm{%s}" % a)

    for a in latex:
        l = "{[\\mathrm{%s}]}" % latex[a]
        x = x.replace("[\\mathrm{" + a + "}]", l)

    return x


def latex_reaction_fmt(s):
    res = []
    seen = {}
    
    for i in s.split():
        x = i

        if seen.has_key(x):
            continue

        if "+" == x:
            continue

        seen[x] = True
        n = 0

        for j in s.split():
            if x == j:
                n += 1

        if x not in latex:
            if n > 1:
                res.append(str(n) + " " + x)
            else:
                res.append(x)                    
            continue
                
        for j in latex:
            if i == j:
                x = "$\\mathrm{%s}$" % latex[j]
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
    # FIXME add contants to notcounted
    n = 0

    for i in a.split():
        if "+" != i and i not in notcounted:
            n += 1
    
    x = re.sub(r'E([\+\-.]*\d{1,})', r' $\\times 10^{\1}$ ', s)

    if 0 == n:
        # FIXME Better croak here?
        print("WARNING: Suspicious values latex_rc(%s, %s)" % (s, a))
        return x
    elif 1 == n:
        return "%s $s^{-1}$" % x
    
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
    r.sort(cmp_rate)
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
            t.sort(cmp_dx)
            r.append("*".join(t).replace("^2", "**2"))

    return " ".join(r)


def latex_vrate(v):
    exec("r = str(%s)" % symbols(" * ".join(v.split("*"))))
    dbg("VRATE = %s" % r)
    r = rate_sorted(r)
    r = latex_sub2(r)
    return r


def latex_output(fbase, src):
    fname = "%s.tex" % fbase
    n = len(rctns)
    neq = len(xdot_raw)
    
    try:
        fp = open(fname, "w")
    except:
        raise

    fp.write("%%%% START SRC %s\n" % os.path.basename(src))
    dbg("LaTeX FILE %s" % fname)

    with open(src) as f:
        for l in f:
            fp.write("%%%% %s" % l)

    fp.write("%%%% END SRC %s\n" % os.path.basename(src))

    fp.write("\n")
    fp.write("\\documentclass[12pt,epsf,a4]{article}\n"
             "\\usepackage{tabularx}\n"
             "\\usepackage{adjustbox}\n"
             "\\usepackage{graphicx}\n"
             "\\usepackage[margin=0.75in]{geometry}\n"
             "\\begin{document}\n"
             "\\title{%s}\n"
             "\\date{\\today}\n"
             "\\author{%s}\n"
             "\\maketitle\n\n" % (config["title"], config["author"]))
             
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
             "\n}\n\n" % (n, neq, lsode_rtol, lsode_atol,
                          simulation["MAXIMUM_STEP_SIZE"]))
    
    sx = 0
    for r in rctns:
        if len(r.text) > sx:
            sx = len(r.text)

    # Elementary reaction table
    fp.write("\\begin{table}\n"
             "\\begin{adjustbox}{width=\\textwidth}\n"
             "\\begin{tabular}{lrclll}\n")

    fp.write(" & \multicolumn{3}{c}{Reactions} & $k_{forward}$ & $k_{reverse}$\\\\\n\\hline\n")
    sx += 3
    fmt = "  %%-%ds (R%%d)\n" % sx
    nr = 0
    ne = 0
    nekinet = []

    for r in rctns:
        nr += 1
        fp.write("%%")
        fp.write(fmt % (r.text, r.i))
        fp.write("(R%d) & " % r.i)
        fp.write("%s & " % latex_reaction_fmt(r.reactants))

        if [1, 0] == r.rates:
            fp.write("%s & " % chemuni(r.i))
        else:
            fp.write("%s & " % chembi(r.i))
        
        fp.write("%s " % latex_reaction_fmt(r.products))

        # Forward reaction rates (or exceptions)
        fp.write(" & ")
        if r.rates[0] and not(r.frate):
            fp.write("%s" % latex_rc(r.kf, r.reactants))
        elif r.rates[0] and r.frate:
            ne += 1
            nekinet.append("$v_{%d} = %s$" % (r.i, latex_vrate(r.frate)))
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

    if len(constant):
        c = []
        for a in constant:
            c.append("%s = %s" % (latex_sub(a, False), constant[a]))
        fp.write("Constants: %s." % ", ".join(c))
        
    fp.write("}\n"
             "\\label{tab1e:reactions}\n"
             "\\end{table}\n\n")

    fp.write("\\begin{eqnarray}\n")

    i = 0
    for dx in xdot:
        dx = dx_sorted(dx)
        dx = latex_sub2(dx)
        dx = eqnarray_rhs(dx, x[i])
        fp.write("%% /* %s */\n" % x[i])
        fp.write("\\frac{d [\\mathrm{%s}]}{dt} & = & %s" %
                 (latex_sub(x[i]), dx))
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
    if len(config["keywords"]):
        fp.write(", ".join(config["keywords"]))
    fp.write("\n\n")

    fp.write("\\begin{thebibliography}{99}\n"
             "\\bibitem{lsoda} A.C. Hindmarsh, {\\em ODEPACK, A Systematized Collection of ODE Solvers}, in {\\em Scientific Computing}, R.S. Stepleman et al. (Eds.), North--Holland, Amsterdam, {\\bf 1983}, pp. 55-64.\n")
    for i in config["bibitem"]:
        fp.write("\\bibitem{%s} %s\n"
                 % (i.split()[0], " ".join(i.split()[1:])))
    fp.write("\\end{thebibliography}\n\n")

    fp.write("\\end{document}\n")
    
    i = 0
    fp.write("\n")
    for v in x:
        fp.write("%%")
        fp.write("  x_%d = %s\n" % (1 + i, v))
        i += 1

    fp.write("\n")

    if len(excess):
        fp.write("%%%% Species in excess: %s\n" % " ".join(excess))
        for a in excess:
            if initial.has_key(a):
                fp.write("%%%% #define %s (%s)\n" % (a, initial[a]))
        fp.write("\n")

    if False:
        fp.write("\n#define NEQ %d\n" % neq)
        fp.write("#define N_REACTIONS %d\n" % n)
        fp.write("#define T_END %s\n" % simulation["T_END"])
        fp.write("#define T_DELTA (1/ (double) %s)\n" % simulation["T_POINTS"])
        fp.write("#define LSODE_ATOL %s\n" % lsode_atol)
        fp.write("#define LSODE_RTOL %s\n" % lsode_rtol)
        if "MAXIMUM_STEP_SIZE" in simulation:
            fp.write("#define LSODE_HMAX %s\n" %
                     simulation["MAXIMUM_STEP_SIZE"])
        if "INITIAL_STEP_SIZE" in simulation:
            fp.write("#define LSODE_H0 %s\n" %
                     simulation["INITIAL_STEP_SIZE"])

        fp.write("\nstatic double kf[N_REACTIONS+1], kr[N_REACTIONS+1];\n")
        
        fp.write("\nstatic void erhelper_init(double *x, double *abstol, double *reltol)\n{\n")

        i = 0
        while i <= neq:
            fp.write("\tabstol[%d] = LSODE_ATOL;\n" % i)
            fp.write("\treltol[%d] = LSODE_RTOL;\n" % i)
            fp.write("\tx0(%d) = 0;\n" % i)
            i += 1
            
        fp.write("\n\t/* initial conditions */\n")
        for a in x:
            if a in initial:
                i = 1 + x.index(a)
                fp.write("\t/* %s */\n\tx0(%d) = %s ;\n" % (a, i, initial[a]))
                
        fp.write("\n}\n")
        
        fp.write("\nstatic void fex(double t, double *x, double *xdot, void *data)\n{\n")

        for i in kinet_keys:
            fp.write("\tdouble %s = %s ;\n" % (i, kinet[i]))
            fp.write("\n")

        # FIXME zeros
        i = 0
        for dx in xdot_raw:
            fp.write("\t/* %s */\n" % x[i])
            fp.write("\txdot(%d) = %s ; \n" % (i + 1, xdot_raw[i]))
            i += 1

        fp.write("\n}\n")

        fp.write("\n#endif\n")

    fp.close()

    try:
        os.chmod(fname, 0644)
    except:
        pass


def lsoda_c_output(fbase):
    fname = "%s.h" % fbase
    mname = "%s.mat" % fbase
    n = len(rctns)
    neq = len(xdot_raw)
    
    try:
        fp = open(fname, "w")

        fp.write("#ifndef __ERHELPER_H__\n#define __ERHELPER_H__\n\n/*\n")
        
        sx = 0
        for r in rctns:
            if len(r.text) > sx:
                sx = len(r.text)

        sx += 3
        fmt = "  %%-%ds (R%%d)\n" % sx
        for r in rctns:
            fp.write(fmt % (r.text, r.i))

        i = 0
        fp.write("\n")
        for v in x:
            fp.write("  x(%d) %s\n" % (1 + i, v))
            i += 1

        fp.write("\n  Plot command:\n\n xplot.sh -t '%s' FILE"
                 % " ".join(x))
            
        fp.write("\n*/\n")

        defs = ("\n#define x(i) (x[i-1])\n"
                "#define x0(i) (x[i])\n"
                "#define xdot(i) (xdot[i-1])\n"
                "#define kf(i) kf[i]\n"
                "#define kr(i) kr[i]\n")
        
        fp.write("%s" % defs)

        if len(excess):
            fp.write("\n/* Species in excess %s */\n" % " ".join(excess))
            for a in excess:
                if initial.has_key(a):
                    fp.write("#define %s (%s)\n" % (a, initial[a]))

        if len(constant):
            fp.write("\n/* Constants %s */\n" % " ".join(constant))
            for a in constant:
                fp.write("#define %s (%s)\n" % (a, constant[a]))

        fp.write("\n#define NEQ %d\n" % neq)
        fp.write("#define N_REACTIONS %d\n" % n)
        fp.write("#define T_END %s\n" % simulation["T_END"])
        fp.write("#define T_DELTA (1/ (double) %s)\n" % simulation["T_POINTS"])
        fp.write("#define LSODE_ATOL %s\n" % lsode_atol)
        fp.write("#define LSODE_RTOL %s\n" % lsode_rtol)
        if "MAXIMUM_STEP_SIZE" in simulation:
            fp.write("#define LSODE_HMAX %s\n" %
                     simulation["MAXIMUM_STEP_SIZE"])
        if "INITIAL_STEP_SIZE" in simulation:
            fp.write("#define LSODE_H0 %s\n" %
                     simulation["INITIAL_STEP_SIZE"])

        fp.write("\nstatic double kf[N_REACTIONS+1], kr[N_REACTIONS+1];\n")
        
        fp.write("\nstatic void erhelper_init(double *x, double *abstol, double *reltol)\n{\n")

        i = 0
        while i <= neq:
            fp.write("\tabstol[%d] = LSODE_ATOL;\n" % i)
            fp.write("\treltol[%d] = LSODE_RTOL;\n" % i)
            fp.write("\tx0(%d) = 0;\n" % i)
            i += 1
        
        fp.write("\n\t/* initial conditions */\n")
        for a in x:
            if a in initial:
                i = 1 + x.index(a)
                fp.write("\t/* %s */\n\tx0(%d) = %s ;\n" % (a, i, initial[a]))

        fp.write("\n\t/* forward reaction rates */\n")
        for r in rctns:
            if r.rates[0]:
                fp.write("\tkf(%d) = %s ;\n" % (r.i, r.kf))
        fp.write("\n\t/* reverse reaction rates */\n")
        for r in rctns:
            if r.rates[1]:
                fp.write("\tkr(%d) = %s ;\n" % (r.i, r.kr))
                
        fp.write("\n}\n")
                
        fp.write("\nstatic void fex(double t, double *x, double *xdot, void *data)\n{\n")

        for i in kinet_keys:
            fp.write("\tdouble %s = %s ;\n" % (i, kinet[i]))
        fp.write("\n")

        # FIXME zeros
        i = 0
        for dx in xdot_raw:
            fp.write("\t/* %s */\n" % x[i])
            fp.write("\txdot(%d) = %s ; \n" % (i + 1, xdot_raw[i]))
            i += 1

        fp.write("\n}\n")

        fp.write("\n#endif\n")

    except:
        raise
    finally:
        fp.close()

    try:
        os.chmod(fname, 0644)
    except:
        pass


def octave_output(fbase):
    fname = "%s.m" % fbase
    mname = "%s.mat" % fbase
    n = len(rctns)

    try:
        fp = open(fname, "w")

        fp.write("#!/usr/bin/octave -qf\n")
        fp.write("# -*-octave-*-\n")

        sx = 0
        for r in rctns:
            if len(r.text) > sx:
                sx = len(r.text)

        sx += 3
        fmt = "# %%-%ds (R%%d)\n" % sx
        for r in rctns:
            fp.write(fmt % (r.text, r.i))

        fp.write("more off\n")
        fp.write("global kf kr ;\n")

        if len(excess):
            fp.write("\n# Species in excess\nglobal %s ;\n" % " ".join(excess))
            for a in excess:
                if initial.has_key(a):
                    fp.write("%s = %s ;\n" % (a, initial[a]))

        if len(constant):
            fp.write("\n# Constants\nglobal %s ;\n" % " ".join(constant))
            for a in constant:
                fp.write("%s = %s ;\n" % (a, constant[a]))
            
        fp.write("\n# forward reaction rates\nkf = zeros(%d, 1) ;\n" % n)
        for r in rctns:
            if r.rates[0]:
                fp.write("kf(%d) = %s ;\n" % (r.i, r.kf))

        fp.write("\n# reverse reaction rates\nkr = zeros(%d, 1) ;\n" % n)
        for r in rctns:
            if r.rates[1]:
                fp.write("kr(%d) = %s ;\n" % (r.i, r.kr))

        fp.write("\n# initial conditions\nx0 = zeros(%d, 1) ;\n" % len(xdot))
        for a in x:
            if a in initial:
                i = 1 + x.index(a)
                fp.write("# %s\nx0(%d) = %s ;\n" % (a, i, initial[a]))
        fp.write("\nfunction xdot = f (x, t)\n")
        fp.write("    global kf kr ;\n")

        if len(excess):
            fp.write("    global %s ;\n" % " ".join(excess))
        if len(constant):
            fp.write("    global %s ;\n" % " ".join(constant))

        fp.write("    xdot = zeros(%d, 1) ;\n" % len(xdot))
        i = 0
        for dx in xdot:
            fp.write("    xdot(%d) = %s ; \n" % (i + 1, xdot[i]))
            i += 1
        
        fp.write("endfunction\n")

        fp.write("\nlsode_options(\"integration method\", \"stiff\") ;\n")

        if "MAXIMUM_STEP_SIZE" in simulation:
            fp.write("lsode_options(\"maximum step size\", %s) ;\n" %
                     simulation["MAXIMUM_STEP_SIZE"])

        fp.write("lsode_options(\"absolute tolerance\", %s) ;\n" % lsode_atol)
        fp.write("lsode_options(\"relative tolerance\", %s) ;\n" % lsode_rtol)
        fp.write("\nt_end = %s ;\n" % simulation["T_END"])
        fp.write("t_points = %s ;\n" % simulation["T_POINTS"])
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
        os.chmod(fname, 0755)
    except:
        pass

def validate_input():
    for a in excess:
        if a not in initial:
            raise Exception("Missing initial value for: %s" % a)

    for p in ["T_END", "T_POINTS"]:
        if p not in simulation:
            raise Exception("Missing SIMULATIONS parameter: %s" % p)


def usage():
    if len(sys.argv) > 0:
        me = sys.argv[0]
    else:
        me = "erhelper.py"

    print("Usage: %s [OPTIONS] FILE" % me)


def main(argv):

    try:
        opts, args = getopt.getopt(argv, "hvl",
                                   ["help", "verbose", "latex", "octave"])
    except getopt.GetoptError as err:
        usage()
        sys.exit(2)

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o in ("-v", "--verbose"):
            config["verbose"] += 1
        elif o in ("-l", "--latex"):
            if config["has_sympy"]:
                config["opt_latex"] = True
            else:
                print("ERROR: sympy is not available")
                sys.exit(1)
        elif o in ("--octave"):
            if config["has_sympy"]:
                config["opt_octave"] = True
            else:
                print("ERROR: sympy is not available")
                sys.exit(1)
        else:
            raise("Unhandled option")

    if len(args) > 0:
        fname = args[0]
    else:
        usage()
        sys.exit(1)
        
    read_r(fname)

    validate_input()

    if config["verbose"] > 1:
        for r in rctns:
            print(r.i)
            print(r.reactants)
            print(r.products)
            print(r.rates)
            print(r.kinet())
            print(r.kinet(False))
            print(r.species())
        
    proc_rspcs()

    proc_r()

    for k in kinet:
        dbg("KINET %s = %s" % (k, kinet[k]))

    i = 0
    for dx in xdot_raw:
        if config["opt_latex"] or config["opt_octave"]:
            dbg("xdot(%d) = %s ; " % (i + 1, xdot[i]))
        dbg("xdot_raw(%d) = %s ; " % (i + 1, xdot_raw[i]))
        i += 1

    lsoda_c_output("erhelper")

    if config["opt_octave"]:
        octave_output("erhelper")

    if config["opt_latex"]:
        latex_output("erhelper", fname)

    dbg("X %s" % x)
    dbg("EXCESS %s" % excess)
    dbg("INITIAL %s" % initial)
    dbg("CONSTANT %s" % constant)

    for ltx in latex:
        dbg("LATEX %s = %s" % (ltx, latex[ltx]))


if "__main__" == __name__:
    main(sys.argv[1:])
