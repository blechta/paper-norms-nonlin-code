"""This script parses computation logs and writes out LaTeX tabular
with results to stdout.
"""
from __future__ import print_function
import glob, os, warnings


# Directory with logs
prefix = ""

# Snippets for LaTeX tabular
prolog = r"""\begin{table}[b]
\centering
\begin{tabular}{lrrrrrrr}
%
Case
& \#cells & $C_{\mathrm{cont}, \mathrm{PF}}$
& $\norm{\epsilon_\mathrm{glob}}_q$
& $\norm{\epsilon_\mathrm{loc}}_q$
& $\norm{\epsilon_\mathrm{flux}}_q$
& $\mathrm{Eff_{\eqref{eq_loc_dual_gal_1}}}$
& $\mathrm{Eff_{\eqref{eq_loc_dual_gal_impr_1}}}$
& $\mathrm{Eff_{\eqref{eq_loc_dual_gal_2}}}$
& $\mathrm{Eff_{\eqref{eq_est_flux}}}$
\\\hline\hline
%
"""
mrow_begin = r"\multirow{%s}{*}{\parbox{3cm}{\centering %s $p=%s$, $N=3$}}" + os.linesep
row = r"& %s & %s & %s & %s & %s & %s & %s & %s & %s \\" + os.linesep
mrow_end = r"\hline" + os.linesep
epilog = r"""\end{tabular}
\caption{Computed quantities of localization inequalities~\eqref{eq_loc_dual_gal_1},
         \eqref{eq_loc_dual_gal_impr_1}, and \eqref{eq_loc_dual_gal_2}
         for the chosen model problems.}
\label{tab_loc}
\end{table}"""

# Read line tagged with 'RESULT' from logs
results = {}
logs = glob.glob(os.path.join(prefix, "*.log"))
for log in logs:
    f = open(log, 'r')
    lines = f.readlines()
    l = [l for l in lines if l[:6]=="RESULT"]
    assert(len(l) in [0, 1])
    if len(l) == 0:
        warnings.warn("There is no 'RESULT' line in file '%s'!" % log)
        continue
    l = l[0]
    l = l[6:]
    l = l.split()
    key = (l[0], l[1])
    val = l[2:]
    if results.has_key(key):
        vals = results[key]
    else:
        vals = results[key] = []
    vals.append(tuple(val))

# Compile tabular code using results
output = prolog
for key in results.keys():
    vals = results[key]
    vals.sort(lambda x,y: 1 if int(x[0]) > int(y[0]) else -1)
    output += mrow_begin % (len(vals), key[0], key[1])
    for val in vals:
        output += row % val
    output += mrow_end
output += epilog

# Write out to stdout
print(output)
