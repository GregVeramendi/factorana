#!/bin/bash
# T10: the legacy C++ simulation (stage 4, TMinLkhd::Simulate) against simulate_factor_model(layout = "legacy").
# Mini model with T = 3 types, N = 2000 base rows, NSIM simulated units in each code, parameters = the legacy typed
# estimates of T9 (results/legacy_s2t3A_2000/fullmodel_par.txt). The C++ problem file is generated from
# prob_hjv_mini.cc: every model gets detail output, the endogenous regressors are declared as in bin/prob_hjv.cc
# (advmath, adveng and the HS_track dummies), and advmath, adveng and HS_track are simulated first.
# Usage: run_t10.sh [NSIM]. Needs the T9 results and data (run_t9.sh). Serial, about 1 minute at NSIM = 100000.
cd "$(dirname "$0")"; N=2000; T=3; NSIM="${1:-100000}"
python3 - <<'EOF' || exit 1
import io
lines = io.open("prob_hjv_mini.cc").read().split("\n")
endo = {"g_swe10": ["HS_track"], "g_sports10": ["HS_track"], "g_gpahs": ["HS_track"],
        "HS_track": ["advmath", "adveng"], "wage1": ["HS_track"], "wage2": ["HS_track"]}
for nm in ["apply", "application", "Dswesat", "swesat", "enroll", "switch_major", "grad1", "grad2", "grad3"]:
    endo[nm] = ["advmath", "adveng", "HS_track"]
out = []; hit = set()
for ln in lines:
    for nm, regs in endo.items():
        if ('mymodel.AddModel("%s",' % nm) in ln:
            code, sep, comment = ln.partition("  //")
            ln = code.rstrip() + "".join(' mymodel.LastModel_SetEndogenousReg("%s");' % r for r in regs) + (sep + comment if sep else "")
            hit.add(nm)
    if "if (mpirank==0) mymodel.PrintModels();" in ln:
        out += ["  mymodel.AllModel_Detailsim();", "  mymodel.SetSimNobs(NSIM);",
                '  { vector<TString> so = V({"advmath","adveng","HS_track"}); mymodel.SetSimModelOrder(so); }']
    out.append(ln)
assert hit == set(endo), set(endo) - hit
io.open("prob_t10_sim.cc", "w").write("\n".join(out))
EOF
OUT=prob_t10_s4t$T legacy_build/build_prob.sh prob_t10_sim.cc serial -DSTAGE=4 -DNOBS=$N -DNTYPES=$T -DNSIM=$NSIM || exit 1
d=results/legacy_s4t${T}_$N; rm -rf $d; mkdir -p $d results/t10
cp results/legacy_s2t${T}A_$N/fullmodel_par.txt $d/ || exit 1
s=$(date +%s); ./prob_t10_s4t$T > $d/run.log 2>&1; echo "legacy simulation: exit $? secs $(( $(date +%s) - s ))"
Rscript t10_sim_compare.R $N $T $NSIM 2>&1 | grep -vE "^Test passed|^$"
