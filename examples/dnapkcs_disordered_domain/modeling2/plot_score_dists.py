import glob
import numpy as np



rex_inp = []
import json
import glob

conn = []
sspr = []
score = []
'''
"restraints": {"SEConnectivityRestraint0": 25.429256684263237, "SecondaryStructureParsimonyRestraint0": 1.7765926597632022, "SEConnectivityRestraint3": 0.0, "SEConnectivityRestraint1": 0.0, "SEConnectivityRestraint2": 0.0, "LoopPairDistanceRestraint9": 10.949064045085471, "LoopPairDistanceRestraint8": 1.5148036716404023, "LoopPairDistanceRestraint1": 0.0, "LoopPairDistanceRestraint0": 2.985286502971212, "LoopPairDistanceRestraint3": 0.0, "LoopPairDistanceRestraint2": 0.0, "LoopPairDistanceRestraint5": 0.024465020751813428, "LoopPairDistanceRestraint4": 0.024465020751813428, "LoopPairDistanceRestraint7": 2.307622516527772, "LoopPairDistanceRestraint6": 0.0, "SingletonRestraint 1": 0.0, "SingletonRestraint 0": 0.0, "SingletonRestraint 2": 0.0, "LoopPairDistanceRestraint18": 15.138280475512147, "LoopPairDistanceRestraint15": 8.22889427369555, "LoopPairDistanceRestraint14": 33.64133963180357, "LoopPairDistanceRestraint17": 10.888531286549188, "LoopPairDistanceRestraint16": 2.646961661266687, "LoopPairDistanceRestraint11": 0.21097489360872715, "LoopPairDistanceRestraint10": 25.09555182052427, "LoopPairDistanceRestraint13": 9.611394489184022, "LoopPairDistanceRestraint12": 6.296125013661004
'''

for outfile in glob.glob("enumerate_multi_*.dat"):
    print outfile
    with open(outfile) as f:
        try:
            lst = json.load(f)
        
            for m in lst['models']:
                rex_inp.append([m['model'], m['score']+(m['restraints']['SecondaryStructureParsimonyRestraint0']-1.73)*500,m['restraints']])
                conn.append(m['restraints']["SEConnectivityRestraint0"] + m['restraints']["SEConnectivityRestraint1"] + m['restraints']["SEConnectivityRestraint2"] + m['restraints']["SEConnectivityRestraint3"])
                score.append(m['score'])
                sspr.append(m['restraints']['SecondaryStructureParsimonyRestraint0'])
        except:
            pass

from matplotlib import pyplot as plt


print min(conn), max(conn)
print min(score), max(score)
print min(sspr), max(sspr)


plt.hist(conn, range(0,200,5))
plt.savefig("Connectivity.png", dpi=300)

plt.clf()

plt.hist(sspr)
plt.savefig("sspr.png", dpi=300)


sc_arr = np.array(score)
sc = sc_arr[~np.isnan(sc_arr)]

print sc

plt.clf()
plt.hist(sc, range(0,300,5))
plt.savefig("Totscore.png", dpi=300)

plt.clf()





