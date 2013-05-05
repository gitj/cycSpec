import glob
import procosf
from multiprocessing import Pool

all = set(glob.glob('/lakitu/data/P2721/gpu0*/rt/rtcs_*.ar'))
fixed = set([x[:-4] for x in glob.glob('/lakitu/data/P2721/gpu0*/rt/rtcs_*.ar.fix')])
tofix = list(all.difference(fixed))

p = Pool(8)

p.map(procosf.fixCyclicDedisp,tofix)