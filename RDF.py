from scipy import spatial
from numpy import sqrt, pi, zeros

def RDF_simple(atoms, r, dr):
    Gr = zeros(len(r))
    
    Maxdist = r[-1] + dr/2.
    
    ParticleTree = spatial.KDTree(atoms)
    
    k = 0
    for CentralP in atoms:
        