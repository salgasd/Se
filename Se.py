import numpy as np
import matplotlib.pyplot as plt


class Se:
    def __init__(self, a, b, c, alpha, beta, gamma, file=False):
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.centers = np.array([0.0, 0.0, 0.0])
        self.atoms = None
        self.metr = np.array([[a*a, a*b*np.cos(gamma),a*c*np.cos(beta)], 
                                     [a*b*np.cos(gamma), b*b, b*c*np.cos(alpha)], 
                                     [a*c*np.cos(beta), b*c*np.cos(alpha), c*c]])
        if file:
            self.file = np.loadtxt('Se6car3.txt', dtype='float')
            self.atoms = []
            for atom in self.file:
                zero = np.array([0.0, 0.0, 0.0])
                product = product = np.sqrt(np.dot(np.dot(zero - atom, self.metr), zero - atom))
                if product and product <= 12.0:
                    self.atoms.append(atom)
            self.atoms = np.array(self.atoms)
        self.dist = np.array([[0.04203, 0.2027, 0.120467],
                              [-0.04203, -0.2027, -0.120467],
                              [0.160233, -0.04203, 0.120467],
                              [-0.160233, 0.04203, -0.120467],
                              [-0.202267, -0.160233, 0.120467],
                              [0.202267,0.160233, -0.120467]
                              ])
        self.volume = a*a*np.sin(np.pi/3.0)
        self.neigh = []
        self.distances = None
    def add_atoms(self):
        self.centers = np.vstack((self.centers, [2./3., 1./3., 1./3.]))
        self.centers = np.vstack((self.centers, [1./3., 2./3., 2./3.]))
        for cen in self.centers:
            self.atoms = np.vstack(cen + self.dist)
    def pr(self):
        #print(self.centers)
        #print(self.atoms)
        print(self.neigh)
        print(self.distances)
    def calc_neighb(self):
        zero = np.array([0., 0., 0.])
        for atom in self.atoms:
            product = np.sqrt(np.dot(np.dot(zero - atom, self.metr), zero - atom))
            if product and product <= 12.0:
                self.neigh.append(product)
        unique, counts = np.unique(self.neigh, return_counts=True)
        self.distances = dict(zip(unique, counts))
    def pl_t(self):
        fig, ax = plt.subplots(1, 1, figsize = (8, 8))
        ax.scatter(self.atoms[:,0], self.atoms[:, 1])
        ax.grid()
    def calc_rdf(self, r, sigma):
        dist = []
        nneig = []
        for key, value in self.distances:
            dist.append(key)
            nneig.append(value)
        pdf = sum(nneig*np.exp(-(np.power(r-dist,2.0)/(2.0*sigma)))/(np.power(2.0*np.pi*sigma,0.5)))
        return pdf/r-4.0*np.pi*r*3*2/(a*b*c*np.sqrt(3))

a = 11.362
b = 11.362
c = 4.429
alpha = np.pi/2.0
beta = np.pi/2.0
gamma = 2.0*np.pi/3.0
sigma_inf = 0.055
gammas = np.ones(105)
gammas[:11] = [0.13, 0.70, 0.55, 0.35, 0.75, 0.85, 0.65, 0.45, 0.90, 0.75, 0.95]
sigmas= sigma_inf*gammas
se = Se(a, b, c, alpha, beta, gamma, True)
#se.add_atoms()
se.calc_neighb()
se.pr()
se.pl_t()
print(len(se.distances))
#print(len(se.atoms))
x = np.linspace(0.0001, 12.0, 500)
fig, ax = plt.subplots(1, 1, figsize = (8, 8))
ax.plot(x, se.calc_rdf(x, sigmas))