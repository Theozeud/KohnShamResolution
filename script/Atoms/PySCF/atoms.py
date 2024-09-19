from pyscf import gto, dft

mol = gto.M(
    atom = 'He 0 0 0',  # Coordonnées de l'atome
    basis = 'sto-3g',   # Base de fonctions (STO-3G pour des calculs de base)
    charge = 0,         # Charge nette de la molécule
    spin = 0,           # Nombre d'électrons non appariés (0 pour un atome d'hydrogène)
    symmetry = True,    # Activer la symétrie
)


mf_hf = dft.RKS(mol)

mf_hf.xc = '' 
energy = mf_hf.kernel()

# Afficher l'énergie fondamentale
print(f"Energy of the ground state: {energy} Hartree")