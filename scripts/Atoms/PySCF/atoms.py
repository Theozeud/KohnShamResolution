from pyscf import gto, dft, scf

mol = gto.M(
    atom = 'He 0 0 0',  # Coordonnées de l'atome
    basis = 'sto-3g',   # Base de fonctions (STO-3G pour des calculs de base)
    charge = 0,         # Charge nette de la molécule
    spin = 0,           # Nombre d'électrons non appariés (0 pour un atome d'hydrogène)
    symmetry = False,    # Activer la symétrie
    #symmetry_subgroup = 'C1' #'D2h'
)

# Configurer le calcul
mf = scf.RHF(mol)

# Exécuter le calcul
energy = mf.kernel()

# Afficher l'énergie fondamentale
print(f"Energy of the ground state: {energy:.6f} Hartree")