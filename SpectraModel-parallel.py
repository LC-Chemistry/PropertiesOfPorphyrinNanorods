import argparse
import math
import multiprocessing
import numpy as np
import scipy as sp

# ---------------------------------------------------------------------------
# Number of CPU cores to use for parallel absorption spectrum calculation.
# Set to None to use all available cores (multiprocessing.cpu_count()).
# Override at the command line with --Workers <int>.
# ---------------------------------------------------------------------------
NUM_WORKERS = None


def compute_cabs_monomer(args):
    """Compute absorption cross section for a single wavelength (monomer case)."""
    v, vmon, GammaSubm, constant = args
    V = round(1 / (v / (10**7)))
    ImAlpha = constant * ((V * vmon * GammaSubm) / ((vmon**2 - V**2)**2 + (V**2) * (GammaSubm**2)))
    return 2 * math.pi * V * ImAlpha


def compute_cabs_aggregate(args):
    """Compute absorption cross section for a single wavelength (aggregate case)."""
    v, vmon, VONE, GammaSubm, constant, umon, uOne = args
    V = round(1 / (v / (10**7)))
    ImAlpha = constant * (
        ((V * vmon * GammaSubm) / ((vmon**2 - V**2)**2 + (V**2) * (GammaSubm**2)) * umon**2) +
        ((V * VONE * GammaSubm) / ((VONE**2 - V**2)**2 + (V**2) * (GammaSubm**2)) * uOne**2)
    )
    return 2 * math.pi * V * ImAlpha


def main():
    # 1. setup the parser
    parser = argparse.ArgumentParser(description="Run a simulation on SLURM")

    # 2. add the arguments you want
    parser.add_argument("--N", type=int, default=24, help="The number of monomers in the aggregate")
    parser.add_argument("--Shape", type=str, help="Aggregate type")
    parser.add_argument("--M", type=int, default=6, help="The number of per ring of the cylindrical aggregate")
    parser.add_argument("--Umon", type=float, default=3.0, help="The transition dipole of the monemer in Debye")
    parser.add_argument("--Vmon", type=float, default=2.3*(10**4), help="Monomer transition energy in per cm")
    parser.add_argument("--Beta", type=float, default=-474.3, help="Beta potential energy in per cm")
    parser.add_argument("--Gamma", type=float, default=-670.8, help="Gamma potenital energy in per cm")
    parser.add_argument("--Zeta", type=float, default=-474.3, help="Zeta potential energy in per cm")
    parser.add_argument("--Xi", type=float, default=-474.3, help="Xi potential energy in per cm")
    parser.add_argument("--Workers", type=int, default=None,
                        help="Number of CPU cores to use (default: value of NUM_WORKERS, or all available cores)")

    # 3. parse them
    args = parser.parse_args()

    # 4. use them in your logic
    n = args.N
    shape = args.Shape
    m = args.M
    umon = args.Umon
    vmon = args.Vmon
    beta = args.Beta
    gamma = args.Gamma
    zeta = args.Zeta
    xi = args.Xi

    # Resolve number of workers: CLI flag > module-level constant > all cores
    num_workers = args.Workers or NUM_WORKERS or multiprocessing.cpu_count()
    print(f"Using {num_workers} worker(s) for parallel spectrum calculation.")

    print("There are ", n, " monomers in the aggregate.")

    #lists to store row indices, column indices, and values of non-zero elements
    row_indices = []
    col_indices = []
    values = []

    #LINEAR AGGREGATE
    if shape == "linear":

        #linear aggregate
        i = 0
        while i < n:
            row_indices.append(i)
            col_indices.append(i)
            values.append(vmon)

            if i - 1 >= 0:
                row_indices.append(i)
                col_indices.append(i - 1)
                values.append(beta)

            if i + 1 < n:
                row_indices.append(i)
                col_indices.append(i + 1)
                values.append(beta)

            i += 1

    #RING AGGREGATE
    if shape == "ring":

        #ring aggregate
        i = 0
        while i < n:
            row_indices.append(i)
            col_indices.append(i)
            values.append(vmon)

            if i - 1 >= 0:
                row_indices.append(i)
                col_indices.append(i - 1)
                values.append(beta)

            if i + 1 < n:
                row_indices.append(i)
                col_indices.append(i + 1)
                values.append(beta)

            i += 1

        row_indices.append(0)
        col_indices.append(n - 1)
        values.append(beta)

        row_indices.append(n - 1)
        col_indices.append(0)
        values.append(beta)

    #DOUBLE-STRANDED AGGREGATE
    if shape == "double-stranded":

        #lists to store row indices, column indices, and values of non-zero elements
        row_indices = []
        col_indices = []
        values = []

        #double-stranded aggregate
        i = 0
        while i < n:
            row_indices.append(i)
            col_indices.append(i)
            values.append(vmon)

            if i > 0:
                row_indices.append(i)
                col_indices.append(i - 1)
                values.append(gamma)

            if i - 1 > 0:
                row_indices.append(i)
                col_indices.append(i - 2)
                values.append(beta)

            if i + 1 < n:
                row_indices.append(i)
                col_indices.append(i + 1)
                values.append(gamma)

            if i + 2 < n:
                row_indices.append(i)
                col_indices.append(i + 2)
                values.append(beta)

            i += 1

    #CYLINDRICAL AGGREGATE VERTICAL-DIPOLE
    if shape == "cylindrical vertical-dipole":

        print("There are ", m, " monomers per cross-sectional ring in the aggregate")

        l = math.ceil(n / m)  #number of rings in cylinder

        #cylinder aggregate - verticle dipole
        i = 0
        k = 0

        while i < n:
            row_indices.append(i)
            col_indices.append(i)
            values.append(vmon)

            if i + 1 + m <= n:
                row_indices.append(i)
                col_indices.append(i + m)
                values.append(beta)
            if i - m >= 0:
                row_indices.append(i)
                col_indices.append(i - m)
                values.append(beta)

        #adds all gamma values
            if i + 1 + 2 * m <= n:
                row_indices.append(i)
                col_indices.append(i + 2 * m)
                values.append(gamma)
            if i - 2 * m >= 0:
                row_indices.append(i)
                col_indices.append(i - 2 * m)
                values.append(gamma)

            #now i need to make the code for the hellish diagnol inputs
            # i need to find a way to determine if the ring the the monomer is on is even or odd

            k = math.ceil((i + 1) / m)

            #even rows
            if k % 2 == 0:
                #these are the 'wrap arounds'
                if (i + 1) % m == 1:
                    row_indices.append(i)
                    col_indices.append(i - 1)
                    values.append(beta)

                    if i + (2 * m) <= n:
                        row_indices.append(i)
                        col_indices.append(i + (2 * m) - 1)
                        values.append(beta)

                #these are the 'standard'
                else:
                    if i - 1 - m >= 0:
                        row_indices.append(i)
                        col_indices.append(i - 1 - m)
                        values.append(beta)

                    if i + m <= n:
                        row_indices.append(i)
                        col_indices.append(i - 1 + m)
                        values.append(beta)

            #odd rows
            else:
                #these are the 'wrap arounds'
                if (i + 1) % m == 0:
                    row_indices.append(i)
                    col_indices.append(i + 1)
                    values.append(beta)

                    if i - m > 0:
                        row_indices.append(i)
                        col_indices.append(i + 1 - (2 * m))
                        values.append(beta)

                # these are the 'standard'
                else:
                    if i + m + 1 < n:
                        row_indices.append(i)
                        col_indices.append(i + m + 1)
                        values.append(beta)

                    if i - m + 1 >= 0:
                        row_indices.append(i)
                        col_indices.append(i - m + 1)
                        values.append(beta)

        i += 1

    #CYLINDRICAL AGGREGATE - HORIZONTAL DIPOLE
    if shape == "cylindrical horizontal-dipole":

        print("There are ", m, " monomers per cross-sectional ring in the aggregate")

        l = math.ceil(n / m)  #number of rings in cylinder

        #cylinder aggregate - verticle dipole
        i = 0

        while i < n:
            row_indices.append(i)
            col_indices.append(i)
            values.append(vmon)

            if i + 1 + m < n:
                row_indices.append(i)
                col_indices.append(i + m)
                values.append(beta)
            if i - m >= 0:
                row_indices.append(i)
                col_indices.append(i - m)
                values.append(beta)

        #adds all gamma values
            if i + 1 + 2 * m <= n:
                row_indices.append(i)
                col_indices.append(i + 2 * m)
                values.append(gamma)
            if i - 2 * m >= 0:
                row_indices.append(i)
                col_indices.append(i - 2 * m)
                values.append(gamma)

            #now i need to make the code for the hellish diagnol inputs
            # i need to find a way to determine if the ring the the monomer is on is even or odd

            k = math.ceil((i + 1) / m)

            #even rows
            if k % 2 == 0:
                #these are the 'wrap arounds'
                if (i + 1) % m == 1:
                    row_indices.append(i)
                    col_indices.append(i - 1)
                    values.append(beta)

                    if i + (2 * m) <= n:
                        row_indices.append(i)
                        col_indices.append(i + (2 * m) - 1)
                        values.append(beta)

                #these are the 'standard'
                else:
                    if i - 1 - m >= 0:
                        row_indices.append(i)
                        col_indices.append(i - 1 - m)
                        values.append(beta)

                    if i + m <= n:
                        row_indices.append(i)
                        col_indices.append(i - 1 + m)
                        values.append(beta)

                if (i + 1) % m == 1:
                    row_indices.append(i)
                    col_indices.append(i - 1 + m)
                    values.append(zeta)

                    row_indices.append(i)
                    col_indices.append(i + 1)
                    values.append(zeta)

                elif (i + 1) % m == 0:
                    row_indices.append(i)
                    col_indices.append(i + 1 - m)
                    values.append(zeta)

                    row_indices.append(i)
                    col_indices.append(i - 1)
                    values.append(zeta)

                else:
                    row_indices.append(i)
                    col_indices.append(i + 1)
                    values.append(zeta)

                    row_indices.append(i)
                    col_indices.append(i - 1)
                    values.append(zeta)

            #odd rows
            else:
                #these are the 'wrap arounds'
                if (i + 1) % m == 0:
                    row_indices.append(i)
                    col_indices.append(i + 1)
                    values.append(beta)

                    if i - m > 0:
                        row_indices.append(i)
                        col_indices.append(i + 1 - (2 * m))
                        values.append(beta)

                # these are the 'standard'
                else:
                    if i + m + 1 < n:
                        row_indices.append(i)
                        col_indices.append(i + m + 1)
                        values.append(beta)

                    if i - m + 1 >= 0:
                        row_indices.append(i)
                        col_indices.append(i - m + 1)
                        values.append(beta)

                if (i + 1) % m == 1:
                    row_indices.append(i)
                    col_indices.append(i - 1 + m)
                    values.append(xi)

                    row_indices.append(i)
                    col_indices.append(i + 1)
                    values.append(xi)

                elif (i + 1) % m == 0:
                    row_indices.append(i)
                    col_indices.append(i + 1 - m)
                    values.append(xi)

                    row_indices.append(i)
                    col_indices.append(i - 1)
                    values.append(xi)

                else:
                    row_indices.append(i)
                    col_indices.append(i + 1)
                    values.append(xi)

                    row_indices.append(i)
                    col_indices.append(i - 1)
                    values.append(xi)

        i += 1

    #create the sparse matrix in the coo format
    #the shape is passed as a tuple (rows, cols)
    sparse_matrix = sp.sparse.coo_matrix((values, (row_indices, col_indices)), shape=(n, n))

    x = np.array(sparse_matrix.toarray())
    eigenvalues, eigenvectors = sp.linalg.eig(sparse_matrix.toarray())  #The i-th eigenvalue corresponds to the i-th eigenvector

    print("Here is the Hamiltionian: ", x)
    print("Eiganvalues: ", eigenvalues)
    print("Eiganvectors: ", eigenvectors)

    # ground state to 1st excited state transition
    vOne = min(eigenvalues)
    position = np.where(eigenvalues == vOne)[0][0]
    Cj = sum(eigenvectors[position])
    uOne = umon * Cj
    VONE = vOne.real
    transitionEnergy = (1 / VONE) * (10**7)

    print("The transition energy from the ground to first excited state is ", transitionEnergy, " nm.")

    h = 6.62607015 * 10**(-34)       # J s   -- Planck's constant
    c = 299792458                     # m/s   -- speed of light
    epsilonNot = 8.854 * 10**(-12)   # F/m   -- vacuum permittivity
    GammaSubm = 1 / (1 / (10**7))    # cm-1  -- spectral linewidth (~1 nm)
    constant = 1 / (2 * math.pi * h * c * epsilonNot)

    v_start = round((1 / vmon) * (10**7))  # units of nm
    wavelengths = list(range(v_start, 581))

    # --- Parallelized absorption spectrum calculation ---
    # Each wavelength is independent, so we distribute them across CPU cores.
    if shape == "monomer":
        task_args = [(v, vmon, GammaSubm, constant) for v in wavelengths]
        with multiprocessing.Pool(processes=num_workers) as pool:
            cabs = pool.map(compute_cabs_monomer, task_args)
    else:
        task_args = [(v, vmon, VONE, GammaSubm, constant, umon, uOne) for v in wavelengths]
        with multiprocessing.Pool(processes=num_workers) as pool:
            cabs = pool.map(compute_cabs_aggregate, task_args)

    Cabs = min(cabs)
    print("The absorbance cross section is ", Cabs, " a.u.")


#calls main
if __name__ == "__main__":
    main()
