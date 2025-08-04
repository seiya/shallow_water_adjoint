import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def main():
    build_dir = Path(__file__).resolve().parents[1] / 'build'
    data = np.load(build_dir / 'taylor_test1_results.npz')
    eps = data['eps']
    fwd = data['fwd']
    rev = data['rev']

    plt.loglog(eps, fwd, 'o-', label='forward AD')
    plt.loglog(eps, rev, 's-', label='reverse AD')
    plt.xlabel('epsilon')
    plt.ylabel('error')
    plt.legend()
    plt.grid(True, which='both')
    plt.savefig(build_dir / 'taylor_test1_convergence.png')

if __name__ == '__main__':
    main()
