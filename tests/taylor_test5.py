import numpy as np
import subprocess
from pathlib import Path


def save_field(path, arr):
    path = Path(path)
    arr.ravel(order='F').astype(np.float64).tofile(path)


def read_energy_cost(path):
    path = Path(path)
    with path.open() as f:
        for line in f:
            if line.startswith('EnergyResidual'):
                return float(line.split()[1])
    raise RuntimeError('EnergyResidual not found')


def init_height(nlon, nlat):
    radius = 6371220.0
    h0 = 10000.0
    h1 = 2000.0
    Lx = 2.0 * np.pi * radius
    Ly = np.pi * radius
    dx = Lx / nlon
    dy = Ly / nlat
    x = (np.arange(nlon) + 0.5) * dx
    y = (np.arange(nlat) + 0.5) * dy
    x0 = 0.5 * Lx
    y0 = 0.75 * Ly
    r0 = radius / 4.0
    b = np.zeros((nlon, nlat))
    for j in range(nlat):
        for i in range(nlon):
            dist = np.hypot(x[i] - x0, y[j] - y0)
            if dist < r0:
                b[i, j] = h1 * (1.0 - dist / r0)
    return h0 - b


def main():
    build_dir = Path(__file__).resolve().parents[1] / 'build'
    subprocess.run([
        'make',
        'shallow_water_test5.out',
        'shallow_water_test5_forward.out',
        'shallow_water_test5_reverse.out',
    ], check=True, cwd=build_dir, capture_output=True)
    exe_base = build_dir / 'shallow_water_test5.out'
    exe_fwd = build_dir / 'shallow_water_test5_forward.out'
    exe_rev = build_dir / 'shallow_water_test5_reverse.out'

    nlon, nlat = 128, 64
    rng = np.random.default_rng(0)
    x = init_height(nlon, nlat)
    d = rng.standard_normal((nlon, nlat))

    x_file = build_dir / 'x5_tay.bin'
    d_file = build_dir / 'd5_tay.bin'
    save_field(x_file, x)
    save_field(d_file, d)

    subprocess.run([str(exe_base), '0', str(x_file)], check=True, cwd=build_dir)
    F0 = read_energy_cost(build_dir / 'cost.log')

    eps_list = np.logspace(-1, -5, num=5)
    diffs = []
    for i, eps in enumerate(eps_list):
        x_eps = x + eps * d
        x_eps_file = build_dir / f'x5_eps_{i}.bin'
        save_field(x_eps_file, x_eps)
        subprocess.run([str(exe_base), '0', str(x_eps_file)], check=True, cwd=build_dir)
        Fe = read_energy_cost(build_dir / 'cost.log')
        diffs.append((Fe - F0) / eps)
    diffs = np.array(diffs)

    res = subprocess.run([str(exe_fwd), '0', str(x_file), str(d_file)],
                         check=True, cwd=build_dir, capture_output=True, text=True)
    energy_ad = float(res.stdout.strip().split()[0])

    res = subprocess.run([str(exe_rev), '0', str(x_file), str(d_file)],
                         check=True, cwd=build_dir, capture_output=True, text=True)
    lines = res.stdout.strip().splitlines()
    grad_dot_d = float(lines[-1].split()[0])

    fwd_err = np.abs(energy_ad - diffs)
    rev_err = np.abs(grad_dot_d - diffs)

    for eps, diff, fe, re in zip(eps_list, diffs, fwd_err, rev_err):
        print(f"eps={eps:.1e} diff={diff:.6e} fwd_err={fe:.6e} rev_err={re:.6e}")

    np.savez(build_dir / 'taylor_test5_results.npz', eps=eps_list, diff=diffs,
             fwd=fwd_err, rev=rev_err)

    assert fwd_err[2] < fwd_err[0]
    assert rev_err[2] < rev_err[0]


if __name__ == '__main__':
    main()
