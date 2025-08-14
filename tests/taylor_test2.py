import numpy as np
import subprocess
from pathlib import Path


def save_field(path, arr):
    path = Path(path)
    arr.ravel(order='F').astype(np.float64).tofile(path)


def read_cost(path):
    path = Path(path)
    with path.open() as f:
        for line in f:
            if line.startswith('MSE'):
                return float(line.split()[1])
    raise RuntimeError('MSE not found')


def geostrophic_height(nx, ny):
    pi = np.pi
    radius = 6371220.0
    g = 9.80616
    day = 86400.0
    omega = 2.0 * pi / (12.0 * day)
    h0 = 10000.0
    u0 = 20.0
    dlat = pi / ny
    lat = -pi/2 + (np.arange(ny) + 0.5) * dlat
    coeff = radius * omega * u0 / g
    hlat = h0 + coeff * np.sin(lat) ** 2
    return np.repeat(hlat[np.newaxis, :], nx, axis=0)


def main():
    build_dir = Path(__file__).resolve().parents[1] / 'build'
    exe_base = build_dir / 'shallow_water_test2.out'
    exe_fwd = build_dir / 'shallow_water_test2_forward.out'
    exe_rev = build_dir / 'shallow_water_test2_reverse.out'

    nx, ny = 128, 64
    rng = np.random.default_rng(0)
    x = geostrophic_height(nx, ny)
    d = rng.standard_normal((nx, ny))

    x_file = build_dir / 'x2.bin'
    d_file = build_dir / 'd2.bin'
    save_field(x_file, x)
    save_field(d_file, d)

    subprocess.run(['mpirun', '--allow-run-as-root', '-n', '4', str(exe_base), '0', str(x_file)],
                   check=True, cwd=build_dir)
    F0 = read_cost(build_dir / 'cost.log')
    assert F0 < 100.0 # 10 m

    eps_list = np.logspace(-1, -5, num=5)
    diffs = []
    for i, eps in enumerate(eps_list):
        x_eps = x + eps * d
        x_eps_file = build_dir / f'x2_eps_{i}.bin'
        save_field(x_eps_file, x_eps)
        subprocess.run(['mpirun', '--allow-run-as-root', '-n', '4', str(exe_base), '0', str(x_eps_file)],
                       check=True, cwd=build_dir)
        Fe = read_cost(build_dir / 'cost.log')
        diffs.append((Fe - F0) / eps)
    diffs = np.array(diffs)

    res = subprocess.run(['mpirun', '--allow-run-as-root', '-n', '4', str(exe_fwd), '0', str(x_file), str(d_file)],
                         check=True, cwd=build_dir, capture_output=True, text=True)
    mse_ad = float(res.stdout.strip().split()[0])

    res = subprocess.run(['mpirun', '--allow-run-as-root', '-n', '4', str(exe_rev), '0', str(x_file), str(d_file)],
                         check=True, cwd=build_dir, capture_output=True, text=True)
    lines = res.stdout.strip().splitlines()
    grad_dot_d = float(lines[-1].split()[0])

    fwd_err = np.abs(mse_ad - diffs)
    rev_err = np.abs(grad_dot_d - diffs)

    for eps, diff, fe, re in zip(eps_list, diffs, fwd_err, rev_err):
        print(f"eps={eps:.1e} diff={diff:.6e} fwd_err={fe:.6e} rev_err={re:.6e}")

    np.savez(build_dir / 'taylor_test2_results.npz', eps=eps_list, diff=diffs,
             fwd=fwd_err, rev=rev_err)

    assert fwd_err[2] < fwd_err[0]
    assert rev_err[2] < rev_err[0]


if __name__ == '__main__':
    main()
