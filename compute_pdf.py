import argparse

import numpy as np
from oct2py import Oct2Py

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-x", default=100)
    parser.add_argument("--n-D", default=8)

    args = parser.parse_args()

    oc = Oct2Py()
    oc.eval("pkg load symbolic")

    n_x = int(args.n_x)
    eps_D = 0.05
    n_D = int(args.n_D)
    x = np.linspace(-3, 3, n_x, endpoint=True)
    D = np.linspace(eps_D, 0.5-eps_D, n_D, endpoint=True)

    print(f"Running with {len(x)} x between {x.min()} and {x.max()}")
    print(f"Running with {len(D)} D between {D.min()} and {D.max()}")
    P = np.zeros((n_x, n_D))

    for i in range(n_x):
        for j in range(n_D):
            P[i,j] = oc.rosenblatt_dist(x[i], D[j], "pdf", 50, 5)

    np.savetxt("pdf.txt", P)
    np.savetxt("x.txt", x)
    np.savetxt("D.txt", D)

