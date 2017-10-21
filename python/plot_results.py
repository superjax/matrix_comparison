#!/usr/env/python

import matplotlib.pyplot as plt
import csv

size = []
TNT = []
GSL = []
Eigen = []

def plot_results(filename):
    csvfile = open(filename, 'rb')

    reader = csv.reader(csvfile, delimiter
    =',')

    for i, row in enumerate(reader):
        if i == 0:
            continue
        size.append(float(row[0]))
        Eigen.append(float(row[1]))
        TNT.append(float(row[2]))
        GSL.append(float(row[3]))

    plt.figure(1)
    plt.title("Symmetric Eigenvalue Decomposition")
    ax = plt.subplot(111)
    plt.plot(size, TNT, label="TNT")
    plt.plot(size, GSL, label="GSL")
    plt.plot(size, Eigen, label="Eigen")
    plt.legend()
    ax.set_yscale("log")
    plt.ylabel('$\mu$s')
    plt.xlabel('matrix size')
    plt.show()

if __name__ == '__main__':
    plot_results('../build/csv_symm_results.csv')
    plot_results('../build/csv_asymm_results.csv')
