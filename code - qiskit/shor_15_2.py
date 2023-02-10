from tkinter import W
from qiskit import *
import numpy as np
from qiskit.visualization import plot_histogram
import qiskit.tools.jupyter
import matplotlib.pyplot as plt
import matplotlib 
import os
import pandas
import operator

sim = Aer.get_backend('aer_simulator')
shots = 1000

## Create 7mod15 gate
N = 15
m = int(np.ceil(np.log2(N)))

U_qc = QuantumCircuit(m)
U_qc.x(range(m))
U_qc.swap(1, 2)
U_qc.swap(2, 3)
U_qc.swap(0, 3)
state = U_qc.save_state
U = U_qc.to_gate()
U.name ='{}Mod{}'.format(7, N)

def cU_multi(k):
    circ = QuantumCircuit(m)
    for _ in range(2**k):
        circ.append(U, range(m))
    
    U_multi = circ.to_gate()
    U_multi.name = '7Mod15_[2^{}]'.format(k)
    
    cU_multi = U_multi.control()
    return cU_multi

def qft(n):
    """Creates an n-qubit QFT circuit"""
    circuit = QuantumCircuit(n)
    def swap_registers(circuit, n):
        for qubit in range(n//2):
            circuit.swap(qubit, n-qubit-1)
        return circuit
    def qft_rotations(circuit, n):
        """Performs qft on the first n qubits in circuit (without swaps)"""
        if n == 0:
            return circuit
        n -= 1
        circuit.h(n)
        for qubit in range(n):
            circuit.cp(np.pi/2**(n-qubit), qubit, n)
        qft_rotations(circuit, n)
    
    qft_rotations(circuit, n)
    swap_registers(circuit, n)
    return circuit

'''
# QPE circuit for Shor shortend
t = 3 
shor_QPE = QuantumCircuit(t+m, t)
shor_QPE.h(range(t))

shor_QPE.x(t)
for idx in range(t-1):
    shor_QPE.append(cU_multi(idx), [idx]+ list(range(t,t+m)))

qft_dag = qft(t).inverse()
qft_dag.name = 'QFT+'

shor_QPE.append(qft_dag, range(t))
shor_QPE.measure(range(t), range(t))

shor_QPE.draw()

shor_QPE_trans = transpile(shor_QPE, sim)
count_QPE = sim.run(shor_QPE_trans, shots=shots).result().get_counts()
key_new = [str(int(key,2)/2**3) for key in count_QPE.keys()]
count_new_QPE = dict(zip(key_new, count_QPE.values()))
plot_histogram(count_new_QPE)
'''
t = 2*m

shor_Orig = QuantumCircuit(t+m, t)
shor_Orig.h(range(t))

shor_Orig.x(t+m-1)
for idx in range(t):
    shor_Orig.append(cU_multi(idx), [idx]+ list(range(t,t+m)))

qft_dag = qft(t).inverse()
qft_dag.name = 'QFT+'

shor_Orig.append(qft_dag, range(t))
shor_Orig.measure(range(t), range(t))
    
#shor_Orig.draw('mpl', scale=0.3, fold=-1)      #disegna il circuito quantistico

shor_Orig_trans = transpile(shor_Orig, sim)
count_Orig = sim.run(shor_Orig_trans, shots=shots).result().get_counts()    #dizionario con i conteggi delle simulazioni
count_Orig_sorted = sorted(count_Orig.items(), key=operator.itemgetter(0))  #lista di tuple ordinata secondo gli output in binario (nelle keys)

key_new = [str(int(key, 2)/2**t) for key in count_Orig.keys()] #lista con le fasi uscite dal QPE (come stringhe) non ordinate
key_y = [str(int(count_Orig_sorted[i][0], 2)) for i in range(len(count_Orig_sorted))] #lista con i valori ordinati (|y>) della misura del primo registro (come stringa) 
key_numbers = [int(count_Orig_sorted[i][0], 2) for i in range(len(count_Orig_sorted))] #lista con i valori ordinati (|y>) della misura del primo registro (come integer)
vals = [ count_Orig_sorted[i][1] for i in range(len(count_Orig_sorted)) ] #conteggi ordinati delle simulazioni per ogni fase
vals_str = [ str(count_Orig_sorted[i][1]) for i in range(len(count_Orig_sorted)) ]
#matplotlib figures
plt.style.use('seaborn-whitegrid')
fig, ax = plt.subplots()
ax.bar(key_numbers, vals, width=10, edgecolor="white", linewidth=0.7)
plt.show()

#gnuplot files
file1 = open("file1.txt", "w")
for i in range(len(key_y)):
    file1.write(key_y[i] + "    " + vals_str[i] + "\n")
file1.close()

#qiskit histogram
count_new_Orig = dict(zip(key_new, count_Orig.values()))
count_Orig_y = dict(zip(key_y, count_Orig.values()))
plot_histogram( count_new_Orig, title='circuit simulation result No noise')

#-----------------------------------------WITH NOISE-------------------------------------------------
from qiskit.test.mock import FakeMelbourne
from qiskit.providers.aer import AerSimulator

backend = FakeMelbourne()
sim_Melborne = AerSimulator.from_backend(backend)

shots=1000

shorOrig_trans = transpile(shor_Orig, backend, optimization_level=3)
count_shorOrig_noise = sim_Melborne.run(shor_Orig_trans, shots=shots).result().get_counts()
key_new = [str(np.round(int(key,2)/2**t,3)) for key in count_shorOrig_noise.keys()]
count_new_Orig_noise = dict(zip(key_new, count_shorOrig_noise.values()))
fig, ax = plt.subplots(2,1, figsize=(30,13))
fig.suptitle('Simulation results for the order finding circuit of $7^{r} mod 15 = 1$', fontsize=23)
plot_histogram(count_new_Orig, ax=ax[0])
plot_histogram(count_new_Orig_noise, ax=ax[1])
ax[0].set_title('sim No noise', fontsize=16)
ax[1].set_title('sim on Melbourne', fontsize=16)
plt.show()
