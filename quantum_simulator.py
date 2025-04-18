import numpy as np

"""
This takes as input an OPENQASM 2.0 file and integer number of shots.
It computes the quantum sate before measurement as a complex vector
and produces the number of counts for each possible measurement outcome.
"""

# Quantum Gates
I = np.array([[1, 0], [0, 1]])
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]])
T = np.array([[1, 0], [0, np.exp(1j / 4 * np.pi)]])
S = np.array([[1, 0], [0, 1j]])
T_dag = np.array([[1, 0], [0, np.exp(-1j / 4 * np.pi)]])
S_dag = np.array([[1, 0], [0, -1j]])
H = 1 / np.sqrt(2) * np.array([[1, 1], [1, -1]])

CNOT = np.array([[1, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 0, 1],
                 [0, 0, 1, 0]])


def parse_openqasm(file_path):
    pass

def simulate_quantum_circuit(file_path, shots):
    pass

def compute_statevector(circuit):
    pass


qregs = dict()
cregs = dict()

gate_to_unitary = {'h': H, 
                   'x': X, 
                   'y': Y, 
                   'z': Z, 
                   'i': I, 
                   's': S, 
                   'sdg': S_dag, 
                   't': T, 
                   'tdg': T_dag,
                   'cx': CNOT,
                   'ccx': None}

def interpret_lines(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    operations = []
    for line in lines:
        tokens = line.strip().split()

        # Ignore empty lines / comments
        if (len(tokens) == 0 or tokens[0].startswith('//') # comments
            or tokens[0] == 'include' # include "qelib1.inc"
            or tokens[0] == 'OPENQASM' # OPENQASM 2.0
            or tokens[0] == 'barrier'):
            continue

        # Create quantum register
        if tokens[0] == 'qreg' or tokens[0] == 'creg':
            reg_name_size = tokens[1]
            reg_name = reg_name_size.split('[')[0]
            # Gets the number in between the [ and ]
            reg_size = int(reg_name_size[reg_name_size.index('[') + 1 : reg_name_size.index(']')])
            if tokens[0] == 'qreg':
                qregs[reg_name] = [0] * reg_size
            else:
                cregs[reg_name] = [0] * reg_size

        
        # Apply gates
        if tokens[0] in gate_to_unitary:
            gate = tokens[0]
            

# # This is kind of a naive approach
# def get_one_qubit_gate(gate_name, qubit_num, total_num_qubits):
#     unitary = 1
#     for i in range(total_num_qubits):
#         if i == qubit_num:
#             unitary = np.kron(unitary, gate_to_unitary[gate_name])
#         else:
#             unitary = np.kron(unitary, I)
#     return unitary

def get_one_qubit_gate(gate_name, qubit_num, num_qubits):
    unitary = 1
    for i in reversed(range(num_qubits)):
        # Apply the identity for each qubit except the one we are operating on
        if i == qubit_num:
            unitary = np.kron(unitary, gate_to_unitary[gate_name])
        else:
            unitary = np.kron(unitary, I)

def get_cnot(num_ctrl, num_target, num_qubits):
    unitary1 = 1
    unitary2 = 1
    # Projector based decomposition of control gates
    for i in reversed(range(num_qubits)):
        if i == num_target:
            unitary1 = np.kron(unitary1, np.array([1, 0], [0, 0]))
            unitary2 = np.kron(unitary2, np.array([0, 0], [0, 1]))
        else:
            unitary1 = np.kron(unitary1, I)
            if i == num_ctrl:
                unitary2 = np.kron(unitary2, X)
            else:
                unitary2 = np.kron(unitary2, I)

    return unitary1 + unitary2
     
