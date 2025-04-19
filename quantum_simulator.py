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
                qregs[reg_name] = np.array([0] * (int(2 ** reg_size)))
                qregs[reg_name][0] = 1
                qregs[reg_name] = qregs[reg_name].reshape(-1, 1)
                print("Quantum register: ", qregs[reg_name])
            else:
                cregs[reg_name] = np.array([0] * reg_size)

        if (tokens[0] == 'measure'):
            # Measure qubit
            qreg_name = tokens[1].split('[')[0]
            qubit_num = int(tokens[1].split('[')[1].split(']')[0])
            creg_name = tokens[3].split('[')[0]
            creg_num = int(tokens[3].split('[')[1].split(']')[0])

            # Measure the qubit, store result in classical register, and update statevector
            prob_measure_1 = 0
            for i in range(len(qregs[qreg_name])):
                if i & (1 << qubit_num):
                    prob_measure_1 += abs(qregs[qreg_name][i][0]) ** 2
            prob_measure_0 = 1 - prob_measure_1
            print("prob_measure_0: ", prob_measure_0)
            print("prob_measure_1: ", prob_measure_1)
            # Update classical register
            result = cregs[creg_name][creg_num] = np.random.choice([0, 1], p=[prob_measure_0, prob_measure_1])
            print("result: ", result)

            normalization_factor = sum(abs(qregs[qreg_name]) ** 2)
            # Update quantum register
            if result == 0:
                for i in range(len(qregs[qreg_name])):
                    if i & (1 << qubit_num):
                        qregs[qreg_name][i][0] = 0
                    # else:
                    #     qregs[qreg_name][i] /= np.sqrt(prob_measure_0)
            else:
                for i in range(len(qregs[qreg_name])):
                    if not (i & (1 << qubit_num)):
                        qregs[qreg_name][i][0] = 0
                    # else:
                    #     qregs[qreg_name][i] /= np.sqrt(prob_measure_1)

            normalization_factor = 0
            for i in range(len(qregs[qreg_name])):
                normalization_factor += abs(qregs[qreg_name][i][0]) ** 2
            print("normalization factor = ", normalization_factor)

            # Normalize the quantum register
            if normalization_factor != 0:
                for i in range(len(qregs[qreg_name])):
                    qregs[qreg_name][i][0] /= np.sqrt(normalization_factor)
                # qregs[qreg_name] /= np.sqrt(normalization_factor)

        # Apply gates
        if tokens[0] in gate_to_unitary:
            gate = tokens[0]
            reg_name = tokens[1].split('[')[0]

            n_qubits = int(np.log2(len(qregs[reg_name])))
            if gate == 'cx':
                control_qubit = int(tokens[1].split('[')[1].split(']')[0])
                target_qubit = int(tokens[2].split('[')[1].split(']')[0])
                cnot_gate = get_cnot(control_qubit, target_qubit, n_qubits)
                qregs[reg_name] = cnot_gate @ qregs[reg_name]
            else:
                qubit_num = int(tokens[1].split('[')[1].split(']')[0])
                unitary = get_one_qubit_gate(gate, qubit_num, n_qubits)
                print("Unitary: ", unitary)
                print("qregs[reg_name]: ", qregs[reg_name])
                qregs[reg_name] = unitary @ qregs[reg_name]
                # pass
                # CNOT gate
                # qubits = tokens[1].split(',')
                # qubit_num = []
                # for q in qubits:
                #     if q.startswith('['):
                #         q = q[q.index('[') + 1 : q.index(']')]
                #     qubit_num.append(int(q))
                # operations.append((gate, qubit_num[0], qubit_num[1]))

    print("Classical registers: ", cregs)
    print("Quantum registers: ", qregs)        

def get_one_qubit_gate(gate_name, qubit_num, num_qubits):
    unitary = 1
    for i in reversed(range(num_qubits)):
        # Apply the identity for each qubit except the one we are operating on
        if i == qubit_num:
            unitary = np.kron(unitary, gate_to_unitary[gate_name])
        else:
            unitary = np.kron(unitary, I)

    return unitary

def get_cnot(num_ctrl, num_target, num_qubits):
    print("in get_cnot")
    unitary1 = 1
    unitary2 = 1
    # Projector based decomposition of control gates
    for i in reversed(range(num_qubits)):
        if i == num_target:
            unitary1 = np.kron(unitary1, np.array([[1, 0], [0, 0]]))
            unitary2 = np.kron(unitary2, np.array([[0, 0], [0, 1]]))
        else:
            unitary1 = np.kron(unitary1, I)
            if i == num_ctrl:
                unitary2 = np.kron(unitary2, X)
            else:
                unitary2 = np.kron(unitary2, I)
    print("ans = ", unitary1 + unitary2)
    return unitary1 + unitary2
     
if __name__ == '__main__':
    # Example usage
    file_path = 'bell_state.qasm'
    # shots = 1024
    interpret_lines(file_path)
    # simulate_quantum_circuit(file_path, shots)