import numpy as np
import matplotlib.pyplot as plt

"""
This takes as input an OPENQASM 2.0 file and integer number of shots.
It computes the quantum state before measurement as a complex vector
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

qregs = dict()
qregs_sv = dict() # Statevector form
cregs = dict()

# p for each gate
p_x = p_y = p_z = p_t = p_s = p_tdag = p_sdag = p_h = p_cnot = 1

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

gate_to_p = dict()

def interpret_lines(file_path, noise, print_state=False):
    have_measured = False

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
            add_register(tokens)

        if tokens[0] == 'measure':
            # Print out the qubit state
            if not have_measured and print_state:
                qreg_name = tokens[1].split('[')[0]
                print(f"qregs {qreg_name} statevector before measurement: \n{qregs_sv[qreg_name]}")
                have_measured = True
            measure(tokens)

        # Apply gates
        if tokens[0] in gate_to_unitary:
            apply_quantum_gate(tokens, noise)


def apply_quantum_gate(tokens, noise):
    gate = tokens[0]
    reg_name = tokens[1].split('[')[0]

    n_qubits = int(np.log2(len(qregs[reg_name])))
    if gate == 'cx':
        control_qubit = int(tokens[1].split('[')[1].split(']')[0])
        target_qubit = int(tokens[2].split('[')[1].split(']')[0])
        cnot_gate = get_cnot(control_qubit, target_qubit, n_qubits)
        result = cnot_gate @ qregs[reg_name] @ cnot_gate.conj().T
        if noise:
            result = gate_to_p[gate] * result + (1 - gate_to_p[gate]) * 1 / 2 * get_depolarization_matrix(control_qubit, n_qubits)
        qregs[reg_name] = result

        qregs_sv[reg_name] = cnot_gate @ qregs_sv[reg_name]
    else:
        qubit_num = int(tokens[1].split('[')[1].split(']')[0])
        unitary = get_one_qubit_gate(gate, qubit_num, n_qubits)
        result = unitary @ qregs[reg_name] @ unitary.conj().T
        if noise:
            result = gate_to_p[gate] * result + (1 - gate_to_p[gate]) * 1 / 2 * get_depolarization_matrix(qubit_num, n_qubits)
        qregs[reg_name] = result

        qregs_sv[reg_name] = unitary @ qregs_sv[reg_name]

def measure(tokens):
    # Parse qreg, creg
    qreg_name = tokens[1].split('[')[0]
    qubit_num = int(tokens[1].split('[')[1].split(']')[0])
    num_qubits = int(np.log2(len(qregs[qreg_name])))
    creg_name = tokens[3].split('[')[0]
    creg_num = int(tokens[3].split('[')[1].split(']')[0])

    P0, P1 = measurement_operators(qubit_num, num_qubits)
    rho = qregs[qreg_name]

    # Compute probabilities
    prob_measure_1 = min(max(0, np.trace(P1 @ rho).real), 1)
    prob_measure_0 = 1 - prob_measure_1

    # Update classical register
    result = cregs[creg_name][creg_num] = np.random.choice([0, 1], p=[prob_measure_0, prob_measure_1])

    # Zero out remaining entries based on what the qubit collapses to
    if result == 0:
        new_rho = P0 @ rho @ P0
        qregs[qreg_name] = new_rho / prob_measure_0 if prob_measure_0 != 0 else new_rho

    else:
        new_rho = P1 @ rho @ P1
        qregs[qreg_name] = new_rho / prob_measure_1 if prob_measure_1 != 0 else new_rho


def measurement_operators(qubit_num, num_qubits):
    # Projector based decomposition of measurement operators
    P0 = np.array([[1, 0], [0, 0]])
    P1 = np.array([[0, 0], [0, 1]])
    I = np.eye(2)

    M0 = 1
    M1 = 1
    for i in reversed(range(num_qubits)):
        if i == qubit_num:
            M0 = np.kron(M0, P0)
            M1 = np.kron(M1, P1)
        else:
            M0 = np.kron(M0, I)
            M1 = np.kron(M1, I)

    return M0, M1


def add_register(tokens):
    reg_name_size = tokens[1]
    reg_name = reg_name_size.split('[')[0]
            # Gets the number in between the [ and ]
    reg_size = int(reg_name_size[reg_name_size.index('[') + 1 : reg_name_size.index(']')])
    if tokens[0] == 'qreg':
        dim = int(2 ** reg_size)
        psi = np.zeros((dim, 1), dtype=complex)
        psi[0] = 1
        rho = psi @ psi.conj().T
        qregs[reg_name] = rho
        qregs_sv[reg_name] = psi
    else:
        cregs[reg_name] = np.array([0] * reg_size)       


def get_one_qubit_gate(gate_name, qubit_num, num_qubits):
    unitary = 1
    for i in reversed(range(num_qubits)):
        # Apply the identity for each qubit except the one we are operating on
        if i == qubit_num:
            unitary = np.kron(unitary, gate_to_unitary[gate_name])
        else:
            unitary = np.kron(unitary, I)

    return unitary

def get_depolarization_matrix(qubit_num, num_qubits):
    unitary = 1
    for i in reversed(range(num_qubits)):
        # Apply the identity for each qubit except the one we are operating on
        if i == qubit_num:
            unitary = np.kron(unitary, I / 2)
        else:
            unitary = np.kron(unitary, I)

    return unitary


def get_cnot(num_ctrl, num_target, num_qubits):
    unitary1 = 1
    unitary2 = 1
    P0 = np.array([[1, 0], [0, 0]])
    P1 = np.array([[0, 0], [0, 1]])
    # Projector based decomposition of control gates
    for i in reversed(range(num_qubits)):
        if i == num_ctrl:
            unitary1 = np.kron(unitary1, P0)
            unitary2 = np.kron(unitary2, P1)
        elif i == num_target:
            unitary1 = np.kron(unitary1, I)
            unitary2 = np.kron(unitary2, X)
        else:
            unitary1 = np.kron(unitary1, I)
            unitary2 = np.kron(unitary2, I)
    return unitary1 + unitary2
     

def simulate_quantum_circuit(file_path, shots, noise):
    # Get frequency for each classical output
    measurement_outcome_counts = dict()

    for i in range(shots):
        # Reset the quantum registers
        for reg_name in qregs:
            dim = qregs[reg_name].shape[0]
            psi = np.zeros((dim, 1), dtype=complex)
            psi[0] = 1
            rho = psi @ psi.conj().T
            qregs[reg_name] = rho
            qregs_sv[reg_name] = psi

        # Reset the classical registers
        for reg_name in cregs:
            cregs[reg_name] = np.zeros(len(cregs[reg_name]))
        # Interpret the lines again to apply gates and measurements
        interpret_lines(file_path, noise, print_state=(True if i == 0 else False))
        
        measurement_outcome = tuple(tuple(cregs[reg_name]) for reg_name in cregs)
        if measurement_outcome in measurement_outcome_counts:
            measurement_outcome_counts[measurement_outcome] += 1
        else:
            measurement_outcome_counts[measurement_outcome] = 1            

    print("Measurement outcomes:")
    for outcome, count in measurement_outcome_counts.items():
        print(f"{outcome}: {count} counts")

    plot_measurement_outcomes(measurement_outcome_counts)

def plot_measurement_outcomes(counts):
    labels = []
    for outcome in counts:
        label = ''
        for bits in outcome:
            for bit in reversed(bits):
                label += str(int(bit))
        labels.append(label)
    
    frequencies = [counts[outcome] for outcome in counts]
    label_frequencies_pairs = list(zip(labels, frequencies))
    label_frequencies_pairs.sort(key = lambda p : int(p[0], 2))
    labels, frequencies = zip(*label_frequencies_pairs)

    plt.figure(figsize=(10, 5))
    plt.bar(labels, frequencies, color='blue', edgecolor='black', width = 0.3)
    plt.xlabel('Measurement outcome')
    plt.ylabel('Counts')
    plt.title('Quantum Measurement Outcomes, Shots = ' + str(sum(frequencies)))
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    # file_path = 'test_4.qasm'
    file_path = input("Enter QASM file path: ")
    shots = int(input("Enter number of shots: "))
    noise = input("Model noisy gates? (Y/N): ").upper() == "Y"
    print("noise inputted: ", noise)
    if noise:
        print("Do you want a single p for all gates? (Y/N): ")
        single_p = input().upper() == "Y"
        if single_p:
            print("Set p: ")
            p = float(input())
            p_x = p_y = p_z = p_t = p_s = p_tdag = p_sdag = p_h = p_cnot = p

            gate_to_p = {'h': p_h,
              'x': p_x,
              'y': p_y,
              'z': p_z,
              'i': 1,
              's': p_s,
              'sdg': p_sdag,
              't': p_t,
              'tdg': p_tdag,
              'cx': p_cnot,
              'ccx': None}
        else:
            print("Set p for each gate: ")
            print("p for X: ")
            p_x = float(input())
            print("p for Y:")
            p_y = float(input())
            print("p for Z: ")
            p_z = float(input())
            print("p for T: ")
            p_t = float(input())
            print("p for S: ")
            p_s = float(input())
            print("p for T_dag: ")
            p_tdag = float(input())
            print("p for S_dag: ")
            p_sdag = float(input())
            print("p for H: ")
            p_h = float(input())
            print("p for CNOT: ")
            p_cnot = float(input())

            gate_to_p = {'h': p_h,
              'x': p_x,
              'y': p_y,
              'z': p_z,
              'i': 1,
              's': p_s,
              'sdg': p_sdag,
              't': p_t,
              'tdg': p_tdag,
              'cx': p_cnot,
              'ccx': None}

    interpret_lines(file_path, noise)
    simulate_quantum_circuit(file_path, shots, noise)