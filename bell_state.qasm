// Create a quantum circuit for a Bell state

// Declare quantum registers and classical registers
qreg q[3];  // 2-qubit quantum register
creg c[3];  // 2-bit classical register for measurement

// Apply a Hadamard gate to the second qubit
h q[2];

// Apply a CNOT gate with q[0] as the control and q[1] as the target
cx q[2], q[0];

// Measure both qubits
measure q[0] -> c[0];
measure q[2] -> c[2];
