// Create a quantum circuit for a Bell state

// Declare quantum registers and classical registers
qreg q[2];  // 2-qubit quantum register
creg c[2];  // 2-bit classical register for measurement

// Apply a Hadamard gate to the first qubit
h q[1];

// Apply a CNOT gate with q[0] as the control and q[1] as the target
cx q[1], q[0];

// Measure both qubits
measure q[0] -> c[0];
measure q[1] -> c[1];
