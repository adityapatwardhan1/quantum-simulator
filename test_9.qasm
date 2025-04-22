OPENQASM 2.0;
include "qelib1.inc";

qreg q[2];
creg c[2];

// Step 1: Put qubit 0 in superposition
h q[0];          // (|0⟩ + |1⟩)/√2

// Step 2: Add global phase to |1⟩ with T gate
// Step 3: Add complex phase to qubit 1
y q[1];          // Y has off-diagonal ±i
t q[0];          // Now amplitudes are 1 and e^{iπ/4}

// Step 4: Entangle them
cx q[0], q[1];   // Entanglement with complex interference
s q[1];          // Multiply |1⟩ by i


// Measurement
measure q[0] -> c[0];
measure q[1] -> c[1];
