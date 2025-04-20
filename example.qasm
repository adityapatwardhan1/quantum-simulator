OPENQASM 2.0;
include "qelib1.inc";

qreg q[2];
creg c[2];

x q[0];
x q[1];
x q[0];
x q[0];

// h q[0];    // put qubit into superposition
// h q[0];    // apply Hadamard again, returns to |0⟩
measure q[0] -> c[0];
measure q[1] -> c[1];