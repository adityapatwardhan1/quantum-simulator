OPENQASM 2.0;
include "qelib1.inc";
qreg q[2]; // 2 qubits
creg c[2]; // 2 classical bits

h q[1]; // hadamard on qubit 1 (leftmost bit)
cx q[1], q[0]; // controlled not with qubit 1 as the control and qubit 0 as the target
measure q[0] -> c[0]; // measure q0
measure q[1] -> c[1]; // measure q1

// @columns [0,1,2,3]