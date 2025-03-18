# Circuit analyser
The task required IO of test files to create a circuit model as part of coursework for [structured programming](https://www.bath.ac.uk/catalogues/2019-2020/ee/EE20084.html). The assignment required the development of an executable that would read in circuit data and compute electrical properties at the output based on a input file. See example input below.
# Circuit Description

The following circuit consists of alternating series of inductors (L) and capacitors (C) connected in a specific sequence.

## Circuit Components

- **Inductors (L)**: 
  - All inductors have a value of \(1.000 \times 10^{-4} \, H\)
  
- **Capacitors (C)**: 
  - All capacitors have a value of \(1.000 \times 10^{-9} \, F\)

### Node Connections:

- \( n1 = 1, n2 = 2, L = 1.000 \times 10^{-4} \, H \)
- \( n1 = 2, n2 = 0, C = 1.000 \times 10^{-9} \, F \)
- \( n1 = 2, n2 = 3, L = 1.000 \times 10^{-4} \, H \)
- \( n1 = 3, n2 = 0, C = 1.000 \times 10^{-9} \, F \)
- \( n1 = 3, n2 = 4, L = 1.000 \times 10^{-4} \, H \)
- \( n1 = 4, n2 = 0, C = 1.000 \times 10^{-9} \, F \)
- \( n1 = 4, n2 = 5, L = 1.000 \times 10^{-4} \, H \)
- \( n1 = 5, n2 = 0, C = 1.000 \times 10^{-9} \, F \)
- ... (continues with similar connections until n1 = 101)

## Terms

- **Input Voltage (VT)**: 5 V
- **Source Resistance (RS)**: 316 Ω
- **Load Resistance (RL)**: 316 Ω

## Frequency Range

- **Start Frequency (Fstart)**: 10.0 Hz
- **End Frequency (Fend)**: \(1.0 \times 10^6\) Hz
- **Number of Frequencies (Nfreqs)**: 150

## Output Parameters

The following output parameters are provided by the simulation:

- **Input Voltage (Vin)**: V
- **Output Voltage (Vout)**: V
- **Input Current (Iin)**: A
- **Output Current (Iout)**: A
- **Input Power (Pin)**: W
- **Output Power (Pout)**: W
- **Output Impedance (Zout)**: Ohms
- **Input Impedance (Zin)**: Ohms
- **Voltage Gain (Av)**: 
- **Current Gain (Ai)**: 
