# MSDTOOLBOX
Multi-Stage Decimation Toolbox
The Multi-Stage Decimation (MSD) toolbox consists of four packages, as follow:
1. Design
2. Implementation [not-included yet]
3. Estimation
4. Troubleshooting
The Design package is used to design and tune the overall decimation parameters, such as number of decimation stages (k), decimation factor at each stage (M), and filter length at each stage (N). However, the main advantage of the MSD is in the controllability and observability capabilities. The user can observe the
in-band noise (IBN) as well as the computation effort (RT) at each stage.
Further, the designer can choose between several filter architectures such as cascaded integrator comb (CIC), half-band (HB), and multi-band (MB) either narrow-band or wide-band.
The Implementation package is used to consider fixed point representation for VHDL implementation throughout quantization. The quantization bit-width (Q) is tuned throughout iteration process within acceptable penalty in the IBN. Further, it exports a VHDL top level model for RTL implementation.
The Estimation package is used to estimate power dissipation in decimation filter based on polyphase decimation (PPD) filter architecture only.
The Troubleshooting package is used in correlation with the VHDL models for
troubleshooting simulation problems.
P.S. A follow update is planned for bandpass decimation filters and multi-stage IIR filters.
