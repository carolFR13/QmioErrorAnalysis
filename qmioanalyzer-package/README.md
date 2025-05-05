### **Analysis of qmio measurements.**

To analyze the measurements taken with the qmio quantum computer a set of functions and classes have been declared in ```qmio_utils.py```. 

Let's explain a bit the ```ErrorAnalyzer``` class defined in this module. The idea of this class it to compute the number of errors from a given array of measurements. The measurements would consist in strings of bits with 0 and 1 according to the measured eigenvalue for each qubit. 

Two scenarios were considered here in order to compute this errors:
1. **Repetition_period $\geq 500 \ \mu s$**
   
In this case, the states after executing the corresponding circuit are expected to go back to $\ket{0}$ state because this time is comparable to $3-5 T_1$ order of magnitude. Considering this, the errors in each measurement are defined only by considering the instantiation of the states:
- If an X-gate is applied : the states are initialized to $\ket{1} \ \rightarrow$ a measurement of $0$ is recorded as an error.
- If nothing is applied : the states are initialized to $\ket{0} \ \rightarrow$ a measurement of $1$ is recorded as an error.

1. **Repetition_period $< 500 \ \mu s$**

The overall circuit can be understood as multiple X-gate applied consecutively without enough time between them for the states to go back to $\ket{0}$ state. In this case, the number of errors considered in each measurement has to be defined from the previous measurement, by simulating the application of an X-gate and comparing the simulated result to the measured one.
As an example, we can consider that a measurement in the iteration $i$ gives $1101$, for $i+1$ we expect $0010$ (which is simply the result of applying an X-gate to the previous measuremet), if instead $0000$ is measured, we register $1$ error. Here we consider convinient to define two types of errors:

- True errors: If a 0 measurement is expected and instead we measure a 1.
- False errors:  If a 1 measurement is expected and instead we measure a 0.

The nomenclature is chosen here, besides for the relation between the $0-1$ false and true nomenclature, has to due with the fact that we are expecting from previous analysis of this type in the literature an **asymmetry** between the transitions $\ket{0} \ \rightarrow \ket{1}$ and $\ket{1} \ \rightarrow \ket{0}$, being the former the one not found in that analysis. The idea here is to verify the certainty of this hypothesis. We set the **null hypothesis** to be the one where we **do not find the asymmetry**, which would be true if true errors are measured, and false in the other case.