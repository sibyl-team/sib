## Sibilla


Belief Propagation for inference in epidemics. This code implements and expands the method described in

https://www.nature.com/articles/srep27538

* Requirements:

- A C++11 compiler
- python3
- pybind11
- header only boost libraries

* Compilation:

Adjust makefile and type make. You'll obtain a standalone CLI ./sib executable, plus a dynamic library containing a python module `sib`

* In brief:

```python3
> import sib
> # a contact between 1 and 2 at time 0, with transmission probability 0.5
> # a contact between 2 and 3 at time 1, with transmission probability 0.4
> contacts = [(1,2,0,0.5), (2,3,1,0.4)]
> # node 2 was observed in state I (states are S=0,I=1,R=2) at time 4. 
> observations = [(2,1,4)]
> f = sib.FactorGraph(contacts, observations)
> sib.converge(f)
> sib.bt
```
