# Prototype 1 Linear Motor Optimization

### Note:
The simulation of the prototype 1 linear motor requires `pyfea`, a custom tool I built. Hence I would recommend not running it as the tool is quite anti-user. In the nature that if you didn't develop it you, its not worth using. Do note however if your reading this after `pyfea` release than it is worth running this script.

### Assumptions
- Secant phase inductance is used throughout the PointToPoint simulation due to numerical instability from 
calculating inductance via flux linkage over current and than back-emf from delta of flux linkage over time step

- Duty cycle is approximation for the dynamic nature of the motor. It scales down the maximal power loss the motor sees when calculating the asymptotic temperatures of both the slot and poles
due to the fact that a 3d printer never sees true 100% power. Steady state != printing basically.

## Use prototype 1 simulation   
It is as simple as installing `pyfea` via either pip install or setup-tools and than just running simulation.py or optimization.py in this folder. To change parameters use configuration.uiv and follow the .uiv standard notation if your changing element beyond their numerical values.