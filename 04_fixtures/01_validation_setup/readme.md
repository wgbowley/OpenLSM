### 01_validation_setup

The validation setup consists of a main board designed to collect position data from the magnetic linear encoder, as well as thermal and acceleration data from ICs on the board. It then sends that data over RS-485 to a Teensy 4.1, which controls the motor using a SimpleFOC shield while also logging the data.

<div align="center">
  <img src="../../05_media/04_fixtures/01_validation_setup/rev_1/stand.png" width="700" alt="Test-stand">
</div>

---
