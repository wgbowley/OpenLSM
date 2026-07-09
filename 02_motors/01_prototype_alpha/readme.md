# Prototype Alpha

> [!WARNING]
> This is a conceptual design intended to explore design ideas and engineering trade-offs.
> It has not been manufactured, experimentally validated, or verified for performance.
> Do not assume the design is suitable for fabrication without further analysis.
>
> Revisions 1 and 2 share the same electromagnetic configuration, only the thermal management and mechanical systems have been modified.

# Overview

Prototype Alpha is an `ironless planar linear motor` with a polylactic acid (PLA) armature featuring 6 slots, hand wound using 0.2 mm diameter enameled copper wire and 5 mm wide Kapton tape, with 2 slots in-series per phase (WYE). The stator, similar to the armature, was printed in PLA and had 4 pole pairs per armature length and 10 pole pairs total. 

A linear encoder `(AS5311)` was used to measure the motor position and a SimpleFOC shield was used to control the motor's 3 phases using closed-loop control.

## Mechanical & Thermal

The pole length was `10 mm`, pole width was `40 mm` and the pole-to-pole pitch was `14 mm`. The armature had slot lengths of `9 mm` with slot thickness of `3 mm` and a slot holder thickness of `3 mm`. The slot heights were `6 mm` and the slot-to-slot pitch was `18.66 mm`. This design had no passive or active thermal design, instead relying on convection from the coil surfaces to the atmosphere to dissipate heat.

## Electromagnetic

As stated above, the pole pitch was `14 mm` and the slot pitch was `18.66 mm`, with `2 slots` in series per phase and 4 pole pairs per armature length, achieving a pole-slot ratio of `6 slots - 8 poles` and hence a synchronous frequency:

$$ v = f \cdot 2 \cdot p_{\text{pitch}} \;\rightarrow\; f = \frac{v}{2 \cdot p_{\text{pitch}} } = \frac{v}{0.028} $$

Due to the `6 slots - 8 poles` configuration, this motor is a fractional-slot concentrated motor, which may have helped decrease cogging torque because the poles and slots do not line up evenly and hence have no preferred configuration.

## Construction

The armature slots were wound individually with approximately the same number of turns. The first slot was then connected in series with the 4th slot, the second slot with the 5th slot, etc. Those `2 slot` phases were then connected in a WYE configuration using a proto-board, and a 3-pin JST connector was used as the motor input. The neutral phase tap was left exposed as a common reference point for measuring.

The stator was constructed using through-thickness magnetized N52 poles, with each pole rotated `180 degrees` with respect to the last to produce an effective stator field. The poles were attached using generic superglue, and an `MG-series` linear rail was mounted to the back of the stator, which connected the armature and stator together via U-shaped mounting plates.

## Results

The motor moved and was able to move `1 kg` load with ease, but the amount of precision and force is unknown. The coil forms collapsed in on themselves due to thermal stress, which caused the slots to detach from the armature and increase friction against the stator. 

The control scheme seemed to work well, but the exact precision was not quantified due to thermal stress. It seemed to have issues with deriving `PID` control parameters, possibly due to the fractional-slot configuration or inconsistent resistance and inductance per phase.

## Conclusions 

The `planar` linear motor did work and produce force, but the amount could not be quantified, nor could the control scheme be fully tested. Future iterations should focus on improving thermal characteristics by managing phase resistance, as `ironless` motors lose a lot of energy to `copper losses`. Beyond thermal considerations, future iterations should aim for similar inductance and resistance for each of the motor's phases to potentially improve control characteristics.