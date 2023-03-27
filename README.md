# FEMM-Coilgun
Easy &amp; quick multistage coilgun design utility for FEMM.

https://user-images.githubusercontent.com/80536083/228027307-fd278073-ac77-4a2f-aae5-a9051f97a105.mp4

## Requirements
 - The FEMM software (https://www.femm.info/wiki/HomePage)
 - pyFEMM
 - matplotlib
 - moviepy
 
## How to Use
See **example_coilgun.py**. All lengths are in mm and currents are in A.

A coilgun is defined as follows:
```python
# a coilgun starts with a loading block, where the projectile sits before the first coil is energized
# loader parameterss: index, inner radius, thickness, length, material, axial position
loader = LoadingBlock(0, 16, 75, 100, "1010 Steel", 100)

# then, the coil stages are defined
# stage parameters: index, number of turns, length, inner radius, wire diameter,
#                   barrel rear thickness, barrel fore thickness, barrel material, axial position
stage1 = Stage(1, 400, 100, 16, 2.5, 50, 35, "1010 Steel", 0)
stage2 = Stage(2, 400, 100, 16, 2.5, 35, 35, "1010 Steel", -100)
stage3 = Stage(3, 400, 100, 16, 2.5, 35, 25, "1010 Steel", -200)
stage4 = Stage(4, 400, 100, 16, 2.5, 25, 25, "1010 Steel", -300)

# all stages, starting with the loader, are put in a list
stages = [loader, stage1, stage2, stage3, stage4]

# the gun is defined
gun = Coilgun(stages)
```

A projectile is defined as follows:
```python
# projectile parameters: radius, length, material, initial position, initial velocity, material density
# (initial position and velocity will likely be overridden by analysis function and mostly exists for debugging)
proj = Projectile(15, 100, "1010 Steel", 50, 0, 7870)
```

Then, a function must be defined to control the coil activation sequence. The function will take the coilgun stage and the projectile instance as arguments and return the current. An example activation function is provided below:
```python
def activate_one_forward(proj, stage):
    z_proj = proj.z_pos
    z_stage = stage.z_pos
    maxdist = stage.L * 1.5

    if z_proj - z_stage > 0 and z_stage + maxdist > z_proj:
        return 60
    else:
        return 0
```
The function above provides 60 Amps to a stage when the projectile approaches it, and turns it off when the projectile gets to the middle of the stage.

With all components and the activation function defined, we are ready to perform the analysis. Currently, there are two analysis methods implemented.

**Discrete position** method iterates projectile position along the barrel, calculating force and acceleration on the projectile.

**Discrete time** method starts calculating forces and accelerations from the initial position and performs Euler integration until the projectile exits the barrel.

These analyses can be performed as follows:
```python
# the last 3 arguments are image view padding, plot max. magnetic density, discrete position interval in mm
PerformAnalysisDiscretePosition(gun, proj, activate_one_forward, "mycoilgun", "mycoilgunfile", 30, 3, 20)

# the last 3 arguments are image view padding, plot max. magnetic density, discrete time interval in seconds
PerformAnalysisDiscreteTime(gun, proj, activate_one_forward, "mycoilgun", "mycoilgunfile", 30, 3, 0.0025)

# ...or both can be performed at once
PerformAnalysisFull(gun, proj, activate_one_forward, "mycoilgun", "mycoilgunfile", 30, 3, 0.0025, 20)
```

## Results

See **example_results** directory.

A video file displaying the magnetic flux density animation is exported, with individual frames available as bitmaps. The muzzle velocity is printed. Force, acceleration and velocity data are plotted.

For more data, the .fem files that will be generated when you run the script can be opened in FEMM.
