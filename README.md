### Asymetry Program and File System
*Ivan Zhenchuk - U of R - 2025-08-11*
---

Asymetry program is made to calculate beem spin asymetry for specific particle, plot it, and extract sigma LTprime from 
the asymetry histogram. This program works with its own file system (hardcoded) to store results.

## How to use it?

To prepare the program for execution run the `make` command in asym_program directory.

# Making asymetry calculation per phi bin

To start processing data, drop a root file into the -input directory. To perform iterations (or single run) you need to 
make or fill an instructional csv iterator_inp.csv. In the iterator_inp.csv you need to specify run criteria such as
filename, particle, shms angle, t min, tmax. iterator_inp.csv must contain a header line and be in the order that example
criteria are provided in. Each row is a new run.

```
iterator_inp.csv Example:
filename,particle,shms angle,t min,tmax
Q0p5W2p40right_highe_kaon,kaonL,right,0.06225,0.09975
Q0p5W2p40right_highe_kaon,kaonL,right,0.09975,0.12375
...
```

It is important to consider the naming convention of your root files. CheckNaming() assigns the branch names that are
used for cuts and histogram creation. CheckNaming() should be modified or reviewed to match your root tree naming convention.

Before running the asymmetry calculation it is important to fill out the cut_config.json file in the -config directory.
Cut config provides easy editing of cut ranges for Coin Time, HMS cherenkov and calorimeter, PID, and yield range
for final asymetry calculation.

List of files
1. ROOT data file with positive and negative helicities for main data and dummy data.
2. iterator_inp.csv file to specify run criteria and number of runs
3. cut_config.json file to specify cuts 

Once these files are provided you can tun the iterator.py with python and use the "Asym" option.

# Kinematic calculations

To calculate kinematics for a specific shms angle and data, you need to have file 1 and 3 from previous section prepared.
The difference is you also need a different csv criteria file. kin_wa_inp.csv is used to provide particle,shms pos,tmin,tmax
for each run of kinematic calculation. 

```
kin_wa_inp.csv Example:
particle,shms pos,tmin,tmax
kaonL,left,0.06225,0.09975
kaonL,left,0.09975,0.12375
...

```

To create the kinematic per t-bin plots run the iterator.py with python and use the "Kin" option

# Weighted avereage for SHMS angles

This program is designed to combine the center, right, and left SHMS angles with an error weighted average.
The weighted average is calculated using created Asym plots for center, right, and left (CRL) SHMS angles.
This option also requires a csv inout. rcl_wa_inp.csv is much smaller and depends on the name of Asym plots
created by Asym option. 

rcl_wa_inp.csv Example:
particle,files tbin
kaonL,-t0.062250-0.099750
kaonL,-t0.099750-0.123750
kaonL,-t0.123750-0.191250

WA option can be used as long as the Asym files are present in -output directory. WA is outputted as a plot
of average kinematic value per average Mandel t value for selected t bin.

## To Consider

This program is designed to calculate asymetry for kaon Lambda on "MMK" missing mass histogram. Also, it asumes the user
has:

1. 3 SHMS setting
2. Dummy data and set dummy scaling factor
3. set cut values and uses specific cut variables
4. Desides t-bins in a separate program
5. Uses 15 phi bins for asymmetry plot.

These asumptions are semi-hardcoded into the code and require switching of most parameters when switching datasets.
The dataset that perfectly fits most setting is Q2=0.5 and W=2.40. For different settings it is important to investigate
MMK/MMPi and other branches to deside necessary cut.

For PID cuts, the program uses corrected RF distribution. asym_main.cpp has heavy gas and aerogel cherenkov cuts
coded in but not utilized. In case of a setting having no RF data, hgcer and aero can be used for PID.