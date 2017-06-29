# README #
This Simple.h, Simple.cpp are LAMMPS scripts for the simplified lubrication force. The theory (or rather the algebra) behind the simplification is given in the pdf file.

### What is this repository for? ###

* Lubrication force between two particles.
* 1.0

### How do I get set up? ###

* Copy Simple.h and Simple.cpp in LAMMPS installation's src folder and recreate an executable.
* pair_lubricate.h
* Modify the two particle collison code: collision.lmp, and compare it using the matlab scripts provided.

### Contribution guidelines ###

* Please email Ranga for updates or features needed. You are welcome to modify and use this code to your specs.

### Who do I talk to? ###

* Ranga (r.radhakrishnan@ed.ac.uk)


### List of things to work on for future updates ###

1. Make the cutoffs with arithmetic mean by default.

### Acknowledgements ###
Tim Najuch for leading the discussion on this topic and helping me develop this code.
Dr. Chris Ness for verifying the theory behind the forces.
Prof. Jin Sun for his kind support and encouragement.