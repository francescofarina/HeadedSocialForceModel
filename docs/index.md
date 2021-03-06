# The Headed Social Force Model

The Headed Social Force Model (HSFM) is an enhancement of the traditional Social Force Model (SFM) which explicitly accounts for the pedestrians’ heading. The motion of each individual is described by means of a dynamic model. Forces and torques driving the dynamics of each pedestrian are generated with the purpose of maximizing the realism of the resulting trajectories. In doing so, several conflicting objectives have to be taken into account. In low density scenarios, the pedestrians’ motion should be as smooth as possible, consistently with what is empirically observed. In these circumstances, lateral motions should be avoided because individuals walk ahead most of the time. On the contrary, in crowded or cluttered environments, the interaction among pedestrians, as well as between pedestrians and the environment, is stronger and determines most of the pedestrians’ trajectories. The solution proposed in the HSFM consists in computing the model inputs as suitable functions of the force terms adopted in the traditional SFM.

![Image](analysis.jpg)

The Headed Social Force Model has been introduced in our recent work:

>[Farina, F., Fontanelli, D., Garulli, A., Giannitrapani, A., & Prattichizzo, D. (2017). Walking Ahead: The Headed Social Force Model. PloS one, 12(1), e0169734.](http://dx.doi.org/10.1371/journal.pone.0169734). 

A preliminary version of the model has been presented in:
>[Farina, F., Fontanelli, D., Garulli, A., Giannitrapani, A., & Prattichizzo, D. (2016, December). When Helbing meets Laumond: The Headed Social Force Model. In Decision and Control (CDC), 2016 IEEE 55th Conference on (pp. 3548-3553). IEEE.](10.1109/CDC.2016.7798802)

You can cite these article if you use this code and our model for your research.

Both a MATLAB and a Python 3 implementation are available.

## Guideline to the Python 3 code
Create a python virtual environment through `virtualenv`. If you do not have `virtualenv` installed, install it by writing in your terminal

    pip3 install virtualenv

Create a new directory *dir* and clone the repository inside it.
In the same directory create a new environment

    virtualenv ENV
    
Then activate it

    source ENV/bin/activate
    
and install the required packages

    pip install numpy matplotlib scipy
    
Now, go to the python scripts directory

    cd HeadedSocialForceModel/python_scripts
    
Here you can run the file launcher.py.

    python launcher.py

If you want to do other simulations or change the model parameter follow the instructions in launcher.py and modify the parameters in aux_functions.py.

## Guideline to the MATLAB code
Open the function map_def.m and follow the instructions to define a map.

Open launcher.m and for each group to be created, insert:
- the number of individuals; 
- the starting point;
- the waypoints.

If needed, specify different model parameters.
