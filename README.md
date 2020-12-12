# infection_risk_calculator
In this project, We propose an infection risk algorithm which accepts building data and a set of parameters regarding occupants and infection rates of the surrounding community. Code and assumptions made in the algorithm will be clearly explained to users for transparency. The resulting algorithm will be referenced in our quarter 2 project. 




## Opening up the project

### Responsibilities

Source codes are located in /src/. 

To see the notebook containing the example EDA process, open Checkpoint_2_EDA in /notebooks.

To use run.py, create a notebook at the root directory and import it, then input main(target), 
where target is a list of targets. Example run is: run.main(["eda_co2", "eda_humidity", "eda_light"]). If you want to run it in command line, input python run.py [targets].  
We currently have the following targets available: "eda_co2", "eda_humidity", "eda_light" and "data". Those "eda_(sensors_types)" targets gives visualization of correlation between PIR sensors and input sensors, and returns a dataframes containing the correlation in each rooms. The "data" simply print string "Data is already cleaned in checkpoint 2" because it's already cleaned by the author. More information about EDA can be found in the EDA notebooks mentioned above. 



### Responsibilities:

* Etienne Doidic worked on the report and project structure.
* Nicholas Kho cleaned data and created graphs and notebooks, helped with EDA. 
* Zhexu Li developed correlation calculations, developed source codes and updated most of the notebook. 
