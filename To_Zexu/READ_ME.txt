This code deals with B2 to R martensitic transformation in NiTi system under the effect of Ni4Ti3 precipitate.The major file includes:

time_evolution.c:
This file calculate the chemical potential and calculate the governing equation.

inhom_v3.c:
This file calculate the local stress field and elastic driving force. The boundary condition, external load, SFTS and elastic modulus

constants.h:
Here we regulated the system size (128 cube), dt, gradient coefficient keta, magnitude of random noise and time.

global_variables.c:
This file regulates the total step and step for printing an output file.

initialization.c:
Here we set the initial condition of the simulation. You can decide if you want to start from a file in the input_files or manual input.


To run the code, you can upload the entire folder to your osc account. Type "make" to compile it and use the following command to let it run in the compute node:
qsub run_job.sh

While running, you can check the status through qstat -u YOUR_USER_NAME

