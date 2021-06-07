# My population infection
A simple model for simulating virus spreading, written in C and MPI. The simulation relies on a world split into a grid of equally-sized countries. 
Individuals follow a linear motion, and each country is assigned to a separate MPI process. The spreading distance of the virus, the exposure time to get infected, the duration of the infection and of the immunity can be configured.
At the end of each simulated day, the program produces a summary with the count of susceptible, infected and immune individuals for each country.

Read the [project report](/releases/latest/download/mpi_report) for more detailed information.

## Getting started
First clone this repo on your local machine:
```
git clone https://github.com/fuljo/my-population-infection.git
cd ./my-population-infection
```

You can choose to compile and run the project in the local MPI environment, or run it with Docker.

### Compiling and running locally
1. Install the required dependencies (assuming you are on Ubuntu or Debian):
    ```
    apt-get install build-essential gcc gfortran binutils \
                    openssh libopenmpi-dev
    ```
2. Compile the source code:
    ```
    cd ./src
    make
    ```
3. Will produce the executable `my-population-infection`
4. Run the program (see below or run with `--help` for the full list of parameters):
    ```
    mpirun -np 4 --oversubscribe ./my-population-infection \
        -N 1000 \
        -I 100 \
        -W 2000 \
        -L 2000 \
        -w 1000 \
        -l 1000 \
        -v 1.4 \
        -d 30 \
        --t-infection=10 \
        --t-recovery=120 \
        --t-immunity=120 \
        --sim-step=1 \
        --sim-length=1 \
        --log-level INFO
    ```
    Alternatively you can use the `run` target in `src/Makefile`.
    
    See the OpenMPI documentation to run the simulation over different computing nodes.

5. Print the summary
    ```
    cat ./results/summary.csv | column -t -s,
    ```

### Using Docker
Alternatively you can use Docker and Docker Compose to run the simulation, even across different containers, without the need of configuring MPI on your machine.

1. Build the docker image
    ```
    make docker
    ```
2. Edit the `docker-compose.yml` to set the parameters of the program. Remember that the total number of processes must match the number of countries.
3. Start the containers:
    ```
    make compose
    ```
4. Print the summary
    ```
    cat ./results/summary.csv | column -t -s,
    ```

### Parameters
The simulation parameters are given as command line options to the main executable.
Here is the output of `my-population-infection --help`:
```
Usage: my-population-infection [OPTION...]
            -N int -I int -W int -L int -w int -l int -d float -v float
            --sim-step seconds --sim-length days
A simple model for virus spreading.

 Population options
  -I, --inf-individuals=INT  Number of initially infected individuals
  -N, --num-individuals=INT  Number of individuals

 World options (lengths in meters)
  -l, --country-length=INT   Length of a single country, must divide L
  -L, --world-length=INT     Length of the world rectangle
  -w, --country-width=INT    Width of a single country, must divide W
  -W, --world-width=INT      Width of the world rectangle

 Individual options (times in seconds)
  -d, --spreading-distance=FLOAT   Maximum spreading distance in meters
      --t-immunity=INT       Duration of immunity period after recovering
                             (default 90 days)
      --t-infection=INT      Minimum continuous exposure time for getting
                             infected (default 10 min)
      --t-recovery=INT       Time needed for recovering (default 10 days)
  -v, --velocity=FLOAT       Moving speed for and individual in m/s

 Simulation options
      --rand-seed=INT        Seed for PRNG. (default time(NULL))
      --sim-length=INT       Length of the simulation in days
      --sim-step=INT         Simulation step in seconds

 Logging options
      --log-level=[TRACE|DEBUG|INFO|WARN|ERROR|FATAL]
                             Logging level (default INFO)
      --write-trace          Write the file results/trace.csv with details
                             about each individual at each time step

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version

Mandatory or optional arguments to long options are also mandatory or optional
for any corresponding short options.

Must be run in an MPI environment where the total number of processes is equal
to the number of countries (W/w * L/l).
Produces a summary in ./results/summary.csv with the number of susceptible,
infected and immune individuals at each time step, and a file
./results/trace_{country}.csv for each country if the --write-trace flag is
given.
```

## Tools

Aside from the main program, we provide some complementary tools in the `tools` directory.

### Profiler
Measures the execution time with a varying number of countries or individuals. Results are written to a csv file.
```
./profile_countries.sh min increment max 
./profile_individuals.sh min increment max 
```
After this you can produce a `png` plot for each of the csv files by running
```
python3 ./profiler_plot.py
```

![Profile countries](/assets/profile_countries_1_20.png) ![Profile individuals](/assets/profile_individuals_10000_60000.png)

### Animation
1. Run the simulation with the `--write-trace` flag, so each node will produce a `./results/trace_{country}.csv` on its local filesystem, for each of its countries.
2. Gather these files together in a single `results` directory (this is done by default if you use our Docker compose setup).
3. Change the parameters at the end of `trace_animation.py` to match those of your simulation. Setting `t_target=None` will produce a complete animation, but you can use a value in seconds to cut it to the desired (simulated) time.
4. Produce the animation as a `mp4` video:
    ```
    python3 ./trace_animation.py
    ```

![Animation](/assets/anim_1000.gif)

## Report
A PDF project report describing the program, its design principles and including a performance analysis can be downloaded [here](/releases/latest/download/mpi_report.pdf).
Alternatively it can be locally compiled with XeLaTeX:
```
cd report
latexmk
```

## Credits
This project was created as part of the exam for the 2020 edition of the *Middleware Technologies for Distributed Systems* at Politecnico di Milano.

It contains the 3rd party library [log.c](https://github.com/rxi/log.c), distributed under MIT license.

## License
This program is distributed under the MIT license.