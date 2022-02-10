# ERSimu - Elementary Reaction Simulation

## About

A tool for exploring chemical kinetics and dynamics in elementary 
reaction systems.

## Design

![ersimu process](./doc/process.svg)

## Input file syntax

Lines starting with the hash (#) is a single line comment.

    # Comment

Species in excess, initial conditions, and numerical constants
are defined as follows.

    EXCESS H2O
    INITIAL H 1.29
    CONSTANT PI 3.14159

Reactions are defined in familiar-looking form which must be a single line.

    A ==> B

Above reaction implies law of mass action with unity reaction rate.
The reaction rate can be given explicitly by appending the line with
the vertical bar (|) and the reaction rate, either a numeric or a
constant.

    A ==> B  | 1E-5

Exceptions to mass action kinetics is also supported in the following
form by using the at sign (@).

    H2Q + BrO2 + BrO2 ==> HBrO2 + HBrO2 + Q @ K6*H2Q*BrO2

Numerical parameters for the simulation are defined by lines starting
with the `SIMULATION` keyword.

    SIMULATION NAME example_simulation
    SIMULATION PLOT X
    SIMULATION RUN 1
    SIMULATION T_END 55
    SIMULATION T_POINTS 10
	SIMULATION ATOL 1.1E-13
	SIMULATION RTOL 1.1E-13
    SIMULATION MAXIMUM_STEP_SIZE 0.01

In above example, the time range is equivalent to Octave syntax
`linspace(0, 55, 55*10)`.

The `LATEX` keyword specifies formatting mapping.

    LATEX H2O H_{2}O

Bibliography items can be specified as follows.

    BIBITEM TestSet \url{https://www.dm.uniba.it/~testset}

Refer to the files in [examples](examples) directory for complete
simulations.

## Contributing

All contributions are welcome. Bug reports, suggestions and feature
requests can be reported by creating a new
[issue](https://github.com/ptrktn/ersimu/issues). Code and
documentation contributions should be provided by creating a [pull
request](https://github.com/ptrktn/ersimu/pulls) (here is a good
[tutorial](https://www.dataschool.io/how-to-contribute-on-github/)).
Use imperative in commit messages.

## License

Licensed under the GNU General Public License Version 3, refer to the
file [LICENSE](LICENSE) for more information.
