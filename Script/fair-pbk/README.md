# Running the python scripts

## Installation

To install the required packages, type:

```
pip install -r requirements.txt
```

## Compile model

Create and annotate SBML model from antimony file and annotations csv:

```
python ./Script/fair-pbk/create_sbml.py
```

## Run simulation scenarios

Run validation scenarios defined in `.yaml` files in the `scenarios` folder:

```
python ./Script/fair-pbk/run_validation.py
```

