# parameters-class
A simple, flexible Python class to manage configurable parameters in the RFM and the SSM1D.


## Example Usage
```python
# Define the case and the resolution
# specify what gases and CIA-pairs to include (see generate_case)
gases     = ['H2O']
ciapairs  = []

# dynamically create an argument dictionary
def generate_args(exp, gases, ciapairs, **thermo):
    return {'exp': exp, **{f'gas{i+1}': gas for i, gas in enumerate(gases)}, 'valid_ciapairs': ciapairs, **thermo}
args = generate_args('earth',gases,ciapairs,RHs=0.75,RHmid=0.99,RHtrp=0.75,uniform=1,Ts=315,Tmid=250,Ttrp=200)
        
# create a class instance and generate an RFM or SSM1D case from argument dictionary
par = parameters()
par.generate_case(**args)

# to switch between SSM1D version with and without CTM absorption
par.update_use_ssm1d_ctm(True) # or False
```

## Installation
```bash
pip install .
```
