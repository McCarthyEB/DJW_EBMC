This will need a virtual environment to work, set up with python 3.9.2, and the following lines appended to virtualenv/bin/activate:

```
module purge
unset PYTHONPATH
unset LD_LIBRARY_PATH
module load python/3.9.2
export PATH=~/bin/ase/bin:$PATH
export PYTHONPATH=~/bin/ase/lib/pythonX.X/site-packages:$PYTHONPATH
export PYTHONPATH=$HOME/ase_env/lib/python3.9/site-packages:$PYTHONPATH
module load compiler/gnu/9/2.0
```

ase.3.23.0 must be present in bin for this to work!

```
git clone --branch 3.23.0 https://gitlab.com/ase/ase.git
cd ase
python setup.py install --user

```

