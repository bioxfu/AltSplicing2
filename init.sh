source activate rMATS
tar -xzf rMATS.4.0.2.tgz
chmod +x $PWD/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py
export PATH=$PWD/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/:$PATH

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/miniconda2/envs/gmatic/lib/
export LD_LIBRARY_PATH
