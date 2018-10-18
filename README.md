# AltSplicing2 Workflow Tutorial

### 1. Install rMATS and PAIRADISE (optional)
```
## create conda env
conda create -n rMATS numpy scipy pysam samtools star rmats2sashimiplot

## which version to use
python -c 'import sys; print sys.maxunicode'
## 1114111: UCS4
## 65535: UCS2

## add rMATS to PATH
tar -xzf rMATS.4.0.2.tgz
chmod +x $PWD/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py
export PATH=$PWD/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/:$PATH

## The PAIRADISE statistical model is applicable to many forms of allele-specific isoform variation (e.g. RNA editing), and can be used as a generic statistical model for RNA-seq studies involving **paired replicates**.
PAIRADISE can be used as a stand-alone statistical model in R by downloading the and installing the file 'PAIRADISE_RPackage_1.0.tar.gz'. The following command will install the package within R:

install.packages('PAIRADISE_RPackage_1.0.tar.gz', repos = NULL, type = 'source')
```


### 2. Login the fat node
```
ssh node2
```

### 3. Clone the repository
```
cd ~/Project/${GROUP}_${DATE}
git clone https://github.com/bioxfu/AltSplicing2
cd AltSplicing2
```

### 4. Initiate the project
```
. init.sh
```


### 5. Create *workflow.sh* based on the example
```
cp example/workflow.sh workflow.sh

# edit workflow.sh
```


### 6. Run the workflow
```
nohup ./workflow.sh &
```

### Tips
rMATSexe: error while loading shared libraries: libgsl.so.0: cannot open shared object file: No such file or directory
```
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/miniconda2/envs/gmatic/lib/
export LD_LIBRARY_PATH
```

