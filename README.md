# AltSplicing2 Workflow Tutorial

### 1. Install rMATS (optional)
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


### 7. Run the workflow
```
nohup ./workflow.sh &
```

