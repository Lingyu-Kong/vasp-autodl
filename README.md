## 先安装环境依赖

```bash
sudo apt update
sudo apt install -y \
  make \
  build-essential \
  g++ \
  gfortran \
  libopenblas-dev \
  mpich \
  libfftw3-dev \
  libhdf5-mpich-dev
```

检查环境
```bash
mpichversion
h5pcc -show
```

## 确保移除所有openmpi

一般来讲运行```sudo apt remove --purge -y openmpi-bin openmpi-common libopenmpi-dev```

十分建议运行```apt list --installed | grep -i openmpi```来检查还有没有其他openmpi相关包，一并删除

然后清理```sudo apt autoremove -y```

## 手动编译安装scalapack

```bash
git clone https://github.com/Reference-ScaLAPACK/scalapack.git
cd scalapack
mkdir build && cd build
cmake .. \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_COMPILER=mpicc \
  -DBLAS_LIBRARIES=/usr/lib/x86_64-linux-gnu/libopenblas.so \ 
  -DLAPACK_LIBRARIES=/usr/lib/x86_64-linux-gnu/libopenblas.so \
  -DCMAKE_INSTALL_PREFIX=/usr/local/scalapack-mpich
make -j4
sudo make install
```

这里需要确认libopenblas的位置，一般在指令所示位置，如果不放心可以先运行```apt list --installed | grep -i openblas```看一眼

执行完毕后就手动编译安装了scalapack到```/usr/local/scalapack-mpich```。然后添加如下到```.bashrc```里面: ```export LD_LIBRARY_PATH=/usr/local/scalapack-mpich/lib:$LD_LIBRARY_PATH```。添加完后重启shell

## 修改makefile.include修改链接改成mpich

目前这里给的makefile.include已经修改完了

## 编译vasp

```bash
cd vasp.6.3.0
mkdir bin
make DEPS=1 -j
```

## 修改路径配置和并行计算设置

在```.bashrc```里面添加

```bash
export OMP_NUM_THREADS=1
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
export PMG_VASP_PSP_DIR=~/autodl-tmp/vasp_autodl/VASP_PP_pymatgen
export PATH=/root/autodl-tmp/vasp_autodl/vasp.6.3.0/bin/:$PATH
```

注意把里面的路径换成正确的
