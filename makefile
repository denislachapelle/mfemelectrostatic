#
# DL230604
# this file compile all programs developped while experimenting with mfem.
#

#MFEM_DIR ?= /mnt/c/mfem-4.5.2
#MFEM_BUILD_DIR ?= /mnt/c/mfem-4.5.2

MFEM_DIR ?= /home/denislachapelle2003/fem/mfem-4.6
MFEM_BUILD_DIR ?= /home/denislachapelle2003/fem/mfem-4.6

#COMMON_LIB = -L$(MFEM_BUILD_DIR)/miniapps/common -lmfem-common

all: electrostatic

electrostatic: electrostatic.cpp
	g++ -g -o electrostatic  -std=c++11 -I$(MFEM_DIR) electrostatic.cpp  -L$(MFEM_BUILD_DIR) -lmfem -L$(MFEM_BUILD_DIR)/miniapps/common -lmfem-common -lrt
#	g++ -o electrostatic -O3 -std=c++11 -I$(MFEM_DIR) electrostatic.cpp  -L$(MFEM_BUILD_DIR) -lmfem -L$(MFEM_BUILD_DIR)/miniapps/common -lmfem-common -lrt

