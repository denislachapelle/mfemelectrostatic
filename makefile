#
# DL230604
# this file compile all programs developped while experimenting with mfem.
#

MFEM_DIR ?= /mnt/c/mfem-4.6
MFEM_BUILD_DIR ?= /mnt/c/mfem-4.6

#MFEM_DIR ?= /home/denislachapelle2003/fem/mfem-4.6
#MFEM_BUILD_DIR ?= /home/denislachapelle2003/fem/mfem-4.6

#COMMON_LIB = -L$(MFEM_BUILD_DIR)/miniapps/common -lmfem-common

all: electrostatic magnetostatic magnetostatic2d skineffect2d

electrostatic: electrostatic.cpp
	g++ -g -o electrostatic  -std=c++11 -I$(MFEM_DIR) electrostatic.cpp  -L$(MFEM_BUILD_DIR) -lmfem -L$(MFEM_BUILD_DIR)/miniapps/common -lrt
#	g++ -o electrostatic -O3 -std=c++11 -I$(MFEM_DIR) electrostatic.cpp  -L$(MFEM_BUILD_DIR) -lmfem -L$(MFEM_BUILD_DIR)/miniapps/common -lmfem-common -lrt

magnetostatic: magnetostatic.cpp
	g++ -g -o magnetostatic  -std=c++11 -I$(MFEM_DIR) magnetostatic.cpp  -L$(MFEM_BUILD_DIR) -lmfem -L$(MFEM_BUILD_DIR)/miniapps/common -lrt

magnetostatic2d: magnetostatic2d.cpp
	g++ -g -o magnetostatic2d  -std=c++11 -I$(MFEM_DIR) magnetostatic2d.cpp  -L$(MFEM_BUILD_DIR) -lmfem -L$(MFEM_BUILD_DIR)/miniapps/common -lrt

skineffect2d: skineffect2d.cpp
	g++ -g -o skineffect2d  -std=c++11 -I$(MFEM_DIR) skineffect2d.cpp  -L$(MFEM_BUILD_DIR) -lmfem -L$(MFEM_BUILD_DIR)/miniapps/common -lrt

