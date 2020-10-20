################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CMD_SRCS += \
../DSP2834x_Headers_BIOS.cmd \
../TMS320C28346.cmd 

CFG_SRCS += \
../clock.cfg 

ASM_SRCS += \
../DSP2834x_CodeStartBranch.asm \
../DSP2834x_usDelay.asm 

C_SRCS += \
../DSP2834x_DefaultIsr.c \
../DSP2834x_GlobalVariableDefs.c \
../DSP2834x_PieCtrl.c \
../DSP2834x_PieVect.c \
../DSP2834x_Sci.c \
../DSP2834x_Spi.c \
../DSP2834x_SysCtrl.c \
../function.c \
../main.c 

OBJS += \
./DSP2834x_CodeStartBranch.obj \
./DSP2834x_DefaultIsr.obj \
./DSP2834x_GlobalVariableDefs.obj \
./DSP2834x_PieCtrl.obj \
./DSP2834x_PieVect.obj \
./DSP2834x_Sci.obj \
./DSP2834x_Spi.obj \
./DSP2834x_SysCtrl.obj \
./DSP2834x_usDelay.obj \
./function.obj \
./main.obj 

ASM_DEPS += \
./DSP2834x_CodeStartBranch.pp \
./DSP2834x_usDelay.pp 

C_DEPS += \
./DSP2834x_DefaultIsr.pp \
./DSP2834x_GlobalVariableDefs.pp \
./DSP2834x_PieCtrl.pp \
./DSP2834x_PieVect.pp \
./DSP2834x_Sci.pp \
./DSP2834x_Spi.pp \
./DSP2834x_SysCtrl.pp \
./function.pp \
./main.pp 

GEN_MISC_DIRS += \
./configPkg/ 

GEN_CMDS += \
./configPkg/linker.cmd 

GEN_OPTS += \
./configPkg/compiler.opt 

GEN_FILES += \
./configPkg/linker.cmd \
./configPkg/compiler.opt 

GEN_FILES__QUOTED += \
"configPkg\linker.cmd" \
"configPkg\compiler.opt" 

GEN_MISC_DIRS__QUOTED += \
"configPkg\" 

C_DEPS__QUOTED += \
"DSP2834x_DefaultIsr.pp" \
"DSP2834x_GlobalVariableDefs.pp" \
"DSP2834x_PieCtrl.pp" \
"DSP2834x_PieVect.pp" \
"DSP2834x_Sci.pp" \
"DSP2834x_Spi.pp" \
"DSP2834x_SysCtrl.pp" \
"function.pp" \
"main.pp" 

OBJS__QUOTED += \
"DSP2834x_CodeStartBranch.obj" \
"DSP2834x_DefaultIsr.obj" \
"DSP2834x_GlobalVariableDefs.obj" \
"DSP2834x_PieCtrl.obj" \
"DSP2834x_PieVect.obj" \
"DSP2834x_Sci.obj" \
"DSP2834x_Spi.obj" \
"DSP2834x_SysCtrl.obj" \
"DSP2834x_usDelay.obj" \
"function.obj" \
"main.obj" 

ASM_DEPS__QUOTED += \
"DSP2834x_CodeStartBranch.pp" \
"DSP2834x_usDelay.pp" 

ASM_SRCS__QUOTED += \
"../DSP2834x_CodeStartBranch.asm" \
"../DSP2834x_usDelay.asm" 

C_SRCS__QUOTED += \
"../DSP2834x_GlobalVariableDefs.c" \
"../DSP2834x_PieCtrl.c" \
"../DSP2834x_PieVect.c" \
"../DSP2834x_Sci.c" \
"../DSP2834x_Spi.c" \
"../DSP2834x_SysCtrl.c" \
"../function.c" \
"../main.c" 


