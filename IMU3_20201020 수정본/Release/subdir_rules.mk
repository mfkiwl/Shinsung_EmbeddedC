################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Each subdirectory must supply rules for building sources it contributes
DSP2834x_CodeStartBranch.obj: ../DSP2834x_CodeStartBranch.asm $(GEN_OPTS) $(GEN_HDRS)
	@echo 'Building file: $<'
	@echo 'Invoking: C2000 Compiler'
	"C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/bin/cl2000" -v28 -ml -mt --float_support=fpu32 -O2 --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --include_path="C:/ti/controlSUITE/libs/math/IQmath/v15c/include" --include_path="C:/ti/controlSUITE/libs/math/FPUfastRTS/V100/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_common/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_headers/include" --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --advice:performance=all --diag_warning=225 --diag_warning=255 --display_error_number --diag_wrap=off --gen_func_subsections=on --preproc_with_compile --preproc_dependency="DSP2834x_CodeStartBranch.pp" $(GEN_OPTS__FLAG) "$<"
	@echo 'Finished building: $<'
	@echo ' '

DSP2834x_DefaultIsr.obj: ../DSP2834x_DefaultIsr.c $(GEN_OPTS) $(GEN_HDRS)
	@echo 'Building file: $<'
	@echo 'Invoking: C2000 Compiler'
	"C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/bin/cl2000" -v28 -ml -mt --float_support=fpu32 -O2 --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_headers/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_common/include" --include_path="C:/ti/controlSUITE/libs/math/FPUfastRTS/V100/include" --include_path="C:/ti/controlSUITE/libs/math/IQmath/v15c/include" --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --advice:performance=all --diag_warning=225 --diag_warning=255 --display_error_number --diag_wrap=off --gen_func_subsections=on --preproc_with_compile --preproc_dependency="DSP2834x_DefaultIsr.pp" $(GEN_OPTS__FLAG) "$<"
	@echo 'Finished building: $<'
	@echo ' '

DSP2834x_GlobalVariableDefs.obj: ../DSP2834x_GlobalVariableDefs.c $(GEN_OPTS) $(GEN_HDRS)
	@echo 'Building file: $<'
	@echo 'Invoking: C2000 Compiler'
	"C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/bin/cl2000" -v28 -ml -mt --float_support=fpu32 -O2 --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --include_path="C:/ti/controlSUITE/libs/math/IQmath/v15c/include" --include_path="C:/ti/controlSUITE/libs/math/FPUfastRTS/V100/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_common/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_headers/include" --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --advice:performance=all --diag_warning=225 --diag_warning=255 --display_error_number --diag_wrap=off --gen_func_subsections=on --preproc_with_compile --preproc_dependency="DSP2834x_GlobalVariableDefs.pp" $(GEN_OPTS__FLAG) "$<"
	@echo 'Finished building: $<'
	@echo ' '

DSP2834x_PieCtrl.obj: ../DSP2834x_PieCtrl.c $(GEN_OPTS) $(GEN_HDRS)
	@echo 'Building file: $<'
	@echo 'Invoking: C2000 Compiler'
	"C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/bin/cl2000" -v28 -ml -mt --float_support=fpu32 -O2 --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --include_path="C:/ti/controlSUITE/libs/math/IQmath/v15c/include" --include_path="C:/ti/controlSUITE/libs/math/FPUfastRTS/V100/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_common/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_headers/include" --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --advice:performance=all --diag_warning=225 --diag_warning=255 --display_error_number --diag_wrap=off --gen_func_subsections=on --preproc_with_compile --preproc_dependency="DSP2834x_PieCtrl.pp" $(GEN_OPTS__FLAG) "$<"
	@echo 'Finished building: $<'
	@echo ' '

DSP2834x_PieVect.obj: ../DSP2834x_PieVect.c $(GEN_OPTS) $(GEN_HDRS)
	@echo 'Building file: $<'
	@echo 'Invoking: C2000 Compiler'
	"C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/bin/cl2000" -v28 -ml -mt --float_support=fpu32 -O2 --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --include_path="C:/ti/controlSUITE/libs/math/IQmath/v15c/include" --include_path="C:/ti/controlSUITE/libs/math/FPUfastRTS/V100/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_common/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_headers/include" --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --advice:performance=all --diag_warning=225 --diag_warning=255 --display_error_number --diag_wrap=off --gen_func_subsections=on --preproc_with_compile --preproc_dependency="DSP2834x_PieVect.pp" $(GEN_OPTS__FLAG) "$<"
	@echo 'Finished building: $<'
	@echo ' '

DSP2834x_Sci.obj: ../DSP2834x_Sci.c $(GEN_OPTS) $(GEN_HDRS)
	@echo 'Building file: $<'
	@echo 'Invoking: C2000 Compiler'
	"C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/bin/cl2000" -v28 -ml -mt --float_support=fpu32 -O2 --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --include_path="C:/ti/controlSUITE/libs/math/IQmath/v15c/include" --include_path="C:/ti/controlSUITE/libs/math/FPUfastRTS/V100/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_common/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_headers/include" --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --advice:performance=all --diag_warning=225 --diag_warning=255 --display_error_number --diag_wrap=off --gen_func_subsections=on --preproc_with_compile --preproc_dependency="DSP2834x_Sci.pp" $(GEN_OPTS__FLAG) "$<"
	@echo 'Finished building: $<'
	@echo ' '

DSP2834x_Spi.obj: ../DSP2834x_Spi.c $(GEN_OPTS) $(GEN_HDRS)
	@echo 'Building file: $<'
	@echo 'Invoking: C2000 Compiler'
	"C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/bin/cl2000" -v28 -ml -mt --float_support=fpu32 -O2 --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --include_path="C:/ti/controlSUITE/libs/math/IQmath/v15c/include" --include_path="C:/ti/controlSUITE/libs/math/FPUfastRTS/V100/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_common/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_headers/include" --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --advice:performance=all --diag_warning=225 --diag_warning=255 --display_error_number --diag_wrap=off --gen_func_subsections=on --preproc_with_compile --preproc_dependency="DSP2834x_Spi.pp" $(GEN_OPTS__FLAG) "$<"
	@echo 'Finished building: $<'
	@echo ' '

DSP2834x_SysCtrl.obj: ../DSP2834x_SysCtrl.c $(GEN_OPTS) $(GEN_HDRS)
	@echo 'Building file: $<'
	@echo 'Invoking: C2000 Compiler'
	"C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/bin/cl2000" -v28 -ml -mt --float_support=fpu32 -O2 --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --include_path="C:/ti/controlSUITE/libs/math/IQmath/v15c/include" --include_path="C:/ti/controlSUITE/libs/math/FPUfastRTS/V100/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_common/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_headers/include" --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --advice:performance=all --diag_warning=225 --diag_warning=255 --display_error_number --diag_wrap=off --gen_func_subsections=on --preproc_with_compile --preproc_dependency="DSP2834x_SysCtrl.pp" $(GEN_OPTS__FLAG) "$<"
	@echo 'Finished building: $<'
	@echo ' '

DSP2834x_usDelay.obj: ../DSP2834x_usDelay.asm $(GEN_OPTS) $(GEN_HDRS)
	@echo 'Building file: $<'
	@echo 'Invoking: C2000 Compiler'
	"C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/bin/cl2000" -v28 -ml -mt --float_support=fpu32 -O2 --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --include_path="C:/ti/controlSUITE/libs/math/IQmath/v15c/include" --include_path="C:/ti/controlSUITE/libs/math/FPUfastRTS/V100/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_common/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_headers/include" --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --advice:performance=all --diag_warning=225 --diag_warning=255 --display_error_number --diag_wrap=off --gen_func_subsections=on --preproc_with_compile --preproc_dependency="DSP2834x_usDelay.pp" $(GEN_OPTS__FLAG) "$<"
	@echo 'Finished building: $<'
	@echo ' '

configPkg/linker.cmd: ../clock.cfg
	@echo 'Building file: $<'
	@echo 'Invoking: XDCtools'
	"C:/ti/xdctools_3_31_00_24_core/xs" --xdcpath="/packages;C:/ti/ccsv6/ccs_base;" xdc.tools.configuro -o configPkg -t ti.targets.C28_float -p ti.platforms.tms320x28:TMS320C28346 -r release -c "C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6" --compileOptions "-g --optimize_with_debug" "$<"
	@echo 'Finished building: $<'
	@echo ' '

configPkg/compiler.opt: | configPkg/linker.cmd
configPkg/: | configPkg/linker.cmd

function.obj: ../function.c $(GEN_OPTS) $(GEN_HDRS)
	@echo 'Building file: $<'
	@echo 'Invoking: C2000 Compiler'
	"C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/bin/cl2000" -v28 -ml -mt --float_support=fpu32 -O2 --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --include_path="C:/ti/controlSUITE/libs/math/IQmath/v15c/include" --include_path="C:/ti/controlSUITE/libs/math/FPUfastRTS/V100/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_common/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_headers/include" --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --advice:performance=all --diag_warning=225 --diag_warning=255 --display_error_number --diag_wrap=off --gen_func_subsections=on --preproc_with_compile --preproc_dependency="function.pp" $(GEN_OPTS__FLAG) "$<"
	@echo 'Finished building: $<'
	@echo ' '

main.obj: ../main.c $(GEN_OPTS) $(GEN_HDRS)
	@echo 'Building file: $<'
	@echo 'Invoking: C2000 Compiler'
	"C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/bin/cl2000" -v28 -ml -mt --float_support=fpu32 -O2 --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --include_path="C:/ti/controlSUITE/libs/math/IQmath/v15c/include" --include_path="C:/ti/controlSUITE/libs/math/FPUfastRTS/V100/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_common/include" --include_path="C:/ti/controlSUITE/device_support/c2834x/v112/DSP2834x_headers/include" --include_path="C:/ti/ccsv6/tools/compiler/ti-cgt-c2000_6.4.6/include" --advice:performance=all --diag_warning=225 --diag_warning=255 --display_error_number --diag_wrap=off --gen_func_subsections=on --preproc_with_compile --preproc_dependency="main.pp" $(GEN_OPTS__FLAG) "$<"
	@echo 'Finished building: $<'
	@echo ' '


