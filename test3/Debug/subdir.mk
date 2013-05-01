################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../asd.cu 

CU_DEPS += \
./asd.d 

OBJS += \
./asd.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	nvcc -O3 -maxrregcount 30 --use_fast_math -ftz true -v -Xcompiler -O3 -Xptxas -v -gencode arch=compute_30,code=sm_30 -odir "" -M -o "$(@:%.o=%.d)" "$<"
	nvcc --compile --use_fast_math -O3 -Xcompiler -O3 -Xptxas -v -ftz true -gencode arch=compute_30,code=compute_30 -gencode arch=compute_30,code=sm_30 -maxrregcount 30 -v -Xcompiler -fopenmp -lgomp  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


