// Program to display the properties of all the NVIDIA-capable GPUs present
// on this machine
#include <stdio.h>
#include "cuda.h"
#include "cuda_runtime.h"

const char* GetComputeModeStr(int computeMode);
const char* GetAsyncEngineDescription(int value);


int main(int argc, char **argv) {
    // Retrieve the number of CUDA devices in this machine
    int numDevices;
    cudaGetDeviceCount(&numDevices);
    printf("There are %d CUDA devices\n", numDevices);
    if (numDevices == 0) {
        return 0;
    }

    // Show the properties of each device found
    for (int i=0; i < numDevices; i++) {
        cudaDeviceProp properties;
        cudaGetDeviceProperties(&properties, i);
        printf("\nDevice %d name '%s'\n", i, properties.name);
        printf("   Major.Minor: %d.%d\n", properties.major, properties.minor);
        printf("   Maximum global memory size: %lu bytes\n", properties.totalGlobalMem);
        printf("   Maximum constant memory size: %lu bytes\n", properties.totalConstMem);
        printf("   Maximum shared memory size per block: %lu bytes\n", properties.sharedMemPerBlock);
        printf("   Maximum thread block dimensions: %d x %d x %d\n",
            properties.maxThreadsDim[0], properties.maxThreadsDim[1], properties.maxThreadsDim[2]);
        printf("   Maximum grid dimensions: %d x %d x %d\n",
            properties.maxGridSize[0], properties.maxGridSize[1], properties.maxGridSize[2]);
        printf("   32-bit registers per block: %d\n", properties.regsPerBlock);
        printf("   Maximum number of threads per block: %d\n", properties.maxThreadsPerBlock);
        printf("   Clock rate: %d kHz\n", properties.clockRate);
        printf("   Peak memory clock rate: %d kHz\n", properties.memoryClockRate);
        printf("   Number of multi-processors: %d\n", properties.multiProcessorCount);
        printf("   Maximum number of threads per multi-processor: %d\n", properties.maxThreadsPerMultiProcessor);
        printf("   Kernel execution timeout: %s\n", properties.kernelExecTimeoutEnabled ? "ENABLED" : "DISABLED");
        printf("   Is GPU integrated in the motherboard: %s\n", properties.integrated ? "YES" : "NO");
        printf("   Can map host memory: %s\n", properties.canMapHostMemory ? "YES" : "NO");
        printf("   GPU compute mode: %s\n", GetComputeModeStr(properties.computeMode));
        printf("   Support multiple kernel execution concurrently: %s\n", properties.concurrentKernels ? "YES" : "FALSE");
        printf("   ECC status: %s\n", properties.ECCEnabled ? "ENABLED" : "DISABLED");
        // printf("   PCI bus id: %d\n", properties.pciBusId);
        // printf("   PCI device id: %d\n", properties.pciDeviceId);
        // printf("   PCI domain id: %d\n", properties.pciDomainId);
        printf("   TCC driver: %s\n", properties.tccDriver ? "YES" : "NO");
        printf("   Async engine capability: %s\n", GetAsyncEngineDescription(properties.asyncEngineCount));
        printf("   Shares address space with host: %s\n", properties.unifiedAddressing ? "YES" : "NO");
        printf("   Memory bus width: %d bits\n", properties.memoryBusWidth);
        printf("   L2 cache size: %d bytes\n", properties.l2CacheSize);
        printf("   Warp size: %d\n", properties.warpSize);
        printf("   Device overlap: %s\n", properties.deviceOverlap ? "YES" : "NO");
    }
}


const char* GetComputeModeStr(int computeMode)
{
    switch (computeMode) {
        case cudaComputeModeDefault:
            return "Default";

        case cudaComputeModeExclusive:
            return "Exclusive";

        case cudaComputeModeProhibited:
            return "Prohibited";

        case cudaComputeModeExclusiveProcess:
            return "ExclusiveProcess";

        default:
            break;
    }
    return "Unknown";
}


const char* GetAsyncEngineDescription(int value)
{
    switch (value) {
        case 0:
            return "concurrent memory copy and kernel execution is not supported";

        case 1:
            return "GPU can copy memory (single direction) while executing kernel";

        case 2:
            return "GPU can copy memory (both directions) while executing kernel";
    }
    return "Unknown";
}

