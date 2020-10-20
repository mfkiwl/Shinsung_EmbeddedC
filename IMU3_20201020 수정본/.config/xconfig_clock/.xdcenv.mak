#
_XDCBUILDCOUNT = 0
ifneq (,$(findstring path,$(_USEXDCENV_)))
override XDCPATH = C:/ti/bios_6_46_05_55/packages;C:/ti/ccsv6/ccs_base;D:/Source/4.WORK/11.IMU3/.config
override XDCROOT = C:/ti/xdctools_3_32_00_06_core
override XDCBUILDCFG = ./config.bld
endif
ifneq (,$(findstring args,$(_USEXDCENV_)))
override XDCARGS = 
override XDCTARGETS = 
endif
#
ifeq (0,1)
PKGPATH = C:/ti/bios_6_46_05_55/packages;C:/ti/ccsv6/ccs_base;D:/Source/4.WORK/11.IMU3/.config;C:/ti/xdctools_3_32_00_06_core/packages;..
HOSTOS = Windows
endif
