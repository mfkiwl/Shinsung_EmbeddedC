#
_XDCBUILDCOUNT = 
ifneq (,$(findstring path,$(_USEXDCENV_)))
override XDCPATH = C:/ti/bios_6_46_05_55/packages;C:/ti/ccsv6/ccs_base
override XDCROOT = C:/ti/xdctools_3_31_00_24_core
override XDCBUILDCFG = ./config.bld
endif
ifneq (,$(findstring args,$(_USEXDCENV_)))
override XDCARGS = 
override XDCTARGETS = 
endif
#
ifeq (0,1)
PKGPATH = C:/ti/bios_6_46_05_55/packages;C:/ti/ccsv6/ccs_base;C:/ti/xdctools_3_31_00_24_core/packages;..
HOSTOS = Windows
endif
