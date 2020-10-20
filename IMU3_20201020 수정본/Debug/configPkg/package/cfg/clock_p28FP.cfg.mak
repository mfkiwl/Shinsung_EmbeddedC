# invoke SourceDir generated makefile for clock.p28FP
clock.p28FP: .libraries,clock.p28FP
.libraries,clock.p28FP: package/cfg/clock_p28FP.xdl
	$(MAKE) -f D:\TI_project\IMU3/src/makefile.libs

clean::
	$(MAKE) -f D:\TI_project\IMU3/src/makefile.libs clean

