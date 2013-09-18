
MY_OS = $(shell uname)
ifeq ($(origin local_makefile), undefined)
local_makefile = make/Makefile.$(MY_OS)
endif
include $(shell echo ${local_makefile})

