# expect iRRAM installation root directory to be $(IRRAM)

# settings, individual section

#IRRAM           = /home/info04/brausse/IVP-tests/iRRAM_2013_01/installed
#PKG_CONFIG_PATH = ~/bin/installed/lib64/pkgconfig
include Makefile.paths

# g++ <= 4.7 cannot handle the syntax properly
#CXX = clang++

PKG_CONFIG      = PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config

override CFLAGS   += -O2 -Wall -DNDEBUG -pg #--coverage #-Wextra -pedantic
override CFLAGS   += -Wextra -Wno-unused-parameter -Wno-unused-function -Wno-switch
override CXXFLAGS := $(CXXFLAGS) -std=c++11 $(CFLAGS) -Wno-tautological-compare
override CFLAGS   := -std=c11 $(CFLAGS)

EXES      = \
	ivp \
	pendulum \
	nbody \

ivp_OBJS     = ivp.o
ivp_LDFLAGS  = -L $(IRRAM)/lib -Wl,-rpath -Wl,$(IRRAM)/lib -pg
ivp_LDLIBS   = -lstdc++ -lm -lmpfr -lgmp -liRRAM # -lgcov
ivp_CXXFLAGS = -I $(IRRAM)/include #-pg

ifneq ($(PICARD),)
ivp_CXXFLAGS += -DMETHOD_PICARD=$(PICARD)
endif

pendulum_OBJS    = pendulum-vis.o
pendulum_PKGS    = cairo sdl
pendulum_CFLAGS  =
pendulum_LDFLAGS =
pendulum_LDLIBS  = -lm

nbody_OBJS    = nbody-vis.o ring-buf.o
nbody_PKGS    = cairo sdl
nbody_CFLAGS  = #-pg
nbody_LDFLAGS = #-pg
nbody_LDLIBS  = -lm

# rules, generic section

define EXE_template
ifneq "$$($(1)_PKGS)" ""
$(1)_PKG_CONFIG := $$(PKG_CONFIG) $$($(1)_PKGS)
else
$(1)_PKG_CONFIG := \#
endif
$$($(1)_OBJS): override CFLAGS   += $$($(1)_CFLAGS) `$$($(1)_PKG_CONFIG) --cflags`
$$($(1)_OBJS): override CXXFLAGS += $$($(1)_CXXFLAGS) `$$($(1)_PKG_CONFIG) --cflags`
$(1): override LDFLAGS += $$($(1)_LDFLAGS) `$$($(1)_PKG_CONFIG) --libs-only-L --libs-only-other`
$(1): override LDLIBS  += $$($(1)_LDLIBS) `$$($(1)_PKG_CONFIG) --libs-only-l`
$(1): $$($(1)_OBJS)
	$$(CC) $$(LDFLAGS) $$^ $$(LDLIBS) $$(OUTPUT_OPTION)
.PHONY: $(1)-clean $(1)-prof
$(1)-clean:
	$(RM) $(1) $$($(1)_OBJS)
$(1)-prof: override $(1)_CFLAGS += -pg
$(1)-prof: override $(1)_CXXFLAGS += -pg
$(1)-prof: override $(1)_LDFLAGS += -pg
$(1)-prof: $(1)
OBJS += $$($(1)_OBJS)
endef

.PHONY: all clean

all: ivp $(EXES)

$(foreach exe,$(EXES),$(eval $(call EXE_template,$(exe))))

$(OBJS): %.o: $(wildcard *.h)

clean:
	$(RM) $(OBJS) $(EXES)
